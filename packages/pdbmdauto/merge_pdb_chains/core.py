"""
merge-pdb-chains core — pure Python logic, no BoCoFlow dependencies.

Extracts selected chains from a PDB structure and writes them to a single
merged PDB file. Replaces the legacy PyMOL-based merge with BioPython's
Bio.PDB module.

For DNA/RNA chains, atoms are written as HETATM records (MODELLER/ProMod3
convention for non-protein chains in multi-chain homology modeling).

A modified residue that the file's MODRES records declare (selenomethionine, MSE, above
all) is kept, written as its standard parent, because Generate Alignment rebuilds it as
that parent and ProMod3 requires the alignment and this structure to agree residue for
residue. Every other HETATM residue (ligands, ions, waters) is left out.
"""

import json
import os
from dataclasses import dataclass, field

from Bio.PDB import PDBIO, PDBParser, Select


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class MergeResult:
    """Result of PDB chain merge."""

    output_pdb: str  # Path to merged PDB file
    selected_chains: list = field(default_factory=list)
    chain_types: dict = field(default_factory=dict)
    chain_type_file: str = ""
    chain_name_file: str = ""


# ---------------------------------------------------------------------------
# BioPython chain selector
# ---------------------------------------------------------------------------

# The heavy atoms of each standard amino acid, by PDB atom name. A declared modified
# residue keeps the atoms its parent also has; the rest (a phosphate, a methyl, ...) go.
_AA_HEAVY_ATOMS = {
    "ALA": {"N", "CA", "C", "O", "CB"},
    "ARG": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"},
    "ASP": {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"},
    "CYS": {"N", "CA", "C", "O", "CB", "SG"},
    "GLN": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"},
    "GLU": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
    "GLY": {"N", "CA", "C", "O"},
    "HIS": {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"},
    "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
    "LYS": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
    "MET": {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"},
    "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"N", "CA", "C", "O", "CB", "CG", "CD"},
    "SER": {"N", "CA", "C", "O", "CB", "OG"},
    "THR": {"N", "CA", "C", "O", "CB", "OG1", "CG2"},
    "TRP": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3",
            "CZ2", "CZ3", "CH2"},
    "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "VAL": {"N", "CA", "C", "O", "CB", "CG1", "CG2"},
}

# Atoms renamed on the way to the parent: (modified residue, atom) -> (name, element).
_RENAMED_ATOMS = {("MSE", "SE"): ("SD", "S")}


def read_modres(pdb_path: str) -> dict:
    """The modified residues a PDB file's MODRES records declare, with their standard parents.

    The same reading as Generate Alignment's (gen_ali/core.py), which rebuilds these residues
    as their parents; the two must agree.

    Returns:
        dict of (chain_id, resid) -> (modified name, standard parent name),
        e.g. ("A", 113) -> ("MSE", "MET").
    """
    modres = {}
    with open(pdb_path, "r", errors="replace") as f:
        for line in f:
            if not line.startswith("MODRES"):
                continue
            try:
                resid = int(line[18:22])
            except ValueError:
                continue
            modres[(line[16:17], resid)] = (line[12:15].strip(), line[24:27].strip())
    return modres


def _kept_as_parent(modres: dict, chain_id: str, resid: int, resname: str):
    """The standard parent a declared modified residue is written as, or None."""
    declared = modres.get((chain_id, resid))
    if not declared or declared[0] != resname or declared[1] not in _AA_HEAVY_ATOMS:
        return None
    return declared[1]


class ChainSelector(Select):
    """Select specific chains from a PDB structure."""

    def __init__(self, chain_ids: list, dna_chains: set = None, modres: dict = None):
        self.chain_ids = set(chain_ids)
        self.dna_chains = dna_chains or set()
        self.modres = modres or {}

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

    def accept_residue(self, residue):
        # A modified residue the file's MODRES records declare stays: Generate Alignment
        # rebuilds it as its parent, and merge_chains() writes it as that parent. Left out,
        # the structure was one residue short of the alignment wherever one occurred, and
        # ProMod3 stopped with "Alignment-structure mismatch" (3I2V chain A, MSE 113).
        hetflag, resid, _icode = residue.id
        if hetflag.startswith("H_") and _kept_as_parent(
            self.modres, residue.get_parent().id, resid, residue.get_resname()
        ):
            return True
        # Skip water and other hetero residues (ligands, ions, undeclared modified residues).
        #
        # BioPython hetflag conventions (residue.id[0]):
        #   ' '   — standard amino acid or nucleotide (from PDB ATOM record)
        #   'W'   — water (HOH, DOD, etc.)
        #   'H_*' — other hetero (from HETATM: ligands, ions, cofactors,
        #           modified residues)
        #
        # pdbmdauto uses upstream pdb_fasta_biopython to extract sequences
        # from standard residues only; if we let hetero residues through
        # here, the merged PDB's chain iteration order would include
        # waters / ligands after the protein, and per-chain alignments
        # (from gen_ali) would not match the structure's residue count.
        # The downstream fix_residues_promod3 → pka_gmx_em → gmx_solv_ion
        # chain expects a protein-only structure and adds water + ions
        # fresh, so stripping here is correct for the MD-prep pipeline.
        return residue.id[0] == " "

    def accept_atom(self, atom):
        return True


class HetatomDnaWriter(PDBIO):
    """Custom PDBIO that writes DNA/RNA chain atoms as HETATM records.

    MODELLER/ProMod3 convention: non-protein chains use HETATM record type
    in multi-chain homology modeling alignment files.
    """

    def __init__(self, dna_chains: set = None):
        super().__init__()
        self._dna_chains = dna_chains or set()

    def _get_atom_line(self, atom, hetfield, segid, atom_number, resname,
                       resseq, icode, chain_id, charge="  "):
        """Override to convert ATOM→HETATM for DNA/RNA chains."""
        # For DNA chains, force HETATM record type
        if chain_id in self._dna_chains and hetfield == " ":
            hetfield = "H"
        return super()._get_atom_line(
            atom, hetfield, segid, atom_number, resname,
            resseq, icode, chain_id, charge
        )


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def select_chains(available_chains: list, selection: str = "all") -> list:
    """Parse chain selection string into a sorted list of chain IDs."""
    if selection.strip().lower() == "all":
        return sorted(available_chains)
    selected = [c.strip() for c in selection.split(",") if c.strip()]
    for chain_id in selected:
        if chain_id not in available_chains:
            raise ValueError(f"Chain '{chain_id}' not in available: {available_chains}")
    return sorted(selected)


def merge_chains(
    pdb_path: str,
    output_path: str,
    selected_chains: list,
    chain_types: dict = None,
    case_name: str = "structure",
) -> str:
    """Extract selected chains from a PDB and write to a merged file.

    Args:
        pdb_path: Path to input PDB file.
        output_path: Path for output merged PDB file.
        selected_chains: List of chain IDs to include.
        chain_types: Dict of chain_id → "P1" or "DL". DNA chains get HETATM.
        case_name: Structure identifier.

    Returns:
        Path to the merged PDB file.
    """
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(case_name, pdb_path)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Identify DNA chains for HETATM conversion
    dna_chains = set()
    if chain_types:
        for chain_id, ctype in chain_types.items():
            if ctype == "DL":
                dna_chains.add(chain_id)

    # Use custom writer if DNA chains need HETATM conversion
    if dna_chains:
        io = HetatomDnaWriter(dna_chains=dna_chains)
    else:
        io = PDBIO()

    modres = read_modres(pdb_path)
    io.set_structure(structure)
    selector = ChainSelector(selected_chains, dna_chains, modres)
    io.save(output_path, selector)
    if modres:
        write_modified_as_parents(output_path, modres)

    return output_path


def _atom_name_field(name: str, element: str) -> str:
    """A PDB atom-name field (columns 13-16): a one-letter element's name starts in 14."""
    if len(name) < 4 and len(element) == 1:
        return f" {name:<3}"
    return f"{name:<4}"


def write_modified_as_parents(pdb_path: str, modres: dict) -> int:
    """Rewrite, in place, each declared modified residue as its standard parent.

    Its HETATM records become ATOM records under the parent's name. Atoms the parent also
    has are kept, renamed where the parent calls them otherwise (selenomethionine's SE is
    methionine's SD, a sulfur); the rest, and any hydrogens, are dropped. Returns the number
    of atom records rewritten.
    """
    with open(pdb_path, "r") as f:
        lines = f.readlines()
    out, rewritten = [], 0
    for line in lines:
        if not line.startswith("HETATM"):
            out.append(line)
            continue
        try:
            resid = int(line[22:26])
        except ValueError:
            out.append(line)
            continue
        resname = line[17:20].strip()
        parent = _kept_as_parent(modres, line[21:22], resid, resname)
        if parent is None:
            out.append(line)
            continue
        name = line[12:16].strip()
        element = line[76:78].strip() if len(line) >= 78 else ""
        name, element = _RENAMED_ATOMS.get((resname, name), (name, element))
        if name not in _AA_HEAVY_ATOMS[parent] and name != "OXT":
            continue
        body = line.rstrip("\n").ljust(80)
        out.append(
            "ATOM  " + body[6:12] + _atom_name_field(name, element) + body[16:17]
            + f"{parent:>3}" + body[20:76] + f"{element:>2}" + body[78:].rstrip() + "\n"
        )
        rewritten += 1
    with open(pdb_path, "w") as f:
        f.writelines(out)
    return rewritten


def write_chain_metadata(
    output_dir: str,
    case_name: str,
    selected_chains: list,
    chain_types: dict,
) -> tuple:
    """Write chain type and chain name JSON files.

    Returns:
        (chain_type_file_path, chain_name_file_path)
    """
    os.makedirs(output_dir, exist_ok=True)

    chain_type_file = os.path.join(output_dir, "chain_type.json")
    with open(chain_type_file, "w") as f:
        json.dump(chain_types, f, indent=2)

    chain_name_file = os.path.join(output_dir, "chain_name.json")
    with open(chain_name_file, "w") as f:
        json.dump(selected_chains, f, indent=2)

    return chain_type_file, chain_name_file


def process_merge(
    pdb_path: str,
    output_dir: str,
    case_name: str,
    selected_chains_str: str = "all",
    chain_types: dict = None,
    merge_folder_name: str = "Merge",
    merge_file_name: str = "merge.pdb",
) -> MergeResult:
    """High-level: merge selected chains and write metadata.

    Args:
        pdb_path: Path to source PDB file.
        output_dir: Base output directory.
        case_name: Case identifier.
        selected_chains_str: "all" or comma-separated chain IDs.
        chain_types: Dict of chain_id → "P1"/"DL" (from gen_multi_chain_ali).
        merge_folder_name: Subfolder for merged output.
        merge_file_name: Output filename.

    Returns:
        MergeResult with file paths and chain metadata.
    """
    # Get available chains from PDB
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(case_name, pdb_path)
    available = []
    for model in structure:
        for chain in model:
            available.append(chain.id)
        break

    # Select chains. "all" means every chain the alignment covers when the alignment step said
    # which those are: a chain it set aside (a short peptide of residues it cannot build, such
    # as 6LU7's inhibitor, chain C) has no alignment, and merged here it would reach ProMod3's
    # template anyway, as 1.2.6's free tripeptide did.
    if chain_types and selected_chains_str.strip().lower() == "all":
        available = [c for c in available if c in chain_types]
    selected = select_chains(available, selected_chains_str)

    # Filter chain_types to selected chains only
    if chain_types:
        filtered_types = {c: t for c, t in chain_types.items() if c in selected}
    else:
        filtered_types = {c: "P1" for c in selected}

    # Merge
    merge_dir = os.path.join(output_dir, merge_folder_name)
    merge_path = os.path.join(merge_dir, merge_file_name)
    merge_chains(pdb_path, merge_path, selected, filtered_types, case_name)

    # Write metadata
    ct_file, cn_file = write_chain_metadata(merge_dir, case_name, selected, filtered_types)

    return MergeResult(
        output_pdb=merge_path,
        selected_chains=selected,
        chain_types=filtered_types,
        chain_type_file=ct_file,
        chain_name_file=cn_file,
    )
