"""
merge-pdb-chains core — pure Python logic, no BoCoFlow dependencies.

Extracts selected chains from a PDB structure and writes them to a single
merged PDB file. Replaces the legacy PyMOL-based merge with BioPython's
Bio.PDB module.

For DNA/RNA chains, atoms are written as HETATM records (MODELLER/ProMod3
convention for non-protein chains in multi-chain homology modeling).
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

class ChainSelector(Select):
    """Select specific chains from a PDB structure."""

    def __init__(self, chain_ids: list, dna_chains: set = None):
        self.chain_ids = set(chain_ids)
        self.dna_chains = dna_chains or set()

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

    def accept_residue(self, residue):
        # Skip water and other hetero residues (ligands, ions, modified residues).
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

    io.set_structure(structure)
    selector = ChainSelector(selected_chains, dna_chains)
    io.save(output_path, selector)

    return output_path


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

    # Select chains
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
