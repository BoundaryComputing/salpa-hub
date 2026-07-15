"""
gen-ali core — pure Python logic, no BoCoFlow dependencies.

Generates per-chain alignment files (.seq and .ali) for homology modeling
tools (ProMod3 / MODELLER). Combines present residues from a PDB structure
with missing residues (from PDB header) to produce template/target alignment
pairs.

Alignment format (.ali):
  >P1;case_chainA
  structure:Merge:start_resid:chain:+count:chain::::
  MVLM---------ALRM*

  >P1;case_chainA_full
  sequence:::::::::
  MVLMAAAALRMALRMMM*

Where '-' marks positions of residues missing in the template structure.
"""

import csv
import os
import textwrap
from dataclasses import dataclass, field

from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import PDBParser


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class PresentResidue:
    """A residue present in the PDB structure."""

    chain: str
    resid: int
    res_name_three: str
    res_name_one: str


@dataclass
class MissingResidue:
    """A residue missing from the PDB structure (from header or CSV)."""

    chain: str
    ssseq: int
    res_name: str
    one_letter: str


@dataclass
class ChainAlignmentResult:
    """Result of alignment generation for a single chain."""

    chain_id: str
    chain_type: str  # "P1" (protein) or "DL" (DNA/RNA)
    template_seq: str  # Present residues + '-' for missing
    full_seq: str  # Complete sequence (no gaps)
    num_present: int
    num_missing: int
    start_resid: int
    seq_agree: bool  # Whether reconstructed seq matches FASTA
    seq_file: str = ""  # Path to .seq file
    ali_file: str = ""  # Path to .ali file


# ---------------------------------------------------------------------------
# PDB residue extraction
# ---------------------------------------------------------------------------

# DNA/RNA nucleotide codes (2-letter PDB residue names)
_NUCLEOTIDE_MAP = {
    "DA": "A", "DC": "C", "DG": "G", "DT": "T",  # DNA
    "A": "A", "C": "C", "G": "G", "U": "U",  # RNA
}


def extract_present_residues(pdb_path: str, case_name: str = "structure") -> dict:
    """Extract residues present in a PDB structure, per chain.

    Args:
        pdb_path: Path to PDB file.
        case_name: Identifier for the structure.

    Returns:
        dict of chain_id -> list of PresentResidue, sorted by resid.
    """
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(case_name, pdb_path)

    chains = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            residues = []

            for residue in chain.get_residues():
                hetflag = residue.get_full_id()[3][0]
                resname = residue.get_resname().strip()
                resid = residue.get_full_id()[3][1]

                # Skip water and non-standard HETATM
                if hetflag not in (" ", ""):
                    continue

                # DNA/RNA: 2-letter codes
                if len(resname) <= 2 and resname in _NUCLEOTIDE_MAP:
                    one_letter = _NUCLEOTIDE_MAP[resname]
                elif resname in protein_letters_3to1:
                    one_letter = protein_letters_3to1[resname]
                else:
                    continue  # Skip unknown residues

                residues.append(PresentResidue(
                    chain=chain_id,
                    resid=resid,
                    res_name_three=resname,
                    res_name_one=one_letter,
                ))

            if residues:
                residues.sort(key=lambda r: r.resid)
                chains[chain_id] = residues
        break  # First model only

    return chains


def extract_missing_from_pdb(pdb_path: str, case_name: str = "structure") -> dict:
    """Extract missing residues from PDB header.

    Args:
        pdb_path: Path to PDB file.
        case_name: Identifier for the structure.

    Returns:
        dict of chain_id -> list of MissingResidue, sorted by ssseq.
    """
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(case_name, pdb_path)

    missing = {}
    header = structure.header
    if not header.get("has_missing_residues"):
        return missing

    for entry in header.get("missing_residues", []):
        chain_id = entry.get("chain", "?")
        res_name = entry.get("res_name", "UNK").strip()
        ssseq = entry.get("ssseq", 0)

        # DNA/RNA check
        if len(res_name) <= 2 and res_name in _NUCLEOTIDE_MAP:
            one_letter = _NUCLEOTIDE_MAP[res_name]
        else:
            one_letter = protein_letters_3to1.get(res_name, "X")

        mr = MissingResidue(
            chain=chain_id,
            ssseq=ssseq,
            res_name=res_name,
            one_letter=one_letter,
        )
        missing.setdefault(chain_id, []).append(mr)

    for chain_id in missing:
        missing[chain_id].sort(key=lambda r: r.ssseq)

    return missing


# ---------------------------------------------------------------------------
# CSV reading (from pdb_fasta_biopython output)
# ---------------------------------------------------------------------------

def read_missing_residues_csv(csv_path: str) -> list:
    """Read missing residues from a CSV file (pdb_fasta_biopython output).

    Expected columns: chain, res_name, one_letter, ssseq, model, insertion

    Returns:
        list of MissingResidue objects, sorted by ssseq.
    """
    residues = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            residues.append(MissingResidue(
                chain=row["chain"],
                ssseq=int(row["ssseq"]),
                res_name=row["res_name"],
                one_letter=row["one_letter"],
            ))
    residues.sort(key=lambda r: r.ssseq)
    return residues


# ---------------------------------------------------------------------------
# FASTA reading (for validation)
# ---------------------------------------------------------------------------

def read_fasta_sequence(fasta_path: str) -> str:
    """Read the first sequence from a FASTA file.

    Returns:
        Sequence string (no header, no newlines).
    """
    sequence_lines = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence_lines:
                    break  # Only first record
                continue
            sequence_lines.append(line)
    return "".join(sequence_lines)


# ---------------------------------------------------------------------------
# Alignment generation
# ---------------------------------------------------------------------------

def generate_chain_alignment(
    present: list,
    missing: list,
    fasta_sequence: str = None,
    append_end: bool = True,
) -> ChainAlignmentResult:
    """Generate alignment for a single chain.

    Merges present and missing residues sorted by residue number to produce:
    - template_seq: present residues at their positions, '-' for missing
    - full_seq: complete sequence (all residues)

    Args:
        present: list of PresentResidue for this chain.
        missing: list of MissingResidue for this chain.
        fasta_sequence: Optional full FASTA sequence for validation.
        append_end: Whether to append '*' terminator.

    Returns:
        ChainAlignmentResult with template and full sequences.
    """
    if not present:
        raise ValueError("No present residues provided")

    chain_id = present[0].chain

    # Determine chain type from present residues
    is_nucleic = any(len(r.res_name_three) <= 2 for r in present)
    chain_type = "DL" if is_nucleic else "P1"

    # Build position map: resid -> (one_letter, is_present)
    positions = {}
    for r in present:
        positions[r.resid] = (r.res_name_one, True)
    for r in missing:
        if r.ssseq not in positions:  # Don't overwrite present
            positions[r.ssseq] = (r.one_letter, False)

    # Sort by residue number
    sorted_resids = sorted(positions.keys())

    full_seq_chars = []
    template_seq_chars = []
    for resid in sorted_resids:
        one_letter, is_present = positions[resid]
        full_seq_chars.append(one_letter)
        template_seq_chars.append(one_letter if is_present else "-")

    full_seq = "".join(full_seq_chars)
    template_seq = "".join(template_seq_chars)

    # Validate against FASTA if provided
    seq_agree = True
    if fasta_sequence:
        seq_agree = full_seq == fasta_sequence

    # Append end marker
    if append_end:
        full_seq += "*"
        template_seq += "*"

    start_resid = sorted_resids[0] if sorted_resids else 0
    num_present = sum(1 for _, is_present in positions.values() if is_present)
    num_missing = sum(1 for _, is_present in positions.values() if not is_present)

    return ChainAlignmentResult(
        chain_id=chain_id,
        chain_type=chain_type,
        template_seq=template_seq,
        full_seq=full_seq,
        num_present=num_present,
        num_missing=num_missing,
        start_resid=start_resid,
        seq_agree=seq_agree,
    )


# ---------------------------------------------------------------------------
# File writers
# ---------------------------------------------------------------------------

def write_seq_file(
    output_path: str,
    case_name: str,
    chain_id: str,
    full_sequence: str,
    chain_type: str = "P1",
) -> str:
    """Write raw FASTA alignment .seq file.

    Format:
        >P1;case_chainA
        sequence:::::::::
        MVLMAALRM...*

    Args:
        output_path: Directory to write the file.
        case_name: Case identifier.
        chain_id: Chain identifier.
        full_sequence: Complete sequence (with '*' terminator if desired).
        chain_type: "P1" for protein, "DL" for DNA.

    Returns:
        Path to written file.
    """
    os.makedirs(output_path, exist_ok=True)
    file_path = os.path.join(output_path, "raw_fasta_record.seq")

    first_line = f">{chain_type};{case_name}{chain_id}"
    second_line = "sequence:::::::::"
    third_line = textwrap.fill(full_sequence, width=60)

    with open(file_path, "w") as f:
        f.write(f"{first_line}\n{second_line}\n{third_line}\n")

    return file_path


def write_ali_file(
    output_path: str,
    case_name: str,
    chain_id: str,
    template_seq: str,
    full_seq: str,
    start_resid: int,
    num_present: int,
    chain_type: str = "P1",
) -> str:
    """Write MODELLER/ProMod3 alignment .ali file.

    Format:
        >P1;case_chainA
        structure:Merge:start:chain:+count:chain::::
        MVLM---------ALRM*

        >P1;case_chainA_full
        sequence:::::::::
        MVLMAAAALRMALRM*

    Args:
        output_path: Directory to write the file.
        case_name: Case identifier.
        chain_id: Chain identifier.
        template_seq: Sequence with '-' for missing residues.
        full_seq: Complete sequence (no gaps).
        start_resid: First residue number in PDB.
        num_present: Number of present residues.
        chain_type: "P1" for protein, "DL" for DNA.

    Returns:
        Path to written file.
    """
    os.makedirs(output_path, exist_ok=True)
    file_path = os.path.join(output_path, "homology.ali")

    # Template entry (structure with gaps)
    first_line = f">{chain_type};{case_name}{chain_id}"
    second_line = ":".join([
        "structure",
        "Merge",
        str(start_resid),
        chain_id,
        f"+{num_present}",
        chain_id,
        "", "", "", "",
    ])
    # break_on_hyphens=False: '-' is a gap character, not a word break
    third_line = textwrap.fill(template_seq, width=60, break_on_hyphens=False)

    # Target entry (full sequence)
    first_line_2nd = f">{chain_type};{case_name}{chain_id}_full"
    second_line_2nd = "sequence:::::::::"
    third_line_2nd = textwrap.fill(full_seq, width=60, break_on_hyphens=False)

    with open(file_path, "w") as f:
        f.write(f"{first_line}\n{second_line}\n{third_line}\n\n")
        f.write(f"{first_line_2nd}\n{second_line_2nd}\n{third_line_2nd}\n")

    return file_path


# ---------------------------------------------------------------------------
# High-level processing
# ---------------------------------------------------------------------------

def process_chain(
    case_name: str,
    chain_id: str,
    output_dir: str,
    present_residues: list,
    missing_residues: list,
    fasta_sequence: str = None,
    append_end: bool = True,
) -> ChainAlignmentResult:
    """Process a single chain: generate alignment and write files.

    Args:
        case_name: Case identifier.
        chain_id: Chain identifier.
        output_dir: Base output directory (chain subfolder will be created).
        present_residues: list of PresentResidue.
        missing_residues: list of MissingResidue.
        fasta_sequence: Optional FASTA sequence for validation.
        append_end: Whether to append '*' terminator.

    Returns:
        ChainAlignmentResult with file paths populated.
    """
    chain_dir = os.path.join(output_dir, chain_id)

    alignment = generate_chain_alignment(
        present_residues,
        missing_residues,
        fasta_sequence=fasta_sequence,
        append_end=append_end,
    )

    seq_path = write_seq_file(
        chain_dir, case_name, chain_id, alignment.full_seq, alignment.chain_type,
    )
    ali_path = write_ali_file(
        chain_dir, case_name, chain_id,
        alignment.template_seq, alignment.full_seq,
        alignment.start_resid, alignment.num_present, alignment.chain_type,
    )

    alignment.seq_file = seq_path
    alignment.ali_file = ali_path

    return alignment


def process_all_chains(
    pdb_path: str,
    output_dir: str,
    case_name: str,
    chain_ids: list = None,
    missing_csv_paths: dict = None,
    fasta_paths: dict = None,
    append_end: bool = True,
) -> dict:
    """Process all chains in a PDB: extract data, generate alignments, write files.

    Args:
        pdb_path: Path to PDB file.
        output_dir: Base output directory.
        case_name: Case identifier.
        chain_ids: Optional list of chain IDs to process (default: all).
        missing_csv_paths: Optional dict of chain_id -> CSV path (from predecessor).
        fasta_paths: Optional dict of chain_id -> FASTA path (for validation).
        append_end: Whether to append '*' terminator.

    Returns:
        dict with keys:
            - chain_results: dict of chain_id -> ChainAlignmentResult
            - seq_agree_all: bool (True if all chains agree)
    """
    # Extract present residues from PDB
    present_by_chain = extract_present_residues(pdb_path, case_name)

    # Get missing residues: from CSVs if provided, otherwise from PDB header
    missing_by_chain = {}
    if missing_csv_paths:
        for chain_id, csv_path in missing_csv_paths.items():
            if os.path.exists(csv_path):
                missing_by_chain[chain_id] = read_missing_residues_csv(csv_path)
    else:
        raw_missing = extract_missing_from_pdb(pdb_path, case_name)
        for chain_id, residues in raw_missing.items():
            missing_by_chain[chain_id] = residues

    # Determine which chains to process
    if chain_ids is None:
        chain_ids = sorted(present_by_chain.keys())

    # Read FASTA sequences for validation
    fasta_seqs = {}
    if fasta_paths:
        for chain_id, fasta_path in fasta_paths.items():
            if os.path.exists(fasta_path):
                fasta_seqs[chain_id] = read_fasta_sequence(fasta_path)

    # Process each chain
    chain_results = {}
    for chain_id in chain_ids:
        if chain_id not in present_by_chain:
            continue

        result = process_chain(
            case_name=case_name,
            chain_id=chain_id,
            output_dir=output_dir,
            present_residues=present_by_chain[chain_id],
            missing_residues=missing_by_chain.get(chain_id, []),
            fasta_sequence=fasta_seqs.get(chain_id),
            append_end=append_end,
        )
        chain_results[chain_id] = result

    seq_agree_all = all(r.seq_agree for r in chain_results.values())

    return {
        "chain_results": chain_results,
        "seq_agree_all": seq_agree_all,
    }
