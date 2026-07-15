"""
pdb-fasta-biopython core — pure Python logic, no BoCoFlow dependencies.

Provides:
- RCSB PDB API fetching (PDB + FASTA by ID)
- PDB structure parsing via BioPython
- Sequence extraction with chain info
- Missing residues extraction from PDB header
- FASTA file writing (combined + per-chain)
- Missing residues CSV output
"""

import csv
import os
from dataclasses import dataclass, field

import requests
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqIO import write as seqio_write
from Bio.SeqRecord import SeqRecord

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}/download"


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ChainInfo:
    """Sequence and metadata for a single chain."""

    chain_id: str
    sequence: str
    length: int
    protein_residues: int = 0
    dna_rna_residues: int = 0
    chain_type: str = "protein"  # "protein" or "nucleic_acid"


@dataclass
class MissingResidue:
    """A single missing residue entry from the PDB header."""

    chain: str
    res_name: str  # 3-letter code
    one_letter: str  # 1-letter code
    ssseq: int  # sequence number
    model: str = ""
    insertion: str = ""


@dataclass
class ParseResult:
    """Result of parsing a PDB structure."""

    chains: dict = field(default_factory=dict)  # chain_id -> ChainInfo
    missing_residues: dict = field(default_factory=dict)  # chain_id -> [MissingResidue]


# ---------------------------------------------------------------------------
# RCSB API
# ---------------------------------------------------------------------------

def fetch_pdb_from_rcsb(pdb_id: str) -> str:
    """Fetch PDB file content from RCSB by PDB ID.

    Args:
        pdb_id: 4-character PDB identifier (e.g. "3LZ0").

    Returns:
        PDB file content as text.

    Raises:
        ValueError: If the PDB ID is not found or the request fails.
    """
    url = RCSB_PDB_URL.format(pdb_id=pdb_id.upper())
    resp = requests.get(url, timeout=30)
    if resp.status_code != 200:
        raise ValueError(
            f"Failed to fetch PDB {pdb_id} from RCSB (HTTP {resp.status_code})"
        )
    return resp.text


def fetch_fasta_from_rcsb(pdb_id: str) -> str:
    """Fetch FASTA sequence from RCSB by PDB ID.

    Args:
        pdb_id: 4-character PDB identifier.

    Returns:
        FASTA file content as text.

    Raises:
        ValueError: If the PDB ID is not found or the request fails.
    """
    url = RCSB_FASTA_URL.format(pdb_id=pdb_id.upper())
    resp = requests.get(url, timeout=30)
    if resp.status_code != 200:
        raise ValueError(
            f"Failed to fetch FASTA for {pdb_id} from RCSB (HTTP {resp.status_code})"
        )
    return resp.text


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

def parse_pdb_structure(pdb_source: str, case_name: str, is_file: bool = True):
    """Parse a PDB file into a BioPython Structure object.

    Args:
        pdb_source: File path (if is_file=True) or PDB text content.
        case_name: Identifier for the structure.
        is_file: If True, pdb_source is a file path; otherwise raw text.

    Returns:
        BioPython Structure object.
    """
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    if is_file:
        return parser.get_structure(case_name, pdb_source)
    else:
        # Write to temp file for parser
        import tempfile

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False
        ) as tmp:
            tmp.write(pdb_source)
            tmp_path = tmp.name
        try:
            return parser.get_structure(case_name, tmp_path)
        finally:
            os.unlink(tmp_path)


def extract_sequences(structure, include_hetatm: bool = False) -> dict:
    """Extract amino acid / nucleic acid sequences from a BioPython Structure.

    Args:
        structure: BioPython Structure object.
        include_hetatm: If True, include HETATM records (marked as 'X').

    Returns:
        dict of chain_id -> ChainInfo.
    """
    chains = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            sequence = []
            protein_residues = 0
            dna_rna_residues = 0

            for residue in chain.get_residues():
                hetflag = residue.get_full_id()[3][0]
                resname = residue.get_resname()

                # Skip HETATM unless requested
                if not include_hetatm and hetflag not in (" ", ""):
                    continue

                # DNA/RNA: 2-letter codes like "DA", "DC", "DG", "DT"
                if len(resname.strip()) == 2:
                    sequence.append(resname.strip()[1])
                    dna_rna_residues += 1
                # Standard amino acids
                elif resname in protein_letters_3to1:
                    sequence.append(protein_letters_3to1[resname])
                    protein_residues += 1
                # Non-standard residues
                elif include_hetatm:
                    sequence.append("X")

            if sequence:
                seq_str = "".join(sequence)
                chains[chain_id] = ChainInfo(
                    chain_id=chain_id,
                    sequence=seq_str,
                    length=len(seq_str),
                    protein_residues=protein_residues,
                    dna_rna_residues=dna_rna_residues,
                    chain_type=(
                        "protein" if protein_residues > dna_rna_residues
                        else "nucleic_acid"
                    ),
                )
    return chains


def extract_missing_residues(structure) -> dict:
    """Extract missing residues from PDB header.

    Args:
        structure: BioPython Structure object (must have header info).

    Returns:
        dict of chain_id -> list of MissingResidue.
    """
    missing = {}
    header = structure.header
    if not header.get("has_missing_residues"):
        return missing

    for entry in header.get("missing_residues", []):
        chain_id = entry.get("chain", "?")
        res_name = entry.get("res_name", "UNK")
        one_letter = protein_letters_3to1.get(res_name, "X")

        mr = MissingResidue(
            chain=chain_id,
            res_name=res_name,
            one_letter=one_letter,
            ssseq=entry.get("ssseq", 0),
            model=str(entry.get("model", "")),
            insertion=entry.get("insertion", ""),
        )
        missing.setdefault(chain_id, []).append(mr)

    return missing


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_fasta_files(
    chains: dict,
    output_dir: str,
    case_name: str,
    split_chains: bool = True,
    pdb_basename: str = "",
) -> dict:
    """Write FASTA files from extracted chain sequences.

    Args:
        chains: dict of chain_id -> ChainInfo.
        output_dir: Directory to write FASTA files.
        case_name: Case identifier for file naming.
        split_chains: If True, write one FASTA per chain.
        pdb_basename: Original PDB filename for descriptions.

    Returns:
        dict of output file labels -> file paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_files = {}

    if split_chains:
        for chain_id, info in chains.items():
            fasta_path = os.path.join(
                output_dir, f"chain_{chain_id}.fasta"
            )
            record = SeqRecord(
                Seq(info.sequence),
                id=f"{case_name}_chain_{chain_id}",
                description=(
                    f"Chain {chain_id} from {pdb_basename or case_name} - "
                    f"{info.chain_type} - {info.length} residues"
                ),
            )
            seqio_write(record, fasta_path, "fasta")
            output_files[f"chain_{chain_id}"] = fasta_path
    else:
        combined_path = os.path.join(output_dir, "combined.fasta")
        records = []
        for chain_id, info in chains.items():
            records.append(
                SeqRecord(
                    Seq(info.sequence),
                    id=f"{case_name}_chain_{chain_id}",
                    description=(
                        f"Chain {chain_id} - {info.chain_type} - "
                        f"{info.length} residues"
                    ),
                )
            )
        seqio_write(records, combined_path, "fasta")
        output_files["combined"] = combined_path

    return output_files


def write_missing_residues_csv(missing: dict, output_dir: str) -> dict:
    """Write missing residues to CSV files (one per chain).

    Args:
        missing: dict of chain_id -> list of MissingResidue.
        output_dir: Directory to write CSV files.

    Returns:
        dict of chain_id -> CSV file path.
    """
    os.makedirs(output_dir, exist_ok=True)
    csv_files = {}

    for chain_id, residues in missing.items():
        csv_path = os.path.join(
            output_dir, f"missing_residues_chain_{chain_id}.csv"
        )
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                ["chain", "res_name", "one_letter", "ssseq", "model", "insertion"]
            )
            for mr in residues:
                writer.writerow(
                    [mr.chain, mr.res_name, mr.one_letter, mr.ssseq,
                     mr.model, mr.insertion]
                )
        csv_files[chain_id] = csv_path

    return csv_files
