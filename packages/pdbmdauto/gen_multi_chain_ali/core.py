"""
gen-multi-chain-ali core — pure Python logic, no BoCoFlow dependencies.

Merges per-chain .ali alignment files into a single multi-chain alignment
file for ProMod3 / MODELLER. Each chain's template and target sequences are
concatenated with '/' chain separators. DNA/RNA chains have their residue
letters replaced with '.' (MODELLER convention).

Multi-chain .ali format:
  >P1;case_name
  structure:merge::first_protein_chain::last_protein_chain::::
  MVLM---ALRM/........../GGCC---TTAA*

  >P1;case_name_full
  sequence:::::::::
  MVLMAAAALRM/........../GGCCAATTAA*
"""

import os
import textwrap
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ChainAliEntry:
    """Parsed content of one chain's .ali file."""

    chain_id: str
    chain_type: str  # "P1" (protein) or "DL" (DNA/RNA)
    template_seq: str  # Present residues + '-' for missing, ends with '*'
    full_seq: str  # Complete sequence, ends with '*'
    template_header: str  # e.g. "structure:Merge:1:A:+97:A::::"
    full_header: str  # e.g. "sequence:::::::::"


@dataclass
class MultiChainResult:
    """Result of multi-chain alignment merge."""

    ali_file: str  # Path to merged .ali file
    chain_types: dict = field(default_factory=dict)  # chain_id -> "P1" or "DL"
    protein_chains: list = field(default_factory=list)
    dna_chains: list = field(default_factory=list)
    selected_chains: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# .ali file parser
# ---------------------------------------------------------------------------

def parse_ali_file(ali_path: str, chain_id: str = "") -> ChainAliEntry:
    """Parse a per-chain .ali file into structured data.

    Expected format (two entries separated by blank line):
        >P1;case_chainA
        structure:Merge:1:A:+97:A::::
        MVLM---------ALRM*

        >P1;case_chainA_full
        sequence:::::::::
        MVLMAAAALRMALRM*

    Args:
        ali_path: Path to .ali file.
        chain_id: Chain identifier (used if not derivable from file).

    Returns:
        ChainAliEntry with parsed data.
    """
    with open(ali_path, "r") as f:
        content = f.read()

    # Split into two entries by blank line
    blocks = [b.strip() for b in content.split("\n\n") if b.strip()]
    if len(blocks) != 2:
        raise ValueError(
            f"Expected 2 alignment blocks in {ali_path}, got {len(blocks)}"
        )

    # Parse each block
    template_entry = _parse_ali_block(blocks[0])
    full_entry = _parse_ali_block(blocks[1])

    # Determine chain type from first header line
    chain_type = template_entry["type_code"]  # "P1" or "DL"

    return ChainAliEntry(
        chain_id=chain_id,
        chain_type=chain_type,
        template_seq=template_entry["sequence"],
        full_seq=full_entry["sequence"],
        template_header=template_entry["descriptor"],
        full_header=full_entry["descriptor"],
    )


def _parse_ali_block(block_text: str) -> dict:
    """Parse a single alignment block (header + descriptor + sequence lines).

    Returns:
        dict with keys: type_code, name, descriptor, sequence
    """
    lines = block_text.strip().split("\n")
    if len(lines) < 3:
        raise ValueError(f"Alignment block needs at least 3 lines, got {len(lines)}")

    # Line 1: >TYPE;name
    header = lines[0]
    if not header.startswith(">"):
        raise ValueError(f"Expected '>' header, got: {header}")
    type_and_name = header[1:]  # Remove '>'
    type_code, name = type_and_name.split(";", 1) if ";" in type_and_name else (type_and_name, "")

    # Line 2: descriptor (e.g. "structure:Merge:1:A:..." or "sequence:::::::::")
    descriptor = lines[1]

    # Lines 3+: sequence (may span multiple lines, wrapped at 60 chars)
    sequence = "".join(lines[2:])

    return {
        "type_code": type_code,
        "name": name,
        "descriptor": descriptor,
        "sequence": sequence,
    }


# ---------------------------------------------------------------------------
# Chain selection
# ---------------------------------------------------------------------------

def select_chains(
    available_chains: list,
    selection: str = "all",
) -> list:
    """Parse chain selection string into a list of chain IDs.

    Args:
        available_chains: All available chain IDs.
        selection: "all" or comma-separated chain IDs (e.g. "A,B,C").

    Returns:
        Sorted list of selected chain IDs.

    Raises:
        ValueError: If a requested chain is not available.
    """
    if selection.strip().lower() == "all":
        return sorted(available_chains)

    selected = [c.strip() for c in selection.split(",") if c.strip()]
    for chain_id in selected:
        if chain_id not in available_chains:
            raise ValueError(
                f"Chain '{chain_id}' not in available chains: {available_chains}"
            )
    return sorted(selected)


# ---------------------------------------------------------------------------
# Multi-chain merge
# ---------------------------------------------------------------------------

def _mask_dna_sequence(seq: str) -> str:
    """Replace DNA/RNA residue letters with '.' (MODELLER convention).

    Preserves '*' (terminator), '/' (chain separator), and '-' (gap).
    """
    return "".join(
        c if c in ("*", "/", "-") else "."
        for c in seq
    )


def merge_chain_alignments(
    chain_entries: list,
    selected_chains: list = None,
) -> tuple:
    """Merge per-chain alignment entries into multi-chain sequences.

    For non-last chains, the terminal '*' is replaced with '/' (chain separator).
    For DNA/RNA chains, residue characters are replaced with '.'.

    Args:
        chain_entries: list of ChainAliEntry, one per chain.
        selected_chains: Optional ordered list of chain IDs to include.
            If None, uses all entries in order.

    Returns:
        (merged_template_seq, merged_full_seq, chain_types) where
        chain_types is a dict of chain_id -> "P1" or "DL".
    """
    if selected_chains:
        by_id = {e.chain_id: e for e in chain_entries}
        entries = [by_id[c] for c in selected_chains if c in by_id]
    else:
        entries = chain_entries

    if not entries:
        raise ValueError("No chain entries to merge")

    merged_template = []
    merged_full = []
    chain_types = {}

    for i, entry in enumerate(entries):
        is_last = (i == len(entries) - 1)

        template = entry.template_seq
        full = entry.full_seq

        # Replace terminal '*' with '/' for non-last chains
        if not is_last:
            if template.endswith("*"):
                template = template[:-1] + "/"
            if full.endswith("*"):
                full = full[:-1] + "/"

        # Mask DNA/RNA residues with '.'
        if entry.chain_type == "DL":
            template = _mask_dna_sequence(template)
            full = _mask_dna_sequence(full)

        merged_template.append(template)
        merged_full.append(full)
        chain_types[entry.chain_id] = entry.chain_type

    return "".join(merged_template), "".join(merged_full), chain_types


# ---------------------------------------------------------------------------
# File writer
# ---------------------------------------------------------------------------

def write_multi_chain_ali(
    output_dir: str,
    case_name: str,
    merged_template: str,
    merged_full: str,
    chain_types: dict,
    filename: str = "homology.ali",
) -> str:
    """Write the merged multi-chain .ali file.

    Format:
        >P1;case_name
        structure:merge::first_protein_chain::last_protein_chain::::
        MERGED_TEMPLATE_SEQ*

        >P1;case_name_full
        sequence:::::::::
        MERGED_FULL_SEQ*

    Args:
        output_dir: Directory for the output file.
        case_name: Case identifier.
        merged_template: Merged template sequence (chains joined with '/').
        merged_full: Merged full sequence (chains joined with '/').
        chain_types: dict of chain_id -> "P1" or "DL".
        filename: Output filename (default: "homology.ali").

    Returns:
        Path to written file.
    """
    os.makedirs(output_dir, exist_ok=True)
    file_path = os.path.join(output_dir, filename)

    # Find first and last protein chains for header
    protein_chains = [c for c, t in chain_types.items() if t == "P1"]
    first_protein = protein_chains[0] if protein_chains else ""
    last_protein = protein_chains[-1] if protein_chains else ""

    # Validate raw sequence lengths before wrapping
    if len(merged_template) != len(merged_full):
        raise ValueError(
            f"Sequence length mismatch: "
            f"template={len(merged_template)}, full={len(merged_full)}"
        )

    # Template entry
    ali_code = "P1"
    first_line = f">{ali_code};{case_name}"
    second_line = ":".join([
        "structure", "merge", "", first_protein, "", last_protein,
        "", "", "", "",
    ])
    # Use break_on_hyphens=False — '-' is a gap character, not a word break
    third_line = textwrap.fill(merged_template, width=60, break_on_hyphens=False)

    # Target entry
    first_line_2nd = f">{ali_code};{case_name}_full"
    second_line_2nd = "sequence:::::::::"
    third_line_2nd = textwrap.fill(merged_full, width=60, break_on_hyphens=False)

    with open(file_path, "w") as f:
        f.write(f"{first_line}\n{second_line}\n{third_line}\n\n")
        f.write(f"{first_line_2nd}\n{second_line_2nd}\n{third_line_2nd}\n")

    return file_path


# ---------------------------------------------------------------------------
# High-level processing
# ---------------------------------------------------------------------------

def process_multi_chain_ali(
    ali_file_paths: dict,
    output_dir: str,
    case_name: str,
    selected_chains_str: str = "all",
    merge_folder_name: str = "Merge",
) -> MultiChainResult:
    """Process per-chain .ali files into a single multi-chain alignment.

    Args:
        ali_file_paths: dict of chain_id -> path to per-chain .ali file.
        output_dir: Base output directory.
        case_name: Case identifier.
        selected_chains_str: "all" or comma-separated chain IDs.
        merge_folder_name: Subfolder name for merged output.

    Returns:
        MultiChainResult with merged file path and chain metadata.
    """
    # Parse all .ali files
    chain_entries = []
    for chain_id in sorted(ali_file_paths.keys()):
        ali_path = ali_file_paths[chain_id]
        if not os.path.exists(ali_path):
            raise FileNotFoundError(f"ALI file not found: {ali_path}")
        entry = parse_ali_file(ali_path, chain_id=chain_id)
        chain_entries.append(entry)

    # Select chains
    available = [e.chain_id for e in chain_entries]
    selected = select_chains(available, selected_chains_str)

    # Merge
    merged_template, merged_full, chain_types = merge_chain_alignments(
        chain_entries, selected_chains=selected,
    )

    # Write output
    merge_dir = os.path.join(output_dir, merge_folder_name)
    ali_file = write_multi_chain_ali(
        merge_dir, case_name, merged_template, merged_full, chain_types,
    )

    protein_chains = [c for c, t in chain_types.items() if t == "P1"]
    dna_chains = [c for c, t in chain_types.items() if t == "DL"]

    return MultiChainResult(
        ali_file=ali_file,
        chain_types=chain_types,
        protein_chains=protein_chains,
        dna_chains=dna_chains,
        selected_chains=selected,
    )
