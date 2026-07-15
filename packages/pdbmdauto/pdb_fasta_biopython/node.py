"""
pdb-fasta-biopython — BoCoFlow node wrapper.

Parses PDB structures and extracts FASTA sequences using BioPython.
Supports two input modes:
  - pdb_id: Fetch PDB + FASTA from RCSB API by PDB ID
  - local_file: Process a local PDB file

Features: chain splitting, HETATM handling, missing residues analysis.
"""

import os
from datetime import datetime

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter,
    FileParameterEdit,
    FolderParameter,
    SelectParameter,
    StringParameter,
)

try:
    from .core import (
        extract_missing_residues,
        extract_sequences,
        fetch_fasta_from_rcsb,
        fetch_pdb_from_rcsb,
        parse_pdb_structure,
        write_fasta_files,
        write_missing_residues_csv,
    )
except ImportError:
    # Server environment: BioPython not installed.
    # OPTIONS still work — functions only called in PIXI_SUBPROCESS.
    extract_missing_residues = extract_sequences = fetch_fasta_from_rcsb = None
    fetch_pdb_from_rcsb = parse_pdb_structure = write_fasta_files = None
    write_missing_residues_csv = None


class PdbFastaBiopython(Node):
    """
    Parse PDB structures and extract FASTA sequences using BioPython.

    Supports two modes:
    - **PDB ID mode**: Fetch PDB + FASTA from RCSB by ID (e.g. "3LZ0")
    - **Local file mode**: Process a local PDB file

    Features:
    - Automatic chain splitting into individual FASTA files
    - HETATM record handling (ligands, modified residues)
    - Missing residues extraction from PDB header
    - Protein and DNA/RNA residue support

    Input: Optional predecessor data (case_name, working_path)
    Output: FASTA file(s), optional missing residues CSV
    """

    # NOTE: Metadata (name, hashtags, num_in, num_out) comes from meta.toml.
    # NOTE: force_to_run is inherited from Node.BASE_OPTIONS — do NOT add it here.

    category = "io"
    tags = ["protein", "pdb", "fasta", "file-conversion", "sequence-extraction"]

    OPTIONS = {
        "input_mode": SelectParameter(
            "Input Mode",
            options=["pdb_id", "local_file"],
            default="pdb_id",
            docstring="'pdb_id' to fetch from RCSB API, 'local_file' for a local PDB.",
        ),
        "pdb_id": StringParameter(
            "PDB ID",
            docstring="4-character RCSB PDB identifier (e.g. '3LZ0'). Used in pdb_id mode.",
        ),
        "pdb_file": FileParameterEdit(
            "PDB File",
            docstring="Path to a local PDB file. Used in local_file mode.",
        ),
        "case_name": StringParameter(
            "Case Name",
            docstring="Identifier for this case. Falls back to PDB ID, predecessor data, or 'protein'.",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory where output files will be written.",
        ),
        "split_chains": BooleanParameter(
            "Split by Chain",
            default=True,
            docstring="Create separate FASTA files for each chain.",
        ),
        "include_hetatm": BooleanParameter(
            "Include HETATM",
            default=False,
            docstring="Include HETATM records (ligands, modified residues) in sequence extraction.",
        ),
        "check_pdb_header": BooleanParameter(
            "Check PDB Header",
            default=False,
            docstring="Extract missing residues information from the PDB header.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute PDB parsing and FASTA extraction."""
        log_message("Starting execution of PdbFastaBiopython")

        try:
            result = NodeResult()
            result.metadata["execution_time"] = datetime.now().isoformat()

            # --- Read parameters ---
            input_mode = flow_vars["input_mode"].get_value() or "pdb_id"
            pdb_id = (flow_vars["pdb_id"].get_value() or "").strip().upper()
            pdb_file = flow_vars["pdb_file"].get_value() or ""
            case_name = flow_vars["case_name"].get_value() or ""
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            split_chains = flow_vars["split_chains"].get_value()
            include_hetatm = flow_vars["include_hetatm"].get_value()
            check_pdb_header = flow_vars["check_pdb_header"].get_value()

            # --- Resolve case_name ---
            if not case_name:
                if input_mode == "pdb_id" and pdb_id:
                    case_name = pdb_id.lower()
                elif predecessor_data and predecessor_data[0]:
                    input_data = predecessor_data[0]
                    case_name = input_data.get("case_name", "protein")
                else:
                    case_name = "protein"

            result.metadata["case_name"] = case_name

            # Create case_name subfolder — all downstream nodes work inside this
            case_dir = os.path.join(output_dir, case_name)
            os.makedirs(case_dir, exist_ok=True)
            output_dir = case_dir  # All outputs go into the case folder

            log_message(f"Mode: {input_mode}, Case: {case_name}")
            log_message(f"Output: {output_dir}")

            # --- Acquire PDB ---
            pdb_path = None

            if input_mode == "pdb_id":
                if not pdb_id:
                    raise NodeException(
                        "parameter", "PDB ID is required in pdb_id mode."
                    )
                log_message(f"Fetching PDB {pdb_id} from RCSB...")
                pdb_text = fetch_pdb_from_rcsb(pdb_id)
                pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
                with open(pdb_path, "w") as f:
                    f.write(pdb_text)
                log_message(f"Saved PDB to {pdb_path}")

                # Also fetch the RCSB FASTA for reference
                log_message(f"Fetching RCSB FASTA for {pdb_id}...")
                try:
                    rcsb_fasta = fetch_fasta_from_rcsb(pdb_id)
                    rcsb_fasta_path = os.path.join(output_dir, f"{pdb_id}_rcsb.fasta")
                    with open(rcsb_fasta_path, "w") as f:
                        f.write(rcsb_fasta)
                    log_message(f"Saved RCSB FASTA to {rcsb_fasta_path}")
                    result.files["output"]["rcsb_fasta"] = self.format_output_path(
                        rcsb_fasta_path
                    )
                except Exception as e:
                    log_message(f"Warning: Could not fetch RCSB FASTA: {e}")

            elif input_mode == "local_file":
                pdb_path = self.resolve_path(pdb_file)
                if not pdb_path or not os.path.exists(pdb_path):
                    raise NodeException("parameter", f"PDB file not found: {pdb_path}")
            else:
                raise NodeException(
                    "parameter",
                    f"Unknown input_mode: {input_mode}. Use 'pdb_id' or 'local_file'.",
                )

            result.files["input"]["pdb_file"] = self.format_output_path(pdb_path)

            # --- Parse PDB ---
            log_message("Parsing PDB structure...")
            structure = parse_pdb_structure(pdb_path, case_name, is_file=True)

            # --- Extract sequences ---
            chains = extract_sequences(structure, include_hetatm=include_hetatm)
            if not chains:
                raise NodeException(
                    "execution",
                    "No sequences extracted from PDB file. "
                    "Check if the file contains valid residues.",
                )

            log_message(
                f"Extracted sequences from {len(chains)} chain(s): "
                f"{', '.join(chains.keys())}"
            )

            # --- Write FASTA files ---
            pdb_basename = os.path.basename(pdb_path)
            fasta_files = write_fasta_files(
                chains, output_dir, case_name, split_chains, pdb_basename
            )
            for label, path in fasta_files.items():
                result.files["output"][label] = self.format_output_path(path)
                log_message(f"Wrote FASTA: {path}")

            # --- Missing residues (optional) ---
            missing_csv_files = {}
            if check_pdb_header:
                log_message("Checking PDB header for missing residues...")
                missing = extract_missing_residues(structure)
                if missing:
                    missing_csv_files = write_missing_residues_csv(missing, output_dir)
                    for chain_id, csv_path in missing_csv_files.items():
                        result.files["output"][
                            f"missing_residues_{chain_id}"
                        ] = self.format_output_path(csv_path)
                        log_message(
                            f"Wrote missing residues for chain {chain_id}: {csv_path}"
                        )
                else:
                    log_message("No missing residues found in PDB header.")

            # --- Build chain_info dict for result data ---
            chain_info = {}
            for cid, info in chains.items():
                chain_info[cid] = {
                    "length": info.length,
                    "protein_residues": info.protein_residues,
                    "dna_rna_residues": info.dna_rna_residues,
                    "type": info.chain_type,
                }

            result.data.update(
                {
                    "case_name": case_name,
                    "input_mode": input_mode,
                    "working_path": self.format_output_path(output_dir),
                    "num_chains": len(chains),
                    "chain_info": chain_info,
                }
            )

            result.success = True
            result.message = (
                f"Successfully extracted sequences from {len(chains)} chain(s) "
                f"in {case_name}"
            )

            log_message(f"PdbFastaBiopython completed for case: {case_name}")
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            log_message(f"Error in PdbFastaBiopython: {str(e)}")
            raise NodeException("pdb to fasta conversion", str(e))
