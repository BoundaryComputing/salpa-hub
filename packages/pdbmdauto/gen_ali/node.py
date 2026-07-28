"""
gen-ali — BoCoFlow node wrapper.

Generates per-chain alignment files (.seq and .ali) for homology modeling.
Reads PDB structure + missing residues from predecessor (pdb_fasta_biopython)
and produces template/target alignment pairs.

Output .ali files are consumed by downstream nodes:
  - gen_multi_chain_ali (merges per-chain alignments)
  - fix_residues_promod3 (fills missing residues using alignment)
"""

import os
from datetime import datetime

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter,
    FolderParameter,
    StringParameter,
)

try:
    from .core import (
        extract_missing_from_pdb,
        extract_present_residues,
        process_all_chains,
        read_fasta_sequence,
    )
except ImportError:
    # Stage 2. node_runner puts the node's directory on sys.path and imports
    # node.py as a TOP-LEVEL module, so there is no package for `.core` to be
    # relative to. Without this the next stage ran instead and every symbol
    # below was None by the time execute() called it.
    try:
        from core import (
            extract_missing_from_pdb,
            extract_present_residues,
            process_all_chains,
            read_fasta_sequence,
        )
    except ImportError:
        extract_missing_from_pdb = extract_present_residues = None
        process_all_chains = read_fasta_sequence = None

class GenAli(Node):
    """
    Generate per-chain alignment files for homology modeling.

    Reads a PDB structure and missing residues information (from predecessor
    pdb_fasta_biopython node or PDB header) and generates:
    - .seq files: raw FASTA alignment per chain
    - .ali files: template/target alignment pairs for ProMod3 or MODELLER

    The template entry shows the existing structure with '-' at missing
    positions. The target entry shows the complete sequence.

    Input: predecessor data from pdb_fasta_biopython (PDB file, chain info,
           optional missing residues CSVs, optional FASTA files)
    Output: Per-chain .seq and .ali files in chain-specific subdirectories
    """

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name",
            default="",
            docstring="Case identifier. Leave empty to use predecessor data.",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            default="",
            docstring="Output directory. Leave empty to use predecessor working_path.",
        ),
        "append_end_in_seq": BooleanParameter(
            "Append End Marker",
            default=True,
            docstring="Append '*' terminator at the end of sequences.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute alignment file generation."""
        try:
            from bocoflow_core.stream_logger import stream_log
        except ImportError:
            stream_log = lambda msg, **kw: log_message(msg)

        try:
            result = NodeResult()
            stream_log("Starting alignment generation", node_id=self.node_id, progress=0)

            # --- Parse predecessor data ---
            if not predecessor_data or not predecessor_data[0]:
                raise NodeException(
                    "gen_ali", "No predecessor data. Connect pdb_fasta_biopython upstream."
                )

            # In BF2 pipeline, predecessor_data[0] is the flat result.data dict
            # (NOT a nested {data, files, metadata} structure)
            input_data = predecessor_data[0]

            # --- Resolve parameters (node options > predecessor) ---
            case_name = (
                flow_vars["case_name"].get_value()
                or input_data.get("case_name")
                or "protein"
            )

            output_dir_param = flow_vars["output_dir"].get_value()
            if output_dir_param:
                output_dir = self.resolve_path(output_dir_param)
            else:
                # Predecessor's working_path (rel: or abs: prefixed)
                working_path = input_data.get("working_path", "")
                output_dir = self.resolve_path(working_path) if working_path else ""

            if not output_dir:
                raise NodeException("gen_ali", "No output directory specified.")

            append_end = flow_vars["append_end_in_seq"].get_value()

            os.makedirs(output_dir, exist_ok=True)
            log_message(f"Case: {case_name}, Output: {output_dir}")

            # --- Find PDB file in output directory ---
            # Scan for .pdb files (pdb_fasta_biopython writes {PDB_ID}.pdb)
            pdb_path = None
            for f in os.listdir(output_dir):
                if f.endswith(".pdb"):
                    pdb_path = os.path.join(output_dir, f)
                    break
            if not pdb_path:
                raise NodeException("gen_ali", f"No PDB file found in {output_dir}")

            log_message(f"PDB file: {pdb_path}")
            stream_log("Parsing PDB structure", node_id=self.node_id, progress=20)

            # --- Collect chain info and find files by scanning directory ---
            chain_info = input_data.get("chain_info", {})
            chain_ids = sorted(chain_info.keys()) if chain_info else None

            # Scan for missing residues CSVs and FASTA files
            missing_csv_paths = {}
            fasta_paths = {}
            for f in os.listdir(output_dir):
                if f.startswith("missing_residues_chain_") and f.endswith(".csv"):
                    chain_id = f.replace("missing_residues_chain_", "").replace(".csv", "")
                    missing_csv_paths[chain_id] = os.path.join(output_dir, f)
                elif f.endswith(".fasta") and "_chain_" in f:
                    chain_id = f.split("_chain_")[-1].replace(".fasta", "")
                    fasta_paths[chain_id] = os.path.join(output_dir, f)

            if missing_csv_paths:
                log_message(f"Found missing residues CSVs for chains: {sorted(missing_csv_paths.keys())}")
            else:
                log_message("No missing residues CSVs found; will extract from PDB header")

            stream_log("Generating alignments", node_id=self.node_id, progress=40)

            # --- Process all chains ---
            processing_result = process_all_chains(
                pdb_path=pdb_path,
                output_dir=output_dir,
                case_name=case_name,
                chain_ids=chain_ids,
                missing_csv_paths=missing_csv_paths if missing_csv_paths else None,
                fasta_paths=fasta_paths if fasta_paths else None,
                append_end=append_end,
            )

            chain_results = processing_result["chain_results"]
            seq_agree_all = processing_result["seq_agree_all"]

            stream_log("Writing output files", node_id=self.node_id, progress=80)

            # --- Build result ---
            chain_alignment_results = {}
            for chain_id, cr in chain_results.items():
                # Record output files
                result.files["output"][f"ali_{chain_id}"] = self.format_output_path(cr.ali_file)
                result.files["output"][f"seq_{chain_id}"] = self.format_output_path(cr.seq_file)

                chain_alignment_results[chain_id] = {
                    "seq_agree": cr.seq_agree,
                    "chain_type": cr.chain_type,
                    "num_present": cr.num_present,
                    "num_missing": cr.num_missing,
                    "start_resid": cr.start_resid,
                }

                log_message(
                    f"Chain {chain_id}: {cr.num_present} present, "
                    f"{cr.num_missing} missing, type={cr.chain_type}, "
                    f"agree={cr.seq_agree}"
                )

            # Preserve input PDB file reference
            result.files["input"]["pdb_file"] = self.format_output_path(pdb_path)

            # Populate result data (preserve upstream data for downstream nodes)
            formatted_output_dir = self.format_output_path(output_dir)
            result.data.update({
                "case_name": case_name,
                "working_path": formatted_output_dir,
                "pdb_chain_list": sorted(chain_results.keys()),
                "chain_info": chain_info,  # Preserved from predecessor
                "chain_alignment_results": chain_alignment_results,
                "seq_agree_pdb_fasta": seq_agree_all,
            })

            result.metadata.update({
                "case_name": case_name,
                "output_dir": formatted_output_dir,
                "execution_time": datetime.now().isoformat(),
            })

            result.success = True
            result.message = (
                f"Generated alignment files for {len(chain_results)} chain(s). "
                f"Sequence agreement: {'all match' if seq_agree_all else 'MISMATCH detected'}"
            )

            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            log_message(f"Error in GenAli: {str(e)}")
            raise NodeException("gen_ali", str(e))
