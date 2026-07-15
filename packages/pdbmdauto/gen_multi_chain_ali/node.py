"""
gen-multi-chain-ali — BoCoFlow node wrapper.

Merges per-chain .ali alignment files into a single multi-chain alignment.
Reads per-chain .ali files from predecessor (gen_ali) and produces a
merged alignment file for multi-chain homology modeling.

Output is consumed by:
  - merge_pdb_chains (uses chain type information)
  - fix_residues_promod3 (uses the merged alignment for gap-filling)
"""

import os
from datetime import datetime

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FolderParameter,
    StringParameter,
)

try:
    from .core import process_multi_chain_ali
except ImportError:
    process_multi_chain_ali = None


class GenMultiChainAli(Node):
    """
    Merge per-chain alignment files into a single multi-chain alignment.

    Reads per-chain .ali files produced by gen_ali and merges them into
    a single alignment file for multi-chain homology modeling with
    ProMod3 or MODELLER.

    Features:
    - Chain selection (all or specific chains)
    - DNA/RNA masking with '.' characters (MODELLER convention)
    - '/' chain separators between chains
    - Proper protein chain range in header

    Input: predecessor data from gen_ali (per-chain .ali files, chain info)
    Output: Single merged .ali file in Merge/ subdirectory
    """

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name",
            default="",
            docstring="Case identifier. Leave empty to use predecessor data.",
        ),
        "merge_folder_name": StringParameter(
            "Merge Folder Name",
            default="Merge",
            docstring="Subfolder name for the merged alignment file.",
        ),
        "selected_chain_ids": StringParameter(
            "Selected Chain IDs",
            default="all",
            docstring="Chains to include: 'all' or comma-separated (e.g. 'A,B,C').",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute multi-chain alignment merge."""
        try:
            from bocoflow_core.stream_logger import stream_log
        except ImportError:
            stream_log = lambda msg, **kw: log_message(msg)

        try:
            result = NodeResult()
            stream_log("Starting multi-chain alignment merge", node_id=self.node_id, progress=0)

            # --- Parse predecessor data ---
            if not predecessor_data or not predecessor_data[0]:
                raise NodeException(
                    "gen_multi_chain_ali",
                    "No predecessor data. Connect gen_ali upstream.",
                )
            # In BF2 pipeline, predecessor_data[0] is the flat result.data dict
            input_data = predecessor_data[0]

            # --- Resolve parameters ---
            case_name = (
                flow_vars["case_name"].get_value()
                or input_data.get("case_name")
                or "protein"
            )
            merge_folder_name = flow_vars["merge_folder_name"].get_value() or "Merge"
            selected_chains_str = flow_vars["selected_chain_ids"].get_value() or "all"

            # Resolve output directory from predecessor's working_path
            working_path = input_data.get("working_path", "")
            output_dir = self.resolve_path(working_path) if working_path else ""
            if not output_dir:
                raise NodeException("gen_multi_chain_ali", "No output directory from predecessor.")

            log_message(f"Case: {case_name}, Output: {output_dir}, Selection: {selected_chains_str}")

            stream_log("Collecting per-chain .ali files", node_id=self.node_id, progress=20)

            # --- Find per-chain .ali files by scanning chain subdirectories ---
            ali_file_paths = {}
            chain_list = input_data.get("pdb_chain_list", [])
            for chain_id in chain_list:
                ali_path = os.path.join(output_dir, chain_id, "homology.ali")
                if os.path.exists(ali_path):
                    ali_file_paths[chain_id] = ali_path

            # Fallback: scan all subdirs for homology.ali
            if not ali_file_paths:
                for entry in os.listdir(output_dir):
                    ali_path = os.path.join(output_dir, entry, "homology.ali")
                    if os.path.isdir(os.path.join(output_dir, entry)) and os.path.exists(ali_path):
                        ali_file_paths[entry] = ali_path

            if not ali_file_paths:
                raise NodeException(
                    "gen_multi_chain_ali",
                    f"No per-chain .ali files found in {output_dir}",
                )

            log_message(f"Found .ali files for chains: {sorted(ali_file_paths.keys())}")

            stream_log("Merging chain alignments", node_id=self.node_id, progress=50)

            # --- Process ---
            merge_result = process_multi_chain_ali(
                ali_file_paths=ali_file_paths,
                output_dir=output_dir,
                case_name=case_name,
                selected_chains_str=selected_chains_str,
                merge_folder_name=merge_folder_name,
            )

            stream_log("Writing merged alignment", node_id=self.node_id, progress=80)

            # --- Build result ---
            formatted_ali_file = self.format_output_path(merge_result.ali_file)
            formatted_output_dir = self.format_output_path(output_dir)
            merge_dir = os.path.join(output_dir, merge_folder_name)
            formatted_merge_dir = self.format_output_path(merge_dir)

            result.files["output"]["multi_chain_ali"] = formatted_ali_file
            # PDB file reference not available in flat predecessor data

            result.data.update({
                "case_name": case_name,
                "working_path": formatted_output_dir,
                "merge_folder": formatted_merge_dir,
                "pdb_chain_list": input_data.get("pdb_chain_list", []),
                "selected_pdb_chain_list": merge_result.selected_chains,
                "chain_types": merge_result.chain_types,
                "chain_info": input_data.get("chain_info", {}),
                "chain_alignment_results": input_data.get("chain_alignment_results", {}),
            })

            result.metadata.update({
                "case_name": case_name,
                "output_dir": formatted_output_dir,
                "execution_time": datetime.now().isoformat(),
            })

            result.success = True
            result.message = (
                f"Merged alignment for {len(merge_result.selected_chains)} chains "
                f"({len(merge_result.protein_chains)} protein, "
                f"{len(merge_result.dna_chains)} DNA/RNA)"
            )

            log_message(result.message)
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            log_message(f"Error in GenMultiChainAli: {str(e)}")
            raise NodeException("gen_multi_chain_ali", str(e))
