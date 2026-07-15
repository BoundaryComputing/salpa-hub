"""
merge-pdb-chains — BoCoFlow node wrapper.

Extracts selected chains from the original PDB file and merges them into
a single structure file. Uses BioPython's Bio.PDB (replaces legacy PyMOL).

Input: predecessor data from gen_multi_chain_ali (chain types, working_path)
Output: Merged PDB file + chain metadata JSON files in Merge/ subdirectory
"""

import os
from datetime import datetime

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    StringParameter,
)

try:
    from .core import process_merge
except ImportError:
    process_merge = None


class MergePdbChains(Node):
    """
    Merge selected PDB chains into a single structure file.

    Reads chain information from predecessor (gen_multi_chain_ali) and
    extracts the selected chains from the original PDB file into a
    merged PDB. DNA/RNA chains are written with HETATM records
    (MODELLER/ProMod3 convention).

    Output files (in Merge/ subdirectory):
    - merge.pdb: Merged structure with selected chains
    - {case}_chain_type.json: Chain type mapping (P1/DL)
    - {case}_chain_name.json: Selected chain IDs list

    Input: predecessor data (working_path, chain_types, pdb_chain_list)
    Output: Merged PDB file path and chain metadata
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
            docstring="Subfolder name for merged output files.",
        ),
        "merge_file_name": StringParameter(
            "Merge File Name",
            default="merge.pdb",
            docstring="Output filename for the merged PDB.",
        ),
        "selected_chain_ids": StringParameter(
            "Selected Chain IDs",
            default="all",
            docstring="Chains to include: 'all' or comma-separated (e.g. 'A,B,C').",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute PDB chain merge."""
        try:
            from bocoflow_core.stream_logger import stream_log
        except ImportError:
            stream_log = lambda msg, **kw: log_message(msg)

        try:
            result = NodeResult()
            stream_log("Starting PDB chain merge", node_id=self.node_id, progress=0)

            # --- Parse predecessor data (flat dict) ---
            if not predecessor_data or not predecessor_data[0]:
                raise NodeException("merge_pdb_chains", "No predecessor data.")
            input_data = predecessor_data[0]

            # --- Resolve parameters ---
            case_name = (
                flow_vars["case_name"].get_value()
                or input_data.get("case_name")
                or "protein"
            )
            merge_folder_name = flow_vars["merge_folder_name"].get_value() or "Merge"
            merge_file_name = flow_vars["merge_file_name"].get_value() or "merge.pdb"
            selected_chains_str = flow_vars["selected_chain_ids"].get_value() or "all"

            # Resolve output directory
            working_path = input_data.get("working_path", "")
            output_dir = self.resolve_path(working_path) if working_path else ""
            if not output_dir:
                raise NodeException("merge_pdb_chains", "No output directory from predecessor.")

            log_message(f"Case: {case_name}, Output: {output_dir}")

            # --- Find PDB file ---
            pdb_path = None
            for f in os.listdir(output_dir):
                if f.endswith(".pdb"):
                    pdb_path = os.path.join(output_dir, f)
                    break
            if not pdb_path:
                raise NodeException("merge_pdb_chains", f"No PDB file found in {output_dir}")

            log_message(f"PDB file: {pdb_path}")
            stream_log("Merging chains", node_id=self.node_id, progress=30)

            # Get chain types from predecessor
            chain_types = input_data.get("chain_types", {})

            # --- Process ---
            merge_result = process_merge(
                pdb_path=pdb_path,
                output_dir=output_dir,
                case_name=case_name,
                selected_chains_str=selected_chains_str,
                chain_types=chain_types,
                merge_folder_name=merge_folder_name,
                merge_file_name=merge_file_name,
            )

            stream_log("Writing output", node_id=self.node_id, progress=80)

            # --- Build result ---
            formatted_output_dir = self.format_output_path(output_dir)
            merge_dir = os.path.join(output_dir, merge_folder_name)
            formatted_merge_dir = self.format_output_path(merge_dir)

            result.files["output"]["merged_pdb"] = self.format_output_path(merge_result.output_pdb)
            result.files["output"]["chain_type_json"] = self.format_output_path(merge_result.chain_type_file)
            result.files["output"]["chain_name_json"] = self.format_output_path(merge_result.chain_name_file)

            result.data.update({
                "case_name": case_name,
                "working_path": formatted_output_dir,
                "merge_folder": formatted_merge_dir,
                "output_pdb": self.format_output_path(merge_result.output_pdb),
                "selected_pdb_chain_list": merge_result.selected_chains,
                "chain_types": merge_result.chain_types,
                "pdb_chain_list": input_data.get("pdb_chain_list", []),
                "chain_info": input_data.get("chain_info", {}),
            })

            result.metadata.update({
                "case_name": case_name,
                "output_dir": formatted_output_dir,
                "execution_time": datetime.now().isoformat(),
            })

            result.success = True
            result.message = (
                f"Merged {len(merge_result.selected_chains)} chains into {merge_file_name}"
            )

            log_message(result.message)
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            log_message(f"Error in MergePdbChains: {str(e)}")
            raise NodeException("merge_pdb_chains", str(e))
