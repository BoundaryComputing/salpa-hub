"""
gen-gmx-ndx — BoCoFlow node wrapper.

Generates custom GROMACS index groups (OriHeavy, OriBackBone) that distinguish
original crystallographic residues from homology-modeled residues.

Runs BEFORE solvation — reads vacuum pdb2gmx.gro from pka_gmx_em.
The missing_residues CSVs come from the case folder (created by pdb_fasta_biopython).

Input: Vacuum structure (.gro) + missing residues CSVs (from case folder)
Output: NDX file with OriHeavy and OriBackBone groups (in gmx/ folder)
"""

import os

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit,
    StringParameter,
)

try:
    from .core import generate_ori_ndx
except ImportError:
    # Stage 2. node_runner puts the node's directory on sys.path and imports
    # node.py as a TOP-LEVEL module, so there is no package for `.core` to be
    # relative to. Without this the next stage ran instead and every symbol
    # below was None by the time execute() called it.
    try:
        from core import generate_ori_ndx
    except ImportError:
        generate_ori_ndx = None

class GenGmxNdx(Node):
    """
    Generate custom GROMACS index groups for position restraints.

    Creates OriHeavy (all heavy atoms from original residues) and
    OriBackBone (backbone atoms from original residues) groups.
    Used by gmx_md_relax (full_4step) with freezegrps for vacuum relaxation.

    Runs after pka_gmx_em, before gmx_md_relax and solvation.
    """

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name", default="",
            docstring="Leave empty to use predecessor data.",
        ),
        "input_structure": FileParameterEdit(
            "Structure File (.gro/.pdb)", default="",
            docstring="Vacuum structure from pka_gmx_em. Leave empty: auto-discovers pdb2gmx.gro.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        try:
            from bocoflow_core.stream_logger import stream_log
        except ImportError:
            stream_log = lambda msg, **kw: log_message(msg)

        try:
            result = NodeResult()
            stream_log("Generating index groups", node_id=self.node_id, progress=0)

            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}
            case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "protein")

            # gmx/ folder from predecessor (pka_gmx_em)
            working_path = input_data.get("working_path", "")
            gmx_dir = self.resolve_path(working_path) if working_path else ""

            # Structure file: auto-discover from predecessor output_gro or scan gmx/
            structure = self.resolve_path(flow_vars["input_structure"].get_value()) or ""

            if not structure and input_data.get("output_gro"):
                structure = self.resolve_path(input_data["output_gro"])

            if not structure and gmx_dir:
                for fname in ["pdb2gmx.gro", "em_hbonds.gro"]:
                    fp = os.path.join(gmx_dir, fname)
                    if os.path.exists(fp):
                        structure = fp
                        break

            if not structure:
                raise NodeException("gen_gmx_ndx", "Structure file required.")

            # Missing residues CSVs are in the CASE folder (parent of gmx/)
            case_dir = os.path.dirname(gmx_dir) if gmx_dir else os.path.dirname(structure)

            # NDX output goes in gmx/ folder
            ndx_path = os.path.join(gmx_dir or os.path.dirname(structure), "index.ndx")

            stream_log("Building OriHeavy/OriBackBone groups", node_id=self.node_id, progress=30)

            ndx_result = generate_ori_ndx(
                structure_path=structure,
                ndx_path=ndx_path,
                missing_csv_dir=case_dir,
                chain_ids=input_data.get("pdb_chain_list"),
            )

            if not ndx_result.success:
                raise NodeException("gen_gmx_ndx", f"Index generation failed: {ndx_result.log}")

            stream_log("Index groups created", node_id=self.node_id, progress=90)

            # Pass through predecessor data + add NDX
            result.data.update({
                "case_name": case_name,
                "working_path": self.format_output_path(gmx_dir or os.path.dirname(structure)),
                "output_gro": input_data.get("output_gro", self.format_output_path(structure)),
                "output_top": input_data.get("output_top", ""),
                "output_ndx": self.format_output_path(ndx_result.output_ndx),
                "n_ori_heavy": ndx_result.n_ori_heavy,
                "n_ori_backbone": ndx_result.n_ori_backbone,
            })
            result.success = True
            result.message = f"OriHeavy: {ndx_result.n_ori_heavy}, OriBackBone: {ndx_result.n_ori_backbone}"
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("gen_gmx_ndx", str(e))
