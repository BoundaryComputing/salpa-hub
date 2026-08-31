"""
gmx-solv-ion — BoCoFlow node wrapper.

GROMACS solvation and ionization. If no GRO/TOP are provided, automatically
runs pdb2gmx on the PDB from predecessor to generate them.

Input: PDB or GRO + TOP (from predecessor or explicit), MDP (bundled default or explicit)
Output: Solvated/ionized GRO + updated TOP + NDX
"""

import os
import subprocess
from datetime import datetime

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit,
    FloatParameter,
    SelectParameter,
    StringParameter,
)

try:
    from .core import process_solv_ion
except ImportError:
    # Stage 2. node_runner puts the node's directory on sys.path and imports
    # node.py as a TOP-LEVEL module, so there is no package for `.core` to be
    # relative to. Without this the next stage ran instead and every symbol
    # below was None by the time execute() called it.
    try:
        from core import process_solv_ion
    except ImportError:
        process_solv_ion = None

class GmxSolvIon(Node):
    """
    GROMACS solvation and ionization.

    If connected to a predecessor (e.g., fix_residues_promod3 or pdb2pqr),
    auto-discovers the PDB/PQR file and runs pdb2gmx to generate GRO + TOP.
    Then sets up box, adds water, neutralizes with ions.

    Includes a bundled ions.mdp for genion preprocessing. Users can override
    with their own MDP file.
    """

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="", docstring="Leave empty to use predecessor data."),
        "run_label": StringParameter("Run Label", default="md", docstring="Simulation label."),
        "force_field": SelectParameter("Force Field", options=["amber99sb", "charmm27", "oplsaa", "gromos53a6"], default="amber99sb", docstring="Force field for pdb2gmx."),
        "water_model": SelectParameter("Water Model", options=["tip3p", "spc", "spce", "tip4p"], default="tip3p", docstring="Water model for pdb2gmx."),
        "box_size": StringParameter("Box Size (nm)", default="5 5 5", docstring="Box dimensions 'X Y Z'. Use '0 0 0' for auto triclinic."),
        "ion_conc": FloatParameter("Ion Concentration (mol/L)", default=0.15, docstring="Na+/Cl- concentration. Set 0 to skip."),
        "input_top_file": FileParameterEdit("Topology File (.top)", default="", docstring="Topology from pka_gmx_em. Leave empty: auto-discovers gmx_em/pdb2gmx.top."),
        "input_gro_file": FileParameterEdit("Structure File (.gro)", default="", docstring="Structure from pka_gmx_em. Leave empty: auto-discovers gmx_em/pdb2gmx.gro."),
        "input_mdp_file": FileParameterEdit("Parameters File (.mdp)", default="node:demo_data/ions.mdp", docstring="MDP for genion preprocessing."),
    }

    def execute(self, predecessor_data, flow_vars):
        try:
            from bocoflow_core.stream_logger import stream_log
        except ImportError:
            stream_log = lambda msg, **kw: log_message(msg)

        try:
            result = NodeResult()
            stream_log("Starting solvation & ionization", node_id=self.node_id, progress=0)

            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}
            case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "protein")

            working_path = input_data.get("working_path", "")
            case_dir = self.resolve_path(working_path) if working_path else ""

            # Resolve explicit file inputs
            gro = self.resolve_path(flow_vars["input_gro_file"].get_value()) or ""
            top = self.resolve_path(flow_vars["input_top_file"].get_value()) or ""
            mdp = self.resolve_path(flow_vars["input_mdp_file"].get_value()) or ""

            # Auto-discover from predecessor data (pka_gmx_em output)
            if not gro and input_data.get("output_gro"):
                gro = self.resolve_path(input_data["output_gro"])
            if not top and input_data.get("output_top"):
                top = self.resolve_path(input_data["output_top"])

            if not gro or not top:
                raise NodeException("gmx_solv_ion", "GRO and TOP files required. Connect to pka_gmx_em or provide explicit paths.")

            # Use the same gmx/ folder as pka_gmx_em (all GROMACS ops in one dir)
            output_dir = case_dir if case_dir else os.path.dirname(gro)

            # Use bundled default MDP if none provided
            if not mdp:
                node_dir = getattr(self, "_node_dir", None) or os.path.dirname(__file__)
                default_mdp = os.path.join(node_dir, "demo_data", "ions.mdp")
                if os.path.exists(default_mdp):
                    mdp = default_mdp
                    log_message(f"Using bundled default MDP: {mdp}")
                else:
                    raise NodeException("gmx_solv_ion", "No MDP file provided and no bundled default found.")

            # Ensure output dir exists
            os.makedirs(output_dir, exist_ok=True)

            # Use existing NDX from predecessor (gen_gmx_ndx created it with OriHeavy/OriBackBone)
            ndx = ""
            if input_data.get("output_ndx"):
                ndx = self.resolve_path(input_data["output_ndx"])
            if not ndx:
                ndx = os.path.join(output_dir, "index.ndx")
            # Only generate if NDX doesn't exist at all (standalone mode without gen_gmx_ndx)
            if not os.path.exists(ndx):
                # "q" accepts make_ndx's default groups. stdin, not a shell
                # pipe, so a spaced path stays one argument (bocoflow#104).
                subprocess.run(
                    ["gmx", "make_ndx", "-f", gro, "-o", ndx],
                    input="q\n", capture_output=True, text=True,
                    cwd=output_dir, timeout=30,
                )

            if not all([gro, top, mdp]):
                raise NodeException("gmx_solv_ion", "GRO, TOP, and MDP files are required.")

            stream_log("Running GROMACS solvation pipeline", node_id=self.node_id, progress=30)

            solv_result = process_solv_ion(
                gro_file=gro, top_file=top, mdp_file=mdp, ndx_file=ndx,
                output_dir=output_dir, case_name=case_name,
                run_label=flow_vars["run_label"].get_value() or "md",
                box_size=flow_vars["box_size"].get_value() or "5 5 5",
                ion_conc=flow_vars["ion_conc"].get_value(),
            )

            if not solv_result.success:
                log_message(f"GROMACS log:\n{solv_result.log}")
                raise NodeException("gmx_solv_ion", f"Solvation/ionization failed: {solv_result.log[-1200:]}")

            stream_log("Solvation complete", node_id=self.node_id, progress=90)

            result.data.update({
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "output_gro": self.format_output_path(solv_result.output_gro),
                "output_top": self.format_output_path(solv_result.output_top),
                "output_ndx": self.format_output_path(solv_result.output_ndx),
                "run_label": flow_vars["run_label"].get_value(),
            })
            result.success = True
            result.message = "Solvation and ionization complete"
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("gmx_solv_ion", str(e))
