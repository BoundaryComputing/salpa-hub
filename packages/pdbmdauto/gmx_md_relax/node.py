"""
gmx-md-relax — BoCoFlow node wrapper.

GROMACS relaxation with protocol selector:
- full_4step: Progressive unfreezing (nvt_fixOri → nvt_fixOriBackbone → mm1 → mm2)
- em_only: Single steepest-descent EM (quick testing)
- custom: Single step with user-provided MDP

The full_4step protocol uses OriHeavy/OriBackBone groups from the NDX file
(created by gen_gmx_ndx) via GROMACS freezegrps mechanism.
"""

import os
from datetime import datetime

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit,
    IntegerParameter,
    SelectParameter,
    StringParameter,
)

try:
    from .core import process_full_4step, process_single_step
except ImportError:
    from core import process_full_4step, process_single_step


class GmxMdRelax(Node):
    """
    GROMACS relaxation/equilibration with protocol selector.

    Protocols:
    - full_4step (default): 4-step progressive unfreezing protocol from legacy pdbmdauto.
      Uses OriHeavy/OriBackBone freezegrps to first equilibrate water/ions around the
      rigid protein, then progressively release constraints. Essential for systems with
      homology-modeled missing residues.
    - em_only: Single energy minimization (steepest descent). Quick for testing.
    - custom: Single grompp+mdrun with user-provided MDP file.
    """

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name", default="",
            docstring="Leave empty to use predecessor data.",
        ),
        "protocol": SelectParameter(
            "Relaxation Protocol",
            options=["full_4step", "em_only", "custom"],
            default="full_4step",
            docstring=(
                "full_4step: 4-step progressive unfreezing (NVT frozen→NVT backbone→CG EM×2). "
                "em_only: single steepest-descent EM. "
                "custom: single step with user MDP."
            ),
        ),
        "input_top_file": FileParameterEdit(
            "Topology File (.top)", default="",
            docstring="Topology. Leave empty: auto-discovers from predecessor.",
        ),
        "input_gro_file": FileParameterEdit(
            "Structure File (.gro)", default="",
            docstring="Solvated structure. Leave empty: auto-discovers from predecessor.",
        ),
        "input_mdp_file": FileParameterEdit(
            "Parameters File (.mdp)", default="",
            docstring="Custom MDP (only for 'custom' protocol). Leave empty for bundled defaults.",
        ),
        "input_ndx_file": FileParameterEdit(
            "Index File (.ndx)", default="",
            docstring="Index file with OriHeavy/OriBackBone groups. Leave empty: auto-discovers.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        try:
            from bocoflow_core.stream_logger import stream_log
        except ImportError:
            stream_log = lambda msg, **kw: log_message(msg)

        try:
            result = NodeResult()
            protocol = flow_vars["protocol"].get_value() or "full_4step"
            stream_log(
                f"Starting relaxation: {protocol}",
                node_id=self.node_id, progress=0,
            )

            input_data = (
                predecessor_data[0]
                if predecessor_data and predecessor_data[0]
                else {}
            )
            case_name = (
                flow_vars["case_name"].get_value()
                or input_data.get("case_name", "protein")
            )

            working_path = input_data.get("working_path", "")
            gmx_dir = self.resolve_path(working_path) if working_path else ""

            # Resolve file inputs
            gro = self.resolve_path(flow_vars["input_gro_file"].get_value()) or ""
            top = self.resolve_path(flow_vars["input_top_file"].get_value()) or ""
            mdp = self.resolve_path(flow_vars["input_mdp_file"].get_value()) or ""
            ndx = self.resolve_path(flow_vars["input_ndx_file"].get_value()) or ""

            # Auto-discover from predecessor data
            if not gro and input_data.get("output_gro"):
                gro = self.resolve_path(input_data["output_gro"])
            if not top and input_data.get("output_top"):
                top = self.resolve_path(input_data["output_top"])
            if not ndx and input_data.get("output_ndx"):
                ndx = self.resolve_path(input_data["output_ndx"])

            # Scan gmx/ folder for files if still missing
            if gmx_dir and (not gro or not top or not ndx):
                if os.path.isdir(gmx_dir):
                    for f in os.listdir(gmx_dir):
                        fp = os.path.join(gmx_dir, f)
                        if not gro and f == "ion.gro":
                            gro = fp
                        elif not top and f == "topol.top":
                            top = fp
                        elif not ndx and f == "index.ndx":
                            ndx = fp

            # Same gmx/ folder for all GROMACS operations
            output_dir = gmx_dir if gmx_dir else (
                os.path.dirname(gro) if gro else "."
            )

            if not gro or not top or not ndx:
                missing = [
                    n for n, v in [("GRO", gro), ("TOP", top), ("NDX", ndx)]
                    if not v
                ]
                raise NodeException(
                    "gmx_md_relax",
                    f"Missing files: {', '.join(missing)}",
                )

            # ── Execute protocol ──────────────────────────────────────────
            if protocol == "full_4step":
                stream_log(
                    "Running 4-step protocol: nvt_fixOri → nvt_fixOriBackbone → mm1 → mm2",
                    node_id=self.node_id, progress=10,
                )

                # MDP files are bundled in demo_data/
                node_dir = getattr(self, "_node_dir", None) or os.path.dirname(__file__)
                mdp_dir = os.path.join(node_dir, "demo_data")

                relax_result = process_full_4step(
                    gro_file=gro,
                    top_file=top,
                    ndx_file=ndx,
                    output_dir=output_dir,
                    mdp_dir=mdp_dir,
                )

            elif protocol == "em_only":
                stream_log(
                    "Running single EM step",
                    node_id=self.node_id, progress=10,
                )

                if not mdp:
                    node_dir = getattr(self, "_node_dir", None) or os.path.dirname(__file__)
                    mdp = os.path.join(node_dir, "demo_data", "em.mdp")

                relax_result = process_single_step(
                    gro_file=gro,
                    top_file=top,
                    mdp_file=mdp,
                    ndx_file=ndx,
                    output_dir=output_dir,
                    run_label="em",
                )

            elif protocol == "custom":
                if not mdp:
                    raise NodeException(
                        "gmx_md_relax",
                        "Custom protocol requires an MDP file.",
                    )

                stream_log(
                    f"Running custom step: {os.path.basename(mdp)}",
                    node_id=self.node_id, progress=10,
                )

                run_label = os.path.splitext(os.path.basename(mdp))[0]
                relax_result = process_single_step(
                    gro_file=gro,
                    top_file=top,
                    mdp_file=mdp,
                    ndx_file=ndx,
                    output_dir=output_dir,
                    run_label=run_label,
                )

            else:
                raise NodeException(
                    "gmx_md_relax",
                    f"Unknown protocol: {protocol}",
                )

            if not relax_result.success:
                log_message(f"Relaxation log:\n{relax_result.log}")
                raise NodeException(
                    "gmx_md_relax",
                    f"Relaxation failed:\n{relax_result.log[:500]}",
                )

            stream_log(
                f"Relaxation complete: {', '.join(relax_result.steps_completed)}"
                f" (max force: {relax_result.max_force:.1f} kJ/mol/nm"
                f"{', SAFE' if relax_result.em_safe else ', HIGH'})",
                node_id=self.node_id, progress=90,
            )

            result.data.update({
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "output_gro": self.format_output_path(relax_result.output_gro),
                "output_top": self.format_output_path(top),
                "output_ndx": self.format_output_path(ndx),
                "protocol": protocol,
                "steps_completed": relax_result.steps_completed,
                "max_force": relax_result.max_force,
                "em_safe": relax_result.em_safe,
            })

            result.success = True
            result.message = (
                f"Relaxation '{protocol}' complete: "
                f"{len(relax_result.steps_completed)} steps, "
                f"max force {relax_result.max_force:.1f} kJ/mol/nm"
            )
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("gmx_md_relax", str(e))
