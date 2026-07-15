"""MD Analysis: Trajectory Center — PBC-correct an MD trajectory so the
metallopeptide is whole and centred before downstream analysis.

GROMACS writes coordinates wrapped into the periodic box, so a solute
that drifts across a box face is split across the boundary. DSSP (the
md_analysis_helix node) then mis-assigns secondary structure because
backbone H-bond geometry is wrong. This node makes the solute whole
(unwrap via the bond graph), centres it in the box, and wraps the
solvent/membrane back in — the ``gmx trjconv -pbc whole -center
-pbc mol`` recipe, done with MDAnalysis so the result is identical to
what the (also MDAnalysis-backed) analysis nodes see.

    gmx_mdrun(_local) → md_traj_center → md_analysis_helix / _distance

Output keys: ``output_trajectory`` (+ alias ``trajectory``),
``output_gro`` / ``tpr`` (topology pass-through for the analysis
nodes), ``case_name``, ``working_path``.
"""
from __future__ import annotations

import os

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FolderParameter, IntegerParameter, SelectParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import run_traj_center
except ImportError:  # script-mode / server introspection fallback
    try:
        from core import run_traj_center  # type: ignore
    except ImportError:
        run_traj_center = None


def _from_predecessors(predecessor_data, *keys):
    """First present value for any of ``keys`` across the predecessors'
    data dicts (supports nested ``output_files``)."""
    for pred in (predecessor_data or []):
        if not pred:
            continue
        for k in keys:
            if k in pred and pred[k]:
                return pred[k]
        of = pred.get("output_files") or {}
        for k in keys:
            if k in of and of[k]:
                return of[k]
    return None


class MdTrajCenter(Node):
    """Make the solute whole, centre it, wrap solvent — PBC-correct."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default=""),
        "input_trajectory": FileParameterEdit(
            "Trajectory (.xtc/.trr)",
            default="",
            docstring=(
                "MD trajectory to PBC-correct. Leave empty to "
                "auto-discover the trajectory output of a predecessor "
                "MD-run node."
            ),
            optional=True,
        ),
        "input_topology": FileParameterEdit(
            "Topology (.tpr)",
            default="",
            docstring=(
                "Topology for the trajectory. A .tpr is required — "
                "unwrapping needs the bond graph, which a .tpr carries "
                "and a bare .gro does not. Auto-discovers a predecessor "
                "MD-run node's .tpr."
            ),
            optional=True,
        ),
        "solute_selection": StringParameter(
            "Solute Selection",
            default="protein",
            docstring=(
                "MDAnalysis selection seeding the solute to keep whole "
                "and centre. It is expanded to whole connected molecules "
                "('same fragment as'), so the default 'protein' pulls in "
                "the bonded SnP fragment too — the whole metallopeptide "
                "is centred, solvent/membrane is not."
            ),
        ),
        "stride": IntegerParameter(
            "Frame Stride",
            default=1,
            docstring="Write every Nth frame (>=1).",
        ),
        "extract_first": SelectParameter(
            "Extract Solute First (large systems)",
            default="no",
            options=["no", "yes"],
            docstring=(
                "Two-pass mode: stream once to write a SOLUTE-ONLY temp "
                "xtc, then unwrap + centre on that small subset. Output "
                "is solute-only, ~90× smaller for a DPPC membrane system. "
                "Recommended whenever the downstream analyses only need "
                "the solute (DSSP, distance to metal, RMSD, …)."
            ),
        ),
        "output_prefix": StringParameter(
            "Output Prefix",
            default="centered",
            docstring="Prefix → <prefix>.xtc + <prefix>.gro.",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for the PBC-corrected trajectory.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting trajectory centering...", node_id=self.node_id,
                   progress=0)
        if run_traj_center is None:
            raise NodeException("setup",
                "core.py could not be imported — run this node in a pixi "
                "env with MDAnalysis available.")
        try:
            result = NodeResult()
            input_data = (predecessor_data[0]
                          if predecessor_data and predecessor_data[0] else {})
            case_name = (flow_vars["case_name"].get_value()
                         or input_data.get("case_name", "case"))

            traj = flow_vars["input_trajectory"].get_value() or \
                _from_predecessors(predecessor_data, "trajectory",
                                   "output_trajectory")
            topo = flow_vars["input_topology"].get_value() or \
                _from_predecessors(predecessor_data, "tpr", "output_tpr",
                                   "structure")
            if not traj:
                raise NodeException("setup",
                    "No trajectory — set 'Trajectory' or connect an "
                    "MD-run predecessor.")
            if not topo:
                raise NodeException("setup",
                    "No topology — set 'Topology' (.tpr) or connect an "
                    "MD-run predecessor.")
            traj = self.resolve_path(traj)
            topo = self.resolve_path(topo)
            for label, p in (("trajectory", traj), ("topology", topo)):
                if not p or not os.path.isfile(p):
                    raise NodeException("setup", f"{label} not found: {p}")

            selection = flow_vars["solute_selection"].get_value() or "protein"
            stride = max(1, int(flow_vars["stride"].get_value() or 1))
            extract_first = (flow_vars["extract_first"].get_value() or "no") == "yes"
            prefix = flow_vars["output_prefix"].get_value() or "centered"
            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)
            out_traj = os.path.join(out_dir, f"{prefix}.xtc")
            out_gro = os.path.join(out_dir, f"{prefix}.gro")

            stream_log(f"PBC-correcting {os.path.basename(traj)} "
                       f"(solute {selection!r}, stride {stride}"
                       f"{', extract-first' if extract_first else ''})...",
                       node_id=self.node_id, progress=30)
            res = run_traj_center(topo, traj, selection, out_traj,
                                  out_gro=out_gro, stride=stride,
                                  extract_first=extract_first)

            stream_log(
                f"Centered {res['n_frames_out']} frames "
                f"({res['n_solute_atoms']} solute atoms made whole)",
                node_id=self.node_id, progress=100)

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(out_dir),
                "output_trajectory": self.format_output_path(out_traj),
                "trajectory": self.format_output_path(out_traj),
                "output_gro": self.format_output_path(res["out_gro"]
                                                      or out_gro),
                # the input .tpr stays valid (atom order unchanged) — pass
                # it through so the analysis nodes auto-discover a topology
                "tpr": self.format_output_path(topo),
                "n_frames_in": res["n_frames_in"],
                "n_frames_out": res["n_frames_out"],
                "n_solute_atoms": res["n_solute_atoms"],
            }
            result.files["output"] = {
                "trajectory": self.format_output_path(out_traj),
                "gro": self.format_output_path(res["out_gro"] or out_gro),
            }
            result.success = True
            result.message = (
                f"PBC-corrected: {res['n_frames_out']} frames, solute "
                f"({res['n_solute_atoms']} atoms) made whole + centred")
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("md traj center", str(e))
