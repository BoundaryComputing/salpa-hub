"""MD Analysis: α-Helix Content — per-frame DSSP helix content of a
peptide along an MD trajectory.

Consumes a trajectory (.xtc/.trr) + topology (.tpr/.gro/.pdb) — either
set explicitly or auto-discovered from a predecessor MD-run node
(``gmx_mdrun`` / ``gmx_mdrun_local``) — runs DSSP frame by frame, and
writes a CSV time series of helix content plus a per-residue helix
propensity. This is one half of the SnP-peptide case deliverable:
"α-helix content over the trajectory".

    gmx_mdrun(_local) → md_analysis_helix → <case>_helix.csv

Output keys: ``output_helix_csv``, ``helix_summary`` (mean/std/min/max
helix fraction), ``case_name``, ``working_path``.
"""
from __future__ import annotations

import os

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FolderParameter, IntegerParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import (
        CSV_HEADER, helix_csv_rows, run_helix_analysis, summarize, write_csv,
    )
except ImportError:  # script-mode / server introspection fallback
    try:
        from core import (  # type: ignore
            CSV_HEADER, helix_csv_rows, run_helix_analysis, summarize,
            write_csv,
        )
    except ImportError:
        CSV_HEADER = None
        helix_csv_rows = run_helix_analysis = summarize = write_csv = None


def _from_predecessors(predecessor_data, *keys):
    """Return the first present value for any of ``keys`` across the
    predecessors' data dicts (supports nested ``output_files``)."""
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


class MdAnalysisHelix(Node):
    """Per-frame DSSP α-helix content of the peptide along a trajectory."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default=""),
        "input_trajectory": FileParameterEdit(
            "Trajectory (.xtc/.trr)",
            default="",
            docstring=(
                "MD trajectory. Leave empty to auto-discover the "
                "trajectory output of a predecessor MD-run node."
            ),
            optional=True,
        ),
        "input_topology": FileParameterEdit(
            "Topology (.tpr/.gro/.pdb)",
            default="",
            docstring=(
                "Topology for the trajectory — a .tpr is richest. Leave "
                "empty to auto-discover (.tpr preferred, else .gro) from "
                "a predecessor MD-run node."
            ),
            optional=True,
        ),
        "peptide_selection": StringParameter(
            "Peptide Selection",
            default="protein",
            docstring=(
                "MDAnalysis selection string for the residues to assess "
                "(DSSP is meaningful only on protein). Default 'protein' "
                "catches the peptide and excludes the SnP fragment + "
                "solvent."
            ),
        ),
        "stride": IntegerParameter(
            "Frame Stride",
            default=1,
            docstring="Analyse every Nth frame (>=1).",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for <case>_helix.csv.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting helix analysis...", node_id=self.node_id,
                   progress=0)
        if run_helix_analysis is None:
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
                _from_predecessors(predecessor_data, "tpr", "structure",
                                   "output_tpr", "output_gro")
            if not traj:
                raise NodeException("setup",
                    "No trajectory — set 'Trajectory' or connect an "
                    "MD-run predecessor.")
            if not topo:
                raise NodeException("setup",
                    "No topology — set 'Topology' or connect an MD-run "
                    "predecessor.")
            traj = self.resolve_path(traj)
            topo = self.resolve_path(topo)
            for label, p in (("trajectory", traj), ("topology", topo)):
                if not p or not os.path.isfile(p):
                    raise NodeException("setup", f"{label} not found: {p}")

            selection = flow_vars["peptide_selection"].get_value() or "protein"
            stride = max(1, int(flow_vars["stride"].get_value() or 1))
            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)

            stream_log(f"Running DSSP ({os.path.basename(traj)}, "
                       f"stride {stride})...", node_id=self.node_id,
                       progress=25)
            res = run_helix_analysis(topo, traj, selection, stride)

            csv_path = os.path.join(out_dir, f"{case_name}_helix.csv")
            write_csv(csv_path, CSV_HEADER,
                      helix_csv_rows(res["frames"], res["times"],
                                     res["counts"], res["fracs"],
                                     res["n_residues"]))
            summary = summarize(res["fracs"])

            stream_log(
                f"Helix analysis: {summary.get('n_frames', 0)} frames, "
                f"mean helix {summary.get('mean_frac_helix', 0):.1%}",
                node_id=self.node_id, progress=100)

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(out_dir),
                "output_helix_csv": self.format_output_path(csv_path),
                "helix_summary": summary,
                "n_residues": res["n_residues"],
                "per_residue_propensity": dict(zip(
                    (f"{rn}{ri}" for rn, ri in
                     zip(res["resnames"], res["resids"])),
                    res["per_residue_propensity"])),
            }
            result.files["output"] = {
                "helix_csv": self.format_output_path(csv_path),
            }
            result.success = True
            result.message = (
                f"Helix content: mean "
                f"{summary.get('mean_frac_helix', 0):.1%} over "
                f"{summary.get('n_frames', 0)} frames "
                f"({res['n_residues']} residues)")
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("md analysis helix", str(e))
