"""MD Analysis: Residue–Metal Distance — per-frame distance from a set
of probe atoms (the Tyr/His quenchers) to the metal centre (Sn) along
an MD trajectory.

Consumes a trajectory (.xtc/.trr) + topology (.tpr/.gro/.pdb) — set
explicitly or auto-discovered from a predecessor MD-run node — and
writes a CSV time series of minimum-image (PBC-aware) distances, in
Ångström, one column per probe atom plus the per-frame minimum. This
is the second half of the SnP-peptide case deliverable: the
"Tyr/His-OH ↔ Sn distance time series".

    gmx_mdrun(_local) → md_analysis_distance → <case>_distance.csv

Output keys: ``output_distance_csv``, ``distance_summary`` (closest
approach, per-frame-minimum statistics), ``case_name``, ``working_path``.
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
        distance_csv_header, distance_csv_rows, run_distance_analysis,
        summarize, write_csv,
    )
except ImportError:  # script-mode / server introspection fallback
    try:
        from core import (  # type: ignore
            distance_csv_header, distance_csv_rows, run_distance_analysis,
            summarize, write_csv,
        )
    except ImportError:
        distance_csv_header = distance_csv_rows = run_distance_analysis = None
        summarize = write_csv = None


def _from_predecessors(predecessor_data, *keys):
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


class MdAnalysisDistance(Node):
    """Per-frame distance from the Tyr/His probe atoms to the Sn centre."""

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
        "metal_selection": StringParameter(
            "Metal Selection",
            default="name Sn1 SN Sn",
            docstring=(
                "MDAnalysis selection for the metal centre — must "
                "resolve to exactly one atom. Default matches the SnP "
                "fragment's tin (named `Sn1` in the converted topology). "
                "If a .tpr is the topology, `type Sn` also works; narrow "
                "to e.g. 'resname mol and name Sn1' if ambiguous."
            ),
        ),
        "probe_selection": StringParameter(
            "Probe Selection",
            default=("(resname TYR and name OH) or "
                     "(resname HIS HID HIE HIP and name ND1 NE2)"),
            docstring=(
                "MDAnalysis selection for the probe atoms whose distance "
                "to the metal is tracked. Default: every tyrosine "
                "hydroxyl O and histidine ring N — the PCET quenchers."
            ),
        ),
        "stride": IntegerParameter(
            "Frame Stride",
            default=1,
            docstring="Analyse every Nth frame (>=1).",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for <case>_distance.csv.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting residue–metal distance analysis...",
                   node_id=self.node_id, progress=0)
        if run_distance_analysis is None:
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

            metal_sel = flow_vars["metal_selection"].get_value() or "name SN Sn"
            probe_sel = flow_vars["probe_selection"].get_value()
            if not probe_sel:
                raise NodeException("setup", "Probe Selection is empty.")
            stride = max(1, int(flow_vars["stride"].get_value() or 1))
            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)

            stream_log(f"Measuring distances ({os.path.basename(traj)}, "
                       f"stride {stride})...", node_id=self.node_id,
                       progress=25)
            try:
                res = run_distance_analysis(topo, traj, metal_sel, probe_sel,
                                            stride)
            except ValueError as ex:
                raise NodeException("execution", str(ex))

            csv_path = os.path.join(out_dir, f"{case_name}_distance.csv")
            write_csv(csv_path, distance_csv_header(res["labels"]),
                      distance_csv_rows(res["frames"], res["times"],
                                        res["dist_matrix"]))
            summary = summarize(res["frames"], res["dist_matrix"],
                                res["labels"])

            stream_log(
                f"Distance analysis: {summary.get('n_frames', 0)} frames, "
                f"closest approach {summary.get('closest_approach_A', 0)} Å "
                f"({summary.get('closest_probe', '?')})",
                node_id=self.node_id, progress=100)

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(out_dir),
                "output_distance_csv": self.format_output_path(csv_path),
                "distance_summary": summary,
                "probes": res["labels"],
            }
            result.files["output"] = {
                "distance_csv": self.format_output_path(csv_path),
            }
            result.success = True
            result.message = (
                f"Distance: closest approach "
                f"{summary.get('closest_approach_A', 0)} Å "
                f"({summary.get('closest_probe', '?')} @ frame "
                f"{summary.get('closest_frame', '?')}); "
                f"{summary.get('n_probes', 0)} probes, "
                f"{summary.get('n_frames', 0)} frames")
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("md analysis distance", str(e))
