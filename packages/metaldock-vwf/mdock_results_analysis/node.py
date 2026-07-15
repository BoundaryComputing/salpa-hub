"""MetalDock Results Analysis — wraps metaldock_modules.results_analysis.

Parse binding energies from the DLG, compute ligand efficiencies, find
protein residues contacting each pose, and (optionally) RMSD vs a reference
geometry. Reads the DLG + pose list + cleaned protein PDB + heavy-atom count
forwarded by the upstream nodes. Terminal node of the metaldock pipeline.
"""

import json
import os
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FloatParameter, FolderParameter, IntegerParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log

_NODE_DIR = os.path.dirname(os.path.abspath(__file__))


def _ensure_metaldock_modules():
    try:
        import metaldock_modules  # noqa: F401  installed in env / already on path
        return
    except ImportError:
        pass
    candidates = [
        os.environ.get("METALDOCK_SRC"),
        os.path.join(_NODE_DIR, "scripts"),
        os.path.join(_NODE_DIR, "..", "_vendor"),   # bundled with the package (release)
        os.path.abspath(os.path.join(_NODE_DIR, "..", "..", "..", "src")),
    ]
    for c in candidates:
        if c and os.path.isdir(os.path.join(c, "metaldock_modules")):
            if c not in sys.path:
                sys.path.insert(0, c)
            return c
    raise NodeException("setup", "Cannot locate metaldock_modules. Set METALDOCK_SRC.")


def _merge_predecessors(predecessor_data):
    merged = {}
    for pred in (predecessor_data or []):
        if not isinstance(pred, dict):
            continue
        scope = pred.get("data") if isinstance(pred.get("data"), dict) else pred
        for k, v in scope.items():
            if k not in ("success", "message", "metadata", "files"):
                merged[k] = v
    return merged


class MdockResultsAnalysis(Node):
    """Binding energies, ligand efficiency, interacting residues, RMSD."""

    category = "Metal Docking"
    tags = ["metaldock", "analysis", "binding-energy", "rmsd", "interacting-residues"]

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory", docstring="Working dir for the analysis summary JSON.",
        ),
        "dlg_path": FileParameterEdit(
            "Docking Log (.dlg)", default="", optional=True,
            docstring="Leave empty to auto-discover `dlg_path` from a predecessor.",
        ),
        "protein_pdb": FileParameterEdit(
            "Protein PDB", default="", optional=True,
            docstring="Cleaned protein for residue contacts. Leave empty to inherit "
                      "`cleaned_pdb` from a predecessor.",
        ),
        "reference_xyz": FileParameterEdit(
            "Reference XYZ (optional)", default="", optional=True,
            docstring="Crystal/reference geometry; if set, RMSD is computed per pose.",
        ),
        "n_heavy_atoms": IntegerParameter(
            "Heavy Atom Count", default=0, optional=True,
            docstring="Heavy atoms in the complex (for ligand efficiency). Leave 0 to "
                      "inherit `n_heavy_atoms` from a predecessor.",
        ),
        "cutoff": FloatParameter(
            "Contact Cutoff (Å)", default=4.0,
            docstring="Distance cutoff for interacting-residue detection.",
        ),
        "num_poses": IntegerParameter(
            "Poses to Analyze", default=0, optional=True,
            docstring="Max poses to analyze. Leave 0 for all available poses.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting results analysis...", node_id=self.node_id, progress=0)
        try:
            _ensure_metaldock_modules()
            from pathlib import Path
            from metaldock_modules import results_analysis

            carried = _merge_predecessors(predecessor_data)
            result = NodeResult()

            case_name = (
                flow_vars["case_name"].get_value() or carried.get("case_name") or "complex"
            )
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            dlg_ref = flow_vars["dlg_path"].get_value() or carried.get("dlg_path")
            if not dlg_ref:
                raise NodeException("setup", "No dlg_path (param or predecessor).")
            dlg_path = Path(self.resolve_path(dlg_ref))

            protein_ref = flow_vars["protein_pdb"].get_value() or carried.get("cleaned_pdb")
            if not protein_ref:
                raise NodeException("setup", "No protein PDB (param or predecessor cleaned_pdb).")
            protein_pdb = Path(self.resolve_path(protein_ref))

            pose_refs = carried.get("pose_xyz_paths") or []
            pose_xyz_paths = [Path(self.resolve_path(p)) for p in pose_refs]
            if not pose_xyz_paths:
                raise NodeException("setup", "No pose_xyz_paths found in predecessor data.")

            n_heavy = int(flow_vars["n_heavy_atoms"].get_value() or 0) or int(
                carried.get("n_heavy_atoms") or 0
            )
            cutoff = float(flow_vars["cutoff"].get_value() or 4.0)
            num_poses = int(flow_vars["num_poses"].get_value() or 0) or None

            ref_xyz_ref = flow_vars["reference_xyz"].get_value()
            reference_xyz = Path(self.resolve_path(ref_xyz_ref)) if ref_xyz_ref else None

            stream_log("Parsing energies + residue contacts...",
                       node_id=self.node_id, progress=40)
            analysis = results_analysis.analyze_docking_results(
                dlg_path=dlg_path,
                pose_xyz_paths=pose_xyz_paths,
                protein_pdb=protein_pdb,
                n_heavy_atoms=n_heavy,
                num_poses=num_poses,
                reference_xyz=reference_xyz,
                cutoff=cutoff,
            )

            # interacting_residues is list[list[(name, id)]] — make JSON-clean
            residues = [
                [[str(name), str(rid)] for (name, rid) in pose_res]
                for pose_res in analysis.get("interacting_residues", [])
            ]

            summary = {
                "case_name": case_name,
                "binding_energies": analysis.get("binding_energies", []),
                "binding_efficiencies": analysis.get("binding_efficiencies", []),
                "interacting_residues": residues,
                "n_heavy_atoms": n_heavy,
            }
            if "rmsd_values" in analysis:
                summary["rmsd_values"] = analysis["rmsd_values"]
                summary["rmsd_stats"] = analysis["rmsd_stats"]

            summary_path = Path(output_dir) / f"{case_name}_analysis.json"
            summary_path.write_text(json.dumps(summary, indent=2))

            best_energy = min(summary["binding_energies"]) if summary["binding_energies"] else None

            result.data = {
                **carried,
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "analysis_json": self.format_output_path(str(summary_path)),
                "binding_energies": summary["binding_energies"],
                "binding_efficiencies": summary["binding_efficiencies"],
                "interacting_residues": residues,
                "best_binding_energy": best_energy,
            }
            if "rmsd_values" in summary:
                result.data["rmsd_values"] = summary["rmsd_values"]
                result.data["rmsd_stats"] = summary["rmsd_stats"]
            result.metadata["case_name"] = case_name
            result.files["output"] = {"analysis_json": result.data["analysis_json"]}
            result.success = True
            n_res = len(residues[0]) if residues else 0
            result.message = (
                f"Best ΔG {best_energy:.2f} kcal/mol, {n_res} contacting residues (pose 1)"
                if best_energy is not None
                else "Analysis complete"
            )
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("results analysis", str(e))
