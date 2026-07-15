"""EasyParm MOL2 Generation Node — wraps antechamber + atom type correction

Converts XYZ → PDB, runs AmberTools antechamber to generate MOL2,
then corrects atom types for metal-containing molecules.

Predecessor data flow: reads Bond Detection outputs (distance/angle files)
from predecessor_data when connected, copies to work_dir for correction scripts.
Forwards upstream file references for downstream nodes.
"""

import os
import shutil
import subprocess
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FolderParameter, IntegerParameter, SelectParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

EASYPARM_SCRIPTS = os.environ.get("EASYPARM_SCRIPTS", "")


def _find_scripts_dir(node_dir):
    if EASYPARM_SCRIPTS and os.path.isdir(EASYPARM_SCRIPTS):
        return EASYPARM_SCRIPTS
    for c in [os.path.join(node_dir, "scripts"),
              os.path.join(node_dir, "..", "..", "collect", "easyPARM", "scripts")]:
        if os.path.isdir(c):
            return os.path.abspath(c)
    raise NodeException("setup", "Cannot find easyPARM scripts. Set EASYPARM_SCRIPTS env var.")


def _get_from_predecessors(predecessor_data, key):
    """Search all predecessor outputs for a key."""
    for pred in (predecessor_data or []):
        if pred and key in pred:
            return pred[key]
    return None


def _ensure_in_workdir(resolve_fn, ref, work_dir, filename):
    """Resolve a path reference and copy to work_dir if needed. Returns local path."""
    if not ref:
        return None
    source = resolve_fn(ref)
    if not source or not os.path.isfile(source):
        return None
    dest = os.path.join(work_dir, filename)
    if os.path.abspath(source) != os.path.abspath(dest):
        shutil.copy2(source, dest)
    return dest


class EpMol2Generation(Node):
    """Generate MOL2 via antechamber with metal-aware atom type correction."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "xyz_file": FileParameterEdit(
            "XYZ File",
            docstring="Optimized geometry in XYZ format",
        ),
        "pdb_template": FileParameterEdit(
            "PDB Template (atom names)",
            default="",
            optional=True,
            docstring="Optional PDB whose atom names overlay the XYZ coords "
            "(atom-by-atom, in order). When omitted and a predecessor "
            "exposes 'output_pdb' (e.g. snp_builder), that is used. "
            "When neither, atoms are named <element><counter>.",
        ),
        "charge": IntegerParameter(
            "Total Charge",
            default=0,
            docstring="Total system charge",
        ),
        "multiplicity": IntegerParameter(
            "Spin Multiplicity",
            default=1,
            docstring="Spin multiplicity of the system",
        ),
        "atom_type": SelectParameter(
            "Atom Type Scheme",
            default="gaff",
            options=["gaff", "gaff2", "amber"],
            docstring="AMBER atom type scheme for antechamber",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for MOL2 and PDB output",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting MOL2 generation...", node_id=self.node_id, progress=0)

        try:
            result = NodeResult()

            # --- Read parameters with predecessor fallback ---
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}
            case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "complex")
            xyz_file = self.resolve_path(flow_vars["xyz_file"].get_value())
            pdb_template = self.resolve_path(flow_vars["pdb_template"].get_value() or "")
            if not pdb_template:
                ref = _get_from_predecessors(predecessor_data, "output_pdb")
                if ref:
                    pdb_template = self.resolve_path(ref)
            charge = flow_vars["charge"].get_value()
            mult = flow_vars["multiplicity"].get_value()
            at_type = flow_vars["atom_type"].get_value()
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            # --- Get Bond Detection files from predecessor (for correction scripts) ---
            # These are needed by atomtype_helper.py and atomtype_detector.py
            bond_det_files = {}
            for key, filename in [
                ("output_distance", "distance.dat"),
                ("output_angle", "angle.dat"),
                ("output_distance_type", "distance_type.dat"),
                ("output_dihedral", "dihedral.dat"),
                ("output_metal_number", "metal_number.dat"),
            ]:
                ref = _get_from_predecessors(predecessor_data, key)
                if ref:
                    local = _ensure_in_workdir(self.resolve_path, ref, output_dir, filename)
                    if local:
                        bond_det_files[key] = self.format_output_path(local)

            scripts = _find_scripts_dir(self._node_dir or ".")
            py = sys.executable

            def run_step(argv, label):
                proc = subprocess.run(
                    argv, cwd=output_dir, capture_output=True, text=True,
                )
                if proc.returncode != 0:
                    tail = (proc.stderr or proc.stdout or "").strip()[-2000:]
                    raise NodeException(
                        "mol2 generation",
                        f"{label} failed (exit {proc.returncode}):\n{tail}",
                    )

            # 1. XYZ → PDB (optionally overlaying atom names from a template
            # PDB so meaningful labels survive into the final lib).
            stream_log("Converting XYZ to PDB...", node_id=self.node_id, progress=10)
            pdb_file = os.path.join(output_dir, "COMPLEX.pdb")
            xyz_argv = [py, os.path.join(scripts, "xyz_to_pdb.py"), xyz_file, pdb_file]
            if pdb_template and os.path.isfile(pdb_template):
                xyz_argv += ["--template", pdb_template]
            run_step(xyz_argv, "xyz_to_pdb.py")

            # 2. Antechamber
            stream_log("Running antechamber...", node_id=self.node_id, progress=30)
            run_step(
                ["antechamber",
                 "-i", "COMPLEX.pdb", "-fi", "pdb",
                 "-o", "COMPLEX.mol2", "-fo", "mol2",
                 "-s", "2", "-rn", "mol",
                 "-nc", str(charge), "-m", str(mult),
                 "-at", at_type, "-dr", "no", "-j", "5"],
                "antechamber",
            )

            # 3. Atom type correction (GAFF/GAFF2) — needs distance/angle files from Bond Detection
            if at_type in ("gaff", "gaff2"):
                stream_log("Correcting atom types...", node_id=self.node_id, progress=60)
                run_step(
                    [py, os.path.join(scripts, "03_correct_mol2.py"), output_dir],
                    "03_correct_mol2.py",
                )
                run_step(
                    [py, os.path.join(scripts, "atomtype_helper.py"),
                     "COMPLEX.mol2", "distance_type.dat", "COMREF.mol2"],
                    "atomtype_helper.py",
                )
                run_step(
                    [py, os.path.join(scripts, "atomtype_detector.py"),
                     "COMREF.mol2", "distance.dat", "angle.dat"],
                    "atomtype_detector.py",
                )

            # 4. Revise special cases
            stream_log("Revising atom types...", node_id=self.node_id, progress=80)
            run_step(
                [py, os.path.join(scripts, "Revise_Atom_Type.py")],
                "Revise_Atom_Type.py",
            )

            mol2 = os.path.join(output_dir, "COMPLEX.mol2")
            pdb = os.path.join(output_dir, "COMPLEX.pdb")

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "output_mol2": self.format_output_path(mol2),
                "output_pdb": self.format_output_path(pdb),
                # Forward Bond Detection outputs for downstream nodes (FF Assembly, Library Gen)
                **bond_det_files,
            }
            if pdb_template and os.path.isfile(pdb_template):
                result.data["pdb_template_used"] = self.format_output_path(pdb_template)
            # Also forward any Bond Detection outputs we didn't copy (preserve full set)
            for key in ["output_distance", "output_angle", "output_dihedral",
                        "output_distance_type", "output_metal_number"]:
                if key not in result.data:
                    ref = _get_from_predecessors(predecessor_data, key)
                    if ref:
                        result.data[key] = ref

            result.files["output"] = {
                "mol2": self.format_output_path(mol2),
                "pdb": self.format_output_path(pdb),
            }
            result.success = True
            result.message = "MOL2 generation complete"

            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("mol2 generation", str(e))
