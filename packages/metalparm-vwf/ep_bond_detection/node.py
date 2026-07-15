"""EasyParm Bond Detection Node — wraps 02_get_bond_angle.py

Detects molecular connectivity (bonds, angles, dihedrals) from optimized
XYZ coordinates using covalent radii with metal-aware tolerance.
"""

import os
import subprocess
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import FileParameterEdit, FolderParameter, StringParameter
from bocoflow_core.stream_logger import stream_log

# easyPARM scripts location — set via env var or auto-detect from node package
EASYPARM_SCRIPTS = os.environ.get("EASYPARM_SCRIPTS", "")


def _find_scripts_dir(node_dir):
    """Find easyPARM scripts dir: env var > node-bundled scripts/ > metal-md source tree."""
    if EASYPARM_SCRIPTS and os.path.isdir(EASYPARM_SCRIPTS):
        return EASYPARM_SCRIPTS
    candidates = [
        os.path.join(node_dir, "scripts"),                                     # bundled with node (release)
        os.path.join(node_dir, "..", "..", "collect", "easyPARM", "scripts"),  # metal-md source fallback
    ]
    for c in candidates:
        if os.path.isdir(c):
            return os.path.abspath(c)
    raise NodeException("setup", "Cannot find easyPARM scripts directory. Set EASYPARM_SCRIPTS env var.")


class EpBondDetection(Node):
    """Detect bonds, angles, and dihedrals from XYZ geometry."""

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name",
            default="complex",
            docstring="Name for this parameterization case",
        ),
        "xyz_file": FileParameterEdit(
            "XYZ File",
            docstring="Optimized geometry in XYZ format (from QM calculation)",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for output files (distance.dat, angle.dat, etc.)",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting bond detection...", node_id=self.node_id, progress=0)

        try:
            result = NodeResult()
            case_name = flow_vars["case_name"].get_value() or "complex"
            xyz_file = self.resolve_path(flow_vars["xyz_file"].get_value())
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            scripts_dir = _find_scripts_dir(self._node_dir or ".")
            script = os.path.join(scripts_dir, "02_get_bond_angle.py")

            stream_log("Running 02_get_bond_angle.py...", node_id=self.node_id, progress=20)

            proc = subprocess.run(
                [sys.executable, script, xyz_file],
                cwd=output_dir,
                capture_output=True,
                text=True,
            )

            if proc.returncode != 0:
                raise NodeException("execution", f"02_get_bond_angle.py failed: {proc.stderr}")

            # Count detected features
            dist_path = os.path.join(output_dir, "distance.dat")
            angle_path = os.path.join(output_dir, "angle.dat")
            n_bonds = sum(1 for _ in open(dist_path)) if os.path.exists(dist_path) else 0
            n_angles = sum(1 for _ in open(angle_path)) if os.path.exists(angle_path) else 0

            # Declare outputs with explicit keys for predecessor data flow
            output_paths = {}
            file_map = {
                "output_distance": "distance.dat",
                "output_angle": "angle.dat",
                "output_dihedral": "dihedral.dat",
                "output_distance_type": "distance_type.dat",
                "output_metal_number": "metal_number.dat",
            }
            for key, fname in file_map.items():
                path = os.path.join(output_dir, fname)
                if os.path.exists(path):
                    output_paths[key] = self.format_output_path(path)

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "n_bonds": n_bonds,
                "n_angles": n_angles,
                **output_paths,
            }
            result.files["input"] = {"xyz": self.format_output_path(xyz_file)}
            result.files["output"] = {k: v for k, v in output_paths.items()}
            result.success = True
            result.message = f"Detected {n_bonds} bonds, {n_angles} angles"

            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("bond detection", str(e))
