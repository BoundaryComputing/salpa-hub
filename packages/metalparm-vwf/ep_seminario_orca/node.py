"""EasyParm Seminario Node (ORCA) — wraps Seminario_method_ORCA.py

Extracts bond and angle force constants from QM Hessian matrix using
the Modified Seminario Method. Also identifies equivalent atoms for
RESP charge constraints.

Predecessor data flow: reads distance/angle/dihedral files from Bond Detection
predecessor. Copies to output_dir if they're in a different directory.
"""

import os
import shutil
import subprocess
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import FileParameterEdit, FolderParameter, StringParameter
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
    """Resolve a path reference and copy to work_dir if needed."""
    if not ref:
        return None
    source = resolve_fn(ref)
    if not source or not os.path.isfile(source):
        return None
    dest = os.path.join(work_dir, filename)
    if os.path.abspath(source) != os.path.abspath(dest):
        shutil.copy2(source, dest)
    return dest


class EpSeminarioOrca(Node):
    """Compute force constants from ORCA Hessian via Seminario method."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "hess_file": FileParameterEdit(
            "ORCA Hessian File",
            docstring="ORCA .hess file from frequency calculation",
        ),
        "log_file": FileParameterEdit(
            "ORCA Output File",
            docstring="ORCA .out/.log file with Mulliken/CHELPG charges",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Working directory. Bond Detection files (distance/angle/dihedral) "
                      "are auto-discovered from predecessor or expected in this directory.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting Seminario method...", node_id=self.node_id, progress=0)

        try:
            result = NodeResult()

            # --- Read parameters with predecessor fallback ---
            #
            # 3-tier resolution for hess_file / log_file:
            #   1. Explicit config value in the node panel
            #   2. predecessor_data['output_hess'] / ['output_out']
            #      (ep_orca_run emits these keys in its result.data)
            #   3. None — surfaces a clear error before subprocess.run.
            #
            # Without tier 2, a workflow built as
            #   ep_orca_run → ep_seminario_orca
            # where the user never opens the seminario node panel to
            # type the hess/log paths (which is the natural flow when
            # ORCA hasn't completed yet at workflow-construction time)
            # crashes the seminario subprocess with "expected str,
            # bytes or os.PathLike object, not NoneType". See
            # snp-full-pipeline.spec.ts Phase 2 for the regression
            # scenario.
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}
            case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "complex")

            hess_explicit = flow_vars["hess_file"].get_value() or ""
            if hess_explicit:
                hess_file = self.resolve_path(hess_explicit)
            else:
                hess_ref = _get_from_predecessors(predecessor_data, "output_hess")
                hess_file = self.resolve_path(hess_ref) if hess_ref else None

            log_explicit = flow_vars["log_file"].get_value() or ""
            if log_explicit:
                log_file = self.resolve_path(log_explicit)
            else:
                log_ref = _get_from_predecessors(predecessor_data, "output_out")
                log_file = self.resolve_path(log_ref) if log_ref else None

            if not hess_file or not os.path.isfile(hess_file):
                raise NodeException(
                    "setup",
                    "ORCA Hessian File not resolvable. Either set 'ORCA Hessian File' "
                    "explicitly in the node panel, or connect to an upstream ep_orca_run "
                    "node whose result has output_hess populated (ep_orca_run populates "
                    "this only after the remote SLURM job completes and 'Check Status' "
                    "downloads artifacts).",
                )
            if not log_file or not os.path.isfile(log_file):
                raise NodeException(
                    "setup",
                    "ORCA Output File not resolvable. Either set 'ORCA Output File' "
                    "explicitly in the node panel, or connect to an upstream ep_orca_run "
                    "with output_out populated.",
                )

            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())

            # --- Get prerequisite files from predecessor or work_dir ---
            for key, filename in [
                ("output_distance", "distance.dat"),
                ("output_angle", "angle.dat"),
                ("output_dihedral", "dihedral.dat"),
            ]:
                ref = _get_from_predecessors(predecessor_data, key)
                if ref:
                    _ensure_in_workdir(self.resolve_path, ref, output_dir, filename)

            # Verify prerequisite files exist (from predecessor copy or shared dir)
            for req in ["distance.dat", "angle.dat", "dihedral.dat"]:
                if not os.path.exists(os.path.join(output_dir, req)):
                    raise NodeException("setup",
                        f"Missing {req} in output_dir — connect to Bond Detection "
                        f"node or ensure file exists in working directory")

            scripts = _find_scripts_dir(self._node_dir or ".")

            stream_log("Computing force constants from Hessian...", node_id=self.node_id, progress=20)

            proc = subprocess.run(
                [sys.executable,
                 os.path.join(scripts, "Seminario_method_ORCA.py"),
                 hess_file, log_file],
                cwd=output_dir,
                capture_output=True,
                text=True,
            )

            if proc.returncode != 0:
                raise NodeException("execution", f"Seminario failed: {proc.stderr}")

            # Declare outputs with explicit keys
            fc_path = os.path.join(output_dir, "bond_angle_dihedral_data.dat")
            similar_path = os.path.join(output_dir, "similar.dat")
            atomic_num_path = os.path.join(output_dir, "atomic_number.dat")

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "output_force_constants": self.format_output_path(fc_path) if os.path.exists(fc_path) else "",
                "output_similar": self.format_output_path(similar_path) if os.path.exists(similar_path) else "",
                "output_atomic_number": self.format_output_path(atomic_num_path) if os.path.exists(atomic_num_path) else "",
            }
            result.files["input"] = {
                "hess": self.format_output_path(hess_file),
                "log": self.format_output_path(log_file),
            }
            result.files["output"] = {k: v for k, v in result.data.items() if k.startswith("output_") and v}
            result.success = True
            result.message = "Seminario force constants computed"

            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("seminario", str(e))
