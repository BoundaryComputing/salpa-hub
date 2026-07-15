"""EasyParm Library Generation Node — wraps tleap + 12_generate_lib.py

Generates an AMBER library (.lib) file using tleap from MOL2 + frcmod,
then fixes atomic numbers and bond connectivity for metal atoms.

Predecessor data flow: reads mol2, frcmod, and distance_type from
Force Field Assembly predecessor. Falls back to explicit OPTIONS or
default filenames in work_dir for backward compatibility.
"""

import os
import shutil
import subprocess
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import FileParameterEdit, FolderParameter, SelectParameter, StringParameter
from bocoflow_core.stream_logger import stream_log

EASYPARM_SCRIPTS = os.environ.get("EASYPARM_SCRIPTS", "")
_LEAPRC = {"gaff": "leaprc.gaff", "gaff2": "leaprc.gaff2", "amber": "leaprc.ff19SB"}


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
    """Resolve a path reference and copy to work_dir if needed. Returns local path or None."""
    if not ref:
        return None
    source = resolve_fn(ref)
    if not source or not os.path.isfile(source):
        return None
    dest = os.path.join(work_dir, filename)
    if os.path.abspath(source) != os.path.abspath(dest):
        shutil.copy2(source, dest)
    return dest


class EpLibraryGeneration(Node):
    """Generate AMBER library (.lib) file via tleap and fix atomic numbers
    and bond connectivity for metal atoms."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Working Directory",
            docstring="Working directory. Input files auto-discovered from predecessor or expected here.",
        ),
        "mol2_file": FileParameterEdit(
            "MOL2 File",
            default="",
            docstring="Labeled MOL2 from Force Field Assembly. Leave empty to auto-discover from predecessor.",
        ),
        "frcmod_file": FileParameterEdit(
            "Frcmod File",
            default="",
            docstring="Force field modification file. Leave empty to auto-discover from predecessor.",
        ),
        "distance_type_file": FileParameterEdit(
            "Distance Type File",
            default="",
            docstring="Bond types from Bond Detection. Leave empty to auto-discover from predecessor.",
        ),
        "atom_type": SelectParameter(
            "Atom Type Scheme",
            default="gaff",
            options=["gaff", "gaff2", "amber"],
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting library generation...", node_id=self.node_id, progress=0)

        try:
            result = NodeResult()

            # --- Read parameters with predecessor fallback ---
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}
            case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "complex")
            work_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            at_type = flow_vars["atom_type"].get_value()

            # --- 3-tier file resolution: explicit config → predecessor → default ---

            # MOL2 file
            mol2 = self.resolve_path(flow_vars["mol2_file"].get_value()) or ""
            if not mol2:
                ref = _get_from_predecessors(predecessor_data, "output_mol2")
                if ref:
                    mol2 = _ensure_in_workdir(self.resolve_path, ref, work_dir, "COMPLEX.mol2") or ""
            if not mol2:
                mol2 = "COMPLEX.mol2"

            # Frcmod file
            frcmod = self.resolve_path(flow_vars["frcmod_file"].get_value()) or ""
            if not frcmod:
                ref = _get_from_predecessors(predecessor_data, "output_frcmod")
                if ref:
                    frcmod = _ensure_in_workdir(self.resolve_path, ref, work_dir, "COMPLEX.frcmod") or ""
            if not frcmod:
                frcmod = "COMPLEX.frcmod"

            # Distance type file
            dist_type = self.resolve_path(flow_vars["distance_type_file"].get_value()) or ""
            if not dist_type:
                ref = _get_from_predecessors(predecessor_data, "output_distance_type")
                if ref:
                    dist_type = _ensure_in_workdir(self.resolve_path, ref, work_dir, "distance_type.dat") or ""
            if not dist_type:
                dist_type = "distance_type.dat"

            for name, path in [("MOL2", mol2), ("Frcmod", frcmod), ("Distance type", dist_type)]:
                full = os.path.join(work_dir, path) if not os.path.isabs(path) else path
                if not os.path.exists(full):
                    raise NodeException("setup",
                        f"{name} file not found: {full}\n"
                        f"Connect to Force Field Assembly or set the file path explicitly.")

            scripts = _find_scripts_dir(self._node_dir or ".")
            leaprc = _LEAPRC.get(at_type, "leaprc.gaff")

            # Prefer NEW_COMPLEX.mol2 (ULS-renamed types) + Hybridization_Info.dat
            # so the saved lib's atom types match the frcmod's ULS labels.
            new_mol2_ref = _get_from_predecessors(predecessor_data, "output_new_mol2")
            new_mol2_path = _ensure_in_workdir(self.resolve_path, new_mol2_ref, work_dir, "NEW_COMPLEX.mol2") if new_mol2_ref else None
            if not new_mol2_path:
                candidate = os.path.join(work_dir, "NEW_COMPLEX.mol2")
                if os.path.exists(candidate):
                    new_mol2_path = candidate

            hybrid_ref = _get_from_predecessors(predecessor_data, "output_hybridization_info")
            hybrid_path = _ensure_in_workdir(self.resolve_path, hybrid_ref, work_dir, "Hybridization_Info.dat") if hybrid_ref else None
            if not hybrid_path:
                candidate = os.path.join(work_dir, "Hybridization_Info.dat")
                if os.path.exists(candidate):
                    hybrid_path = candidate

            tleap_mol2 = "NEW_COMPLEX.mol2" if new_mol2_path else mol2
            # Inline the addAtomTypes block only when it actually declares
            # new types. An empty `addAtomTypes { }` is a tleap error.
            add_atom_types_block = ""
            if hybrid_path:
                with open(hybrid_path) as hf:
                    body = hf.read().strip()
                between_braces = ""
                if "{" in body and "}" in body:
                    between_braces = body.split("{", 1)[1].rsplit("}", 1)[0].strip()
                if between_braces:
                    add_atom_types_block = body + "\n"

            stream_log("Writing tleap input...", node_id=self.node_id, progress=10)
            tleap_input = os.path.join(work_dir, "input_library.tleap")
            with open(tleap_input, "w") as f:
                f.write(f"source {leaprc}\n")
                if add_atom_types_block:
                    f.write(add_atom_types_block)
                f.write(f"loadamberparams {frcmod}\n")
                f.write(f'mol = loadmol2 "{tleap_mol2}"\n')
                f.write("check mol\n")
                f.write("charge mol\n")
                f.write("savepdb mol COMPLEX.pdb\n")
                f.write("saveoff mol COMPLEX.lib\n")
                f.write("quit\n")

            stream_log("Running tleap...", node_id=self.node_id, progress=40)
            tl = subprocess.run(
                ["tleap", "-f", tleap_input],
                cwd=work_dir, capture_output=True, text=True,
            )
            # tleap routinely exits nonzero on benign warnings, so success is
            # judged by whether it wrote a non-empty COMPLEX.lib rather than by
            # exit code. Keep the log to diagnose a genuine failure.
            tleap_log = ((tl.stdout or "") + (tl.stderr or "")).strip()
            raw_lib = os.path.join(work_dir, "COMPLEX.lib")
            if not (os.path.exists(raw_lib) and os.path.getsize(raw_lib) > 0):
                raise NodeException(
                    "library generation",
                    f"tleap did not produce a COMPLEX.lib (exit {tl.returncode}). "
                    f"This usually means a missing/invalid frcmod parameter or an "
                    f"unrecognized mol2 atom type. tleap log:\n{tleap_log[-2000:]}",
                )

            stream_log("Fixing library file...", node_id=self.node_id, progress=70)
            fix = subprocess.run(
                [sys.executable, os.path.join(scripts, "12_generate_lib.py")],
                cwd=work_dir, capture_output=True, text=True,
            )
            if fix.returncode != 0:
                tail = (fix.stderr or fix.stdout or "").strip()[-2000:]
                raise NodeException(
                    "library generation",
                    f"12_generate_lib.py failed (exit {fix.returncode}):\n{tail}",
                )

            lib_file = os.path.join(work_dir, "COMPLEX.lib")
            if not (os.path.exists(lib_file) and os.path.getsize(lib_file) > 0):
                raise NodeException(
                    "library generation",
                    "COMPLEX.lib is missing or empty after the tleap + "
                    "12_generate_lib fix step.",
                )

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(work_dir),
                "output_lib": self.format_output_path(lib_file),
            }
            result.files["output"] = {
                "lib": self.format_output_path(lib_file),
            }
            result.success = True
            result.message = "Library generation complete — all files ready"

            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("library generation", str(e))
