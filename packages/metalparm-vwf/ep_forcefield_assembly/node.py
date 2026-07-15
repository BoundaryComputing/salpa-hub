"""EasyParm Force Field Assembly Node — wraps steps 03-13

Runs the full force field assembly chain: MOL2 correction, parmchk2 frcmod,
Unique Labeling Strategy, merge Seminario parameters, remove generic metal
entries, deduplicate, UFF gap-fill, final cleanup.

Predecessor data flow: reads mol2 (from MOL2 Generation), force constants
(from Seminario), and distance files (forwarded from Bond Detection through
MOL2 Generation). Falls back to explicit OPTIONS or default filenames in
work_dir for backward compatibility.
"""

import os
import shutil
import subprocess
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FolderParameter, SelectParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

EASYPARM_SCRIPTS = os.environ.get("EASYPARM_SCRIPTS", "")
_AT_FLAG = {"gaff": "1", "gaff2": "2", "amber": "3"}


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


class EpForcefieldAssembly(Node):
    """Assemble AMBER force field: MOL2 correction, parmchk2, ULS labeling,
    merge Seminario params, UFF fill, cleanup (easyPARM steps 03-13)."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Working Directory",
            docstring="Directory where output files will be written and intermediate files are read from",
        ),
        "mol2_file": FileParameterEdit(
            "MOL2 File",
            default="",
            docstring="Input MOL2 from antechamber. Leave empty to auto-discover from MOL2 Generation predecessor.",
        ),
        "force_constants_file": FileParameterEdit(
            "Force Constants File",
            default="",
            docstring="Seminario force constants. Leave empty to auto-discover from Seminario predecessor.",
        ),
        "distance_file": FileParameterEdit(
            "Distance File",
            default="",
            docstring="Bond distances from bond detection. Leave empty to auto-discover from predecessor.",
        ),
        "distance_type_file": FileParameterEdit(
            "Distance Type File",
            default="",
            docstring="Bond types from bond detection. Leave empty to auto-discover from predecessor.",
        ),
        "atom_type": SelectParameter(
            "Atom Type Scheme",
            default="gaff",
            options=["gaff", "gaff2", "amber"],
        ),
        "multi_metal": BooleanParameter(
            "Multi-Metal System",
            default=False,
            docstring="Enable for systems with connected metal groups",
        ),
        "use_uls": BooleanParameter(
            "Use Unique Labeling Strategy (ULS)",
            default=False,
            docstring="Rename metal-coordinating atom types to unique labels "
            "(n1, o2, …). Off by default — the lib and frcmod use GAFF types "
            "throughout and match by construction. Enable only when combining "
            "multiple parameterized metal complexes in one MD system that "
            "would otherwise share atom-type names.",
        ),
        "scale_fc": BooleanParameter(
            "Apply easyPARM Force-Constant Scaling",
            default=True,
            docstring="Apply the empirical BOND/ANGLE multipliers from "
            "vanilla easyPARM (lines 1269-1299 of 01_easyPARM.sh). Boosts "
            "weak Seminario constants for dative metal bonds: BOND k<20 "
            "× 4.6; ANGLE k<5 × 11.6, k<10 × 7.8, k<20 × 3.6, k<29 × 2.7. "
            "On by default to match upstream easyPARM behavior; turn off "
            "for raw Seminario-averaged constants.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting force field assembly...", node_id=self.node_id, progress=0)

        try:
            result = NodeResult()
            case_name = flow_vars["case_name"].get_value() or "complex"
            work_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            at_type = flow_vars["atom_type"].get_value()
            multi_metal = flow_vars["multi_metal"].get_value()
            use_uls = flow_vars["use_uls"].get_value()
            scale_fc = flow_vars["scale_fc"].get_value()

            # Inherit case_name from predecessor if not set
            if case_name == "complex":
                pred_name = _get_from_predecessors(predecessor_data, "case_name")
                if pred_name:
                    case_name = pred_name

            # --- 3-tier file resolution: explicit config → predecessor → default ---

            # MOL2 file
            mol2 = self.resolve_path(flow_vars["mol2_file"].get_value()) or ""
            if not mol2:
                ref = _get_from_predecessors(predecessor_data, "output_mol2")
                if ref:
                    mol2 = _ensure_in_workdir(self.resolve_path, ref, work_dir, "COMPLEX.mol2") or ""
            if not mol2:
                mol2 = "COMPLEX.mol2"  # default in work_dir

            # Force constants file
            fc_file = self.resolve_path(flow_vars["force_constants_file"].get_value()) or ""
            if not fc_file:
                ref = _get_from_predecessors(predecessor_data, "output_force_constants")
                if ref:
                    fc_file = _ensure_in_workdir(self.resolve_path, ref, work_dir, "bond_angle_dihedral_data.dat") or ""
            if not fc_file:
                fc_file = "bond_angle_dihedral_data.dat"

            # Distance file
            dist_file = self.resolve_path(flow_vars["distance_file"].get_value()) or ""
            if not dist_file:
                ref = _get_from_predecessors(predecessor_data, "output_distance")
                if ref:
                    dist_file = _ensure_in_workdir(self.resolve_path, ref, work_dir, "distance.dat") or ""
            if not dist_file:
                dist_file = "distance.dat"

            # Distance type file
            dist_type = self.resolve_path(flow_vars["distance_type_file"].get_value()) or ""
            if not dist_type:
                ref = _get_from_predecessors(predecessor_data, "output_distance_type")
                if ref:
                    dist_type = _ensure_in_workdir(self.resolve_path, ref, work_dir, "distance_type.dat") or ""
            if not dist_type:
                dist_type = "distance_type.dat"

            # Verify required input files exist
            for name, path in [("MOL2", mol2), ("Force constants", fc_file),
                               ("Distance", dist_file), ("Distance type", dist_type)]:
                full = os.path.join(work_dir, path) if not os.path.isabs(path) else path
                if not os.path.exists(full):
                    raise NodeException("setup",
                        f"{name} file not found: {full}\n"
                        f"Connect upstream nodes or set the file path explicitly.")

            scripts = _find_scripts_dir(self._node_dir or ".")
            py = sys.executable
            s_flag = _AT_FLAG.get(at_type, "1")

            def run_script(name, *args):
                proc = subprocess.run(
                    [py, os.path.join(scripts, name)] + list(args),
                    cwd=work_dir, capture_output=True, text=True,
                )
                if proc.returncode != 0:
                    tail = (proc.stderr or proc.stdout or "").strip()[-2000:]
                    raise NodeException(
                        "forcefield assembly",
                        f"{name} failed (exit {proc.returncode}):\n{tail}",
                    )

            # Step 03: MOL2 correction
            stream_log("[03] MOL2 correction...", node_id=self.node_id, progress=5)
            run_script("03_correct_mol2.py", work_dir)

            # Step 04: parmchk2
            stream_log("[04] parmchk2 frcmod...", node_id=self.node_id, progress=15)
            run_script("generate_preforcefield.py", work_dir)
            pc = subprocess.run(
                ["parmchk2",
                 "-i", "COMPLEX_modified.mol2", "-f", "mol2",
                 "-o", "COMPLEX_modified.frcmod", "-s", s_flag],
                cwd=work_dir, capture_output=True, text=True,
            )
            # Judge parmchk2 by whether it wrote a non-empty frcmod rather than
            # by exit code alone, then surface its log if it didn't. Without
            # this the chain silently continued on a parmchk2 failure and
            # produced a truncated/empty final frcmod while reporting success.
            mod_frcmod = os.path.join(work_dir, "COMPLEX_modified.frcmod")
            if not (os.path.exists(mod_frcmod) and os.path.getsize(mod_frcmod) > 0):
                tail = (pc.stderr or pc.stdout or "").strip()[-2000:]
                raise NodeException(
                    "forcefield assembly",
                    f"parmchk2 produced no COMPLEX_modified.frcmod (exit "
                    f"{pc.returncode}) — check that AmberTools is on PATH and "
                    f"COMPLEX_modified.mol2 is valid:\n{tail}",
                )
            run_script("generate_preforcefield.py", work_dir)

            # Step 05: ULS rename (opt-in) or pass-through (default)
            if use_uls:
                stream_log("[05] Unique Labeling Strategy (ULS rename)...", node_id=self.node_id, progress=25)
                uls_script = "05_prepare_mol2_frcmod_more_atom.py" if multi_metal else "05_prepare_mol2_frcmod.py"
            else:
                stream_log("[05] Pass-through (no-ULS, GAFF types kept)...", node_id=self.node_id, progress=25)
                uls_script = "05_prepare_mol2_frcmod_passthrough.py"
            run_script(uls_script, work_dir)

            # Steps 06-13: force field chain
            chain = [
                ("06", "06_get_atom_type.py", 35),
                ("07", "07_Seminario_forcefield.py", 45),
                ("08", "08_update_forcefield.py", 55),
                ("09", "09_clean_updatedforcefield.py", 65),
                ("10", "10_postclean_updatedforcefield.py", 72),
                ("11", "11_retrieve_uffdata.py", 80),
                ("13", "13_final_clean.py", 88),
            ]
            for step_num, script_name, pct in chain:
                stream_log(f"[{step_num}] {script_name}...", node_id=self.node_id, progress=pct)
                run_script(script_name, work_dir)

            # Rename final frcmod
            frcmod = os.path.join(work_dir, "COMPLEX.frcmod")
            for name in ["filtered_COMPLEX_modified2.frcmod",
                         "updated_updated_COMPLEX_modified2.frcmod"]:
                candidate = os.path.join(work_dir, name)
                if os.path.exists(candidate) and not os.path.exists(frcmod):
                    os.rename(candidate, frcmod)

            # The frcmod chain (steps 03-13) must have produced one of the
            # rename candidates; if neither exists, the chain failed silently
            # upstream (e.g. missing Seminario .dat inputs) and reporting
            # success here would emit a phantom output_frcmod path.
            if not os.path.exists(frcmod):
                raise NodeException(
                    "forcefield assembly",
                    "COMPLEX.frcmod was not produced — the frcmod chain "
                    "(steps 03-13) emitted neither filtered_COMPLEX_modified2.frcmod "
                    "nor updated_updated_COMPLEX_modified2.frcmod. Check the "
                    "Seminario force-constant inputs and the parmchk2 output.",
                )

            # Append UFF NONBON for metal elements (closes the gap that
            # parmchk2 doesn't emit metal NONBON and step 11 only patches
            # existing entries — without this, fresh-tleap saveamberparm
            # fails with "could not find vdW for type (Sn|Ru|Zn|...)").
            stream_log("[metal_vdw] Append UFF NONBON for metals...",
                       node_id=self.node_id, progress=92)
            run_script("metal_nonbon_fill.py", work_dir, "COMPLEX.frcmod")

            # Apply vanilla easyPARM force-constant scaling (the awk block
            # in 01_easyPARM.sh lines 1269-1299). Empirical multipliers
            # boost weak Seminario constants for dative metal bonds.
            if scale_fc:
                stream_log("[fc_scaling] Scale weak BOND/ANGLE force constants...",
                           node_id=self.node_id, progress=94)
                run_script("apply_fc_scaling.py", os.path.join(work_dir, "COMPLEX.frcmod"))

            # Forward both the antechamber-typed mol2 (COMPLEX.mol2) and the
            # ULS-renamed mol2 (NEW_COMPLEX.mol2) plus its addAtomTypes block.
            # Library Generation prefers NEW_COMPLEX.mol2 + Hybridization_Info.dat
            # so the saved lib's atom types match the frcmod's ULS labels.
            new_mol2 = os.path.join(work_dir, "NEW_COMPLEX.mol2")
            mol2_out = os.path.join(work_dir, "COMPLEX.mol2")
            hybrid_info = os.path.join(work_dir, "Hybridization_Info.dat")

            # Resolve distance_type path for forwarding
            dt_full = os.path.join(work_dir, dist_type) if not os.path.isabs(dist_type) else dist_type

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(work_dir),
                "output_frcmod": self.format_output_path(frcmod),
                "output_mol2": self.format_output_path(mol2_out),
                # Forward distance_type for Library Generation
                "output_distance_type": self.format_output_path(dt_full) if os.path.exists(dt_full) else "",
            }
            if os.path.exists(new_mol2):
                result.data["output_new_mol2"] = self.format_output_path(new_mol2)
            if os.path.exists(hybrid_info):
                result.data["output_hybridization_info"] = self.format_output_path(hybrid_info)
            result.files["output"] = {
                "frcmod": self.format_output_path(frcmod),
                "mol2": self.format_output_path(mol2_out),
            }
            if os.path.exists(new_mol2):
                result.files["output"]["new_mol2"] = self.format_output_path(new_mol2)
            result.success = True
            result.message = "Force field assembly complete"

            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("forcefield assembly", str(e))
