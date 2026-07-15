"""EasyParm AMBER → GROMACS Node — wraps ParmEd's prmtop/rst7 → top/gro
conversion.

Mirrors easyPARM's ``amber_converter.py:71-83``: load the AMBER pair
with ParmEd, save as GROMACS top + gro. ParmEd preserves all custom
atom types, cross-FF linkage parameters, and metal MASS+NONBON entries
that landed in the prmtop, so the GROMACS output is parametrically
equivalent.

Predecessor data flow: auto-discovers ``output_prmtop`` and
``output_rst7`` from a fragment-fuse / library-generation predecessor.
Explicit OPTIONS override.

Forwarded data:
  ``output_top``, ``output_gro``, ``case_name``, ``working_path``.
"""

import os
import shutil

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FloatParameter, FolderParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import convert_amber_to_gromacs
except ImportError:
    try:
        from core import convert_amber_to_gromacs  # type: ignore
    except ImportError:  # server-side introspection (no heavy deps yet)
        convert_amber_to_gromacs = None


def _get_from_predecessors(predecessor_data, key):
    for pred in (predecessor_data or []):
        if pred and key in pred:
            return pred[key]
    return None


def _ensure_in_workdir(resolve_fn, ref, work_dir, filename):
    """Stage a referenced file into work_dir; return local path or None."""
    if not ref:
        return None
    source = resolve_fn(ref)
    if not source or not os.path.isfile(source):
        return None
    dest = os.path.join(work_dir, filename)
    if os.path.abspath(source) != os.path.abspath(dest):
        shutil.copy2(source, dest)
    return dest


class EpAmberToGromacs(Node):
    """Convert AMBER prmtop+rst7 → GROMACS top+gro via ParmEd."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for the GROMACS .top and .gro files.",
        ),
        "prmtop_file": FileParameterEdit(
            "AMBER Prmtop",
            default="",
            docstring=(
                "AMBER topology file (.prmtop / .parm7). Leave empty to "
                "auto-discover ``output_prmtop`` from an upstream fuse / "
                "library-generation node."
            ),
            optional=True,
        ),
        "rst7_file": FileParameterEdit(
            "AMBER Restart",
            default="",
            docstring=(
                "AMBER coordinate file (.rst7 / .inpcrd / .ncrst). Leave "
                "empty to auto-discover ``output_rst7`` from a predecessor."
            ),
            optional=True,
        ),
        "output_prefix": StringParameter(
            "Output Prefix",
            default="complex",
            docstring=(
                "Basename for the produced files: ``<prefix>.top`` and "
                "``<prefix>.gro`` (written into the output directory)."
            ),
        ),
        "itp_filename": StringParameter(
            "ITP Filename",
            default="",
            docstring=(
                "Optional. When set, also write a self-contained "
                "``<itp_filename>.itp`` containing the metallopeptide's "
                "[ moleculetype ] block (atoms, bonds, angles, dihedrals, "
                "exclusions). The master .top then ``#include``s this "
                "file instead of inlining the molecule, so you can drop "
                "the .itp into a master topology that mixes the "
                "metallopeptide with water, ions, lipids, or other "
                "molecules. Leave empty to keep the v1.10.0 monolithic "
                ".top (no .itp). The .itp lives next to the .top in the "
                "output directory; basename only (no path separators). "
                "A trailing ``.itp`` is stripped before appending — both "
                "'complex' and 'complex.itp' produce complex.itp."
            ),
            optional=True,
        ),
        "add_box_if_absent": BooleanParameter(
            "Add Cubic Box if Non-Periodic",
            default=True,
            docstring=(
                "If the prmtop has no periodic box (the typical fuse "
                "output), add a cubic box sized to the molecule extent + "
                "2 × box padding. Matches easyPARM's amber_converter. "
                "Disable if you'll solvate later via gmx editconf / gmx "
                "solvate."
            ),
        ),
        "box_padding": FloatParameter(
            "Box Padding (Å)",
            default=10.0,
            docstring=(
                "Half-width of cubic box padding around the molecule. "
                "Only applied when 'Add Cubic Box if Non-Periodic' is on."
            ),
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting AMBER → GROMACS conversion...",
                   node_id=self.node_id, progress=0)

        if convert_amber_to_gromacs is None:
            raise NodeException(
                "setup",
                "ep_amber_to_gromacs core.py could not be imported — run "
                "this node in the metalparm_vwf pixi env (ParmEd required).",
            )

        try:
            result = NodeResult()
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}

            case_name = (flow_vars["case_name"].get_value()
                         or input_data.get("case_name", "complex"))
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            output_prefix = flow_vars["output_prefix"].get_value() or "complex"
            add_box = bool(flow_vars["add_box_if_absent"].get_value())
            padding = float(flow_vars["box_padding"].get_value() or 10.0)
            itp_filename = (flow_vars["itp_filename"].get_value() or "").strip()

            # 3-tier resolution: explicit > predecessor > error
            prmtop = self.resolve_path(flow_vars["prmtop_file"].get_value()) or ""
            if not prmtop:
                ref = _get_from_predecessors(predecessor_data, "output_prmtop")
                if ref:
                    prmtop = _ensure_in_workdir(
                        self.resolve_path, ref, output_dir, "complex.prmtop") or ""
            if not prmtop:
                raise NodeException(
                    "setup",
                    "prmtop not provided and no output_prmtop in predecessor data."
                )

            rst7 = self.resolve_path(flow_vars["rst7_file"].get_value()) or ""
            if not rst7:
                ref = _get_from_predecessors(predecessor_data, "output_rst7")
                if ref:
                    rst7 = _ensure_in_workdir(
                        self.resolve_path, ref, output_dir, "complex.rst7") or ""
            if not rst7:
                raise NodeException(
                    "setup",
                    "rst7 not provided and no output_rst7 in predecessor data."
                )

            stream_log(
                f"Loading {os.path.basename(prmtop)} + {os.path.basename(rst7)}...",
                node_id=self.node_id, progress=20,
            )

            full_prefix = os.path.join(output_dir, output_prefix)
            stats = convert_amber_to_gromacs(
                prmtop_path=prmtop,
                rst7_path=rst7,
                output_prefix=full_prefix,
                add_box_if_absent=add_box,
                box_padding=padding,
                itp_filename=itp_filename or None,
            )

            box_msg = (f"box added ({padding} Å padding)" if stats["box_added"]
                       else "box from prmtop" if stats["has_box"]
                       else "no box (non-periodic)")
            itp_msg = (f"; ITP: {os.path.basename(stats['itp'])}"
                       if stats.get("itp") else "")
            stream_log(
                f"Converted {stats['n_atoms']} atoms; {box_msg}{itp_msg}.",
                node_id=self.node_id, progress=90,
            )

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "output_top": self.format_output_path(stats["top"]),
                "output_gro": self.format_output_path(stats["gro"]),
                "output_itp": (self.format_output_path(stats["itp"])
                               if stats.get("itp") else None),
                "n_atoms": stats["n_atoms"],
            }
            result.files["output"] = {
                "top": self.format_output_path(stats["top"]),
                "gro": self.format_output_path(stats["gro"]),
            }
            if stats.get("itp"):
                result.files["output"]["itp"] = self.format_output_path(
                    stats["itp"])
            result.success = True
            itp_suffix = (f" + {os.path.basename(stats['itp'])}"
                          if stats.get("itp") else "")
            result.message = (
                f"AMBER → GROMACS: {stats['n_atoms']} atoms, "
                f"{output_prefix}.top + {output_prefix}.gro{itp_suffix} ready"
            )
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("amber → gromacs conversion", str(e))
