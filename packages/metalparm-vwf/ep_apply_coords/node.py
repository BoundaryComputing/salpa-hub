"""Apply Coordinates Node — write final complex.rst7 + complex.pdb from
a topology-only complex.prmtop plus aligned source PDBs.

Sister of ``ep_fragment_fuse_topology``. Input wiring:
  - ``output_prmtop`` from ep_fragment_fuse_topology
  - ``output_pdb`` from peptide_builder (peptide structure, NOT moved
    during alignment — peptide_builder writes it with its build coords)
  - ``output_fragment_pdb`` from fragment_align (the aligned fragment)
  - ``peptide_residues`` from ep_fragment_fuse_topology (or any
    upstream node that forwards it)

Predecessor data flow (3-tier explicit > role-specific > generic):
  - ``prmtop_file``:        explicit > ``output_prmtop`` > error
  - ``peptide_pdb``:        explicit > ``output_peptide_pdb`` > role-specific
                              ``peptide_pdb`` predecessor > error
  - ``fragment_pdb``:       explicit > ``output_fragment_pdb`` > error
  - ``peptide_residues``:   explicit > ``peptide_residues`` predecessor > error

Outputs forwarded:
  ``output_rst7``, ``output_pdb``, ``output_prmtop`` (forwarded
  unchanged from predecessor for downstream convenience),
  ``case_name``, ``working_path``.
"""
from __future__ import annotations

import os
import shutil

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FolderParameter, IntegerParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import apply_coords
except ImportError:
    try:
        from core import apply_coords  # type: ignore
    except ImportError:  # server-side introspection (no heavy deps yet)
        apply_coords = None


def _get_from_predecessors(predecessor_data, key):
    for pred in (predecessor_data or []):
        if pred and key in pred:
            return pred[key]
    return None


def _ensure_in_workdir(resolve_fn, ref, work_dir, filename):
    if not ref:
        return None
    source = resolve_fn(ref)
    if not source or not os.path.isfile(source):
        return None
    dest = os.path.join(work_dir, filename)
    if os.path.abspath(source) != os.path.abspath(dest):
        shutil.copy2(source, dest)
    return dest


class EpApplyCoords(Node):
    """Place aligned source-PDB coords onto a topology-only prmtop, write
    final complex.rst7 + complex.pdb."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for the final complex.rst7 + complex.pdb",
        ),
        "prmtop_file": FileParameterEdit(
            "Topology Prmtop",
            default="",
            docstring=(
                "AMBER topology from ep_fragment_fuse_topology. Leave "
                "empty to auto-discover ``output_prmtop`` from a predecessor."
            ),
            optional=True,
        ),
        "peptide_pdb": FileParameterEdit(
            "Peptide PDB",
            default="",
            docstring=(
                "Aligned peptide PDB. Typically the one written by "
                "peptide_builder (fragment_align doesn't move it). "
                "Leave empty to auto-discover ``output_peptide_pdb`` (or "
                "``output_pdb`` from a peptide_builder predecessor)."
            ),
            optional=True,
        ),
        "fragment_pdb": FileParameterEdit(
            "Fragment PDB",
            default="",
            docstring=(
                "Aligned fragment PDB (typically from fragment_align). "
                "Leave empty to auto-discover ``output_fragment_pdb``."
            ),
            optional=True,
        ),
        "peptide_residues": IntegerParameter(
            "Peptide Residue Count",
            default=0,
            docstring=(
                "Residue count of the peptide unit (1-based atom-resid "
                "switches from peptide PDB → fragment PDB at "
                "peptide_residues+1). Leave 0 to auto-inherit from a "
                "predecessor's ``peptide_residues``."
            ),
            optional=True,
        ),
        "output_prefix": StringParameter(
            "Output Prefix",
            default="complex",
            docstring=(
                "Basename for the produced files: ``<prefix>.rst7`` and "
                "``<prefix>.pdb`` (written into the output directory)."
            ),
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting coord-applier...", node_id=self.node_id, progress=0)

        if apply_coords is None:
            raise NodeException(
                "setup",
                "ep_apply_coords core.py could not be imported — run this "
                "node in the metalparm_vwf pixi env (ParmEd required).",
            )

        try:
            result = NodeResult()
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}

            case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "complex")
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)
            output_prefix = flow_vars["output_prefix"].get_value() or "complex"

            def _resolve(option_key: str, pred_keys: list[str], work_filename: str) -> str:
                explicit = self.resolve_path(flow_vars[option_key].get_value()) or ""
                if explicit and os.path.isfile(explicit):
                    return _ensure_in_workdir(self.resolve_path, explicit, output_dir,
                                              os.path.basename(explicit)) or explicit
                for pk in pred_keys:
                    ref = _get_from_predecessors(predecessor_data, pk)
                    if ref:
                        landed = _ensure_in_workdir(self.resolve_path, ref, output_dir,
                                                    work_filename)
                        if landed:
                            return landed
                raise NodeException(
                    "setup",
                    f"{option_key} not provided and no "
                    f"{' / '.join(pred_keys)} in predecessor data.",
                )

            prmtop = _resolve(
                "prmtop_file",
                ["output_prmtop"],
                "complex.prmtop",
            )
            peptide_pdb = _resolve(
                "peptide_pdb",
                ["output_peptide_pdb", "output_pdb"],
                "peptide.pdb",
            )
            fragment_pdb = _resolve(
                "fragment_pdb",
                ["output_fragment_pdb"],
                "fragment.pdb",
            )

            pep_size_raw = flow_vars["peptide_residues"].get_value()
            pep_size = int(pep_size_raw) if pep_size_raw else 0
            if pep_size <= 0:
                ref = _get_from_predecessors(predecessor_data, "peptide_residues")
                if ref:
                    pep_size = int(ref)
            if pep_size <= 0:
                raise NodeException(
                    "setup",
                    "Could not determine peptide residue count. Connect "
                    "ep_fragment_fuse_topology (or any node forwarding "
                    "``peptide_residues``) upstream, or set the option.",
                )

            stream_log(
                f"Mapping coords ({pep_size}-res peptide → fragment frame)...",
                node_id=self.node_id, progress=40,
            )

            full_prefix = os.path.join(output_dir, output_prefix)
            stats = apply_coords(
                prmtop_path=prmtop,
                peptide_pdb_path=peptide_pdb,
                fragment_pdb_path=fragment_pdb,
                output_prefix=full_prefix,
                pep_residues=pep_size,
            )

            stream_log(
                f"Wrote {stats['n_atoms']} atoms "
                f"({stats['pep_atoms']} peptide + {stats['frag_atoms']} fragment).",
                node_id=self.node_id, progress=90,
            )

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                # Forward prmtop unchanged so downstream nodes (e.g.
                # ep_amber_to_gromacs) can keep auto-discovering output_prmtop
                # from this node without needing two predecessors.
                "output_prmtop": self.format_output_path(prmtop),
                "output_rst7": self.format_output_path(stats["rst7"]),
                "output_pdb": self.format_output_path(stats["pdb"]),
                "n_atoms": stats["n_atoms"],
                "peptide_residues": pep_size,
            }
            result.files["output"] = {
                "prmtop": self.format_output_path(prmtop),
                "rst7": self.format_output_path(stats["rst7"]),
                "pdb": self.format_output_path(stats["pdb"]),
            }
            result.success = True
            result.message = (
                f"Applied coords: {stats['n_atoms']} atoms "
                f"({stats['pep_atoms']}+{stats['frag_atoms']}) → "
                f"{output_prefix}.rst7 + {output_prefix}.pdb"
            )
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("apply coords", str(e))
