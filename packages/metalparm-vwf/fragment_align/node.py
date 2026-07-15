"""Fragment Align Node — pure-geometry rigid-body placement before fuse.

What it does (Avogadro / GaussView style — no FF, no QM):

  Given a peptide PDB (from ``peptide_builder``) and a fragment PDB +
  fragment OFF library (from ``snp_builder`` and ``ep_library_generation``),
  compute the rigid-body transformation that places the fragment so its
  bonding atom sits at the standard amide-bond distance and angle from
  the peptide's bonding atom. Apply the transformation via ``tleap``
  (translate + transform + saveoff) and emit an aligned lib + PDB.

Downstream ``ep_fragment_fuse_topology`` then loads the aligned lib via ``loadoff``
and combines it with the peptide topology — the interface bond ends up
at correct geometry without any minimization.

Predecessor data flow (3-tier resolution):
  - ``peptide_pdb``:    explicit > ``output_pdb`` from peptide_builder predecessor
  - ``fragment_pdb``:   explicit > ``output_pdb`` from snp_builder OR
                        ep_library_generation predecessor
  - ``fragment_lib``:   explicit > ``output_lib`` from ep_library_generation
  - ``fragment_frcmod``: explicit > ``output_frcmod`` from ep_library_generation
                        (only loaded so tleap recognises atom types)

The pure-geometry math lives in ``core.py`` so it's unit-testable without
bocoflow_core or a tleap binary.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FolderParameter, StringParameter,
    TextParameter,
)
from bocoflow_core.stream_logger import stream_log

# fuse_helpers ships DEFAULT_INTERFACE_BONDS + parse_interface_bonds; reuse
# them so the schema stays a single source of truth across align + fuse.
try:
    from ep_fragment_fuse_topology.fuse_helpers import (  # type: ignore
        DEFAULT_INTERFACE_BONDS, parse_interface_bonds,
    )
except ImportError:
    # Fallback: when the node is installed standalone the sibling import
    # path may differ. Re-define the default + parser locally as a last resort.
    DEFAULT_INTERFACE_BONDS = [
        {
            "pep_resid": 6, "pep_atom": "CD",
            "frag_resid": 1, "frag_atom": "NH2",
            "pep_remove": ["OE2", "HE2"],
            "frag_remove": ["CM", "HM1", "HM2", "HM3", "CAP", "OAP"],
        }
    ]

    def parse_interface_bonds(raw: str):  # type: ignore[no-redef]
        if not raw or not raw.strip():
            return DEFAULT_INTERFACE_BONDS
        v = json.loads(raw)
        if not isinstance(v, list) or not v:
            raise ValueError("interface_bonds must be a non-empty list")
        return v

try:
    from .core import (
        STANDARD_BOND_LENGTHS, build_align_tleap_script,
        compute_rigid_transformation,
    )
except ImportError:
    try:
        from core import (  # type: ignore
            STANDARD_BOND_LENGTHS, build_align_tleap_script,
            compute_rigid_transformation,
        )
    except ImportError:  # server-side introspection path (no heavy deps)
        STANDARD_BOND_LENGTHS = None
        build_align_tleap_script = None
        compute_rigid_transformation = None


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


class FragmentAlign(Node):
    """Pre-fuse rigid-body placement of the fragment onto the peptide anchor."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for the aligned fragment lib + PDB",
        ),
        "peptide_pdb": FileParameterEdit(
            "Peptide PDB",
            default="",
            docstring=(
                "Peptide structure used for anchor-direction computation. "
                "Leave empty to auto-discover ``output_pdb`` from a "
                "peptide_builder predecessor."
            ),
            optional=True,
        ),
        "fragment_pdb": FileParameterEdit(
            "Fragment PDB",
            default="",
            docstring=(
                "Fragment structure with atom names matching interface_bonds "
                "(e.g. snp_builder's snp_frag.pdb with NH2/CAP/OAP, OR the "
                "antechamber-renamed PDB from ep_library_generation). Leave "
                "empty to auto-discover ``output_pdb`` from a snp_builder OR "
                "ep_library_generation predecessor."
            ),
            optional=True,
        ),
        "fragment_lib": FileParameterEdit(
            "Fragment Library (.lib/.off)",
            default="",
            docstring="AMBER OFF/lib file. Leave empty to auto-discover ``output_lib``.",
            optional=True,
        ),
        "fragment_frcmod": FileParameterEdit(
            "Fragment Frcmod",
            default="",
            docstring=(
                "AMBER frcmod file (loaded so tleap recognises fragment atom "
                "types during transform/saveoff). Leave empty to auto-discover."
            ),
            optional=True,
        ),
        "fragment_resname": StringParameter(
            "Fragment Unit Name",
            default="mol",
            docstring=(
                "Name of the unit inside the fragment lib (the leading "
                "``!!index array str \"NAME\"`` line). antechamber-derived libs "
                "are typically named ``mol``."
            ),
        ),
        "interface_bonds": TextParameter(
            "Interface Bonds (JSON)",
            default=json.dumps(DEFAULT_INTERFACE_BONDS, indent=2),
            docstring=(
                "Same JSON schema as ep_fragment_fuse_topology. Optional extra keys "
                "per bond: ``bond_kind`` (key into STANDARD_BOND_LENGTHS), "
                "``bond_length`` (explicit Å override), ``pep_hybrid`` / "
                "``frag_hybrid`` (override hybridization)."
            ),
        ),
        "bond_lengths": TextParameter(
            "Bond Length Overrides (JSON)",
            default="{}",
            docstring=(
                "Optional per-bond-kind override for STANDARD_BOND_LENGTHS, "
                "e.g. ``{\"C-N_amide\": 1.32}``."
            ),
        ),
        "clash_optimize": BooleanParameter(
            "Clash-free Rotation Scan",
            default=True,
            docstring=(
                "Scan 12 candidate rotations around the new bond axis and "
                "pick the one that maximizes peptide↔fragment atom distance."
            ),
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting fragment align...", node_id=self.node_id, progress=0)

        if compute_rigid_transformation is None:
            raise NodeException("setup",
                "fragment_align core.py could not be imported — run this node "
                "in a pixi env with numpy + biopython available.")

        try:
            result = NodeResult()
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}

            case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "complex")
            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)

            frag_resname = flow_vars["fragment_resname"].get_value() or "mol"
            clash_opt = bool(flow_vars["clash_optimize"].get_value())

            try:
                interface_bonds = parse_interface_bonds(
                    flow_vars["interface_bonds"].get_value() or ""
                )
            except ValueError as ve:
                raise NodeException("setup", str(ve))

            # bond_lengths overrides
            raw_overrides = flow_vars["bond_lengths"].get_value() or "{}"
            try:
                overrides_dict = json.loads(raw_overrides) if raw_overrides.strip() else {}
            except json.JSONDecodeError as e:
                raise NodeException("setup", f"bond_lengths is not valid JSON: {e}")

            # --- 3-tier resolution for inputs ---
            def _resolve_input(option_key: str, predecessor_key: str,
                               local_filename: str, required: bool = True) -> str:
                explicit = self.resolve_path(flow_vars[option_key].get_value()) or ""
                if explicit and os.path.isfile(explicit):
                    return _ensure_in_workdir(self.resolve_path, explicit, out_dir,
                                              os.path.basename(explicit)) or explicit
                ref = _get_from_predecessors(predecessor_data, predecessor_key)
                if ref:
                    landed = _ensure_in_workdir(self.resolve_path, ref, out_dir,
                                                local_filename)
                    if landed:
                        return landed
                if required:
                    raise NodeException("setup",
                        f"{option_key} not provided and no {predecessor_key} "
                        f"in predecessor data.")
                return ""

            peptide_pdb = _resolve_input("peptide_pdb", "output_pdb", "peptide.pdb")
            fragment_pdb = _resolve_input("fragment_pdb", "output_pdb", "fragment.pdb")
            # fragment_pdb may have been resolved to peptide.pdb if both
            # predecessors emit "output_pdb" — disambiguate by checking if
            # the user explicitly provided an upstream-named output
            fragment_lib = _resolve_input("fragment_lib", "output_lib", "fragment.lib")
            fragment_frcmod = _resolve_input(
                "fragment_frcmod", "output_frcmod", "fragment.frcmod", required=False,
            )

            stream_log(
                f"Computing rigid transformation for {len(interface_bonds)} interface bond(s)...",
                node_id=self.node_id, progress=30,
            )

            # Compute transformations for each interface bond. For multi-bond
            # cases (e.g. Zn-finger 2C2H), we currently use the first bond's
            # transformation only; the others are honored for the bond
            # creation in fuse but not for placement.
            ifb = interface_bonds[0]
            try:
                translation, rotation = compute_rigid_transformation(
                    peptide_pdb_path=peptide_pdb,
                    fragment_pdb_path=fragment_pdb,
                    interface_bond=ifb,
                    clash_optimize=clash_opt,
                    bond_length_overrides=overrides_dict,
                )
            except (ValueError, FileNotFoundError) as ex:
                raise NodeException("setup", f"compute_rigid_transformation failed: {ex}")

            # Drive tleap to apply the transformation to the lib
            aligned_lib = os.path.join(out_dir, f"{case_name}_aligned.lib")
            aligned_pdb = os.path.join(out_dir, f"{case_name}_aligned.pdb")
            tleap_script = build_align_tleap_script(
                fragment_lib_basename=os.path.basename(fragment_lib),
                fragment_resname=frag_resname,
                translation=translation,
                rotation=rotation,
                output_lib_basename=os.path.basename(aligned_lib),
                output_pdb_basename=os.path.basename(aligned_pdb),
                fragment_frcmod_basename=os.path.basename(fragment_frcmod) if fragment_frcmod else None,
            )

            tleap_path = os.path.join(out_dir, "align.tleap")
            with open(tleap_path, "w") as f:
                f.write(tleap_script)

            # tleap saveoff APPENDS to existing libs — wipe stale outputs.
            for stale in (aligned_lib, aligned_pdb,
                          os.path.join(out_dir, "leap.log")):
                if os.path.isfile(stale):
                    os.unlink(stale)

            stream_log("Running tleap to apply transformation...",
                       node_id=self.node_id, progress=70)
            proc = subprocess.run(
                ["tleap", "-f", "align.tleap"],
                cwd=out_dir, capture_output=True, text=True,
            )
            leap_log = os.path.join(out_dir, "leap.log")
            if os.path.isfile(leap_log):
                shutil.copy2(leap_log, os.path.join(out_dir, "align.leap.log"))

            if not (os.path.isfile(aligned_lib) and os.path.getsize(aligned_lib) > 0):
                msg = proc.stderr or proc.stdout or "tleap failed"
                raise NodeException("execution",
                    f"tleap did not produce {os.path.basename(aligned_lib)}. "
                    f"stderr tail:\n{msg[-2000:]}")

            stream_log(
                f"Aligned fragment ready: |t|={float(__import__('numpy').linalg.norm(translation)):.3f} Å",
                node_id=self.node_id, progress=100,
            )

            output = {
                "case_name": case_name,
                "working_path": self.format_output_path(out_dir),
                # Role-specific keys for the fragment side (these are what
                # ep_fragment_fuse_topology looks up for fragment_lib/fragment_frcmod
                # — avoids collision with output_lib from a peptide_builder
                # passthrough below).
                "output_fragment_lib": self.format_output_path(aligned_lib),
                "output_fragment_pdb": self.format_output_path(aligned_pdb),
                # Generic keys (kept for back-compat with consumers that look
                # up output_lib without a role hint).
                "output_lib": self.format_output_path(aligned_lib),
                "output_pdb": self.format_output_path(aligned_pdb),
                "interface_bonds": interface_bonds,
                "fragment_resname": frag_resname,
            }
            # Fragment frcmod (own role-specific + generic passthrough)
            if fragment_frcmod and os.path.isfile(fragment_frcmod):
                output["output_fragment_frcmod"] = self.format_output_path(fragment_frcmod)
                output["output_frcmod"] = self.format_output_path(fragment_frcmod)
            # Pass through peptide-side keys from the upstream peptide_builder
            # predecessor so ep_fragment_fuse_topology / ep_fragment_fuse_topology can
            # find peptide.lib + .frcmod, and ep_apply_coords can find the
            # peptide PDB, all via this single fragment_align predecessor.
            for k in ("output_peptide_pdb", "output_peptide_lib",
                      "output_peptide_frcmod", "peptide_residues",
                      "forcefield"):
                v = _get_from_predecessors(predecessor_data, k)
                if v is not None:
                    output[k] = v

            result.data = output
            result.files["output"] = {
                "lib": output["output_lib"],
                "pdb": output["output_pdb"],
                "tleap": self.format_output_path(tleap_path),
            }
            if "output_frcmod" in output:
                result.files["output"]["frcmod"] = output["output_frcmod"]
            result.success = True
            result.message = (
                f"Aligned fragment to peptide ({len(interface_bonds)} interface bond(s))"
            )
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("fragment align", str(e))
