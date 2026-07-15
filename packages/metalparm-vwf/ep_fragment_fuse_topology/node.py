"""Fragment Fuse (Topology) Node — emits only the merged AMBER topology.

Produces ``complex.prmtop`` from a peptide parameterized topology (built
upstream by ``peptide_builder``) and a fragment topology
(``ep_library_generation``). Interface bonds and cap-atom removals are
specified as a JSON list — generalizes to multi-bond cofactors
(Zn-finger 2C2H, heme axial Cys, etc.).

This is half of the v1.12.0 split of the previous monolithic
``ep_fragment_fuse``: this node owns topology fusion only; the sister
``ep_apply_coords`` node consumes this prmtop plus aligned source PDBs
to write the final ``complex.rst7`` + ``complex.pdb``. See
``packages/metalparm-vwf/CHANGELOG.md`` v1.12.0 for the rationale.

Predecessor data flow (3-tier — same as today's fuse):
  - ``peptide_lib``:    explicit > ``output_peptide_lib`` from peptide_builder > error
  - ``peptide_frcmod``: explicit > ``output_peptide_frcmod`` from peptide_builder > error
  - ``fragment_lib``:   explicit > ``output_fragment_lib`` from ep_library_generation > error
                        (fragment_align is NOT required as a predecessor here —
                         alignment only affects coordinates, which this node
                         doesn't write)
  - ``fragment_frcmod``:explicit > ``output_fragment_frcmod`` from same > error
  - ``forcefield``:     explicit > ``forcefield`` from peptide_builder > "ff19SB"

Outputs forwarded:
  ``output_prmtop``, ``output_tleap_script``, ``peptide_residues``,
  ``interface_bonds``, ``case_name``, ``working_path``.

The auto-generated ``<case>_topology.rst7`` (tleap can't suppress it
during ``saveamberparm``) is kept on disk for diagnostics — its coords
reflect whatever the input libs carried. Callers should consume the
final ``complex.rst7`` written by ``ep_apply_coords``, not this one.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FolderParameter, SelectParameter,
    StringParameter, TextParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .fuse_helpers import (
        DEFAULT_INTERFACE_BONDS, build_tleap_script, count_pdb_residues,
        parse_interface_bonds, read_lib_unit_name, rebalance_residue_charges,
    )
except ImportError:
    try:
        from fuse_helpers import (  # type: ignore
            DEFAULT_INTERFACE_BONDS, build_tleap_script, count_pdb_residues,
            parse_interface_bonds, read_lib_unit_name, rebalance_residue_charges,
        )
    except ImportError:  # server-side introspection (no sys.path tweak yet)
        DEFAULT_INTERFACE_BONDS = []
        build_tleap_script = None
        count_pdb_residues = None
        parse_interface_bonds = None
        read_lib_unit_name = None
        rebalance_residue_charges = None


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


class EpFragmentFuseTopology(Node):
    """tleap-based topology fuser: peptide.lib + fragment.lib → complex.prmtop."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for fuse_topology.tleap + complex.prmtop",
        ),
        "peptide_lib": FileParameterEdit(
            "Peptide Library (.lib/.off)",
            default="",
            docstring=(
                "Parameterized peptide AMBER OFF library. Leave empty to "
                "auto-discover ``output_peptide_lib`` (or ``output_lib``) "
                "from a peptide_builder predecessor."
            ),
            optional=True,
        ),
        "peptide_frcmod": FileParameterEdit(
            "Peptide Frcmod",
            default="",
            docstring=(
                "Peptide-side frcmod (typically the empty placeholder from "
                "peptide_builder for pure ff19SB). Leave empty to auto-discover."
            ),
            optional=True,
        ),
        "fragment_lib": FileParameterEdit(
            "Fragment Library (.lib/.off)",
            default="",
            docstring=(
                "Fragment AMBER OFF library — the raw output from "
                "ep_library_generation. Aligning it via fragment_align is "
                "NOT required here (this node writes only topology, not "
                "coordinates) but is harmless. Leave empty to auto-discover "
                "``output_fragment_lib`` (or ``output_lib``)."
            ),
            optional=True,
        ),
        "fragment_frcmod": FileParameterEdit(
            "Fragment Frcmod",
            default="",
            docstring="AMBER frcmod file. Leave empty to auto-discover from predecessor (output_fragment_frcmod / output_frcmod).",
            optional=True,
        ),
        "fragment_resname": StringParameter(
            "Fragment Unit Name",
            default="",
            docstring=(
                "Name of the unit inside the fragment library (e.g. 'mol' "
                "from antechamber, or 'SNP' from a hand-named lib). Leave "
                "empty to auto-detect from the lib's index header."
            ),
        ),
        "peptide_resname": StringParameter(
            "Peptide Unit Name",
            default="",
            docstring=(
                "Name of the unit inside the peptide library (auto-set by "
                "peptide_builder via tleap's saveoff). Leave empty to "
                "auto-detect from the lib's index header."
            ),
        ),
        "interface_bonds": TextParameter(
            "Interface Bonds (JSON)",
            default=json.dumps(DEFAULT_INTERFACE_BONDS, indent=2),
            docstring=(
                "JSON list of interface bonds. Each entry: "
                "{pep_resid, pep_atom, frag_resid, frag_atom, pep_remove, frag_remove}"
            ),
        ),
        "forcefield": SelectParameter(
            "Peptide Force Field",
            default="ff19SB",
            options=["ff19SB", "ff14SB"],
            docstring=(
                "Used to source the protein leaprc. Auto-inherited from "
                "peptide_builder predecessor when present."
            ),
        ),
        "linkage_frcmod": FileParameterEdit(
            "Linkage Frcmod",
            default="",
            docstring=(
                "Optional cross-FF frcmod that fills parameter gaps at the "
                "interface bond. Most GAFF2-ff19SB amide attachments (e.g., "
                "GLU.CD-NH2 in the SnP demo) need a small set of cross-FF "
                "BOND/ANGLE/DIHE/IMPROPER entries that neither FF defines "
                "alone. Bundled examples live under "
                "ep_fragment_fuse_topology/demo_data/linkages/. Auto-discovers "
                "'output_linkage_frcmod' from a predecessor if set."
            ),
            optional=True,
        ),
        "charge_rebalance": BooleanParameter(
            "Rebalance Interface Charges",
            default=True,
            docstring=(
                "After tleap, redistribute the non-integer net charge each "
                "residue is left with when an interface bond's pep_remove / "
                "frag_remove atoms are deleted, so every residue — and the "
                "complex total — is an exact integer. Deleting e.g. GLU's "
                "OE2 (-0.82 e) from a -1 glutamate leaves the residue at "
                "-0.18 e; this spreads that remainder over the residue. A "
                "residue more than 0.4 e from an integer aborts the node "
                "(genuine parameterisation error, not a deletion remainder)."
            ),
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting fragment fuse (topology)...", node_id=self.node_id, progress=0)

        if build_tleap_script is None:
            raise NodeException("setup",
                "ep_fragment_fuse_topology fuse_helpers.py could not be imported "
                "— run this node in a pixi env with ambertools available.")

        try:
            result = NodeResult()
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}

            case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "complex")
            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)

            forcefield = flow_vars["forcefield"].get_value() or ""
            if not forcefield:
                forcefield = _get_from_predecessors(predecessor_data, "forcefield") or "ff19SB"

            raw_bonds = flow_vars["interface_bonds"].get_value() or ""
            try:
                interface_bonds = parse_interface_bonds(raw_bonds)
            except ValueError as ve:
                raise NodeException("setup", str(ve))

            def _resolve(option_key: str, pred_keys: list[str], work_filename: str,
                         required: bool = True) -> str:
                explicit = self.resolve_path(flow_vars[option_key].get_value()) or ""
                if explicit and os.path.isfile(explicit):
                    return _ensure_in_workdir(self.resolve_path, explicit, out_dir,
                                              os.path.basename(explicit)) or explicit
                for pk in pred_keys:
                    ref = _get_from_predecessors(predecessor_data, pk)
                    if ref:
                        landed = _ensure_in_workdir(self.resolve_path, ref, out_dir,
                                                    work_filename)
                        if landed:
                            return landed
                if required:
                    raise NodeException("setup",
                        f"{option_key} not provided and no "
                        f"{' / '.join(pred_keys)} in predecessor data.")
                return ""

            peptide_lib = _resolve(
                "peptide_lib",
                ["output_peptide_lib", "output_lib"],
                "peptide.lib",
            )
            peptide_frcmod = _resolve(
                "peptide_frcmod",
                ["output_peptide_frcmod", "output_frcmod"],
                "peptide.frcmod",
                required=False,
            )
            fragment_lib = _resolve(
                "fragment_lib",
                ["output_fragment_lib", "output_lib"],
                "fragment.lib",
            )
            fragment_frcmod = _resolve(
                "fragment_frcmod",
                ["output_fragment_frcmod", "output_frcmod"],
                "fragment.frcmod",
            )

            if os.path.abspath(peptide_lib) == os.path.abspath(fragment_lib):
                raise NodeException("setup",
                    "peptide_lib and fragment_lib resolved to the same file. "
                    "Wire a peptide_builder upstream so it forwards "
                    "output_peptide_lib (role-specific key), or set both "
                    "options explicitly.")

            frag_resname = flow_vars["fragment_resname"].get_value().strip()
            if not frag_resname:
                frag_resname = read_lib_unit_name(fragment_lib)
            pep_resname = flow_vars["peptide_resname"].get_value().strip()
            if not pep_resname:
                pep_resname = read_lib_unit_name(peptide_lib)

            pep_size = _get_from_predecessors(predecessor_data, "peptide_residues") or 0
            if pep_size <= 0:
                candidate = os.path.join(os.path.dirname(peptide_lib), "peptide.pdb")
                if os.path.isfile(candidate):
                    pep_size = count_pdb_residues(candidate)
            if pep_size <= 0:
                raise NodeException("setup",
                    "Could not determine peptide residue count. Connect "
                    "a peptide_builder upstream so it forwards peptide_residues.")

            linkage_frcmod = self.resolve_path(
                flow_vars["linkage_frcmod"].get_value() or ""
            )
            if not linkage_frcmod:
                ref = _get_from_predecessors(predecessor_data, "output_linkage_frcmod")
                if ref:
                    linkage_frcmod = self.resolve_path(ref)
            linkage_frcmods: list[str] = []
            if linkage_frcmod and os.path.isfile(linkage_frcmod):
                staged = os.path.join(out_dir, os.path.basename(linkage_frcmod))
                if os.path.abspath(linkage_frcmod) != os.path.abspath(staged):
                    shutil.copy2(linkage_frcmod, staged)
                linkage_frcmods.append(staged)

            tleap_script = build_tleap_script(
                forcefield=forcefield,
                fragment_lib=fragment_lib,
                fragment_frcmod=fragment_frcmod,
                peptide_lib=peptide_lib,
                fragment_resname=frag_resname,
                peptide_resname=pep_resname,
                interface_bonds=interface_bonds,
                pep_unit_size=pep_size,
                linkage_frcmods=linkage_frcmods,
            )
            tleap_path = os.path.join(out_dir, "fuse_topology.tleap")
            with open(tleap_path, "w") as f:
                f.write(tleap_script)

            for stale in (os.path.join(out_dir, "complex.prmtop"),
                          os.path.join(out_dir, "complex.rst7"),
                          os.path.join(out_dir, "complex.pdb"),
                          os.path.join(out_dir, "leap.log")):
                if os.path.isfile(stale):
                    os.unlink(stale)

            stream_log(
                f"Running tleap ({pep_size}-res peptide:{pep_resname} + {frag_resname})...",
                node_id=self.node_id, progress=40,
            )
            proc = subprocess.run(
                ["tleap", "-f", "fuse_topology.tleap"],
                cwd=out_dir, capture_output=True, text=True,
            )
            leap_log = os.path.join(out_dir, "leap.log")
            if os.path.isfile(leap_log):
                shutil.copy2(leap_log, os.path.join(out_dir, "fuse_topology.leap.log"))

            prmtop = os.path.join(out_dir, "complex.prmtop")
            if not (os.path.isfile(prmtop) and os.path.getsize(prmtop) > 0):
                msg = proc.stderr or proc.stdout or "tleap failed"
                raise NodeException("execution",
                    f"tleap did not produce complex.prmtop. stderr tail:\n"
                    f"{msg[-2000:]}")

            # Repair the non-integer net charge that interface-atom
            # deletion leaves on the linkage residue(s) — tleap does not.
            # See fuse_helpers.rebalance_residue_charges.
            charge_summary = None
            if bool(flow_vars["charge_rebalance"].get_value()):
                try:
                    charge_summary = rebalance_residue_charges(prmtop)
                except ValueError as ex:
                    raise NodeException("execution",
                        f"interface charge rebalance refused — {ex}")
                if charge_summary["adjusted"]:
                    detail = "; ".join(
                        f"{rn}{ri} {b:+.4f}->{a:+.0f}"
                        for rn, ri, b, a in charge_summary["adjusted"])
                    stream_log(
                        f"Charge rebalance: complex total "
                        f"{charge_summary['total_before']:+.5f} -> "
                        f"{charge_summary['total_after']:+.5f} e "
                        f"(adjusted {detail})",
                        node_id=self.node_id, progress=70,
                    )

            # tleap saveamberparm always writes both prmtop+rst7. We can't
            # suppress the rst7, so rename it to make clear the canonical
            # rst7 comes from ep_apply_coords.
            tleap_rst7 = os.path.join(out_dir, "complex.rst7")
            initial_rst7 = os.path.join(out_dir, f"{case_name}_initial.rst7")
            if os.path.isfile(tleap_rst7):
                shutil.move(tleap_rst7, initial_rst7)
            tleap_pdb = os.path.join(out_dir, "complex.pdb")
            initial_pdb = os.path.join(out_dir, f"{case_name}_initial.pdb")
            if os.path.isfile(tleap_pdb):
                shutil.move(tleap_pdb, initial_pdb)

            stream_log("Topology fuse complete.", node_id=self.node_id, progress=100)
            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(out_dir),
                "output_prmtop": self.format_output_path(prmtop),
                "output_tleap_script": self.format_output_path(tleap_path),
                "peptide_residues": pep_size,
                "interface_bonds": interface_bonds,
                "charge_total": (
                    round(charge_summary["total_after"], 6)
                    if charge_summary else None
                ),
            }
            result.files["output"] = {
                "prmtop": self.format_output_path(prmtop),
                "tleap": self.format_output_path(tleap_path),
            }
            result.success = True
            result.message = (
                f"Topology: {pep_size}-res peptide ({pep_resname}) + {frag_resname} "
                f"via {len(interface_bonds)} interface bond(s)"
            )
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("fragment fuse topology", str(e))
