"""MD Membrane Build — embed a dry metallopeptide GROMACS topology in a
DPPC bilayer (transmembrane), solvate + add ions, and emit a solvated
GROMACS topology for membrane MD.

    ep_amber_to_gromacs → md_membrane_build → (gmx grompp / mdrun)

This is the Case 2 counterpart of md_solvate_gmx (which builds a
mixed-solvent box for Case 1). It follows the same principle — the
solute topology from ep_amber_to_gromacs is the source of truth and is
never re-derived by tleap — see core.py for why that matters.

The peptide is oriented transmembrane by aligning its helix axis to
z geometrically (MEMEMBED, packmol-memgen's default orienter, cannot
handle a metallopeptide — the non-standard fragment residue makes it
fail). The case's long peptide is a YALP-family sequence — a
canonical transmembrane model helix — so this is the intended
embedding. Set ``preoriented`` to skip the auto-orientation if the
input is already oriented (helix axis along z).

Output keys: ``output_top``, ``output_gro``, ``output_itp``,
``case_name``, ``working_path``, plus lipid/water/ion counts.
"""
from __future__ import annotations

import os
import shutil

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FloatParameter, FolderParameter,
    IntegerParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import (
        assemble_membrane_system, compute_solute_charge,
        count_membrane_residues, load_solute_structure,
        orient_peptide_along_z, parametrize_membrane, read_packmol_tolerance,
        run_packmol_memgen, save_gromacs_outputs, split_packed_pdb,
        write_solute_pdb,
    )
except ImportError:  # script-mode / server introspection fallback
    try:
        from core import (  # type: ignore
            assemble_membrane_system, compute_solute_charge,
            count_membrane_residues, load_solute_structure,
            orient_peptide_along_z, parametrize_membrane,
            read_packmol_tolerance, run_packmol_memgen, save_gromacs_outputs,
            split_packed_pdb, write_solute_pdb,
        )
    except ImportError:
        assemble_membrane_system = compute_solute_charge = None
        count_membrane_residues = load_solute_structure = None
        orient_peptide_along_z = read_packmol_tolerance = None
        parametrize_membrane = run_packmol_memgen = None
        save_gromacs_outputs = split_packed_pdb = write_solute_pdb = None


def _from_predecessors(predecessor_data, *keys):
    for pred in (predecessor_data or []):
        if not pred:
            continue
        for k in keys:
            if k in pred and pred[k]:
                return pred[k]
        of = pred.get("output_files") or {}
        for k in keys:
            if k in of and of[k]:
                return of[k]
    return None


class MdMembraneBuild(Node):
    """Embed a dry metallopeptide GROMACS topology in a DPPC bilayer."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default=""),
        "input_top": FileParameterEdit(
            "Topology File (.top)",
            default="",
            docstring=(
                "Dry metallopeptide GROMACS topology from "
                "ep_amber_to_gromacs. Leave empty to auto-discover "
                "'output_top' from a predecessor."
            ),
            optional=True,
        ),
        "input_gro": FileParameterEdit(
            "Structure File (.gro)",
            default="",
            docstring="Dry GROMACS coordinates. Auto-discovers 'output_gro'.",
            optional=True,
        ),
        "input_itp": FileParameterEdit(
            "ITP File (.itp, optional)",
            default="",
            docstring=(
                "Moleculetype .itp the .top #includes, if split. "
                "Auto-discovers 'output_itp'; staged beside the .top so "
                "the include resolves."
            ),
            optional=True,
        ),
        "lipid": StringParameter(
            "Lipid",
            default="DPPC",
            docstring=(
                "Bilayer lipid. Only DPPC is currently supported: the "
                "membrane residue classifier and tleap parametrisation are "
                "keyed to Lipid21 DPPC's split-residue names (PC/PA). Other "
                "lipids (POPC, etc.) route their tail residues into the solute "
                "and fail — the node rejects a non-DPPC value up front."
            ),
        ),
        "box_dist": FloatParameter(
            "Box Distance (Å)",
            default=17.5,
            docstring="packmol-memgen --dist: min solute-to-box padding.",
        ),
        "water_dist": FloatParameter(
            "Water Layer (Å)",
            default=17.5,
            docstring="packmol-memgen --dist_wat: water width each side.",
        ),
        "saltcon": FloatParameter(
            "Salt Concentration (M)",
            default=0.15,
            docstring=(
                "NaCl concentration; 0 = counter-ions only. Default 0.15 M "
                "(physiological)."
            ),
        ),
        "xy_box_A": FloatParameter(
            "XY Box Size (Å, override)",
            default=0.0,
            docstring=(
                "packmol-memgen --distxy_fix: force the XY box to this fixed "
                "size in Å. 0 (default) leaves packmol-memgen to auto-size "
                "XY from the solute extent + Box Distance. Bump this above "
                "the auto-size when the bilayer packing fails to converge "
                "(symptom: packmol's all-together packing loop runs without "
                "lowering its function value). Past Case-2 long-peptide runs "
                "auto-sized to ~80 Å; setting 80 here matches that geometry."
            ),
        ),
        "nloop_all": IntegerParameter(
            "Packmol All-Together Iterations",
            default=0,
            docstring=(
                "packmol-memgen --nloop_all: number of all-together packing "
                "iterations. 0 (default) leaves packmol-memgen's built-in "
                "default. Bump (e.g. 100–200) when the bilayer pack converges "
                "slowly. Independent of --nloop (per-round GENCAN iterations)."
            ),
        ),
        "preoriented": BooleanParameter(
            "Pre-oriented",
            default=False,
            docstring=(
                "If the solute is already oriented with the peptide "
                "helix axis along z, skip the geometric pre-orientation. "
                "Default off — the node aligns the helix axis to z "
                "itself (MEMEMBED cannot orient a metallopeptide)."
            ),
        ),
        "output_prefix": StringParameter(
            "Output Prefix",
            default="complex",
            docstring="Master file prefix → <prefix>.top / <prefix>.gro.",
        ),
        "itp_filename": StringParameter(
            "ITP Filename",
            default="metallopeptide_mem",
            docstring=(
                "Basename for the split moleculetype .itp. Empty = keep "
                "the .top monolithic."
            ),
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory for the solvated membrane topology.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting DPPC membrane build...", node_id=self.node_id,
                   progress=0)
        if run_packmol_memgen is None:
            raise NodeException("setup",
                "core.py could not be imported — run this node in the "
                "metalparm_vwf pixi env (ambertools + parmed).")
        try:
            result = NodeResult()
            input_data = (predecessor_data[0]
                          if predecessor_data and predecessor_data[0] else {})
            case_name = (flow_vars["case_name"].get_value()
                         or input_data.get("case_name", "case"))

            top = flow_vars["input_top"].get_value() or \
                _from_predecessors(predecessor_data, "output_top", "top")
            gro = flow_vars["input_gro"].get_value() or \
                _from_predecessors(predecessor_data, "output_gro", "gro")
            itp = flow_vars["input_itp"].get_value() or \
                _from_predecessors(predecessor_data, "output_itp", "itp")
            if not top or not gro:
                raise NodeException("setup",
                    "Need a GROMACS .top + .gro — set them or connect an "
                    "ep_amber_to_gromacs predecessor.")
            top = self.resolve_path(top)
            gro = self.resolve_path(gro)
            itp = self.resolve_path(itp) if itp else None
            for label, p in (("topology", top), ("structure", gro)):
                if not p or not os.path.isfile(p):
                    raise NodeException("setup", f"{label} not found: {p}")

            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)

            # Stage top (+ itp) into the work dir so the .top's relative
            # #include resolves when ParmEd loads it.
            local_top = os.path.join(out_dir, os.path.basename(top))
            if os.path.abspath(local_top) != os.path.abspath(top):
                shutil.copy2(top, local_top)
            if itp and os.path.isfile(itp):
                local_itp = os.path.join(out_dir, os.path.basename(itp))
                if os.path.abspath(local_itp) != os.path.abspath(itp):
                    shutil.copy2(itp, local_itp)

            lipid = flow_vars["lipid"].get_value() or "DPPC"
            # Only DPPC is supported: split_packed_pdb / count_membrane_residues
            # / parametrize_membrane are hardcoded to Lipid21 DPPC's PC/PA
            # split-residue names. A non-DPPC lipid would silently misclassify
            # its tail residues as solute and crash on the atom-count guard, so
            # reject it here with an actionable message instead.
            if lipid.strip().upper() != "DPPC":
                raise NodeException(
                    "setup",
                    f"Unsupported lipid '{lipid}'. md_membrane_build currently "
                    "only supports DPPC (its residue classification and tleap "
                    "parametrisation are DPPC-specific). Set lipid=DPPC, or "
                    "build the bilayer externally and use md_solvate_gmx.",
                )
            box_dist = float(flow_vars["box_dist"].get_value() or 17.5)
            water_dist = float(flow_vars["water_dist"].get_value() or 17.5)
            saltcon = float(flow_vars["saltcon"].get_value() or 0.0)
            xy_box_A = float(flow_vars["xy_box_A"].get_value() or 0.0)
            nloop_all = int(flow_vars["nloop_all"].get_value() or 0)
            preoriented = bool(flow_vars["preoriented"].get_value())
            prefix = flow_vars["output_prefix"].get_value() or "complex"
            itp_name = flow_vars["itp_filename"].get_value() or None

            stream_log("Loading dry solute topology...",
                       node_id=self.node_id, progress=10)
            solute = load_solute_structure(local_top, gro)
            q = compute_solute_charge(solute)
            if not preoriented:
                stream_log("Orienting peptide helix axis along z...",
                           node_id=self.node_id, progress=18)
                orient_peptide_along_z(solute)
            solute_pdb = write_solute_pdb(
                solute, os.path.join(out_dir, "solute.pdb"))

            packmol_dir = os.path.join(out_dir, "packmol")
            tuning_note = []
            if xy_box_A > 0:
                tuning_note.append(f"--distxy_fix {xy_box_A:g} Å")
            if nloop_all > 0:
                tuning_note.append(f"--nloop_all {nloop_all}")
            stream_log(
                f"Running packmol-memgen ({lipid} bilayer, saltcon {saltcon} M"
                + ("; " + ", ".join(tuning_note) if tuning_note else "")
                + ")...",
                node_id=self.node_id, progress=25)
            packed = run_packmol_memgen(
                solute_pdb, packmol_dir, lipid=lipid, dist=box_dist,
                dist_wat=water_dist, saltcon=saltcon,
                xy_box_A=xy_box_A or None, nloop_all=nloop_all or None)

            stream_log("Splitting packed system (solute / membrane)...",
                       node_id=self.node_id, progress=55)
            packed_solute, packed_membrane, counts = split_packed_pdb(
                packed, packmol_dir)

            stream_log(f"Parametrising membrane "
                       f"({counts['water']} water atoms, "
                       f"{counts['ion']} ion atoms)...",
                       node_id=self.node_id, progress=70)
            mem_prmtop, mem_rst7 = parametrize_membrane(
                packed_membrane, packmol_dir)

            stream_log("Assembling solvated membrane system...",
                       node_id=self.node_id, progress=85)
            # packmol packs non-periodically; the box is derived from the
            # actual packed-coordinate extent + the packmol tolerance so
            # periodic images do not clash (a sub-Å face clash makes the
            # first EM step diverge) — see assemble_membrane_system.
            margin = read_packmol_tolerance(
                os.path.join(packmol_dir, "packmol.inp"))
            system = assemble_membrane_system(
                solute, packed_solute, mem_prmtop, mem_rst7, margin=margin)

            outs = save_gromacs_outputs(
                system, out_dir, prefix=prefix, itp_filename=itp_name)

            # counts from split_packed_pdb are atom counts; report
            # molecule counts from the assembled topology instead.
            res_counts = count_membrane_residues(system)
            total_q = float(sum(a.charge for a in system.atoms))
            stream_log(
                f"Membrane build complete: {res_counts['lipid']} DPPC "
                f"lipids, {res_counts['water']} water, {res_counts['ion']} "
                f"ions; system charge {total_q:+.4f} e",
                node_id=self.node_id, progress=100)

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(out_dir),
                "output_top": self.format_output_path(outs["top"]),
                "output_gro": self.format_output_path(outs["gro"]),
                "output_itp": (self.format_output_path(outs["itp"])
                               if outs["itp"] else None),
                "lipid": lipid,
                "n_lipid": res_counts["lipid"],
                "n_water": res_counts["water"],
                "n_ion": res_counts["ion"],
                "solute_charge": round(q, 5),
                "system_charge": round(total_q, 5),
            }
            result.files["output"] = {
                "top": self.format_output_path(outs["top"]),
                "gro": self.format_output_path(outs["gro"]),
            }
            if outs["itp"]:
                result.files["output"]["itp"] = self.format_output_path(
                    outs["itp"])
            result.success = True
            result.message = (
                f"{lipid} bilayer: {res_counts['lipid']} lipids, "
                f"{res_counts['water']} water, {res_counts['ion']} ions; "
                f"system charge {total_q:+.3f} e")
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("md membrane build", str(e))
