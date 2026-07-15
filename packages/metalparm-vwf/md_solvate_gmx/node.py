"""MD Solvate (GROMACS-side) — wraps raw packmol + ParmEd to solvate a
dry metallopeptide GROMACS topology in a (mixed) solvent box, without
ever round-tripping the solute through tleap re-derivation.

Designed to sit between ``ep_amber_to_gromacs`` and a downstream MD
step:

    ep_amber_to_gromacs → md_solvate_gmx → (gmx grompp / mdrun)

Why this exists alongside ``md_solvate_packmol``:
  packmol-memgen invokes tleap internally to rebuild the AMBER topology
  around the packed solvent. tleap re-types every residue against the
  standard AMBER library; for fragment-fused residues (e.g. a GLU whose
  side chain has lost OE2 and gained a bond to a SnP fragment N),
  tleap auto-completes the "missing" OE2 atom at standard GLU geometry,
  placing it on top of the fragment N. The result is a topology with
  two atoms in the same spot and no covalent bond between them — LJ-SR
  explodes at step 0 of MD.

  md_solvate_gmx avoids this by:
    1. consuming the *already-converted* GROMACS topology from
       ep_amber_to_gromacs — the solute topology is the source of truth;
    2. using raw packmol (not packmol-memgen) for coordinate placement only;
    3. building solvent moleculetypes via tleap on the AmberTools
       libraries (solvents.lib, atomic_ions.lib, frcmod.opc,
       frcmod.ionslm_126_opc), then ParmEd-converting each to GROMACS;
    4. assembling the final Structure entirely in ParmEd's algebra
       (`+`/`*`), then exporting as GROMACS .top + .gro + .itp.

Forwarded data (matches ep_amber_to_gromacs's output keys):
  ``output_top``, ``output_gro``, ``output_itp``, ``case_name``,
  ``working_path``, plus per-solvent counts and box dims.
"""
from __future__ import annotations

import os
from pathlib import Path

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FloatParameter, FolderParameter, IntegerParameter,
    SelectParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import (
        WATER_MODEL_LEAPRC, assemble_solvated_structure,
        build_solvent_unit_structures, compute_box_dimensions_nm,
        compute_ion_counts, compute_solute_charge,
        compute_solvent_counts, estimate_solute_volume_nm3,
        load_solute_structure, parse_solvent_ratio, run_packmol,
        save_gromacs_outputs, write_packmol_input,
        write_single_molecule_pdbs, write_solute_pdb,
    )
    from .helpers import detect_amberhome, get_from_predecessors
except ImportError:  # script-mode fallback
    from core import (  # type: ignore
        WATER_MODEL_LEAPRC, assemble_solvated_structure,
        build_solvent_unit_structures, compute_box_dimensions_nm,
        compute_ion_counts, compute_solute_charge,
        compute_solvent_counts, estimate_solute_volume_nm3,
        load_solute_structure, parse_solvent_ratio, run_packmol,
        save_gromacs_outputs, write_packmol_input,
        write_single_molecule_pdbs, write_solute_pdb,
    )
    from helpers import detect_amberhome, get_from_predecessors  # type: ignore


_VALID_WATER_MODELS = sorted(WATER_MODEL_LEAPRC.keys())  # opc, opc3, tip3p, ...
_VALID_CATIONS = ["K+", "Na+"]
_VALID_ANIONS = ["Cl-"]


class MdSolvateGmx(Node):
    """Solvate a dry GROMACS topology with mixed solvent via raw packmol + ParmEd."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring=(
                "Working dir; node writes complex.{top,gro} + "
                "metallopeptide_solv.itp + packmol/tleap scratch files here."
            ),
        ),

        # ── Inputs (auto-discover from ep_amber_to_gromacs) ────────────
        "input_top": FileParameterEdit(
            "Topology File (.top)",
            default="",
            docstring=(
                "Dry GROMACS topology of the solute. Leave empty to auto-"
                "discover ``output_top`` from an ep_amber_to_gromacs "
                "predecessor."
            ),
            optional=True,
        ),
        "input_gro": FileParameterEdit(
            "Structure File (.gro)",
            default="",
            docstring=(
                "Dry GROMACS coordinates of the solute. Leave empty to "
                "auto-discover ``output_gro``."
            ),
            optional=True,
        ),
        "input_itp": FileParameterEdit(
            "ITP File (.itp, optional)",
            default="",
            docstring=(
                "Optional separate ``metallopeptide.itp``. If provided, "
                "the node ensures it's co-located with the input .top so "
                "ParmEd's ``#include`` directive resolves cleanly. Leave "
                "blank if the input .top is monolithic."
            ),
            optional=True,
        ),

        # ── Solvent + box ──────────────────────────────────────────────
        "solvents": StringParameter(
            "Solvents",
            default="MOH:WAT",
            docstring=(
                "Solvent codes separated by ':'. Currently supported: "
                "'WAT' (water only) or 'MOH:WAT' (methanol + water). "
                "Order matters for ratio matching."
            ),
        ),
        "solvent_ratio": StringParameter(
            "Solvent Ratio (molar)",
            default="2:1",
            docstring=(
                "Numerical ratios separated by ':'. Interpreted as MOLAR "
                "ratio (unlike packmol-memgen's basis-unspecified ratio). "
                "Counts are solved from bulk densities (MeOH 14.8 mol/nm³, "
                "H₂O 33.4 mol/nm³). For single-solvent: '1'."
            ),
        ),
        "padding_A": FloatParameter(
            "Box Padding (Å)",
            default=12.0,
            docstring=(
                "Half-width margin between the solute extent and each "
                "box face. 12 Å is the AMBER convention."
            ),
        ),

        # ── FF pairings ────────────────────────────────────────────────
        "water_model": SelectParameter(
            "Water Model",
            default="opc",
            options=_VALID_WATER_MODELS,
            docstring=(
                "Water model. OPC is recommended for ff19SB (the "
                "force field's design pair). Drives ``leaprc.water.<model>``."
            ),
        ),

        # ── Ions ───────────────────────────────────────────────────────
        "cation": SelectParameter(
            "Cation", default="K+", options=_VALID_CATIONS,
            docstring="Cation species (must exist in atomic_ions.lib).",
        ),
        "anion": SelectParameter(
            "Anion", default="Cl-", options=_VALID_ANIONS,
            docstring="Anion species (must exist in atomic_ions.lib).",
        ),
        "saltcon_M": FloatParameter(
            "Salt Concentration (M)",
            default=0.0,
            docstring=(
                "Salt concentration above neutralisation (mol/L). 0.0 = "
                "neutralise solute charge only. Physiological mimic: 0.15."
            ),
        ),

        # ── Advanced ───────────────────────────────────────────────────
        "random_seed": IntegerParameter(
            "Random Seed",
            default=-1,
            docstring=(
                "packmol RNG seed; -1 = packmol picks its own random "
                "seed each run (non-reproducible). Set a positive int "
                "for reproducible packings."
            ),
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log(
            "Starting md_solvate_gmx (raw packmol + ParmEd)...",
            node_id=self.node_id, progress=0,
        )

        try:
            result = NodeResult()

            case_name = flow_vars["case_name"].get_value() or "complex"
            output_dir = self.resolve_path(
                flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            # ── Resolve inputs (explicit OR predecessor → abs path) ────
            def resolve_input(explicit_ref, *predecessor_keys):
                ref = explicit_ref
                if not ref:
                    for key in predecessor_keys:
                        ref = get_from_predecessors(predecessor_data, key)
                        if ref:
                            break
                if not ref:
                    return ""
                resolved = self.resolve_path(ref) or ""
                return resolved if (resolved and os.path.isfile(resolved)) else ""

            top_path = resolve_input(
                flow_vars["input_top"].get_value(), "output_top")
            gro_path = resolve_input(
                flow_vars["input_gro"].get_value(), "output_gro")
            itp_path = resolve_input(
                flow_vars["input_itp"].get_value(), "output_itp")

            if not top_path:
                raise NodeException(
                    "setup",
                    "input_top not provided and no output_top in predecessor data.",
                )
            if not gro_path:
                raise NodeException(
                    "setup",
                    "input_gro not provided and no output_gro in predecessor data.",
                )

            # If a separate .itp was supplied (the ITP-split convention
            # from ep_amber_to_gromacs v1.11.0+), ParmEd's GROMACS loader
            # resolves its #include directive relative to the .top's
            # directory. Stage the .itp next to the .top in our output
            # dir so the load works regardless of original layout.
            staged_top = Path(output_dir) / "input_complex.top"
            staged_gro = Path(output_dir) / "input_complex.gro"
            # Use shutil.copyfile (not copy) — preserves permissions
            # without trying to copy metadata that may not be available
            # in the pixi env.
            import shutil

            def _stage(src, dst):
                # Skip when src and dst resolve to the same file — happens when
                # this node's output_dir is the SAME working_path that produced
                # the inputs (e.g. ep_amber_to_gromacs → md_solvate_gmx chained
                # in one workflow). Copying a file onto itself raises
                # shutil.SameFileError; the file is already staged correctly.
                src_p = Path(src).resolve()
                dst_p = Path(dst).resolve()
                if src_p != dst_p:
                    shutil.copyfile(src_p, dst_p)

            _stage(top_path, staged_top)
            _stage(gro_path, staged_gro)
            if itp_path:
                staged_itp = Path(output_dir) / Path(itp_path).name
                _stage(itp_path, staged_itp)

            stream_log(
                f"Loading solute: top={Path(top_path).name}, "
                f"gro={Path(gro_path).name}"
                + (f", itp={Path(itp_path).name}" if itp_path else ""),
                node_id=self.node_id, progress=5,
            )

            # ── Load solute Structure ──────────────────────────────────
            solute = load_solute_structure(str(staged_top), str(staged_gro))
            solute_charge = compute_solute_charge(solute)
            stream_log(
                f"Solute: {len(solute.atoms)} atoms, "
                f"net charge {solute_charge:+.3f} e",
                node_id=self.node_id, progress=10,
            )

            # ── Parse solvents + ratio ─────────────────────────────────
            solvents_raw = flow_vars["solvents"].get_value() or "MOH:WAT"
            solvent_codes = [
                c.strip() for c in solvents_raw.split(":") if c.strip()
            ]
            want_moh = "MOH" in solvent_codes
            ratio_str = flow_vars["solvent_ratio"].get_value() or (
                "2:1" if want_moh else "1"
            )
            ratio = parse_solvent_ratio(ratio_str, solvent_codes)

            water_model = flow_vars["water_model"].get_value() or "opc"
            cation = flow_vars["cation"].get_value() or "K+"
            anion = flow_vars["anion"].get_value() or "Cl-"
            saltcon = float(flow_vars["saltcon_M"].get_value() or 0.0)
            padding_A = float(flow_vars["padding_A"].get_value() or 12.0)
            seed = int(flow_vars["random_seed"].get_value() or -1)

            # ── Build solvent unit Structures via tleap ────────────────
            stream_log(
                f"Building solvent units via tleap (water={water_model}, "
                f"moh={want_moh}, cation={cation}, anion={anion})",
                node_id=self.node_id, progress=20,
            )
            env = os.environ.copy()
            env["AMBERHOME"] = detect_amberhome()

            tleap_scratch = Path(output_dir) / "tleap_scratch"
            solvent_units = build_solvent_unit_structures(
                str(tleap_scratch),
                water_model=water_model,
                want_moh=want_moh,
                cation=cation,
                anion=anion,
                env=env,
            )

            # ── Compute box + counts ───────────────────────────────────
            box_nm = compute_box_dimensions_nm(
                solute.coordinates, padding_A=padding_A,
            )
            solute_vol_nm3 = estimate_solute_volume_nm3(solute)
            solvent_counts_map = compute_solvent_counts(
                box_nm=box_nm,
                solute_volume_nm3=solute_vol_nm3,
                ratio=ratio,
            )
            ion_counts = compute_ion_counts(
                solute_charge=solute_charge,
                total_volume_nm3=box_nm[0] * box_nm[1] * box_nm[2],
                saltcon_M=saltcon,
                cation=cation, anion=anion,
            )

            # Merge solvent + ion counts; preserve order (solvents first,
            # then ions) so packmol input and ParmEd assembly agree.
            counts: dict = {}
            for code in solvent_codes:
                if solvent_counts_map.get(code, 0) > 0:
                    counts[code] = solvent_counts_map[code]
            for code, n in ion_counts.items():
                if n > 0:
                    counts[code] = n

            stream_log(
                f"Box: {box_nm[0]:.2f} × {box_nm[1]:.2f} × {box_nm[2]:.2f} nm. "
                f"Counts: {counts}",
                node_id=self.node_id, progress=40,
            )

            # ── Write single-molecule PDB templates + solute PDB ───────
            packmol_work = Path(output_dir) / "packmol_scratch"
            packmol_work.mkdir(parents=True, exist_ok=True)
            solvent_pdb_basenames = write_single_molecule_pdbs(
                solvent_units, str(packmol_work),
            )
            solute_pdb_basename = write_solute_pdb(
                solute, str(packmol_work), center=True,
            )

            # ── Run packmol ────────────────────────────────────────────
            stream_log(
                f"Running packmol (seed={seed})...",
                node_id=self.node_id, progress=55,
            )
            write_packmol_input(
                str(packmol_work),
                solute_pdb=solute_pdb_basename,
                solvent_pdbs=solvent_pdb_basenames,
                solvent_counts=counts,
                box_nm=box_nm,
                seed=seed,
            )
            packed_pdb = run_packmol(
                str(packmol_work), env=env,
            )

            # ── Assemble final Structure in ParmEd ─────────────────────
            stream_log(
                "Assembling solvated Structure in ParmEd (no re-derivation)",
                node_id=self.node_id, progress=75,
            )
            system = assemble_solvated_structure(
                solute_structure=solute,
                solvent_units=solvent_units,
                counts=counts,
                packed_pdb_path=packed_pdb,
                box_nm=box_nm,
            )

            # ── Save GROMACS outputs + split ITP ───────────────────────
            stream_log(
                f"Writing complex.top, complex.gro, metallopeptide_solv.itp...",
                node_id=self.node_id, progress=88,
            )
            saved = save_gromacs_outputs(
                system, output_dir,
                prefix="complex",
                itp_filename="metallopeptide_solv",
            )

            # ── Forward outputs ────────────────────────────────────────
            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "output_top": self.format_output_path(saved["top"]),
                "output_gro": self.format_output_path(saved["gro"]),
                "output_itp": (
                    self.format_output_path(saved["itp"])
                    if saved["itp"] else ""
                ),
                "solvent_counts": counts,
                "box_dimensions_nm": list(box_nm),
                "solute_charge": solute_charge,
            }
            result.files["output"] = {
                "top": self.format_output_path(saved["top"]),
                "gro": self.format_output_path(saved["gro"]),
                "itp": (
                    self.format_output_path(saved["itp"])
                    if saved["itp"] else ""
                ),
                "packmol_log": self.format_output_path(
                    str(packmol_work / "packmol.log")),
            }
            result.success = True
            count_summary = ", ".join(f"{k}={v}" for k, v in counts.items())
            result.message = (
                f"Solvated: {count_summary} in "
                f"{box_nm[0]:.2f} nm cubic box → complex.{{top,gro,itp}}"
            )
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("md_solvate_gmx", str(e))
