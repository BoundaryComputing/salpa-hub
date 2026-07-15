"""MD Solvate (packmol-memgen) — wraps AmberTools' packmol-memgen to
solvate a metallopeptide AMBER topology in a (mixed) solvent box.

The packmol-memgen tool itself orchestrates packmol (for packing) + tleap
(for building the AMBER topology around the packed coordinates). All the
heavy lifting — densities, atom-type compatibility with ff19SB, AMBER-
compatible cosolvent parameter libraries — already lives in the tool. We
are a thin wrapper that translates BoCoFlow node options into the right
command-line flags + handles predecessor auto-discovery + parses the
log to surface molecule counts on the result object.

Designed to sit between ``ep_apply_coords`` (which produces the dry
metallopeptide AMBER topology) and ``ep_amber_to_gromacs`` (which
converts the solvated AMBER topology to GROMACS):

    ep_apply_coords → md_solvate_packmol → ep_amber_to_gromacs

Forwarded data:
  ``output_prmtop``, ``output_rst7``, ``output_pdb``, ``case_name``,
  ``working_path``, plus per-solvent counts.
"""

from __future__ import annotations

import os
import shutil
import subprocess

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FloatParameter, FolderParameter, SelectParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log

# Pure-Python helpers, no bocoflow_core deps — unit-testable in isolation.
try:
    from .helpers import (
        detect_amberhome, get_from_predecessors,
        parse_rst7_box, parse_solvent_counts,
    )
except ImportError:  # script-mode fallback
    from helpers import (  # type: ignore
        detect_amberhome, get_from_predecessors,
        parse_rst7_box, parse_solvent_counts,
    )


# ─── Constants ────────────────────────────────────────────────────────

# packmol-memgen names its outputs by prefixing the input PDB's stem
# with "solvated_" and suffixing with the water-model marker "_wat".
# So `complex.pdb` becomes `solvated_complex_wat.{top,crd,pdb}`.
_OUTPUT_SUFFIX = "_wat"

# packmol-memgen valid water models (must match its --ffwat choices)
_VALID_WATER_MODELS = ["opc", "tip3p", "tip4pew", "tip4pd", "opc3",
                       "spce", "spceb", "fb3"]
# Valid protein force fields (must match its --ffprot choices)
_VALID_FFPROTS = ["ff19SB", "ff14SB", "ff15ipq"]


# ─── Node ─────────────────────────────────────────────────────────────


class MdSolvatePackmol(Node):
    """Solvate a metallopeptide AMBER topology via packmol-memgen."""

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Working dir; node writes solvated_*.{top,crd,pdb} here.",
        ),

        # ── Inputs (auto-discover from ep_apply_coords) ───────────────
        "input_pdb": FileParameterEdit(
            "Input PDB",
            default="",
            docstring=(
                "Metallopeptide PDB (post-coord-transfer). Leave empty to "
                "auto-discover ``output_pdb`` from an ep_apply_coords "
                "predecessor."
            ),
            optional=True,
        ),
        "fragment_lib": FileParameterEdit(
            "Fragment Library (.lib/.off)",
            default="",
            docstring=(
                "AMBER .lib for the non-canonical residue in the input PDB "
                "(e.g. the SnP fragment unit 'mol'). Leave empty to auto-"
                "discover ``output_lib`` from an upstream library-generation "
                "node. packmol-memgen passes this to tleap so it can load "
                "the residue during topology build."
            ),
            optional=True,
        ),
        "fragment_frcmod": FileParameterEdit(
            "Fragment Frcmod",
            default="",
            docstring=(
                "AMBER frcmod paired with fragment_lib. Auto-discovers "
                "``output_frcmod`` from a predecessor when blank."
            ),
            optional=True,
        ),
        "linkage_frcmod": FileParameterEdit(
            "Linkage Frcmod",
            default="",
            docstring=(
                "Optional cross-FF amide linkage frcmod (e.g. the bundled "
                "``amide_glh_gaff2_n.frcmod``). Loaded via packmol-memgen's "
                "--leapline so tleap picks up the cross-FF parameters when "
                "building the solvated topology."
            ),
            optional=True,
        ),

        # ── Solvent + box ─────────────────────────────────────────────
        "solvents": StringParameter(
            "Solvents",
            default="MOH:WAT",
            docstring=(
                "Solvent codes separated by ':'. Pass verbatim to "
                "``packmol-memgen --solvents``. Run ``packmol-memgen "
                "--available_solvents`` to list valid codes (MOH=methanol, "
                "WAT=water, CL3=chloroform, DMS=DMSO, NMA, ACN, ACT, BNZ, "
                "IPH, TFE). For single-solvent water: 'WAT'."
            ),
        ),
        "solvent_ratio": StringParameter(
            "Solvent Ratio",
            default="2:1",
            docstring=(
                "Numerical ratios separated by ':' — pass verbatim to "
                "``packmol-memgen --solvent_ratio``. The basis (molar / "
                "volume / mass) is whatever packmol-memgen uses; inspect "
                "actual molecule counts in the result's ``n_solvent_*`` "
                "fields to verify. For single-solvent: '1'."
            ),
        ),
        "padding_A": FloatParameter(
            "Box Padding (Å)",
            default=12.0,
            docstring=(
                "Distance from the solute to each box face. Passed as "
                "``packmol-memgen --dist``. Cubic box only (the simpler "
                "default). 12 Å is the AMBER convention; 10 Å is common "
                "for small soluble peptides."
            ),
        ),

        # ── FF pairings ───────────────────────────────────────────────
        "ff_protein": SelectParameter(
            "Protein Force Field",
            default="ff19SB",
            options=_VALID_FFPROTS,
            docstring=(
                "Protein FF for tleap. ff19SB is the canonical modern "
                "choice; pair with OPC water for the FF's design "
                "conditions. Passed as ``--ffprot``."
            ),
        ),
        "water_model": SelectParameter(
            "Water Model",
            default="opc",
            options=_VALID_WATER_MODELS,
            docstring=(
                "Water model. **OPC is recommended for ff19SB** (the "
                "force field was parameterised against OPC in the "
                "original paper). TIP3P works but is the design pair for "
                "ff14SB. Passed as ``--ffwat``."
            ),
        ),

        # ── Ions ──────────────────────────────────────────────────────
        "saltcon_M": FloatParameter(
            "Salt Concentration (M)",
            default=0.0,
            docstring=(
                "KCl concentration beyond neutralization (mol/L). 0.0 = "
                "just neutralize the system. Typical physiological-mimic "
                "value: 0.15. Passed as ``--saltcon``."
            ),
        ),

        # ── Advanced tuning ───────────────────────────────────────────
        "extra_leaplines": StringParameter(
            "Extra LEaP Lines",
            default="",
            docstring=(
                "Optional extra ``--leapline`` arguments (one per logical "
                "line, separated by ';;'). For most cases the auto-"
                "generated leap lines (loading linkage_frcmod when set) "
                "are sufficient."
            ),
            optional=True,
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting packmol-memgen solvation...",
                   node_id=self.node_id, progress=0)

        try:
            result = NodeResult()

            case_name = (flow_vars["case_name"].get_value()
                         or "complex")
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            # ── Input resolution: explicit OR predecessor → real abs path ─
            # packmol-memgen accepts absolute paths for --pdb,
            # --ligand_param, and abs paths inside --leapline (verified
            # 2026-05-11). So we just resolve the prefix-tagged BoCoFlow
            # paths to their real filesystem locations and pass them
            # through verbatim — no need to copy/stage inputs into
            # output_dir. The subprocess still runs with cwd=output_dir
            # so the tool's *output* files land where we want them, but
            # *input* files are read from wherever they already live.
            #
            # See dev-notes/node-io-helpers-duplication.md for the
            # rationale (and the earlier defensive stage_local() that
            # this replaces).

            def resolve_input(explicit_ref, *predecessor_keys):
                """Return an absolute filesystem path for the resolved
                explicit option OR the first matching predecessor key, or
                '' if neither yields a real file.
                """
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

            pdb = resolve_input(
                flow_vars["input_pdb"].get_value(), "output_pdb",
            )
            if not pdb:
                raise NodeException(
                    "setup",
                    "input_pdb not provided and no output_pdb in predecessor data.",
                )

            frag_lib = resolve_input(
                flow_vars["fragment_lib"].get_value(), "output_lib",
            )
            frag_frcmod = resolve_input(
                flow_vars["fragment_frcmod"].get_value(), "output_frcmod",
            )
            linkage_frcmod = resolve_input(
                flow_vars["linkage_frcmod"].get_value(),
                "output_linkage_frcmod",
            )

            # ── Build packmol-memgen command ─────────────────────────
            # Output naming is derived from the input PDB's basename
            # by packmol-memgen: complex.pdb → solvated_complex_wat.{top,crd,pdb}
            cmd = [
                shutil.which("packmol-memgen") or "packmol-memgen",
                "--pdb", pdb,                       # ← abs path; no staging
                "--solvate",
                "--solvents", flow_vars["solvents"].get_value() or "MOH:WAT",
                "--solvent_ratio", flow_vars["solvent_ratio"].get_value() or "2:1",
                "--dist", str(float(flow_vars["padding_A"].get_value() or 12.0)),
                "--ffprot", flow_vars["ff_protein"].get_value() or "ff19SB",
                "--ffwat", flow_vars["water_model"].get_value() or "opc",
                "--parametrize",
                "--notprotonate",
                "--log", "packmol_memgen.log",
            ]

            saltcon = float(flow_vars["saltcon_M"].get_value() or 0.0)
            if saltcon > 0:
                cmd += ["--saltcon", str(saltcon)]
            else:
                # We still need --salt_override so packmol-memgen doesn't
                # warn at concentration=0; the neutralizing ions are added
                # regardless.
                cmd += ["--salt_override", "--saltcon", "0.0"]

            if frag_frcmod and frag_lib:
                # --ligand_param expects a colon-joined FRCMOD:LIB pair.
                # Both can be absolute paths.
                cmd += [
                    "--ligand_param", f"{frag_frcmod}:{frag_lib}",
                    "--gaff2",  # our SnP fragment is GAFF2-typed
                ]

            leaplines = []
            if linkage_frcmod:
                # tleap accepts absolute paths in `loadamberparams`.
                leaplines.append(f"loadamberparams {linkage_frcmod}")
            extra = (flow_vars["extra_leaplines"].get_value() or "").strip()
            if extra:
                # Split on ';;' for multi-line input convenience
                leaplines.extend(
                    ln.strip() for ln in extra.split(";;") if ln.strip()
                )
            for ln in leaplines:
                cmd += ["--leapline", ln]

            stream_log(
                f"Solvents: {flow_vars['solvents'].get_value()}  "
                f"ratio: {flow_vars['solvent_ratio'].get_value()}  "
                f"padding: {flow_vars['padding_A'].get_value()} Å  "
                f"FF: {flow_vars['ff_protein'].get_value()} + "
                f"{flow_vars['water_model'].get_value()}",
                node_id=self.node_id, progress=10,
            )
            stream_log(
                f"Command: {' '.join(cmd[:8])} ... "
                f"(see packmol_memgen.log for full output)",
                node_id=self.node_id, progress=15,
            )

            # ── Run packmol-memgen ───────────────────────────────────
            env = os.environ.copy()
            env["AMBERHOME"] = detect_amberhome()

            proc = subprocess.run(
                cmd, cwd=output_dir, env=env,
                capture_output=True, text=True, timeout=1800,
            )
            log_text = (proc.stdout or "") + "\n" + (proc.stderr or "")

            if proc.returncode != 0:
                # Surface the last 30 lines of output to stream_log + raise
                tail = "\n".join(log_text.strip().splitlines()[-30:])
                stream_log(
                    f"packmol-memgen failed (rc={proc.returncode}); log tail:\n{tail}",
                    node_id=self.node_id, progress=50,
                )
                raise NodeException(
                    "execution",
                    f"packmol-memgen exited with code {proc.returncode}. "
                    f"See packmol_memgen.log + leap_.log + packmol.log in "
                    f"{output_dir} for details. Last lines:\n{tail}",
                )

            # ── Resolve output paths ─────────────────────────────────
            # packmol-memgen names outputs as solvated_<stem>_wat.{top,crd,pdb}
            stem = os.path.splitext(os.path.basename(pdb))[0]
            out_prefix = f"solvated_{stem}{_OUTPUT_SUFFIX}"
            out_top = os.path.join(output_dir, f"{out_prefix}.top")
            out_crd = os.path.join(output_dir, f"{out_prefix}.crd")
            out_pdb = os.path.join(output_dir, f"{out_prefix}.pdb")

            for p in (out_top, out_crd, out_pdb):
                if not os.path.isfile(p):
                    raise NodeException(
                        "execution",
                        f"packmol-memgen finished but expected output {p} "
                        f"is missing. Check the log files in {output_dir}.",
                    )

            # ── Parse counts from log ────────────────────────────────
            counts = parse_solvent_counts(log_text)
            box_dims = parse_rst7_box(out_crd)

            stream_log(
                f"Solvated. Counts: "
                f"{', '.join(f'{k}={v}' for k, v in sorted(counts.items()))}",
                node_id=self.node_id, progress=90,
            )

            # ── Forward outputs ──────────────────────────────────────
            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                # Downstream nodes (ep_amber_to_gromacs) expect prmtop/rst7.
                # packmol-memgen emits .top + .crd — these are the same
                # formats with different conventional extensions in the
                # AmberTools toolchain. ParmEd reads either.
                "output_prmtop": self.format_output_path(out_top),
                "output_rst7":   self.format_output_path(out_crd),
                "output_pdb":    self.format_output_path(out_pdb),
                "solvent_counts": counts,
                "box_dimensions_A": box_dims,
            }
            result.files["output"] = {
                "prmtop": self.format_output_path(out_top),
                "rst7":   self.format_output_path(out_crd),
                "pdb":    self.format_output_path(out_pdb),
                "log":    self.format_output_path(
                    os.path.join(output_dir, "packmol_memgen.log")),
            }
            result.success = True
            count_summary = (
                ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
                or "no counts parsed (check log)"
            )
            result.message = (
                f"packmol-memgen solvated: {count_summary} → "
                f"{out_prefix}.{{top,crd,pdb}}"
            )
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("packmol-memgen solvation", str(e))
