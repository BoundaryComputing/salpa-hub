"""
pka-gmx-em — BoCoFlow node wrapper.

Protonation state prediction (PDB2PQR/PROPKA) + GROMACS topology generation
(pdb2gmx) + two-step vacuum energy minimisation.

Replaces the legacy FixPkaGmxEM Docker-based node with direct CLI calls and
a proper protonation bridge (PQR residue names → patched PDB → pdb2gmx).

Output: Energy-minimised GRO + topology TOP, ready for solvation/ionisation.
"""

import os
from datetime import datetime

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter,
    FileParameterEdit,
    FloatParameter,
    FolderParameter,
    IntegerParameter,
    SelectParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import process_pka_gmx_em
except ImportError:
    from core import process_pka_gmx_em


class PkaGmxEm(Node):
    """
    Protonation prediction + GROMACS topology + energy minimisation.

    Takes a PDB file (typically from fix_residues or merge step) and:
    1. Predicts protonation states via PDB2PQR/PROPKA at target pH
    2. Bridges PQR → PDB by renaming titratable residues (HID/HIE/HIP, ASH, GLH, CYX)
    3. Generates GROMACS topology via pdb2gmx -ignh on the patched PDB
    4. Creates triclinic simulation box with padding
    5. Runs two-step energy minimisation (unconstrained → h-bonds)

    Output: Energy-minimised GRO + topology TOP ready for solvation.
    """

    name = "pKa + GROMACS EM"
    node_key = "PkaGmxEm"

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name",
            default="",
            docstring="Leave empty to use predecessor data.",
        ),
        "input_pdb": FileParameterEdit(
            "Input PDB File",
            default="",
            docstring="PDB from fix_residues. Leave empty: auto-discovers Merge/fixed.pdb from predecessor working_path.",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            default="",
            docstring="GROMACS working folder. Leave empty: creates gmx/ in case folder. All downstream nodes share this folder.",
        ),
        "force_field": SelectParameter(
            "Force Field",
            options=["amber99sb-ildn", "amber99sb", "charmm27", "oplsaa"],
            default="amber99sb",
            docstring=(
                "GROMACS force field. PDB2PQR uses matching naming: "
                "AMBER (HID/HIE/HIP) for amber*, CHARMM (HSD/HSE/HSP) for charmm*."
            ),
        ),
        "water_model": SelectParameter(
            "Water Model",
            options=["tip3p", "spc", "spce", "tip4p"],
            default="tip3p",
            docstring="Water model for pdb2gmx.",
        ),
        "box_distance": FloatParameter(
            "Box Distance (nm)",
            default=2.0,
            docstring="Padding distance around solute for simulation box.",
        ),
        "em_steps": IntegerParameter(
            "EM Steps",
            default=1000,
            docstring="Maximum steps per EM stage (2 stages total).",
        ),
        "ph": FloatParameter(
            "pH",
            default=7.0,
            docstring="Target pH for PROPKA protonation state prediction.",
        ),
        "run_pdb2pqr": BooleanParameter(
            "Run PDB2PQR",
            default=True,
            docstring=(
                "Run PDB2PQR/PROPKA for pH-dependent protonation. "
                "If disabled, pdb2gmx uses its own hydrogen-bond analysis."
            ),
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute the protonation + topology + EM pipeline."""
        stream_log(
            "Starting pKa + GROMACS EM...",
            node_id=self.node_id, progress=0,
        )

        try:
            result = NodeResult()
            result.metadata["execution_time"] = datetime.now().isoformat()

            input_data = (
                predecessor_data[0]
                if predecessor_data and predecessor_data[0]
                else {}
            )
            case_name = (
                flow_vars["case_name"].get_value()
                or input_data.get("case_name", "protein")
            )

            # ── Resolve input PDB ─────────────────────────────────────────
            input_pdb = flow_vars["input_pdb"].get_value()
            if not input_pdb:
                # Auto-discover from predecessor working_path
                working_path = input_data.get("working_path", "")
                if working_path:
                    search_dir = self.resolve_path(working_path)
                    input_pdb = self._find_pdb(search_dir)

            if not input_pdb:
                raise NodeException(
                    "pka_gmx_em",
                    "No input PDB file found. Provide one or connect to a predecessor.",
                )

            input_pdb = self.resolve_path(input_pdb)
            if not os.path.exists(input_pdb):
                raise NodeException(
                    "pka_gmx_em", f"Input PDB not found: {input_pdb}"
                )

            log_message(f"Input PDB: {input_pdb}")

            # ── Resolve output directory ──────────────────────────────────
            output_dir = flow_vars["output_dir"].get_value()
            if not output_dir:
                # Single gmx/ folder in the case directory for all GROMACS operations
                working_path = input_data.get("working_path", "")
                if working_path:
                    output_dir = os.path.join(self.resolve_path(working_path), "gmx")
                else:
                    output_dir = os.path.join(os.path.dirname(input_pdb), "gmx")
            output_dir = self.resolve_path(output_dir)

            stream_log(
                "Running PDB2PQR + protonation bridge...",
                node_id=self.node_id, progress=10,
            )

            # ── Run pipeline ──────────────────────────────────────────────
            em_result = process_pka_gmx_em(
                input_pdb=input_pdb,
                output_dir=output_dir,
                case_name=case_name,
                force_field=flow_vars["force_field"].get_value(),
                water_model=flow_vars["water_model"].get_value(),
                box_distance=flow_vars["box_distance"].get_value(),
                em_steps=flow_vars["em_steps"].get_value(),
                ph=flow_vars["ph"].get_value(),
                run_pdb2pqr=flow_vars["run_pdb2pqr"].get_value(),
            )

            if not em_result.success:
                log_message(f"Pipeline log:\n{em_result.log}")
                raise NodeException(
                    "pka_gmx_em",
                    f"Pipeline failed:\n{em_result.log[:500]}",
                )

            stream_log(
                f"EM complete — max force: {em_result.em_max_force:.1f} kJ/mol/nm",
                node_id=self.node_id, progress=90,
            )

            # ── Build result ──────────────────────────────────────────────
            # Pass gmx/ as working_path — all downstream nodes use the same folder
            result.data.update({
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "output_gro": self.format_output_path(em_result.em_gro),
                "output_top": self.format_output_path(em_result.em_top),
                "em_max_force": em_result.em_max_force,
            })

            if em_result.patched_pdb:
                result.data["patched_pdb"] = self.format_output_path(
                    em_result.patched_pdb
                )
            if em_result.pqr_file:
                result.data["pqr_file"] = self.format_output_path(
                    em_result.pqr_file
                )
            if em_result.protonation_changes:
                result.data["protonation_changes"] = {
                    f"{ch or '-'}:{seq}": f"{old}->{new}"
                    for (ch, seq), (old, new)
                    in em_result.protonation_changes.items()
                }

            result.files["input"] = {
                "input_pdb": self.format_output_path(input_pdb),
            }
            result.files["output"] = {
                "em_gro": self.format_output_path(em_result.em_gro),
                "em_top": self.format_output_path(em_result.em_top),
            }

            result.success = True
            n_changes = len(em_result.protonation_changes)
            result.message = (
                f"pKa+GROMACS EM complete for {case_name}: "
                f"{n_changes} protonation change(s), "
                f"max force {em_result.em_max_force:.1f} kJ/mol/nm"
            )

            stream_log(result.message, node_id=self.node_id, progress=100)
            log_message(f"Pipeline log:\n{em_result.log}")
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            log_message(f"Error in PkaGmxEm: {e}")
            raise NodeException("pka_gmx_em", str(e))

    @staticmethod
    def _find_pdb(search_dir):
        """Auto-discover a PDB file in search_dir or its Merge/ subfolder."""
        if not os.path.isdir(search_dir):
            return None

        for subdir in [os.path.join(search_dir, "Merge"), search_dir]:
            if not os.path.isdir(subdir):
                continue
            for fname in sorted(os.listdir(subdir)):
                if fname.endswith("_fixed.pdb") or fname.endswith(".pdb"):
                    return os.path.join(subdir, fname)
        return None
