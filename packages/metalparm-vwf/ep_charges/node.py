"""EasyParm Charge Injection Node — assigns QM-derived partial charges to the
zero-charge MOL2 that antechamber produces.

This is the step the original node refactor dropped: antechamber (in
ep_mol2_generation) emits a MOL2 whose charge column is all 0.0, and EasyParm's
01_easyPARM.sh injects QM charges as a *separate* final step. Without this node
the parameterized fragment carries zero electrostatics.

Sits between MOL2 Generation and Force Field Assembly, and also consumes the
ORCA run:

    ep_mol2_generation ─┐
                        ├─→ ep_charges ─→ ep_forcefield_assembly
    ep_orca_run ────────┘

Auto-discovers ``output_mol2`` from MOL2 Generation and ``output_out`` (ORCA
.out) from ORCA Run. Forwards every predecessor key downstream (so FF Assembly
still sees Bond Detection / Seminario outputs), overriding ``output_mol2`` with
the now-charged copy.

Methods (mirroring EasyParm's ORCA charge menu):
  * orca_resp   — ORCA 6 native !RESP block (recommended, AMBER-standard)
  * orca_chelpg — ORCA CHELPG block (plain ESP fit)
  * resp_vpot   — classic EasyParm path: .vpot grid → AmberTools `resp`
"""

import os
import shutil

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FolderParameter, IntegerParameter, SelectParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from . import core as charges_core
except ImportError:  # script-mode / server introspection fallback
    try:
        import core as charges_core  # type: ignore
    except ImportError:
        charges_core = None


def _get_from_predecessors(predecessor_data, key):
    for pred in (predecessor_data or []):
        if pred and key in pred:
            return pred[key]
    return None


class EpCharges(Node):
    """Inject QM-derived (RESP / CHELPG) partial charges into the MOL2."""

    name = "Charges"
    num_in = 2
    num_out = 1

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "method": SelectParameter(
            "Charge Method",
            default="orca_resp",
            options=["orca_resp", "orca_chelpg", "resp_vpot"],
            docstring="orca_resp: ORCA 6 native !RESP block (recommended). "
            "orca_chelpg: ORCA CHELPG block. "
            "resp_vpot: classic EasyParm path (ORCA .vpot grid → AmberTools resp).",
        ),
        "mol2_file": FileParameterEdit(
            "MOL2 File",
            default="",
            optional=True,
            docstring="Zero-charge MOL2 from MOL2 Generation. "
            "Leave empty to auto-discover (output_mol2).",
        ),
        "orca_out": FileParameterEdit(
            "ORCA .out",
            default="",
            optional=True,
            docstring="ORCA output containing the RESP/CHELPG block. "
            "Leave empty to auto-discover from ORCA Run (output_out).",
        ),
        "charge": IntegerParameter(
            "Total Charge",
            default=0,
            docstring="Net charge (used only by resp_vpot for the resp control file).",
        ),
        "vpot_file": FileParameterEdit(
            "ESP grid (.vpot)", default="", optional=True,
            docstring="resp_vpot only: ORCA ESP grid file.",
        ),
        "similar_file": FileParameterEdit(
            "Atom equivalence (similar.dat)", default="", optional=True,
            docstring="resp_vpot only: EasyParm atom-equivalence file.",
        ),
        "geom_file": FileParameterEdit(
            "Geometry (.hess/$atoms)", default="", optional=True,
            docstring="resp_vpot only: geometry file with an $atoms block (e.g. .hess).",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Where the charged COMPLEX.mol2 is written.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Injecting QM charges...", node_id=self.node_id, progress=0)
        try:
            result = NodeResult()
            method = flow_vars["method"].get_value()
            case_name = (flow_vars["case_name"].get_value()
                         or _get_from_predecessors(predecessor_data, "case_name")
                         or "complex")
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            # --- resolve MOL2 (explicit → predecessor) and copy into work_dir ---
            mol2 = self.resolve_path(flow_vars["mol2_file"].get_value() or "")
            if not mol2:
                ref = _get_from_predecessors(predecessor_data, "output_mol2")
                mol2 = self.resolve_path(ref) if ref else ""
            if not mol2 or not os.path.isfile(mol2):
                raise NodeException("charges", "No MOL2 found (set mol2_file or "
                                    "connect MOL2 Generation).")
            local_mol2 = os.path.join(output_dir, "COMPLEX.mol2")
            if os.path.abspath(mol2) != os.path.abspath(local_mol2):
                shutil.copy2(mol2, local_mol2)

            stream_log(f"Method: {method}", node_id=self.node_id, progress=20)

            kwargs = dict(charge=flow_vars["charge"].get_value(), work_dir=output_dir,
                          scripts_dir=os.path.join(self._node_dir or ".", "scripts"))
            if method in ("orca_resp", "orca_chelpg"):
                orca_out = self.resolve_path(flow_vars["orca_out"].get_value() or "")
                if not orca_out:
                    ref = _get_from_predecessors(predecessor_data, "output_out")
                    orca_out = self.resolve_path(ref) if ref else ""
                if not orca_out or not os.path.isfile(orca_out):
                    raise NodeException("charges", "No ORCA .out found (set orca_out "
                                        "or connect ORCA Run).")
                kwargs["orca_out"] = orca_out
            else:  # resp_vpot
                kwargs.update(
                    vpot_file=self.resolve_path(flow_vars["vpot_file"].get_value() or ""),
                    similar_file=self.resolve_path(flow_vars["similar_file"].get_value() or ""),
                    geom_file=self.resolve_path(flow_vars["geom_file"].get_value() or ""),
                )

            stream_log("Parsing + injecting charges...", node_id=self.node_id, progress=50)
            summary = charges_core.run_charges(local_mol2, method, **kwargs)

            stream_log(
                f"Injected {summary['n_atoms']} charges "
                f"(net {summary['net_charge']:+.4f}, "
                f"max |q| {summary['max_abs_charge']:.3f})",
                node_id=self.node_id, progress=90,
            )

            # Forward all predecessor data, then override the MOL2 with the charged copy.
            merged = {}
            for pred in (predecessor_data or []):
                if pred:
                    merged.update(pred)
            merged.update({
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "output_mol2": self.format_output_path(local_mol2),
                "charge_method": method,
                "charge_net": summary["net_charge"],
            })
            result.data = merged
            result.files["output"] = {"mol2": self.format_output_path(local_mol2)}
            result.success = True
            result.message = (f"Charges injected ({method}): net "
                              f"{summary['net_charge']:+.4f} e over "
                              f"{summary['n_atoms']} atoms")
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("charges", str(e))
