"""xtb_opt — GFN-xTB geometry optimization (semi-empirical, local).

Fast pre-optimizer that sits between ``snp_builder`` and the expensive
``ep_orca_run`` DFT step. Relaxes the rigid geometric placement from the
builder to a near-minimum at GFN2-xTB level in seconds to minutes, cutting
subsequent DFT convergence time substantially.

Predecessor data flow: consumes ``output_xyz`` from an upstream node
(typically ``snp_builder``). Emits ``output_xyz`` that downstream
``ep_bond_detection``, ``ep_mol2_generation``, and ``ep_orca_run`` can all
pick up via the 3-tier resolution pattern.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Dict

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FolderParameter, IntegerParameter, SelectParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import XtbInputConfig, parse_final_energy, run_xtb
except ImportError:
    try:
        from core import XtbInputConfig, parse_final_energy, run_xtb  # type: ignore
    except ImportError:
        XtbInputConfig = parse_final_energy = run_xtb = None  # type: ignore


def _get_from_predecessors(predecessor_data, key):
    for pred in (predecessor_data or []):
        if pred and key in pred:
            return pred[key]
    return None


class XtbOpt(Node):
    """Run GFN-xTB geometry optimization (local, semi-empirical)."""

    name = "xTB Opt"
    node_key = "XtbOpt"
    category = "Force Field Parameterization"
    tags = ["xtb", "gfn", "geometry-optimization", "semi-empirical",
            "metal-complex", "pre-optimization"]

    num_in = 1
    num_out = 1

    ENVIRONMENT = {
        "type": "pixi",
        "name": "metalparm_vwf",
        "pixi_toml": str(Path(__file__).parent.parent / "pixi.toml"),
    }

    OPTIONS = {
        "case_name": StringParameter(
            label="Case Name", default="",
            docstring="Optional tag for logging (inherits from predecessor if empty)",
            optional=True),
        "xyz_file": FileParameterEdit(
            label="XYZ File", default="",
            docstring="Input geometry. Leave empty to auto-read ``output_xyz`` "
                      "from an upstream node (e.g., snp_builder).",
            optional=True),
        "charge": IntegerParameter(
            label="Charge", default=0,
            docstring="Net molecular charge"),
        "multiplicity": IntegerParameter(
            label="Multiplicity", default=1,
            docstring="Spin multiplicity (2S+1). 1 = singlet, 2 = doublet, …"),
        "method": SelectParameter(
            label="Method", default="GFN2-xTB",
            options=["GFN2-xTB", "GFN1-xTB", "GFN0-xTB", "GFN-FF"],
            docstring="Hamiltonian. GFN2-xTB is the default recommended pre-opt."),
        "opt_level": SelectParameter(
            label="Optimization Level", default="normal",
            options=["crude", "loose", "normal", "tight", "verytight", "extreme"],
            docstring="xtb --opt threshold. 'normal' is enough for a DFT seeder."),
        "solvent": StringParameter(
            label="Solvent (ALPB)", default="",
            docstring="Implicit solvent name (e.g. 'methanol', 'water'). "
                      "Blank = vacuum.",
            optional=True),
        "extra_args": StringParameter(
            label="Extra xtb Flags", default="",
            docstring="Free-form CLI flags appended to the xtb command.",
            optional=True),
        "xtb_bin": StringParameter(
            label="xtb Binary (override)", default="",
            docstring="Leave blank to use $XTB_BIN or 'xtb' on PATH.",
            optional=True),
        "output_dir": FolderParameter(
            label="Output Directory", default="",
            docstring="Directory for the xtb run. Defaults to <working_path>/xtb_opt.",
            optional=True),
    }

    def _resolve_xyz(self, predecessor_data, flow_vars) -> str:
        explicit = flow_vars["xyz_file"].get_value()
        if explicit:
            return self.resolve_path(explicit)
        from_pred = _get_from_predecessors(predecessor_data, "output_xyz")
        if from_pred:
            return self.resolve_path(from_pred)
        raise NodeException(
            "xtb_opt",
            "No XYZ: set 'XYZ File' in the node panel or connect to a node that "
            "declares 'output_xyz' (e.g. snp_builder).",
        )

    def execute(self, predecessor_data, flow_vars):
        if XtbInputConfig is None:
            raise NodeException(
                "xtb_opt",
                "core helpers not importable. Run this node in the "
                "metalparm_vwf pixi env (needs xtb + numpy).",
            )

        stream_log("Starting xtb geometry optimization…",
                   node_id=self.node_id, progress=0)

        result = NodeResult()

        case_name = flow_vars["case_name"].get_value() or \
            _get_from_predecessors(predecessor_data, "case_name") or "complex"
        xyz_file = self._resolve_xyz(predecessor_data, flow_vars)

        output_dir = flow_vars["output_dir"].get_value()
        if not output_dir:
            work = _get_from_predecessors(predecessor_data, "working_path")
            output_dir = os.path.join(
                self.resolve_path(work) if work else ".", "xtb_opt")
        else:
            output_dir = self.resolve_path(output_dir)

        cfg = XtbInputConfig(
            xyz_file=xyz_file,
            charge=int(flow_vars["charge"].get_value()),
            multiplicity=int(flow_vars["multiplicity"].get_value()),
            method=flow_vars["method"].get_value(),
            opt_level=flow_vars["opt_level"].get_value(),
            solvent=flow_vars["solvent"].get_value() or None,
            extra_args=flow_vars["extra_args"].get_value() or "",
        )
        stream_log(f"Running {cfg.method} / opt={cfg.opt_level} "
                   f"(charge={cfg.charge}, mult={cfg.multiplicity})",
                   node_id=self.node_id, progress=15)

        xtb_bin = flow_vars["xtb_bin"].get_value() or None
        res = run_xtb(cfg, work_dir=output_dir, xtb_bin=xtb_bin)

        if not res.success:
            log_message(f"xtb_opt failed: {res.message}", level="error")
            raise NodeException("xtb_opt", res.message)

        # Promote xtbopt.xyz to a predictable case-named copy for downstream predecessors
        opt_xyz_named = os.path.join(output_dir, f"{case_name}_xtbopt.xyz")
        if os.path.abspath(res.opt_xyz_file) != os.path.abspath(opt_xyz_named):
            import shutil as _sh
            _sh.copy2(res.opt_xyz_file, opt_xyz_named)

        energy = parse_final_energy(res.log_file)
        stream_log(f"xtb done. Final energy: "
                   f"{energy:.6f} Eh" if energy is not None else "xtb done.",
                   node_id=self.node_id, progress=100)

        result.data = {
            "case_name": case_name,
            "working_path": self.format_output_path(output_dir),
            "output_xyz": self.format_output_path(opt_xyz_named),
            "output_log": self.format_output_path(res.log_file),
            "final_energy_eh": energy,
            "method": cfg.method,
            "opt_level": cfg.opt_level,
        }
        result.files["input"] = {"xyz": self.format_output_path(xyz_file)}
        result.files["output"] = {
            "opt_xyz": result.data["output_xyz"],
            "log":     result.data["output_log"],
        }
        result.success = True
        result.message = (f"xtb {cfg.method} opt complete "
                          f"(E = {energy:.6f} Eh)" if energy is not None
                          else f"xtb {cfg.method} opt complete")
        return result.to_json()
