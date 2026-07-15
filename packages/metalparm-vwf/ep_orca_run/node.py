"""ORCA Run Node — local or HPC/SLURM execution.

Wraps the ORCA quantum-chemistry binary. Inputs: an XYZ geometry (e.g., from
``snp_builder``). Outputs: the artifacts that ``ep_seminario_orca`` expects —
``.out`` (main log), ``.hess`` (Hessian), optimized ``.xyz``, and optional
``.property.txt``. Fits between the geometry builder and the Seminario node.

Execution modes (inherited from HPCNodeBase):
- ``local``  — runs ``orca input.inp`` in a subprocess (needs ``ORCA_BIN`` env
  var or ``orca`` on PATH; ORCA is not available via conda)
- ``remote`` — submits a user-provided SLURM script to an HPC cluster
  via the HPC profile. A formal, Snellius-validated reference template
  ships at ``ep_orca_run/templates/default-slurm.sh`` (matches the
  ``pdbmdauto/gmx_mdrun`` convention). Copy that file's contents into
  the "SLURM Job Script" field in the GUI and edit the cluster-specific
  lines (partition, module tree, memory, scratch). The ``{{VARIABLE}}``
  placeholders are substituted by ``HPCNodeBase`` at submit time. The
  textarea does not auto-prefill yet — see
  ``dev-notes/slurm-script-prefill-not-yet-supported.md`` for the
  underlying GUI-layer gap.

Input modes:
- ``input_mode = "auto"`` (default) — generate an ORCA input file from
  (xyz_file, charge, multiplicity, method, basis, chelpg, ecp, extra_blocks).
  Good default for the easyPARM recipe (``OPT FREQ CHELPG``).
- ``input_mode = "file"`` — user supplies a pre-written ``.inp`` file.

Predecessor data flow: auto-discovers ``output_xyz`` from an upstream builder
(e.g., ``snp_builder``), following the 3-tier pattern used by the rest of
metalparm-vwf.
"""

from __future__ import annotations

import glob
import json
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List

from bocoflow_core.hpc_node import HPCNodeBase
from bocoflow_core.logger import log_message
from bocoflow_core.node import NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit, FolderParameter, IntegerParameter, SelectParameter,
    StringParameter, TextParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import (
        OrcaInputConfig, OrcaRunResult, auto_ecp_block, find_orca_binary,
        heavy_ecp_elements, run_orca, write_orca_input,
    )
except ImportError:  # direct-path import (node_runner)
    try:
        from core import (  # type: ignore
            OrcaInputConfig, OrcaRunResult, auto_ecp_block, find_orca_binary,
            heavy_ecp_elements, run_orca, write_orca_input,
        )
    except ImportError:  # server env (introspect OPTIONS only)
        OrcaInputConfig = OrcaRunResult = None  # type: ignore
        find_orca_binary = run_orca = write_orca_input = None  # type: ignore
        auto_ecp_block = heavy_ecp_elements = None  # type: ignore


def _get_from_predecessors(predecessor_data, key):
    for pred in (predecessor_data or []):
        if pred and key in pred:
            return pred[key]
    return None


def _copy_to_workdir(src: str, dest_dir: str) -> str:
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, os.path.basename(src))
    if os.path.abspath(src) != os.path.abspath(dest):
        shutil.copy2(src, dest)
    return dest


def _resolve_ecp_block(xyz_path: str, user_block: str,
                       node_id: str | None = None) -> str | None:
    """Prefer the user-supplied ecp_block. If empty and the XYZ contains atoms
    with Z >= 37 (def2-ECP threshold), auto-generate a NewECP block for each
    heavy element. Returns None when no ECP block is needed.
    """
    user = (user_block or "").strip()
    if user:
        return user
    if heavy_ecp_elements is None:  # server-env stub
        return None
    try:
        heavy = heavy_ecp_elements(xyz_path)
    except Exception:
        return None
    if not heavy:
        return None
    block = auto_ecp_block(heavy)
    log_message(f"ep_orca_run: auto-detected heavy elements {heavy}; "
                f"injecting def2-ECP block.", level="info")
    try:
        stream_log(f"Auto-ECP applied to {', '.join(heavy)}",
                   node_id=node_id, progress=None)
    except Exception:
        pass
    return block


class EpOrcaRun(HPCNodeBase):
    """Run ORCA (OPT+FREQ+CHELPG recipe by default) locally or on SLURM."""

    name = "ORCA Run"
    node_key = "EpOrcaRun"
    category = "Force Field Parameterization"
    tags = ["orca", "qm", "dft", "hessian", "chelpg", "hpc", "slurm", "metal-complex"]

    num_in = 1
    num_out = 1

    ENVIRONMENT = {
        "type": "pixi",
        "name": "metalparm_vwf",
        "pixi_toml": str(Path(__file__).parent.parent / "pixi.toml"),
    }

    OPTIONS = {
        **HPCNodeBase.HPC_OPTIONS,
        "case_name": StringParameter(
            label="Case Name", default="",
            docstring="Optional tag for logging (inherits from predecessor if empty)",
            optional=True),
        "run_label": StringParameter(
            label="Run Label", default="orca",
            docstring="Prefix for ORCA input/output files (input.inp, input.out, input.hess, …)"),
        "output_dir": FolderParameter(
            label="Working Directory",
            docstring="Directory for the ORCA input + output files"),
        "xyz_file": FileParameterEdit(
            label="XYZ File", default="",
            docstring="Starting geometry in XYZ format. Leave empty to auto-discover "
                      "from predecessor (output_xyz)",
            optional=True),
        "input_mode": SelectParameter(
            label="Input Mode", default="auto",
            options=["auto", "file"],
            docstring="auto: generate ORCA input from OPTIONS; file: use user-supplied .inp"),
        "orca_input_file": FileParameterEdit(
            label="ORCA Input File (.inp)", default="",
            docstring="Used only when input_mode=file",
            optional=True),
        # --- auto-mode ORCA recipe settings ---
        "charge": IntegerParameter(
            label="Total Charge", default=0,
            docstring="Net molecular charge"),
        "multiplicity": IntegerParameter(
            label="Spin Multiplicity", default=1,
            docstring="2S+1 of the electronic ground state"),
        "method": StringParameter(
            label="Method/Functional", default="B3LYP D3BJ",
            docstring="ORCA keywords for the functional/method (e.g., B3LYP D3BJ, PBE0, TPSS)"),
        "basis": StringParameter(
            label="Basis Set", default="def2-SVP",
            docstring="Main basis set keyword (def2-SVP is a safe default for medium-size systems)"),
        "aux_basis": StringParameter(
            label="Auxiliary Basis", default="def2/J",
            docstring="Auxiliary basis for RIJCOSX (def2/J matches def2-SVP/TZVP)"),
        "run_type": StringParameter(
            label="Run Type", default="OPT FREQ",
            docstring="ORCA run-type keywords (OPT FREQ for easyPARM; use OPT, FREQ, SP separately if desired)"),
        "chelpg": SelectParameter(
            label="Compute CHELPG Charges", default="yes",
            options=["yes", "no"],
            docstring="Adds CHELPG analysis — required by ep_seminario_orca downstream"),
        "resp": SelectParameter(
            label="Compute RESP Charges", default="no",
            options=["no", "yes"],
            docstring="Adds ORCA 6 native !RESP (restrained ESP) charges — the "
            "AMBER-standard choice for ep_charges. Parsed from the .out."),
        "ecp_block": TextParameter(
            label="ECP Block (optional)", default="",
            docstring="Full ORCA basis-block for ECPs (e.g., for Sn, I). Example:\n"
                      '%basis\n  NewECP Sn "def2-ECP" end\nend'),
        "extra_blocks": TextParameter(
            label="Extra ORCA Blocks (optional)", default="",
            docstring="Additional ORCA blocks appended to the input (%geom, %scf, %method …)"),
        "nprocs": IntegerParameter(
            label="Processors (local)", default=4,
            docstring="Passed as %pal nprocs for local runs; HPC runs typically set this via SLURM"),
        "memory_mb": IntegerParameter(
            label="Memory per core (MB)", default=4000,
            docstring="Passed as %maxcore"),
        "orca_bin": StringParameter(
            label="ORCA Binary Path (local)", default="",
            docstring="Absolute path to ORCA binary. Falls back to ORCA_BIN env var, then `orca` on PATH",
            optional=True),
    }

    # ------------------------------------------------------------------
    # HPCNodeBase hooks
    # ------------------------------------------------------------------

    def _resolve_xyz(self, predecessor_data, flow_vars, work_dir) -> str:
        """Explicit config → predecessor → error."""
        xyz = self.resolve_path(flow_vars["xyz_file"].get_value()) or ""
        if xyz and os.path.isfile(xyz):
            return _copy_to_workdir(xyz, work_dir)
        ref = _get_from_predecessors(predecessor_data, "output_xyz")
        if ref:
            resolved = self.resolve_path(ref)
            if resolved and os.path.isfile(resolved):
                return _copy_to_workdir(resolved, work_dir)
        raise NodeException("setup",
            "No XYZ file — set xyz_file explicitly or connect an upstream "
            "builder that emits output_xyz.")

    def _prepare_inp_file(self, predecessor_data, flow_vars, work_dir) -> tuple[str, str]:
        """Return (inp_path, run_label). Works for both auto and file modes."""
        run_label = flow_vars["run_label"].get_value() or "orca"
        mode = flow_vars["input_mode"].get_value() or "auto"
        if mode == "file":
            src = self.resolve_path(flow_vars["orca_input_file"].get_value()) or ""
            if not src or not os.path.isfile(src):
                raise NodeException("setup",
                    f"input_mode=file but orca_input_file not found: {src}")
            local = _copy_to_workdir(src, work_dir)
            # Rename to match run_label so all outputs share the prefix
            inp_path = os.path.join(work_dir, f"{run_label}.inp")
            if os.path.abspath(local) != os.path.abspath(inp_path):
                shutil.move(local, inp_path)
            return inp_path, run_label

        # auto mode — need the XYZ
        xyz = self._resolve_xyz(predecessor_data, flow_vars, work_dir)
        chelpg = (flow_vars["chelpg"].get_value() or "yes") == "yes"
        resp = (flow_vars["resp"].get_value() or "no") == "yes"
        cfg = OrcaInputConfig(
            xyz_file=xyz,
            charge=int(flow_vars["charge"].get_value() or 0),
            multiplicity=int(flow_vars["multiplicity"].get_value() or 1),
            method=flow_vars["method"].get_value() or "B3LYP D3BJ",
            basis=flow_vars["basis"].get_value() or "def2-SVP",
            aux_basis=flow_vars["aux_basis"].get_value() or "def2/J",
            run_type=flow_vars["run_type"].get_value() or "OPT FREQ",
            chelpg=chelpg,
            resp=resp,
            nprocs=int(flow_vars["nprocs"].get_value() or 4),
            memory_mb=int(flow_vars["memory_mb"].get_value() or 4000),
            ecp_line=_resolve_ecp_block(
                xyz,
                flow_vars["ecp_block"].get_value() or "",
                node_id=getattr(self, "node_id", None),
            ),
            extra_blocks=(flow_vars["extra_blocks"].get_value() or "").strip() or "",
        )
        inp_path = os.path.join(work_dir, f"{run_label}.inp")
        write_orca_input(cfg, inp_path)
        return inp_path, run_label

    def get_input_files(self, flow_vars: dict) -> List[str]:
        """Files to transfer to HPC.

        The .inp and any auto-generated content are written into the working
        directory at ``run_local()`` time; on HPC we need them BEFORE submit,
        so the caller (HPCNodeBase.execute -> submit) calls this with the
        flow_vars at submission time. We materialize the input + xyz in the
        working dir here and return their absolute paths.
        """
        work_dir = self.resolve_path(flow_vars["output_dir"].get_value())
        os.makedirs(work_dir, exist_ok=True)
        # NB: HPCNodeBase calls this from the worker context; predecessor
        # data is not reachable here. For HPC mode the user must set
        # xyz_file explicitly OR use input_mode=file.
        files: list[str] = []
        mode = flow_vars["input_mode"].get_value() or "auto"
        if mode == "file":
            src = flow_vars["orca_input_file"].get_value()
            if src:
                resolved = self.resolve_path(src)
                if resolved and os.path.isfile(resolved):
                    files.append(_copy_to_workdir(resolved, work_dir))
        else:
            xyz = flow_vars["xyz_file"].get_value()
            if xyz:
                resolved = self.resolve_path(xyz)
                if resolved and os.path.isfile(resolved):
                    local_xyz = _copy_to_workdir(resolved, work_dir)
                    files.append(local_xyz)
                    # also write the .inp next to it so get_template_variables
                    # + SLURM can reference it
                    run_label = flow_vars["run_label"].get_value() or "orca"
                    chelpg = (flow_vars["chelpg"].get_value() or "yes") == "yes"
                    resp = (flow_vars["resp"].get_value() or "no") == "yes"
                    cfg = OrcaInputConfig(
                        xyz_file=local_xyz,
                        charge=int(flow_vars["charge"].get_value() or 0),
                        multiplicity=int(flow_vars["multiplicity"].get_value() or 1),
                        method=flow_vars["method"].get_value() or "B3LYP D3BJ",
                        basis=flow_vars["basis"].get_value() or "def2-SVP",
                        aux_basis=flow_vars["aux_basis"].get_value() or "def2/J",
                        run_type=flow_vars["run_type"].get_value() or "OPT FREQ",
                        chelpg=chelpg,
                        resp=resp,
                        nprocs=int(flow_vars["nprocs"].get_value() or 4),
                        memory_mb=int(flow_vars["memory_mb"].get_value() or 4000),
                        ecp_line=_resolve_ecp_block(
                            local_xyz,
                            flow_vars["ecp_block"].get_value() or "",
                            node_id=getattr(self, "node_id", None),
                        ),
                        extra_blocks=(flow_vars["extra_blocks"].get_value() or "").strip() or "",
                    )
                    inp_path = os.path.join(work_dir, f"{run_label}.inp")
                    write_orca_input(cfg, inp_path)
                    files.append(inp_path)
        return files

    def get_output_files(self, flow_vars: dict) -> List[str]:
        run_label = flow_vars["run_label"].get_value() or "orca"
        return [
            f"{run_label}.out",
            f"{run_label}.hess",
            f"{run_label}.xyz",
            f"{run_label}.property.txt",
            f"{run_label}.chelpg.xyz",
        ]

    def get_output_files_by_category(self, flow_vars: dict) -> Dict[str, List[str]]:
        run_label = flow_vars["run_label"].get_value() or "orca"
        return {
            "essential": [f"{run_label}.out", f"{run_label}.hess", f"{run_label}.xyz"],
            "standard":  [f"{run_label}.property.txt", f"{run_label}.chelpg.xyz"],
            "large":     [f"{run_label}.gbw", f"{run_label}.densities"],
        }

    def get_template_variables(self, flow_vars: dict) -> Dict[str, str]:
        run_label = flow_vars["run_label"].get_value() or "orca"
        mode = flow_vars["input_mode"].get_value() or "auto"
        orca_input = f"{run_label}.inp"
        xyz_file = flow_vars["xyz_file"].get_value() or ""
        return {
            "RUN_LABEL": run_label,
            "ORCA_INPUT_FILE": orca_input,
            "ORCA_OUTPUT_FILE": f"{run_label}.out",
            "INPUT_MODE": mode,
            "XYZ_FILE": os.path.basename(xyz_file) if xyz_file else "",
            "NPROCS": str(flow_vars["nprocs"].get_value() or 4),
        }

    def run_local(self, predecessor_data: list, flow_vars: dict) -> dict:
        stream_log("Starting ORCA run (local)...", node_id=self.node_id, progress=0)
        if run_orca is None:
            raise NodeException("setup", "core.py failed to import (missing deps).")

        result = NodeResult()
        result.metadata["execution_time"] = datetime.now().isoformat()

        input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}
        case_name = flow_vars["case_name"].get_value() or input_data.get("case_name", "orca")
        work_dir = self.resolve_path(flow_vars["output_dir"].get_value())
        os.makedirs(work_dir, exist_ok=True)

        try:
            stream_log("Preparing ORCA input...", node_id=self.node_id, progress=10)
            inp_path, run_label = self._prepare_inp_file(predecessor_data, flow_vars, work_dir)

            binary_override = flow_vars["orca_bin"].get_value() or None
            stream_log(f"Running ORCA ({run_label})...", node_id=self.node_id, progress=30)
            res = run_orca(inp_path, work_dir=work_dir, orca_bin=binary_override)

            if not res.success:
                # Tail the .out if present for a helpful error
                tail = ""
                if res.out_file and os.path.isfile(res.out_file):
                    with open(res.out_file, "rb") as fh:
                        fh.seek(0, os.SEEK_END)
                        size = fh.tell()
                        fh.seek(max(0, size - 2000))
                        tail = fh.read().decode("utf-8", errors="replace")
                raise NodeException("execution",
                    f"{res.message}\n--- {run_label}.out tail ---\n{tail}")

            stream_log("ORCA completed.", node_id=self.node_id, progress=100)
            result.files["input"]["inp"] = self.format_output_path(inp_path)
            result.files["output"]["out"] = self.format_output_path(res.out_file)
            if res.hess_file:
                result.files["output"]["hess"] = self.format_output_path(res.hess_file)
            if res.xyz_file:
                result.files["output"]["xyz"] = self.format_output_path(res.xyz_file)
            if res.property_file:
                result.files["output"]["property"] = self.format_output_path(res.property_file)

            result.data.update({
                "case_name": case_name,
                "run_label": run_label,
                "working_path": self.format_output_path(work_dir),
                "output_out": self.format_output_path(res.out_file),
                "output_hess": self.format_output_path(res.hess_file) if res.hess_file else "",
                "output_xyz_opt": self.format_output_path(res.xyz_file) if res.xyz_file else "",
                "output_property": self.format_output_path(res.property_file) if res.property_file else "",
                "returncode": res.returncode,
            })
            result.success = True
            result.message = f"ORCA OK — {run_label}.{'hess,' if res.hess_file else ''}out"
            return json.loads(result.to_json())
        except NodeException:
            raise
        except Exception as e:
            raise NodeException("orca_run", str(e))

    def build_hpc_result_data(self, flow_vars: dict, local_dir: str) -> dict:
        run_label = flow_vars["run_label"].get_value() or "orca"
        case_name = flow_vars["case_name"].get_value() or "orca"
        produced = {}
        for ext, key in [(".out", "output_out"), (".hess", "output_hess"),
                         (".xyz", "output_xyz_opt"),
                         (".property.txt", "output_property"),
                         (".chelpg.xyz", "output_chelpg")]:
            p = os.path.join(local_dir, f"{run_label}{ext}")
            if os.path.isfile(p):
                produced[key] = self.format_output_path(p)
        return {
            "case_name": case_name,
            "run_label": run_label,
            "working_path": self.format_output_path(local_dir),
            "execution_mode": "remote_hpc",
            **produced,
        }
