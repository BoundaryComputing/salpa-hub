"""
gmx-md-relax core — GROMACS relaxation/equilibration.

Supports three protocols:
- full_4step: Legacy pdbmdauto 4-step relaxation with OriHeavy/OriBackBone freezegrps
- em_only: Single steepest-descent energy minimization (quick testing)
- custom: Single grompp+mdrun with user-provided MDP

Full 4-step protocol (progressive unfreezing):
  1. nvt_fixOri         — NVT, freeze all original heavy atoms (water/ions equilibrate)
  2. nvt_fixOriBackbone — NVT, freeze backbone only (side chains relax)
  3. mm1                — Conjugate gradient EM, no restraints (remove bad contacts)
  4. mm2                — Conjugate gradient EM, no restraints (final cleanup)

The OriHeavy and OriBackBone groups must exist in the NDX file (created by gen_gmx_ndx).
Restraints use GROMACS `freezegrps` mechanism (not posres .itp files).
"""

import os
import re
import subprocess
from dataclasses import dataclass
from shlex import quote as _q  # every path in a shell command goes through this

CHECK_MAX_FORCE = 1000.0  # kJ/mol/nm — threshold for "safe" minimization


@dataclass
class RelaxResult:
    """Result of relaxation protocol."""

    output_gro: str = ""
    output_tpr: str = ""
    run_label: str = ""
    max_force: float = 0.0
    em_safe: bool = False
    steps_completed: list = None
    success: bool = False
    log: str = ""

    def __post_init__(self):
        if self.steps_completed is None:
            self.steps_completed = []


def _run(cmd, cwd=None, timeout=3600):
    """Run shell command. Returns (returncode, combined output)."""
    r = subprocess.run(
        cmd, shell=True, capture_output=True, text=True,
        cwd=cwd, timeout=timeout,
    )
    return r.returncode, (r.stdout or "") + "\n" + (r.stderr or "")


def _extract_max_force(log_path):
    """Extract final 'Maximum force' from GROMACS log."""
    max_force = 0.0
    try:
        with open(log_path) as fh:
            for line in fh:
                m = re.search(r"Maximum force\s*=\s*([\d.eE+\-]+)", line)
                if m:
                    max_force = float(m.group(1))
    except (OSError, ValueError):
        pass
    return max_force


def _purge_temp_files(directory):
    """Remove GROMACS backup files and step PDBs."""
    for f in os.listdir(directory):
        if re.match(r"^#.*|^step.*\.pdb$", f):
            try:
                os.remove(os.path.join(directory, f))
            except OSError:
                pass


def _run_grompp_mdrun(
    mdp_file, gro_file, top_file, ndx_file,
    output_dir, run_label, timeout=3600,
):
    """Run a single grompp + mdrun step.

    Returns (success, output_gro, log_text).
    """
    out_tpr = os.path.join(output_dir, f"{run_label}.tpr")

    mdout = os.path.join(output_dir, f"mdout_{run_label}.mdp")
    cmd = (
        f"gmx grompp -f {_q(mdp_file)} -c {_q(gro_file)} -p {_q(top_file)} "
        f"-o {_q(out_tpr)} -n {_q(ndx_file)} "
        f"-po {_q(mdout)} -maxwarn 10"
    )

    rc, out = _run(cmd, cwd=output_dir)
    if rc != 0:
        return False, "", f"grompp({run_label}) failed (rc={rc}):\n{out}"

    cmd = f"gmx mdrun -v -deffnm {_q(run_label)}"
    rc, out = _run(cmd, cwd=output_dir, timeout=timeout)
    if rc != 0:
        return False, "", f"mdrun({run_label}) failed (rc={rc}):\n{out}"

    output_gro = os.path.join(output_dir, f"{run_label}.gro")
    if os.path.exists(output_gro):
        return True, output_gro, f"{run_label}: OK"
    else:
        return False, "", f"mdrun({run_label}) completed but {run_label}.gro not found"


def process_full_4step(
    gro_file, top_file, ndx_file, output_dir,
    mdp_dir, timeout=3600,
):
    """Run the full 4-step relaxation protocol.

    Steps:
      1. nvt_fixOri         — NVT with freezegrps=OriHeavy
      2. nvt_fixOriBackbone — NVT with freezegrps=OriBackBone
      3. mm1                — CG energy minimization (no restraints)
      4. mm2                — CG energy minimization (no restraints)

    Args:
        gro_file:   Input solvated structure (.gro).
        top_file:   Topology file (.top).
        ndx_file:   Index file with OriHeavy/OriBackBone groups.
        output_dir: Working directory.
        mdp_dir:    Directory with bundled MDP files.
        timeout:    Max seconds per step.

    Returns:
        RelaxResult with final output.
    """
    result = RelaxResult()
    log_lines = []
    os.makedirs(output_dir, exist_ok=True)

    steps = [
        ("nvt_fixOri",         os.path.join(mdp_dir, "nvt_fixOri.mdp")),
        ("nvt_fixOriBackbone", os.path.join(mdp_dir, "nvt_fixOriBackbone.mdp")),
        ("mm1",                os.path.join(mdp_dir, "mm_cg.mdp")),
        ("mm2",                os.path.join(mdp_dir, "mm_cg.mdp")),
    ]

    current_gro = gro_file

    for step_label, mdp_file in steps:
        if not os.path.exists(mdp_file):
            log_lines.append(f"{step_label}: SKIP (MDP not found: {mdp_file})")
            continue

        _purge_temp_files(output_dir)

        ok, out_gro, step_log = _run_grompp_mdrun(
            mdp_file=mdp_file,
            gro_file=current_gro,
            top_file=top_file,
            ndx_file=ndx_file,
            output_dir=output_dir,
            run_label=step_label,
            timeout=timeout,
        )

        log_lines.append(step_log)

        if not ok:
            result.log = "\n".join(log_lines)
            return result

        result.steps_completed.append(step_label)
        current_gro = out_gro

    # Extract max force from final EM step
    mm2_log = os.path.join(output_dir, "mm2.log")
    result.max_force = _extract_max_force(mm2_log)
    result.em_safe = result.max_force <= CHECK_MAX_FORCE

    _purge_temp_files(output_dir)

    result.output_gro = current_gro
    result.run_label = "mm2"
    result.success = True
    result.log = "\n".join(log_lines)

    return result


def process_single_step(
    gro_file, top_file, mdp_file, ndx_file,
    output_dir, run_label="em", timeout=3600,
):
    """Run a single grompp+mdrun step (em_only or custom).

    Returns RelaxResult.
    """
    result = RelaxResult(run_label=run_label)
    os.makedirs(output_dir, exist_ok=True)
    _purge_temp_files(output_dir)

    ok, out_gro, step_log = _run_grompp_mdrun(
        mdp_file=mdp_file,
        gro_file=gro_file,
        top_file=top_file,
        ndx_file=ndx_file,
        output_dir=output_dir,
        run_label=run_label,
        timeout=timeout,
    )

    result.log = step_log

    if ok:
        result.output_gro = out_gro
        result.steps_completed = [run_label]

        # Extract max force if EM
        log_file = os.path.join(output_dir, f"{run_label}.log")
        result.max_force = _extract_max_force(log_file)
        result.em_safe = result.max_force <= CHECK_MAX_FORCE
        result.success = True

    _purge_temp_files(output_dir)
    return result
