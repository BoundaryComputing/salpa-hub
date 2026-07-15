"""Pure-Python helpers for the xtb_opt node — runs GFN-xTB geometry optimization.

Thin shell-out wrapper around the `xtb` CLI. No bocoflow_core imports here so
the helpers can be unit-tested in the plain pipeline env against a fake xtb
binary.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Input configuration
# ---------------------------------------------------------------------------

_GFN_CHOICES = {"GFN2-xTB": "2", "GFN1-xTB": "1", "GFN0-xTB": "0", "GFN-FF": "ff"}
_OPT_LEVELS = {"crude", "loose", "normal", "tight", "verytight", "extreme"}


@dataclass
class XtbInputConfig:
    xyz_file: str
    charge: int = 0
    multiplicity: int = 1         # 2S+1. xtb takes --uhf N where N = mult - 1
    method: str = "GFN2-xTB"      # one of _GFN_CHOICES
    opt_level: str = "normal"     # one of _OPT_LEVELS
    solvent: Optional[str] = None # ALPB implicit solvent (e.g., "methanol")
    extra_args: str = ""          # free-form extra CLI args

    def as_argv(self, xyz_basename: str) -> list[str]:
        """Build the xtb command-line for a run in the directory holding xyz_basename."""
        if self.method not in _GFN_CHOICES:
            raise ValueError(f"Unknown method {self.method!r}. "
                             f"Valid: {sorted(_GFN_CHOICES)}")
        if self.opt_level not in _OPT_LEVELS:
            raise ValueError(f"Unknown opt_level {self.opt_level!r}. "
                             f"Valid: {sorted(_OPT_LEVELS)}")
        argv = [
            "xtb", xyz_basename,
            "--opt", self.opt_level,
            "--gfn", _GFN_CHOICES[self.method],
            "--chrg", str(self.charge),
            "--uhf", str(max(self.multiplicity - 1, 0)),
        ]
        if self.solvent:
            argv += ["--alpb", self.solvent]
        if self.extra_args:
            argv += self.extra_args.split()
        return argv


# ---------------------------------------------------------------------------
# Binary resolution
# ---------------------------------------------------------------------------

def find_xtb_binary(override: Optional[str] = None) -> str:
    """Locate the xtb binary.

    Priority:
      1. Explicit override path (file must exist + be executable)
      2. $XTB_BIN env var
      3. PATH
    """
    if override:
        p = Path(override)
        if not p.is_file() or not os.access(p, os.X_OK):
            raise FileNotFoundError(f"xtb override not found or not executable: {p}")
        return str(p)
    env = os.environ.get("XTB_BIN")
    if env:
        p = Path(env)
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
    path_bin = shutil.which("xtb")
    if path_bin:
        return path_bin
    raise FileNotFoundError(
        "xtb binary not found. Set XTB_BIN or install xtb in the pixi env.")


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

@dataclass
class XtbRunResult:
    success: bool
    returncode: int
    message: str
    work_dir: str = ""
    log_file: str = ""
    opt_xyz_file: str = ""


def run_xtb(config: XtbInputConfig,
            work_dir: str,
            xtb_bin: Optional[str] = None,
            timeout: Optional[int] = None) -> XtbRunResult:
    """Run xtb in ``work_dir``. Copies the input xyz next to the run and returns
    paths to the captured log + the optimized xyz (xtbopt.xyz).
    """
    xyz_src = Path(config.xyz_file).resolve()
    if not xyz_src.is_file():
        return XtbRunResult(False, -1, f"Input XYZ missing: {xyz_src}")

    try:
        binary = xtb_bin or find_xtb_binary()
    except FileNotFoundError as e:
        return XtbRunResult(False, -1, str(e))

    wd = Path(work_dir).resolve()
    wd.mkdir(parents=True, exist_ok=True)
    local_xyz = wd / xyz_src.name
    if xyz_src != local_xyz:
        shutil.copy2(xyz_src, local_xyz)

    argv = config.as_argv(local_xyz.name)
    argv[0] = binary  # replace the 'xtb' stub with the resolved path

    log_path = wd / f"{local_xyz.stem}.xtb.out"
    try:
        with open(log_path, "w") as log:
            proc = subprocess.run(
                argv, cwd=wd, stdout=log, stderr=subprocess.STDOUT,
                timeout=timeout,
            )
    except subprocess.TimeoutExpired:
        return XtbRunResult(False, -1, f"xtb timed out after {timeout}s",
                            work_dir=str(wd), log_file=str(log_path))

    opt_xyz = wd / "xtbopt.xyz"
    if proc.returncode != 0:
        return XtbRunResult(False, proc.returncode,
                            f"xtb failed (rc={proc.returncode}), see {log_path}",
                            work_dir=str(wd), log_file=str(log_path))
    if not opt_xyz.is_file():
        return XtbRunResult(False, proc.returncode,
                            f"xtb exit 0 but xtbopt.xyz missing in {wd}",
                            work_dir=str(wd), log_file=str(log_path))

    return XtbRunResult(True, proc.returncode, "xtb opt completed",
                        work_dir=str(wd), log_file=str(log_path),
                        opt_xyz_file=str(opt_xyz))


def parse_final_energy(log_path: str) -> Optional[float]:
    """Grep the xtb log for the final total energy (Eh). Returns None if absent."""
    try:
        text = Path(log_path).read_text(errors="replace")
    except OSError:
        return None
    # Look for "TOTAL ENERGY   -xxx.yyy Eh"
    for line in reversed(text.splitlines()):
        if "TOTAL ENERGY" in line.upper() and "Eh" in line:
            for tok in line.split():
                try:
                    return float(tok)
                except ValueError:
                    continue
    return None
