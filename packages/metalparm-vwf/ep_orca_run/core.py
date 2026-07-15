"""ORCA Run — core (Level 1) helpers.

Pure-Python helpers that:
- generate an ORCA input file for the common easyPARM recipe
  (OPT + FREQ + CHELPG), with safe defaults for transition/heavy metals
- run ORCA locally via subprocess and collect the artifacts
  (.out, .hess, optimized .xyz, .property.txt) that ep_seminario_orca expects

ORCA cannot be installed from conda; users provide it via:
  * ``ORCA_BIN`` env var (absolute path to the ``orca`` binary), or
  * a system-wide install on ``$PATH``.

For HPC execution the node uses ``HPCNodeBase`` — see node.py. The SLURM
script is user-provided; see ``templates/default-slurm.sh`` for a reference
template with ``{{VARIABLE}}`` placeholders (formal, Snellius-validated).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Set


# ---------------------------------------------------------------------------
# Heavy-atom / ECP detection
# ---------------------------------------------------------------------------

# def2-ECP is applied to atoms with Z >= 37 (Rb and heavier). Below that,
# def2-SVP is already all-electron and no NewECP block is needed.
_Z: dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
    "Rn": 86, "U": 92,
}

# Metals (for auto-ECP we care about Z >= 37) — full set for logging clarity
_METALS = {sym for sym, z in _Z.items() if z >= 21 and sym not in
           {"Ga", "Ge", "As", "Se", "Br", "Kr",
            "In", "Sn", "Sb", "Te", "I", "Xe",
            "Tl", "Pb", "Bi", "Po", "At", "Rn"}} | {"Sn", "Pb"}


def parse_xyz_elements(xyz_file: str | Path) -> List[str]:
    """Return element symbols (in order) from an XYZ file. Case-normalised."""
    elements: List[str] = []
    text = Path(xyz_file).read_text().splitlines()
    if len(text) < 3:
        return elements
    try:
        n = int(text[0].strip())
    except (ValueError, IndexError):
        return elements
    for line in text[2:2 + n]:
        tok = line.split()
        if not tok:
            continue
        sym = tok[0].strip()
        # Handle 2-char symbols like "Sn"
        sym = sym[:1].upper() + sym[1:].lower() if len(sym) > 1 else sym.upper()
        elements.append(sym)
    return elements


def heavy_ecp_elements(xyz_file: str | Path, z_threshold: int = 37) -> List[str]:
    """Unique element symbols in the XYZ with Z >= threshold (default 37 = Rb).
    These are the atoms that need a def2-ECP NewECP block at def2-SVP level.
    """
    elements = parse_xyz_elements(xyz_file)
    heavy: Set[str] = {e for e in elements if _Z.get(e, 0) >= z_threshold}
    return sorted(heavy)


def auto_ecp_block(elements: list[str], ecp_family: str = "def2-ECP") -> str:
    """Render a NewECP block for the given elements. Empty string if none."""
    if not elements:
        return ""
    lines = ["%basis"]
    for e in elements:
        lines.append(f'  NewECP {e} "{ecp_family}" end')
    lines.append("end")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# ORCA input generation
# ---------------------------------------------------------------------------

@dataclass
class OrcaInputConfig:
    """Configuration for auto-generating a single ORCA input file."""

    xyz_file: str                      # absolute or relative path to the XYZ file
    charge: int = 0
    multiplicity: int = 1
    method: str = "B3LYP D3BJ"         # keywords (append def2 basis etc. separately)
    basis: str = "def2-SVP"            # main basis
    aux_basis: str = "def2/J"          # auxiliary basis (RIJ)
    grid: str = "DEFGRID2"
    run_type: str = "OPT FREQ"         # " " separated keywords; user may use OPT or FREQ only
    extra_keywords: str = ""           # raw ORCA keywords appended to the ! line
    chelpg: bool = True                # append CHELPG population analysis
    resp: bool = False                 # append RESP (ORCA 6 native restrained ESP charges)
    nprocs: int = 4                    # passed to %pal
    memory_mb: int = 4000              # %maxcore per core
    ecp_line: Optional[str] = None     # e.g. '%basis NewECP Sn "def2-ECP" end end' — whole block
    extra_blocks: str = ""             # user-appended blocks (e.g., %geom, %scf)

    def render(self) -> str:
        parts = [self.method, self.basis, self.aux_basis, self.grid,
                 self.run_type, self.extra_keywords]
        if self.chelpg:
            parts.append("CHELPG")
        if self.resp:
            parts.append("RESP")
        keywords = " ".join(p.strip() for p in parts if p and p.strip())
        lines = [f"! {keywords}"]
        lines.append(f"%pal nprocs {self.nprocs} end")
        lines.append(f"%maxcore {self.memory_mb}")
        if self.ecp_line:
            lines.append(self.ecp_line)
        if self.extra_blocks:
            lines.append(self.extra_blocks)
        lines.append("")
        lines.append(f"* xyzfile {self.charge} {self.multiplicity} {os.path.basename(self.xyz_file)}")
        lines.append("")
        return "\n".join(lines)


def write_orca_input(config: OrcaInputConfig, out_path: str | Path) -> str:
    """Write the ORCA input file and ensure the XYZ is in the same dir.

    Returns the absolute path to the written .inp file.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(config.render())
    # Mirror the XYZ next to the .inp (ORCA reads it as basename)
    src = Path(config.xyz_file).resolve()
    if src.is_file():
        dest = out_path.parent / src.name
        if src != dest:
            shutil.copy2(src, dest)
    return str(out_path.resolve())


# ---------------------------------------------------------------------------
# Local execution
# ---------------------------------------------------------------------------

def find_orca_binary(override: str | None = None) -> str:
    """Resolve the ``orca`` binary.

    Priority:
      1. ``override`` argument
      2. ``ORCA_BIN`` env var
      3. ``orca`` on $PATH
    Raises ``FileNotFoundError`` otherwise.
    """
    if override:
        p = Path(override).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
        raise FileNotFoundError(f"ORCA binary not found at override: {override}")
    env_bin = os.environ.get("ORCA_BIN", "").strip()
    if env_bin:
        p = Path(env_bin).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
    which = shutil.which("orca")
    if which:
        return which
    raise FileNotFoundError(
        "ORCA binary not found. Set ORCA_BIN env var or add `orca` to PATH.")


@dataclass
class OrcaRunResult:
    success: bool
    returncode: int
    message: str
    out_file: str = ""                 # .out (main output)
    hess_file: str = ""                # .hess (FREQ)
    xyz_file: str = ""                 # .xyz (optimized)
    property_file: str = ""            # .property.txt (if produced)
    chelpg_file: str = ""              # .chelpg.xyz (if present)
    extras: list[str] = field(default_factory=list)   # other produced files


def _collect_outputs(work_dir: Path, base: str) -> dict[str, str]:
    """Collect expected ORCA output artifacts by prefix."""
    result: dict[str, str] = {}
    patterns = {
        "out": f"{base}.out",
        "hess": f"{base}.hess",
        "xyz": f"{base}.xyz",
        "property": f"{base}.property.txt",
        "chelpg": f"{base}.chelpg.xyz",
    }
    for key, fname in patterns.items():
        p = work_dir / fname
        if p.is_file():
            result[key] = str(p.resolve())
    return result


def run_orca(inp_file: str | Path,
             work_dir: str | Path | None = None,
             orca_bin: str | None = None,
             timeout: int | None = None) -> OrcaRunResult:
    """Run ORCA locally: ``orca input.inp > input.out``.

    ORCA writes several artifacts to the working directory using the input
    file's basename as prefix (.out, .hess, .xyz, .property.txt, ...).
    """
    inp_path = Path(inp_file).resolve()
    if not inp_path.is_file():
        return OrcaRunResult(False, -1, f"ORCA input file not found: {inp_path}")

    wd = Path(work_dir).resolve() if work_dir else inp_path.parent
    wd.mkdir(parents=True, exist_ok=True)
    base = inp_path.stem
    out_path = wd / f"{base}.out"

    try:
        binary = find_orca_binary(orca_bin)
    except FileNotFoundError as e:
        return OrcaRunResult(False, -1, str(e))

    # ORCA uses absolute path to its own binary for MPI launching. Per the
    # ORCA manual the invocation is: `{full_path_to_orca} input.inp > output.out`
    with open(out_path, "w") as out_fh:
        proc = subprocess.run(
            [binary, inp_path.name],
            cwd=wd, stdout=out_fh, stderr=subprocess.STDOUT,
            timeout=timeout, check=False,
        )
    artifacts = _collect_outputs(wd, base)
    success = proc.returncode == 0 and out_path.is_file() and out_path.stat().st_size > 0
    msg = "ORCA finished" if success else f"ORCA returned {proc.returncode}"
    return OrcaRunResult(
        success=success, returncode=proc.returncode, message=msg,
        out_file=str(out_path.resolve()),
        hess_file=artifacts.get("hess", ""),
        xyz_file=artifacts.get("xyz", ""),
        property_file=artifacts.get("property", ""),
        chelpg_file=artifacts.get("chelpg", ""),
        extras=[],
    )
