"""
gmx-solv-ion core — GROMACS solvation and ionization.

Pipeline:
1. gmx editconf — define simulation box
2. gmx solvate — add SPC216 water molecules
3. gmx grompp — preprocess for genion
4. gmx genion — add Na+/Cl- ions to neutralize
5. gmx make_ndx — generate/update index groups
"""

import glob
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field


@dataclass
class SolvIonResult:
    """Result of solvation and ionization."""

    output_gro: str = ""
    output_top: str = ""
    output_ndx: str = ""
    box_gro: str = ""
    solv_gro: str = ""
    ion_gro: str = ""
    success: bool = False
    log: str = ""


def _extract_custom_groups(ndx_path: str) -> str:
    """Extract OriHeavy and OriBackBone groups from an NDX file as text.

    Returns the group text (including headers) to be appended after regeneration,
    or empty string if no custom groups found.
    """
    if not os.path.exists(ndx_path):
        return ""
    result_lines = []
    capturing = False
    with open(ndx_path) as f:
        for line in f:
            if line.strip() in ("[ OriHeavy ]", "[ OriBackBone ]"):
                capturing = True
                result_lines.append(line)
                continue
            if line.startswith("[") and capturing:
                # New group that's not ours — check if it's the other custom group
                if line.strip() in ("[ OriHeavy ]", "[ OriBackBone ]"):
                    result_lines.append(line)
                    continue
                capturing = False
                continue
            if capturing:
                result_lines.append(line)
    return "".join(result_lines) if result_lines else ""


def _run_gmx(cmd: str, cwd: str = None) -> tuple:
    """Run a GROMACS shell command. Returns (returncode, stdout+stderr)."""
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, cwd=cwd, timeout=300,
    )
    combined = result.stdout + "\n" + result.stderr
    return result.returncode, combined


def _purge_temp_files(directory: str):
    """Remove GROMACS temporary files (#backups and step*.pdb)."""
    for pattern in [r"^#.*", r"^step.*\.pdb$"]:
        for f in os.listdir(directory):
            if re.match(pattern, f):
                os.remove(os.path.join(directory, f))


def process_solv_ion(
    gro_file: str,
    top_file: str,
    mdp_file: str,
    ndx_file: str,
    output_dir: str,
    case_name: str,
    run_label: str = "md",
    box_size: str = "20 20 20",
    ion_conc: float = 0.15,
    scale_fill: float = 0.57,
) -> SolvIonResult:
    """Run GROMACS solvation and ionization pipeline.

    Args:
        gro_file: Input structure (.gro).
        top_file: Topology file (.top).
        mdp_file: MD parameters file (.mdp).
        ndx_file: Index file (.ndx).
        output_dir: Working directory.
        case_name: Case identifier.
        run_label: Simulation label.
        box_size: Box dimensions "X Y Z" in nm ("0 0 0" for triclinic auto).
        ion_conc: Ion concentration in mol/L (0 to skip).
        scale_fill: Van der Waals scale factor for solvation density.

    Returns:
        SolvIonResult with output file paths.
    """
    result = SolvIonResult()
    log_lines = []
    os.makedirs(output_dir, exist_ok=True)

    # Copy topology to working dir (GROMACS solvate/genion modify it in-place)
    # In single gmx/ folder mode, top_file and output_dir are the same directory
    # so we just copy to a working name
    work_top = os.path.join(output_dir, "topol.top")
    if os.path.abspath(top_file) != os.path.abspath(work_top):
        shutil.copy2(top_file, work_top)

    # Step 1: editconf — define box
    box = [float(x) for x in box_size.split()]
    box_gro = os.path.join(output_dir, "box.gro")

    if box[0] > 0:
        cmd = f"gmx editconf -f {gro_file} -o {box_gro} -box {box[0]} {box[1]} {box[2]}"
    else:
        cmd = f"gmx editconf -f {gro_file} -o {box_gro} -bt triclinic -d 2.0"

    rc, out = _run_gmx(cmd, cwd=output_dir)
    log_lines.append(f"editconf: rc={rc}")
    if rc != 0:
        result.log = "\n".join(log_lines) + "\n" + out
        return result
    result.box_gro = box_gro

    # Step 2: solvate — add water
    solv_gro = os.path.join(output_dir, "solv.gro")
    cmd = f"gmx solvate -cp {box_gro} -cs spc216.gro -p {work_top} -o {solv_gro} -scale {scale_fill}"
    rc, out = _run_gmx(cmd, cwd=output_dir)
    log_lines.append(f"solvate: rc={rc}")
    if rc != 0:
        result.log = "\n".join(log_lines) + "\n" + out
        return result
    result.solv_gro = solv_gro

    # Step 3: grompp — preprocess for genion
    ion_tpr = os.path.join(output_dir, "ion.tpr")
    cmd = f"gmx grompp -f {mdp_file} -c {solv_gro} -p {work_top} -o {ion_tpr} -maxwarn 10"
    rc, out = _run_gmx(cmd, cwd=output_dir)
    log_lines.append(f"grompp: rc={rc}")
    if rc != 0:
        result.log = "\n".join(log_lines) + "\n" + out
        return result

    # Step 4: genion — add ions
    ion_gro = os.path.join(output_dir, "ion.gro")
    if ion_conc > 0:
        cmd = (
            f"echo SOL | gmx genion -s {ion_tpr} -p {work_top} "
            f"-o {ion_gro} -neutral -nname CL -pname NA -conc {ion_conc}"
        )
        rc, out = _run_gmx(cmd, cwd=output_dir)
        log_lines.append(f"genion: rc={rc}")
        if rc != 0:
            result.log = "\n".join(log_lines) + "\n" + out
            return result
    else:
        shutil.copy2(solv_gro, ion_gro)
        log_lines.append("genion: skipped (ion_conc=0)")
    result.ion_gro = ion_gro

    # Regenerate NDX for the solvated structure (adds Water, SOL, Ion, non-Protein groups)
    # but PRESERVE OriHeavy/OriBackBone from gen_gmx_ndx
    out_ndx = ndx_file
    ori_groups = _extract_custom_groups(ndx_file)  # save OriHeavy/OriBackBone
    cmd = f"echo q | gmx make_ndx -f {ion_gro} -o {out_ndx}"
    _run_gmx(cmd, cwd=output_dir)
    if ori_groups:
        with open(out_ndx, "a") as f:
            f.write(ori_groups)
        log_lines.append("make_ndx: regenerated + preserved OriHeavy/OriBackBone")

    # Clean up temp files
    _purge_temp_files(output_dir)

    result.output_gro = ion_gro
    result.output_top = work_top
    result.output_ndx = out_ndx
    result.success = True
    result.log = "\n".join(log_lines)

    return result
