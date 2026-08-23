"""Module 5: AutoDock4 docking execution.

Generate GPF/DPF parameter files → run autogrid4 → run autodock4 →
extract docked poses.

All paths and tool locations are explicit arguments.
No os.chdir() — uses cwd= parameter for subprocess calls.
"""

import logging
import math
import re
import os
import shutil
import subprocess
from pathlib import Path

import networkx as nx
import numpy as np

from .utils import get_lj_params, INTERNAL_PARAM_METALS, resolve_mgltools_interpreter

logger = logging.getLogger(__name__)


# ===================================================================
# Public interface
# ===================================================================

def run_autodock(
    ligand_pdbqt: Path,
    receptor_pdbqt: Path,
    output_dir: Path,
    graph: nx.Graph,
    metal_symbol: str,
    parameter_file: Path,
    num_poses: int = 10,
    box_center: list[float] | None = None,
    box_size: list[float] | None = None,
    scale_factor: float = 0.0,
    random_pos: bool = True,
    lj_params: list[float] | None = None,
    ga_dock: bool = True,
    ga_settings: dict | None = None,
    sa_dock: bool = False,
    sa_settings: dict | None = None,
    autogrid4_path: str = "autogrid4",
    autodock4_path: str = "autodock4",
    prepare_gpf_script: str | Path | None = None,
    python_path: str = "python3",
    mgltools_dir: str | Path | None = None,
    atom_index_mapping: dict | None = None,
    vacant_site: bool = True,
) -> dict:
    """Run the full AutoDock4 docking pipeline.

    Args:
        ligand_pdbqt: Path to the ligand PDBQT file.
        receptor_pdbqt: Path to the receptor PDBQT file.
        output_dir: Directory for docking output files.
        graph: The molecular graph (for box center calculation).
        metal_symbol: Metal atom symbol.
        parameter_file: Path to AutoDock4 parameter file (metal_dock.dat).
        num_poses: Number of docking poses to generate.
        box_center: [x, y, z] center of the docking box. If None, uses metal position.
        box_size: [x, y, z] box dimensions in Angstrom. Default [20, 20, 20].
        scale_factor: If >0 and box_size is all zeros, auto-calculate box size.
        random_pos: Randomize initial ligand position/orientation.
        lj_params: Custom [e_NA, e_OA, e_SA, e_HD] for the metal. None = use standard set.
        ga_dock: Use Genetic Algorithm.
        ga_settings: GA parameter overrides.
        sa_dock: Use Simulated Annealing.
        sa_settings: SA parameter overrides.
        autogrid4_path: Path to autogrid4 binary.
        autodock4_path: Path to autodock4 binary.
        prepare_gpf_script: Path to prepare_gpf4.py (or derived from mgltools_dir).
        python_path: Python interpreter path.
        mgltools_dir: Path to AutoDockTools directory.
        atom_index_mapping: Mapping from graph node → pdbqt_index for PDB output.
        vacant_site: Whether dummy atoms were added.

    Returns:
        Dict with keys:
        - ``dlg_path``: path to the docking log file.
        - ``pose_pdbqt_paths``: list of extracted pose PDBQT files.
        - ``pose_xyz_paths``: list of extracted pose XYZ files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Copy inputs to docking dir
    dock_ligand = output_dir / ligand_pdbqt.name
    dock_receptor = output_dir / receptor_pdbqt.name
    if not dock_ligand.exists():
        shutil.copy2(ligand_pdbqt, dock_ligand)
    if not dock_receptor.exists():
        shutil.copy2(receptor_pdbqt, dock_receptor)

    name_ligand = ligand_pdbqt.stem
    name_protein = receptor_pdbqt.stem

    # Determine box parameters
    if box_center is None:
        box_center = _get_metal_position(graph, metal_symbol)
    if box_size is None:
        box_size = [20.0, 20.0, 20.0]

    npts = _calculate_npts(box_size, scale_factor, graph, metal_symbol)

    # Resolve LJ params
    sym_upper = metal_symbol.strip().upper()
    internal_param = sym_upper in INTERNAL_PARAM_METALS
    if not internal_param and lj_params is None:
        lj_params = get_lj_params(metal_symbol)

    # Reference the parameter library by bare name (see _stage_parameter_file)
    parameter_file = _stage_parameter_file(parameter_file, output_dir)

    # GPF
    gpf_path = output_dir / f"{name_ligand}_{name_protein}.gpf"
    _create_gpf_file(
        dock_ligand, dock_receptor, gpf_path, parameter_file,
        npts, box_center, metal_symbol, lj_params, internal_param,
        prepare_gpf_script=prepare_gpf_script,
        python_path=python_path,
        mgltools_dir=mgltools_dir,
    )

    # DPF
    dpf_path = output_dir / f"{name_ligand}_{name_protein}.dpf"
    _create_dpf_file(
        dpf_path, gpf_path, parameter_file,
        name_ligand, name_protein, num_poses,
        random_pos, ga_dock, ga_settings, sa_dock, sa_settings,
    )

    # Run autogrid + autodock
    dlg_path = output_dir / f"{name_ligand}_{name_protein}.dlg"
    _run_autogrid(output_dir / "autogrid.log", gpf_path, autogrid4_path, cwd=output_dir)
    _run_autodock(output_dir / "autodock.log", dpf_path, autodock4_path, cwd=output_dir)

    # Extract poses
    pose_pdbqt_paths = _write_conformations(dlg_path, output_dir, name_ligand)

    # Clean dummy atoms
    if vacant_site:
        _clean_dummy_atoms_from_pdbqt(pose_pdbqt_paths)

    # Convert to XYZ
    pose_xyz_paths = _pdbqt_poses_to_xyz(pose_pdbqt_paths, name_ligand)

    return {
        "dlg_path": dlg_path,
        "pose_pdbqt_paths": pose_pdbqt_paths,
        "pose_xyz_paths": pose_xyz_paths,
    }


# ===================================================================
# Box parameter helpers
# ===================================================================

def _get_metal_position(graph: nx.Graph, metal_symbol: str) -> list[float]:
    for _, data in graph.nodes(data=True):
        if data.get("element") == metal_symbol:
            return list(data["xyz"])
    raise ValueError(f"Metal '{metal_symbol}' not found in graph")


def _calculate_npts(
    box_size: list[float],
    scale_factor: float,
    graph: nx.Graph,
    metal_symbol: str,
) -> list[int]:
    """Convert box size in Angstrom to grid points."""
    if any(x != 0 for x in box_size) and scale_factor == 0:
        npts_raw = [x * 2.66 for x in box_size]
        return [math.ceil(x) for x in npts_raw]

    if all(x == 0 for x in box_size) and scale_factor != 0:
        return _auto_box_size(graph, metal_symbol, 0.375, scale_factor)

    # Default
    return [math.ceil(20 * 2.66)] * 3


def _auto_box_size(graph: nx.Graph, metal_symbol: str, spacing: float, scale_factor: float) -> list[int]:
    """Auto-calculate box size from ligand extent."""
    x_all, y_all, z_all = [], [], []
    metal_xyz = None

    for _, data in graph.nodes(data=True):
        xyz = data["xyz"]
        x_all.append(xyz[0])
        y_all.append(xyz[1])
        z_all.append(xyz[2])
        if data.get("element") == metal_symbol:
            metal_xyz = np.array(xyz)

    if metal_xyz is None:
        raise ValueError(f"Metal '{metal_symbol}' not found")

    dists = []
    for arr, c in [(x_all, metal_xyz[0]), (y_all, metal_xyz[1]), (z_all, metal_xyz[2])]:
        shifted = np.array(arr) - c
        d = min(abs(shifted.max() - shifted.min()) * scale_factor, 20)
        dists.append(d)

    return [(round(d / spacing)) & (-2) for d in dists]


# ===================================================================
# GPF / DPF file creation
# ===================================================================

def _stage_parameter_file(parameter_file: Path, output_dir: Path) -> Path:
    """Copy the parameter library into output_dir; return its bare filename.

    autogrid4/autodock4 read the GPF and DPF as whitespace-delimited text and
    the format has no quoting syntax, so an absolute ``parameter_file`` path
    containing a space is silently truncated at that space.  The packaged app
    installs nodes under ``~/Library/Application Support/``, which contains
    one.  Every other path in these files is already a bare filename resolved
    against the run directory, so stage this one the same way.
    """
    src = Path(parameter_file)
    dst = Path(output_dir) / src.name
    if not (dst.exists() and src.exists() and src.samefile(dst)):
        shutil.copyfile(src, dst)
    return Path(src.name)


def _create_gpf_file(
    ligand_pdbqt: Path,
    receptor_pdbqt: Path,
    gpf_path: Path,
    parameter_file: Path,
    npts: list[int],
    box_center: list[float],
    metal_symbol: str,
    lj_params: list[float] | None,
    internal_param: bool,
    prepare_gpf_script: str | Path | None = None,
    python_path: str = "python3",
    mgltools_dir: str | Path | None = None,
) -> Path:
    """Generate a Grid Parameter File for autogrid4."""
    if prepare_gpf_script is None:
        if mgltools_dir is None:
            raise ValueError("Provide prepare_gpf_script or mgltools_dir")
        prepare_gpf_script = Path(mgltools_dir) / "prepare_gpf4.py"

    # pythonsh re-splits arguments on whitespace (see resolve_mgltools_interpreter),
    # which destroys any ligand/receptor path containing a space.
    interpreter, env_overrides = resolve_mgltools_interpreter(python_path)
    env = {**os.environ, **env_overrides} if env_overrides else None

    cmd = [
        interpreter,
        str(prepare_gpf_script),
        "-l", str(ligand_pdbqt),
        "-r", str(receptor_pdbqt),
        "-p", f"parameter_file={parameter_file}",
        "-p", f"npts={npts[0]},{npts[1]},{npts[2]}",
        "-p", f"gridcenter={box_center[0]:.6f},{box_center[1]:.6f},{box_center[2]:.6f}",
        "-o", str(gpf_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(gpf_path.parent), env=env)
    if result.returncode != 0:
        raise RuntimeError(
            f"prepare_gpf4 failed (exit {result.returncode}) writing {gpf_path.name}.\n"
            f"stderr: {result.stderr.strip()}"
        )

    # Append metal-specific LJ parameters
    with open(gpf_path, "a") as f:
        # Standard Zn parameters (always included)
        f.write("nbp_r_eps 0.25 23.2135   12 6  NA TZ\n")
        f.write("nbp_r_eps 2.10  3.8453   12 6  OA Zn\n")
        f.write("nbp_r_eps 2.25  7.5914   12 6  SA Zn\n")
        f.write("nbp_r_eps 1.00  0.0000   12 6  HD Zn\n")
        f.write("nbp_r_eps 2.00  0.0060   12 6  NA Zn\n")
        f.write("nbp_r_eps 2.00  0.2966   12 6  N  Zn\n")

        if not internal_param and lj_params is not None:
            f.write(f"nbp_r_eps 2.20  {lj_params[0]:>.4f}   12 10 NA {metal_symbol}\n")
            f.write(f"nbp_r_eps 2.25  {lj_params[1]:>.4f}   12 10 OA {metal_symbol}\n")
            f.write(f"nbp_r_eps 2.30  {lj_params[2]:>.4f}   12 10 SA {metal_symbol}\n")
            f.write(f"nbp_r_eps 1.00  {lj_params[3]:>.4f}   12 6  HD {metal_symbol}\n")

    return gpf_path


def _create_dpf_file(
    dpf_path: Path,
    gpf_path: Path,
    parameter_file: Path,
    name_ligand: str,
    name_protein: str,
    num_poses: int,
    random_pos: bool,
    ga_dock: bool,
    ga_settings: dict | None,
    sa_dock: bool,
    sa_settings: dict | None,
) -> Path:
    """Generate a Docking Parameter File for autodock4."""
    # Parse ligand types from GPF
    with open(gpf_path) as f:
        gpf_lines = [line.split() for line in f]

    for _tokens in gpf_lines:
        if _tokens and _tokens[0] == "ligand_types":
            ligand_type = _tokens[1:]
            break
    else:
        raise RuntimeError(
            f"No 'ligand_types' line in {gpf_path}; the grid parameter file is "
            f"malformed and autogrid4 would produce unusable maps."
        )
    # Remove trailing grid-related entries
    while ligand_type and ligand_type[-1] in ("#", "atom", "types", "in", "ligand"):
        ligand_type.pop()
    ligand_type_str = " ".join(ligand_type)

    # Default GA settings
    ga = {
        "pop_size": 150, "num_evals": 2500000, "num_generations": 27000,
        "elitism": 1, "mutation_rate": 0.02, "crossover_rate": 0.80,
        "window_size": 10,
    }
    if ga_settings:
        ga.update(ga_settings)

    # Default SA settings
    sa = {"temp_reduction_factor": 0.90, "number_of_runs": 50, "max_cycles": 50}
    if sa_settings:
        sa.update(sa_settings)

    with open(dpf_path, "w") as f:
        f.write("autodock_parameter_version 4.2       # used by autodock to validate parameter set\n")
        f.write(f"parameter_file {parameter_file}     # parameter library filename\n")
        f.write("outlev 1                             # diagnostic output level\n")
        f.write("intelec                              # calculate internal electrostatics\n")
        f.write("seed pid time                        # seeds for random generator\n")
        f.write(f"ligand_types {ligand_type_str}      # atoms types in ligand\n")
        f.write(f"fld {name_protein}.maps.fld         # grid_data_file\n")

        for lt in ligand_type:
            f.write(f"map {name_protein}.{lt}.map     # atom-specific affinity map\n")

        f.write(f"elecmap {name_protein}.e.map        # electrostatics map\n")
        f.write(f"desolvmap {name_protein}.d.map      # desolvation map\n\n")
        f.write(f"move {name_ligand}.pdbqt            # small molecule\n")

        if random_pos:
            f.write("tran0 random                         # initial coordinates/A or random\n")
            f.write("quaternion0 random                   # initial orientation\n")
            f.write("dihe0 random                         # initial dihedrals (relative) or random\n")

        if ga_dock and not sa_dock:
            f.write("# GA parameters\n")
            f.write(f"ga_pop_size {ga['pop_size']}          # number of individuals in population\n")
            f.write(f"ga_num_evals {ga['num_evals']}        # maximum number of energy evaluations\n")
            f.write(f"ga_num_generations {ga['num_generations']} # maximum number of generations\n")
            f.write(f"ga_elitism {ga['elitism']}            # top individuals surviving to next generation\n")
            f.write(f"ga_mutation_rate {ga['mutation_rate']} # rate of gene mutation\n")
            f.write(f"ga_crossover_rate {ga['crossover_rate']} # rate of crossover\n")
            f.write(f"ga_window_size {ga['window_size']}    # window size for worst individual\n")
            f.write("ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution\n")
            f.write("ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution\n")
            f.write("# Local Search Parameters\n")
            f.write("sw_max_its 300\n")
            f.write("sw_max_succ 4\n")
            f.write("sw_max_fail 4\n")
            f.write("sw_rho 1.0\n")
            f.write("sw_lb_rho 0.01\n")
            f.write("ls_search_freq 0.06\n")
            f.write("# Activate LGA\n")
            f.write("set_ga\n")
            f.write("set_psw1\n")
            f.write(f"ga_run {num_poses}                   # number of hybrid GA-LS runs\n")

        elif sa_dock and not ga_dock:
            f.write("# SA Parameters\n")
            f.write("tstep 2.0\n")
            f.write("linear_schedule\n")
            f.write("rt0 500\n")
            f.write(f"rtrf {sa['temp_reduction_factor']}\n")
            f.write(f"runs {sa['number_of_runs']}\n")
            f.write(f"cycles {sa['max_cycles']}\n")
            f.write("accs 30000\n")
            f.write("rejs 30000\n")
            f.write("select m\n")
            f.write("trnrf 1.0\n")
            f.write("quarf 1.0\n")
            f.write("dihrf 1.0\n")
            f.write("# Activate SA\n")
            f.write(f"simanneal {num_poses}\n")

        f.write("analysis                             # perform a ranked cluster analysis\n")

    return dpf_path


# ===================================================================
# AutoGrid / AutoDock execution
# ===================================================================

def _run_autogrid(log_path: Path, gpf_path: Path, autogrid4_path: str, cwd: Path) -> None:
    """Execute autogrid4."""
    logger.info("Running autogrid4...")
    with open(log_path, "w") as log:
        subprocess.run(
            [autogrid4_path, "-p", str(gpf_path)],
            stdout=log, stderr=subprocess.STDOUT,
            cwd=str(cwd),
        )


def _run_autodock(log_path: Path, dpf_path: Path, autodock4_path: str, cwd: Path) -> None:
    """Execute autodock4."""
    logger.info("Running autodock4...")
    with open(log_path, "w") as log:
        subprocess.run(
            [autodock4_path, "-p", str(dpf_path)],
            stdout=log, stderr=subprocess.STDOUT,
            cwd=str(cwd),
        )


# ===================================================================
# Pose extraction
# ===================================================================

def _write_conformations(dlg_path: Path, output_dir: Path, name_ligand: str) -> list[Path]:
    """Extract docked poses from the DLG file into individual PDBQT files."""
    poses = []
    mol_id = 1

    with open(dlg_path) as f:
        in_block = False
        block_lines: list[str] = []
        for line in f:
            if "DOCKED: ROOT" in line:
                in_block = True
                block_lines = [line]
            elif in_block and "TER" in line:
                in_block = False
                block_lines.append(line)
                pose_path = output_dir / f"{name_ligand}_{mol_id}.pdbqt"
                with open(pose_path, "w") as out:
                    for bl in block_lines:
                        out.write(bl.replace("DOCKED: ", "", 1))
                poses.append(pose_path)
                mol_id += 1
            elif in_block:
                block_lines.append(line)

    logger.info("Extracted %d poses from DLG", len(poses))
    return poses


def _clean_dummy_atoms_from_pdbqt(pose_paths: list[Path]) -> None:
    """Remove lines containing dummy atoms (DD) from PDBQT files."""
    for path in pose_paths:
        lines = path.read_text().splitlines(keepends=True)
        with open(path, "w") as f:
            for line in lines:
                if "DD" not in line:
                    f.write(line)


def _pdbqt_poses_to_xyz(pose_paths: list[Path], name_ligand: str) -> list[Path]:
    """Convert extracted PDBQT poses to XYZ format."""
    xyz_paths = []
    for pdbqt_path in pose_paths:
        xyz_path = pdbqt_path.with_suffix(".xyz")
        atoms = []
        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    parts = line.strip().split()
                    atoms.append(f"{parts[2]:>2} {float(parts[6]):>8.3f} {float(parts[7]):>8.3f} {float(parts[8]):>8.3f}\n")

        with open(xyz_path, "w") as f:
            f.write(f"{len(atoms)}\n\n")
            for a in atoms:
                f.write(a)

        xyz_paths.append(xyz_path)
    return xyz_paths
