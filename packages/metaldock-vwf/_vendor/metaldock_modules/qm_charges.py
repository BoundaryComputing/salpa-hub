"""Module 3: Quantum-mechanical charge calculation.

Run DFT (via ORCA, Gaussian, or ADF) to obtain CM5 charges and
Mayer bond orders, then enrich the molecular graph.

Each engine is implemented as a set of pure functions.
"""

import json
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from .xyz2graph import build_molecular_graph
from .utils import DATA_DIR

logger = logging.getLogger(__name__)

# GFN2-xTB reports Mulliken charges only — the CM5 column exists in the GFN1
# population analysis and nowhere else. Silently handing back Mulliken charges
# labelled CM5 would be the worst possible failure here, so the parametrisation
# is checked rather than trusted.
XTB_CM5_GFN = 1


# ===================================================================
# Public interface
# ===================================================================

def run_qm_and_enrich_graph(
    graph: nx.Graph,
    xyz_path: Path,
    output_dir: Path,
    engine: str = "orca",
    geom_opt: bool = True,
    charge: int = 0,
    spin: float = 0,
    ncpu: int = 1,
    memory: int = 1500,
    # ORCA-specific
    orca_path: str | Path | None = None,
    orcasimpleinput: str = "PBE def2-TZVP CPCM(Water)",
    orcablocks: str = "",
    # Gaussian-specific
    functional: str = "PBE",
    basis_set: str = "def2-TZVP",
    solvent: str = "",
    dispersion: str = "",
    # ADF-specific
    functional_type: str = "GGA",
    relativity: str = "",
    # xTB-specific
    xtb_path: str | Path | None = None,
    xtb_gfn: int = XTB_CM5_GFN,
    xtb_accuracy: float = 1.0,
    xtb_solvent: str = "",
    xtb_etemp: float = 300.0,
) -> dict:
    """Run QM calculation and add charges + bond orders to the graph.

    Args:
        graph: The molecular graph (from ligand_prep).
        xyz_path: Path to the (canonicalized) XYZ file.
        output_dir: Directory for QM output files.
        engine: One of 'orca', 'gaussian', 'adf', 'xtb'. 'xtb' (GFN1) needs no
            user-supplied binary and is the only engine that runs unattended;
            the rest are DFT and want a licensed/downloaded program.
        geom_opt: If True, run geometry optimization; else single-point.
        charge: Total charge of the system.
        spin: Spin (number of unpaired electrons).
        ncpu: Number of CPU cores.
        memory: Memory in MB (Gaussian).
        orca_path: Path to the ORCA binary directory or executable. If None,
            ORCA must be on PATH or ASE_ORCA_COMMAND must be set.
        orcasimpleinput: ORCA simple input line.
        orcablocks: ORCA blocks.
        functional: DFT functional name.
        basis_set: Basis set name.
        solvent: Solvent name for implicit solvation.
        dispersion: Dispersion correction keyword.
        functional_type: ADF functional type (LDA, GGA, etc.).
        relativity: Relativity setting (e.g. 'ZORA') for ADF.
        xtb_path: xtb binary or its directory. None = resolve from PATH.
        xtb_gfn: GFN parametrisation; only 1 emits CM5 charges.
        xtb_accuracy: xtb '--acc' (lower is tighter).
        xtb_solvent: ALPB implicit solvent name, e.g. 'water'. Empty = gas phase.
        xtb_etemp: Electronic temperature in Kelvin.

    Returns:
        Dict with keys:
        - ``graph``: enriched NetworkX graph (with 'charge' and 'bond_order').
        - ``energy``: QM energy (str).
        - ``charges``: dict mapping 'ELEMENT#' → charge value.
        - ``run_type``: 'geom_opt' or 'single_point'.
        - ``output_xyz``: path to the output XYZ (geometry).
    """
    engine = engine.lower()
    run_type = "geom_opt" if geom_opt else "single_point"
    qm_dir = output_dir / "QM" / run_type
    qm_dir.mkdir(parents=True, exist_ok=True)

    if engine == "orca":
        result = _run_orca(
            xyz_path, qm_dir, geom_opt=geom_opt,
            charge=charge, spin=spin, ncpu=ncpu,
            orcasimpleinput=orcasimpleinput, orcablocks=orcablocks,
            orca_path=orca_path,
        )
    elif engine == "gaussian":
        result = _run_gaussian(
            xyz_path, qm_dir, geom_opt=geom_opt,
            charge=charge, spin=spin, ncpu=ncpu, memory=memory,
            functional=functional, basis_set=basis_set,
            solvent=solvent, dispersion=dispersion,
        )
    elif engine == "adf":
        result = _run_adf(
            xyz_path, qm_dir, geom_opt=geom_opt,
            charge=charge, spin=spin,
            basis_set=basis_set, functional=functional,
            functional_type=functional_type, dispersion=dispersion,
            relativity=relativity, solvent=solvent,
        )
    elif engine == "xtb":
        result = _run_xtb(
            xyz_path, qm_dir, geom_opt=geom_opt,
            charge=charge, spin=spin, ncpu=ncpu,
            xtb_path=xtb_path, gfn=xtb_gfn, accuracy=xtb_accuracy,
            solvent=xtb_solvent, electronic_temperature=xtb_etemp,
        )
    else:
        raise ValueError(f"Unknown QM engine: {engine}")

    energy = result["energy"]
    charges = result["charges"]
    output_xyz = result["output_xyz"]
    log_file = result["log_file"]

    # Rebuild graph from optimized geometry if geom_opt
    graph = build_molecular_graph(output_xyz)

    # Add charges
    for node, data in graph.nodes(data=True):
        element = data.get("element", "")
        key = f"{element.upper()}{node + 1}"
        if key in charges:
            graph.nodes[node]["charge"] = charges[key]

    # Add bond orders
    bond_orders = result.get("bond_orders", {})
    for (a1, a2), order in bond_orders.items():
        if graph.has_edge(a1, a2):
            graph[a1][a2]["bond_order"] = order

    # Remove edges without bond orders (spurious adjacency-only bonds)
    if bond_orders:
        for edge in list(graph.edges()):
            sorted_edge = tuple(sorted(edge))
            if sorted_edge not in bond_orders:
                graph.remove_edge(*edge)

    return {
        "graph": graph,
        "energy": energy,
        "charges": charges,
        "run_type": run_type,
        "output_xyz": output_xyz,
    }


# ===================================================================
# ORCA engine
# ===================================================================

def _run_orca(
    xyz_path: Path,
    qm_dir: Path,
    geom_opt: bool,
    charge: int,
    spin: float,
    ncpu: int,
    orcasimpleinput: str,
    orcablocks: str,
    orca_path: str | Path | None = None,
) -> dict:
    """Run ORCA geometry optimization or single-point via ASE.

    Args:
        orca_path: Path to ORCA binary or directory containing it.
            If a directory, the ``orca`` binary is expected inside it.
            Sets ``ASE_ORCA_COMMAND`` so ASE finds the right binary.
    """
    from ase.io import read as ase_read
    from ase.calculators.orca import ORCA, OrcaProfile
    from ase.optimize.lbfgs import LBFGS

    # Resolve ORCA binary path and set up environment
    orca_bin = None
    if orca_path is not None:
        p = Path(orca_path)
        if p.is_dir():
            orca_bin = str(p / "orca")
            orca_dir = str(p)
        else:
            orca_bin = str(p)
            orca_dir = str(p.parent)
        # ORCA sub-executables need the ORCA dir on PATH
        if orca_dir not in os.environ.get("PATH", ""):
            os.environ["PATH"] = orca_dir + os.pathsep + os.environ.get("PATH", "")
        # ORCA shared libraries (lib/ subfolder) need DYLD_LIBRARY_PATH on macOS
        lib_dir = str(Path(orca_dir) / "lib")
        if Path(lib_dir).is_dir():
            current = os.environ.get("DYLD_LIBRARY_PATH", "")
            if lib_dir not in current:
                os.environ["DYLD_LIBRARY_PATH"] = lib_dir + (os.pathsep + current if current else "")
        logger.info("ORCA binary: %s", orca_bin)

    mult = int(2 * (spin * 0.5) + 1)
    label = "geom" if geom_opt else "single_point"
    out_file = qm_dir / f"{label}.out"
    output_xyz = qm_dir / "output.xyz"

    if not out_file.exists():
        mol = ase_read(str(xyz_path))

        # Build ORCA calculator — ASE 3.28+ uses 'directory' for output location
        orca_kwargs = dict(
            charge=charge,
            mult=mult,
            orcasimpleinput=orcasimpleinput,
            orcablocks=f"%pal nprocs {ncpu} end %output Print[P_hirshfeld] 1 end {orcablocks}",
            directory=str(qm_dir),
        )
        if orca_bin:
            profile = OrcaProfile(command=orca_bin)
            orca_kwargs["profile"] = profile

        mol.calc = ORCA(**orca_kwargs)
        if geom_opt:
            opt = LBFGS(mol)
            opt.run(fmax=0.05)
        mol.get_potential_energy()
        mol.write(str(output_xyz))

        # ASE writes output as 'orca.out' in the directory
        ase_out = qm_dir / "orca.out"
        if ase_out.exists() and not out_file.exists():
            ase_out.rename(out_file)

    # Check convergence
    with open(out_file) as f:
        content = f.read()
    if geom_opt and "SUCCESS" not in content:
        raise RuntimeError(f"ORCA geometry optimization did not converge. See {out_file}")
    if not geom_opt and "ORCA TERMINATED NORMALLY" not in content:
        raise RuntimeError(f"ORCA single-point failed. See {out_file}")

    energy = _orca_extract_energy(out_file)
    charges = _orca_extract_charges(out_file, qm_dir)
    bond_orders = _orca_extract_bond_orders(out_file)

    return {
        "energy": energy,
        "charges": charges,
        "bond_orders": bond_orders,
        "output_xyz": output_xyz,
        "log_file": out_file,
    }


def _orca_extract_energy(log_file: Path) -> str:
    with open(log_file) as f:
        for line in f:
            if line.startswith("FINAL"):
                return line.split()[4]
    return ""


def _orca_extract_charges(log_file: Path, qm_dir: Path) -> dict:
    """Extract Hirshfeld charges from ORCA output and convert to CM5."""
    cm5_model_path = DATA_DIR / "cm5pars.json"
    with open(cm5_model_path) as f:
        cm5_model = json.load(f)

    a0_df = pd.DataFrame.from_dict(cm5_model["A0"])
    rd_df = pd.DataFrame.from_dict(cm5_model["radii"])
    pt_df = pd.DataFrame.from_dict(cm5_model["PeriodicTable"])

    data = _orca_parse_logfile(log_file, pt_df, rd_df)
    qcm5 = _hirshfeld_to_cm5(data, a0_df)

    charges = {}
    for idx, row in qcm5.iterrows():
        charges[f"{row['ATOM'].upper()}{idx + 1}"] = row["QCM5"]

    cm5_file = qm_dir / "CM5_charges.csv"
    qcm5.to_csv(cm5_file, index=False, float_format="%6.4f")

    return charges


def _orca_parse_logfile(log_file: Path, pt_df: pd.DataFrame, rd_df: pd.DataFrame) -> pd.DataFrame:
    pt_df["symbol"] = pt_df["symbol"].map(str.strip)
    sym2num = pt_df.set_index(["symbol"])["atomicNumber"].to_dict()
    num2rad = rd_df.set_index(["RAD_NO"])["VALUE"].to_dict()

    xyz_data = []
    charge_data = []
    id_charges = False
    id_coos = False

    with open(log_file) as f:
        for line in f:
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                id_coos = True
            elif "CARTESIAN COORDINATES (A.U.)" in line:
                id_coos = False
            if "HIRSHFELD ANALYSIS" in line:
                id_charges = True
            elif "TIMINGS" in line:
                id_charges = False
            if id_charges:
                charge_data.append(line.strip().split())
            if id_coos:
                xyz_data.append(line.strip().split())

    hir_charges = pd.DataFrame(charge_data[7:-4], columns=["N", "ATOM", "QHir", "Spin"])
    hir_charges[["N", "QHir", "Spin"]] = hir_charges[["N", "QHir", "Spin"]].apply(pd.to_numeric)
    hir_charges = hir_charges[["N", "QHir"]]

    xyzcoos = pd.DataFrame(xyz_data[2:-2], columns=["ATOM", "X", "Y", "Z"])
    xyzcoos[["X", "Y", "Z"]] = xyzcoos[["X", "Y", "Z"]].apply(pd.to_numeric)

    final_data = pd.concat([xyzcoos, hir_charges], axis=1)
    final_data["AtNum"] = [sym2num[s] for s in final_data.ATOM]
    final_data["RAD"] = [num2rad[s] for s in final_data.AtNum]
    return final_data


def _hirshfeld_to_cm5(df: pd.DataFrame, a0_df: pd.DataFrame) -> pd.DataFrame:
    """Convert Hirshfeld charges to CM5 charges."""
    dvals = _get_avals(a0_df)
    alpha = 2.474
    cm5_charges = []

    for i, r in df.iterrows():
        qcm5 = r.QHir
        for j, p in df.iterrows():
            if r.AtNum != p.AtNum:
                dist = np.sqrt((r.X - p.X)**2 + (r.Y - p.Y)**2 + (r.Z - p.Z)**2)
                factor = np.exp(-1.0 * alpha * (dist - r.RAD - p.RAD))
                qcm5 += factor * dvals[r.AtNum - 1, p.AtNum - 1]
        cm5_charges.append(qcm5)

    df["QCM5"] = np.array(cm5_charges)
    return df


def _get_avals(a0_df: pd.DataFrame) -> np.ndarray:
    num2a = a0_df.set_index(["A0_NO"])["VALUE"].to_dict()
    list_keys = list(num2a.keys())
    n = max(list_keys)
    dvals = np.empty([n, n])
    for i in range(n):
        for j in range(n):
            if i != j:
                dvals[i, j] = num2a[i + 1] - num2a[j + 1]
    # Special coefficients for H-C-N-O interactions
    dvals[0, 5] = 0.0502;  dvals[5, 0] = -0.0502
    dvals[0, 6] = 0.1747;  dvals[6, 0] = -0.1747
    dvals[0, 7] = 0.1671;  dvals[7, 0] = -0.1671
    dvals[5, 6] = 0.0556;  dvals[6, 5] = -0.0556
    dvals[5, 7] = 0.0234;  dvals[7, 5] = -0.0234
    dvals[6, 7] = -0.0346; dvals[7, 6] = 0.0346
    return dvals


def _orca_extract_bond_orders(log_file: Path) -> dict:
    """Extract Mayer bond orders from ORCA output."""
    bond_orders = {}
    with open(log_file) as f:
        content = f.read()

    pattern = r"Mayer bond orders larger than \d+\.\d+\n(.*?)\n\n"
    match = re.search(pattern, content, re.DOTALL)
    if match:
        section = match.group(1)
        bond_pattern = r"B\(\s*(\d+)-\w+\s*,\s*(\d+)-\w+\s*\)\s*:\s*(\d+\.\d+)"
        for a1, a2, order in re.findall(bond_pattern, section):
            atoms = tuple(sorted((int(a1), int(a2))))
            bond_orders[atoms] = float(order)

    return bond_orders


# ===================================================================
# Gaussian engine
# ===================================================================

def _run_gaussian(
    xyz_path: Path,
    qm_dir: Path,
    geom_opt: bool,
    charge: int,
    spin: float,
    ncpu: int,
    memory: int,
    functional: str,
    basis_set: str,
    solvent: str,
    dispersion: str,
) -> dict:
    """Run Gaussian geometry optimization or single-point via ASE."""
    from ase.io import read as ase_read
    from ase.calculators.gaussian import Gaussian, GaussianOptimizer
    import itertools

    mult = int(2 * (spin * 0.5) + 1)
    label = "geom_opt" if geom_opt else "single_point"
    log_file = qm_dir / f"{label}.log"
    chk_file = qm_dir / f"{label}.chk"
    output_xyz = qm_dir / "output.xyz"

    if not chk_file.exists():
        mol = ase_read(str(xyz_path))

        calc_kwargs = dict(
            label=str(qm_dir / label),
            nprocshared=ncpu,
            mem=f"{memory}MB",
            chk=f"{label}.chk",
            xc=functional,
            charge=charge,
            mult=mult,
            basis=basis_set,
            pop="Hirshfeld",
            ioplist=["6/80=1"],
        )
        if solvent:
            calc_kwargs["SCRF"] = f"PCM, solvent={solvent}"
        if dispersion:
            calc_kwargs["EmpiricalDispersion"] = dispersion

        calc = Gaussian(**calc_kwargs)

        if geom_opt:
            opt = GaussianOptimizer(mol, calc)
            opt.run(fmax="tight")
        else:
            mol.calc = calc
            mol.get_potential_energy()
        mol.write(str(output_xyz))

    # Check convergence
    with open(log_file) as f:
        content = f.read()
    if geom_opt and "Optimization completed." not in content:
        raise RuntimeError(f"Gaussian optimization did not converge. See {log_file}")
    if not geom_opt and "SCF Done" not in content:
        raise RuntimeError(f"Gaussian single-point failed. See {log_file}")

    energy = _gaussian_extract_energy(log_file)

    # Extract CM5 charges
    mol_ase = ase_read(str(output_xyz))
    n_atoms = len(mol_ase.positions)
    charges = _gaussian_extract_cm5(log_file, n_atoms)

    # Bond orders
    bond_orders = _gaussian_extract_bond_orders(log_file, n_atoms)

    return {
        "energy": energy,
        "charges": charges,
        "bond_orders": bond_orders,
        "output_xyz": output_xyz,
        "log_file": log_file,
    }


def _gaussian_extract_energy(log_file: Path) -> str:
    with open(log_file) as f:
        for line in f:
            if line.startswith(" SCF Done:"):
                return line.split()[4]
    return ""


def _gaussian_extract_cm5(log_file: Path, n_atoms: int) -> dict:
    """Extract CM5 charges from Gaussian Hirshfeld output."""
    import itertools

    charges = {}
    with open(log_file) as f:
        for line in f:
            if " Hirshfeld charges, spin densities, dipoles, and CM5 charges" in line:
                fin_lines = list(itertools.islice(f, n_atoms + 1))
                fin_lines = [l.strip().split() for l in fin_lines]
                atom_id = 0
                for parts in fin_lines[1:]:
                    atom_id += 1
                    charges[f"{parts[1].upper()}{atom_id}"] = float(parts[7])
                break
    return charges


def _gaussian_extract_bond_orders(log_file: Path, n_atoms: int) -> dict:
    """Extract Mayer bond order matrix from Gaussian output."""
    with open(log_file) as f:
        lines = f.readlines()

    extract_lines = []
    recording = False
    for line in lines:
        if " Atomic Valencies and Mayer Atomic Bond Orders:" in line:
            recording = True
            continue
        elif " Lowdin Atomic Charges:" in line:
            recording = False
            break
        if recording:
            extract_lines.append(line)

    clean_lines = [l.strip() for l in extract_lines if l.strip()]
    if not clean_lines:
        return {}

    matrix = np.zeros((n_atoms, n_atoms))
    n_chunks = n_atoms // 6
    rest_columns = n_atoms % 6

    for chunk in range(n_chunks + 1):
        start_line = 1 + chunk * (n_atoms + 1)
        columns = 6 if chunk < n_chunks else rest_columns
        if columns == 0:
            continue

        for col in range(columns):
            for row in range(n_atoms):
                if start_line + row >= len(clean_lines):
                    break
                parts = clean_lines[start_line + row].split()
                value = float(parts[col + 2])
                matrix[row][chunk * 6 + col] = value

    bond_orders = {}
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if matrix[i][j] > 0.2:
                bond_orders[(i, j)] = matrix[i][j]

    return bond_orders


# ===================================================================
# ADF engine
# ===================================================================

def _run_adf(
    xyz_path: Path,
    qm_dir: Path,
    geom_opt: bool,
    charge: int,
    spin: float,
    basis_set: str,
    functional: str,
    functional_type: str,
    dispersion: str,
    relativity: str,
    solvent: str,
) -> dict:
    """Run ADF geometry optimization or single-point via PLAMS."""
    import scm.plams as scm

    run_type = "geom_opt" if geom_opt else "single_point"
    output_xyz = qm_dir / "output.xyz"
    plams_dir = qm_dir / "plamsjob"
    log_file = plams_dir / "ams.log"

    if not plams_dir.is_dir():
        scm.init(folder=str(qm_dir))
        m = scm.Molecule(str(xyz_path))
        m.properties.charge = str(charge)

        s = scm.Settings()
        s.input.ams.Task = "GeometryOptimization" if geom_opt else "SinglePoint"
        s.input.ams.properties.bondorders = "Yes"
        s.input.adf.bondorders.TypeForAMS = "Mayer"
        s.input.adf.scf.iterations = "500"
        s.input.adf.AtomicChargesTypeForAMS = "CM5"
        s.input.adf.basis.type = basis_set.upper()
        s.input.adf.basis.core = "None"

        # Functional
        ft = functional_type.lower()
        func_map = {"lda": "lda", "gga": "gga", "metagga": "metagga",
                     "hybrid": "hybrid", "metahybrid": "metahybrid"}
        if ft in func_map:
            setattr(s.input.adf.xc, func_map[ft], functional.upper())

        if dispersion:
            s.input.adf.xc.dispersion = dispersion.upper()
        if spin != 0:
            s.input.adf.unrestricted = "yes"
            s.input.adf.spinpolarization = str(spin)
        if relativity == "ZORA":
            s.input.adf.relativity.formalism = "ZORA"
            s.input.adf.relativity.level = "Scalar"
        if solvent:
            s.input.adf.Solvation.Solv = f"Name={solvent}"

        j = scm.AMSJob(molecule=m, settings=s)
        j.run()
        j.results.get_main_molecule().write(str(output_xyz), "xyz")
        scm.finish()

    # Check convergence
    with open(log_file) as f:
        content = f.read()
    if geom_opt and "Geometry optimization converged" not in content:
        raise RuntimeError(f"ADF geometry optimization did not converge. See {log_file}")
    if not geom_opt and "NORMAL TERMINATION" not in content:
        raise RuntimeError(f"ADF single-point did not converge. See {log_file}")

    energy = _adf_extract_energy(log_file)
    charges = _adf_extract_charges(plams_dir)
    bond_orders = _adf_extract_bond_orders(plams_dir / "plamsjob.out")

    return {
        "energy": energy,
        "charges": charges,
        "bond_orders": bond_orders,
        "output_xyz": output_xyz,
        "log_file": log_file,
    }


def _adf_extract_energy(log_file: Path) -> str:
    with open(log_file) as f:
        for line in f:
            if "kcal/mol" in line:
                return line.split()[4]
    return ""


def _adf_extract_charges(plams_dir: Path) -> dict:
    """Extract CM5 charges from ADF output."""
    cm5_file = plams_dir / "CM5_charges"
    # The CM5 charges file is generated by amsreport
    charges = {}
    if cm5_file.exists():
        with open(cm5_file) as f:
            for idx, line in enumerate(f):
                if idx == 0:
                    continue
                parts = line.strip().split()
                if not parts:
                    continue
                key = parts[0].replace("(", "").replace(")", "").upper()
                charges[key] = float(parts[-1])
    return charges


def _adf_extract_bond_orders(out_file: Path) -> dict:
    """Extract Mayer bond orders from ADF output."""
    bond_orders = {}
    if not out_file.exists():
        return bond_orders

    with open(out_file) as f:
        content = f.read()

    pattern = (
        r"Description: Mayer bond orders\n"
        r"Only bonds with bond orders > 0\.200 are printed\.\n\n"
        r" Index  Atom    Index  Atom    BondOrder\n(.*?)\n\n"
    )
    match = re.search(pattern, content, re.DOTALL)
    if match:
        section = match.group(1)
        bond_pattern = r"\s*(\d+)\s+\w+\s+(\d+)\s+\w+\s+([\d.]+)"
        for a1, a2, order in re.findall(bond_pattern, section):
            atoms = tuple(sorted((int(a1) - 1, int(a2) - 1)))
            bond_orders[atoms] = float(order)

    return bond_orders


# ===================================================================
# xTB engine (semi-empirical, open source)
# ===================================================================
#
# GFN1-xTB is the only engine here that needs no user-supplied binary: it is
# LGPL-3.0, ships on conda-forge for every platform, and prints CM5 charges
# and Wiberg bond orders directly — the two quantities this pipeline consumes.
# That makes it the engine a demo workflow can actually run out of the box.
#
# It is semi-empirical, not DFT. On the 1JZI Re reference case it reproduces
# the ORCA metal charge to 0.05 e (Re +0.749 vs +0.704) with a 0.065 e mean
# absolute deviation over all 29 atoms, in under a second rather than minutes.
# Use it to get a pipeline moving and to screen; use ORCA for numbers you
# intend to publish.

def _run_xtb(
    xyz_path: Path,
    qm_dir: Path,
    geom_opt: bool,
    charge: int,
    spin: float,
    ncpu: int,
    xtb_path: str | Path | None = None,
    gfn: int = XTB_CM5_GFN,
    accuracy: float = 1.0,
    solvent: str = "",
    electronic_temperature: float = 300.0,
) -> dict:
    """Run a GFN1-xTB single-point or optimization and read CM5 charges + WBOs.

    Args:
        xyz_path: Canonicalized input geometry.
        qm_dir: Scratch directory. xtb writes ``charges``, ``wbo`` and
            ``xtbopt.xyz`` into its working directory, so it is run with
            ``cwd=qm_dir`` and never pollutes the caller's directory.
        geom_opt: Optimize the geometry (``--opt``) instead of a single point.
        charge: Net molecular charge.
        spin: Unpaired electrons (mapped to ``--uhf``).
        ncpu: Threads for ``--parallel``.
        xtb_path: xtb binary or its containing directory. Falls back to PATH.
        gfn: GFN parametrisation. Must be 1 — see ``XTB_CM5_GFN``.
        accuracy: ``--acc``; lower is tighter.
        solvent: ALPB implicit solvent name (e.g. ``water``). Empty = gas phase.
        electronic_temperature: ``--etemp`` in Kelvin.

    Returns:
        The engine dict consumed by :func:`run_qm_and_enrich_graph`.
    """
    if int(gfn) != XTB_CM5_GFN:
        raise ValueError(
            f"xTB engine requires --gfn {XTB_CM5_GFN} (got {gfn}). Only GFN1-xTB "
            "prints CM5 charges; GFN2-xTB reports Mulliken charges only, which "
            "are not interchangeable with the CM5 charges this pipeline expects."
        )

    xtb_bin = _resolve_xtb_binary(xtb_path)
    qm_dir.mkdir(parents=True, exist_ok=True)

    label = "geom" if geom_opt else "single_point"
    out_file = qm_dir / f"{label}.out"
    output_xyz = qm_dir / "output.xyz"

    if not out_file.exists():
        cmd = [
            xtb_bin,
            str(Path(xyz_path).resolve()),
            "--gfn", str(int(gfn)),
            "--chrg", str(int(charge)),
            "--uhf", str(int(round(spin))),
            "--acc", str(float(accuracy)),
            "--parallel", str(max(1, int(ncpu))),
            "--etemp", str(float(electronic_temperature)),
        ]
        cmd.append("--opt" if geom_opt else "--sp")
        if solvent:
            cmd += ["--alpb", solvent]

        logger.info("Running xtb: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd, cwd=str(qm_dir), capture_output=True, text=True,
            # xtb parallelises through OpenMP; --parallel alone does not cap it.
            env={**os.environ, "OMP_NUM_THREADS": str(max(1, int(ncpu)))},
        )
        # xtb writes its report to stdout; keep stderr too so a failure is legible.
        out_file.write_text(proc.stdout + ("\n" + proc.stderr if proc.stderr else ""))
        if proc.returncode != 0:
            raise RuntimeError(
                f"xtb exited {proc.returncode}. See {out_file}\n"
                f"{proc.stderr.strip()[:2000]}"
            )

    content = out_file.read_text()
    if "normal termination of xtb" not in content:
        raise RuntimeError(f"xtb did not terminate normally. See {out_file}")

    # Geometry: --opt writes xtbopt.xyz; a single point leaves the input geometry.
    if geom_opt:
        opt_xyz = qm_dir / "xtbopt.xyz"
        if not opt_xyz.exists():
            raise RuntimeError(
                f"xtb optimization produced no xtbopt.xyz. See {out_file}"
            )
        shutil.copyfile(opt_xyz, output_xyz)
    elif not output_xyz.exists():
        shutil.copyfile(Path(xyz_path), output_xyz)

    charges = _xtb_extract_charges(out_file, qm_dir)
    bond_orders = _xtb_extract_bond_orders(qm_dir / "wbo")

    return {
        "energy": _xtb_extract_energy(out_file),
        "charges": charges,
        "bond_orders": bond_orders,
        "output_xyz": output_xyz,
        "log_file": out_file,
    }


def _resolve_xtb_binary(xtb_path: str | Path | None) -> str:
    """Return a runnable xtb binary, or explain precisely what is missing."""
    if xtb_path:
        p = Path(xtb_path)
        candidate = p / "xtb" if p.is_dir() else p
        if not candidate.is_file():
            raise FileNotFoundError(f"No xtb binary at {candidate}")
        return str(candidate)

    found = shutil.which("xtb")
    if not found:
        raise FileNotFoundError(
            "xtb was not found on PATH. Install it into the node environment "
            "(conda-forge: `xtb`) or pass its path explicitly."
        )
    return found


def _xtb_extract_energy(log_file: Path) -> str:
    """Total energy in Hartree, as a string (matching the other engines)."""
    match = re.search(
        r"\|\s*TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", log_file.read_text()
    )
    return match.group(1) if match else ""


def _xtb_extract_charges(log_file: Path, qm_dir: Path) -> dict:
    """Read the CM5 column of the GFN1 population analysis.

    The table looks like::

          Mulliken/CM5 charges         n(s)   n(p)   n(d)
             1O    -0.27212 -0.32021   1.724  4.549  0.000
             3Re    0.29185  0.74910   0.596  0.610  5.502

    Index and element are glued together, and the first numeric column is
    Mulliken — the CM5 charge this pipeline wants is the second.
    """
    lines = log_file.read_text().splitlines()
    try:
        start = next(
            i for i, line in enumerate(lines) if "Mulliken/CM5 charges" in line
        )
    except StopIteration:
        raise RuntimeError(
            f"No Mulliken/CM5 charge table in {log_file}. This is what a "
            "GFN2-xTB run looks like; CM5 charges require GFN1-xTB."
        )

    row = re.compile(r"^\s*(\d+)([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
    charges: dict[str, float] = {}
    rows = []
    for line in lines[start + 1:]:
        match = row.match(line)
        if not match:
            break  # the table ends at the first non-matching line
        index, element, _mulliken, cm5 = match.groups()
        charges[f"{element.upper()}{int(index)}"] = float(cm5)
        rows.append((int(index), element, float(_mulliken), float(cm5)))

    if not charges:
        raise RuntimeError(f"Could not parse any CM5 charges from {log_file}")

    # Mirror the ORCA engine, which leaves a readable CM5 table beside its output.
    cm5_file = qm_dir / "CM5_charges.csv"
    pd.DataFrame(rows, columns=["N", "ATOM", "QMulliken", "QCM5"]).to_csv(
        cm5_file, index=False, float_format="%6.4f"
    )
    return charges


def _xtb_extract_bond_orders(wbo_file: Path) -> dict:
    """Read Wiberg bond orders from xtb's ``wbo`` file.

    Format is ``<atom1> <atom2> <order>`` with **1-based** atom indices; the
    molecular graph is 0-based, so they are shifted here — the same conversion
    the ADF reader does, and the opposite of ORCA, which already counts from 0.
    """
    bond_orders: dict[tuple[int, int], float] = {}
    if not wbo_file.exists():
        logger.warning("xtb wrote no wbo file at %s; bond orders unavailable", wbo_file)
        return bond_orders

    for line in wbo_file.read_text().splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            a1, a2, order = int(parts[0]), int(parts[1]), float(parts[2])
        except ValueError:
            continue
        bond_orders[tuple(sorted((a1 - 1, a2 - 1)))] = order

    return bond_orders
