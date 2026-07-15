"""Module 6: Docking results analysis.

Extract binding energies, ligand efficiencies, interacting residues,
and (optionally) RMSD values from AutoDock4 output.

Pure Python — no external tool dependencies (except numpy/scipy).
"""

import logging
import re
import subprocess
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)


# ===================================================================
# Public interface
# ===================================================================

def analyze_docking_results(
    dlg_path: Path,
    pose_xyz_paths: list[Path],
    protein_pdb: Path,
    n_heavy_atoms: int,
    num_poses: int | None = None,
    reference_xyz: Path | None = None,
    cutoff: float = 4.0,
) -> dict:
    """Analyze docking results: energies, efficiencies, contacts, RMSD.

    Args:
        dlg_path: Path to the AutoDock4 docking log (.dlg).
        pose_xyz_paths: List of pose XYZ file paths.
        protein_pdb: Path to the cleaned protein PDB.
        n_heavy_atoms: Number of heavy atoms in the metal complex.
        num_poses: Max number of poses to analyze. Defaults to len(pose_xyz_paths).
        reference_xyz: If given, compute RMSD against this reference.
        cutoff: Distance cutoff (Angstrom) for interacting residues.

    Returns:
        Dict with keys:
        - ``binding_energies``: list of floats (kcal/mol).
        - ``binding_efficiencies``: list of floats (kcal/mol per heavy atom).
        - ``interacting_residues``: list of lists of (residue_name, residue_id) tuples.
        - ``rmsd_values``: list of floats (only if reference_xyz is provided).
        - ``rmsd_stats``: dict with mean, std, var (only if reference_xyz is provided).
    """
    if num_poses is None:
        num_poses = len(pose_xyz_paths)

    # 1. Extract binding energies
    energies = extract_binding_energies(dlg_path, num_poses)
    efficiencies = [e / n_heavy_atoms for e in energies] if n_heavy_atoms > 0 else energies

    # 2. Interacting residues
    residues_per_pose = []
    for pose_path in pose_xyz_paths[:num_poses]:
        res = extract_interacting_residues(pose_path, protein_pdb, cutoff=cutoff)
        residues_per_pose.append(res)

    result: dict = {
        "binding_energies": energies,
        "binding_efficiencies": efficiencies,
        "interacting_residues": residues_per_pose,
    }

    # 3. RMSD (optional)
    if reference_xyz is not None:
        rmsds = []
        for pose_path in pose_xyz_paths[:num_poses]:
            rmsd = calculate_rmsd(reference_xyz, pose_path)
            rmsds.append(rmsd)
        result["rmsd_values"] = rmsds
        result["rmsd_stats"] = {
            "mean": float(np.mean(rmsds)),
            "std": float(np.std(rmsds)),
            "var": float(np.var(rmsds)),
        }

    return result


# ===================================================================
# Binding energy extraction
# ===================================================================

def extract_binding_energies(dlg_path: Path, num_poses: int) -> list[float]:
    """Extract estimated free energies of binding from the DLG file."""
    with open(dlg_path) as f:
        content = f.read()

    matches = re.findall(
        r"Estimated Free Energy of Binding\s*=\s*([-\d.]+)\s*kcal/mol",
        content,
    )
    energies = [float(v) for v in matches]
    return energies[:num_poses]


# ===================================================================
# Interacting residues
# ===================================================================

def extract_interacting_residues(
    pose_xyz: Path,
    protein_pdb: Path,
    cutoff: float = 4.0,
) -> list[tuple[str, str]]:
    """Find protein residues within *cutoff* Angstrom of any pose atom.

    Args:
        pose_xyz: XYZ file of a docked pose.
        protein_pdb: PDB file of the protein.
        cutoff: Distance threshold in Angstrom.

    Returns:
        Sorted list of (residue_name, residue_id) tuples.
    """
    # Read pose coordinates
    pose_coords = []
    with open(pose_xyz) as f:
        for _ in range(2):
            next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                pose_coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    if not pose_coords:
        return []

    # Read protein atom coordinates and residue info
    residue_info = []
    protein_coords = []
    with open(protein_pdb) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.strip().split()
                residue_info.append((parts[3], parts[5]))  # (name, id)
                protein_coords.append([
                    float(parts[6]), float(parts[7]), float(parts[8]),
                ])

    if not protein_coords:
        return []

    # Distance matrix
    dist_matrix = cdist(pose_coords, protein_coords, "euclidean")
    within_cutoff = np.any(dist_matrix <= cutoff, axis=0)

    unique_residues: set[tuple[str, str]] = set()
    for i, close in enumerate(within_cutoff):
        if close:
            unique_residues.add(residue_info[i])

    return sorted(unique_residues, key=lambda r: int(r[1]))


# ===================================================================
# RMSD calculation
# ===================================================================

def calculate_rmsd(
    reference_xyz: Path,
    pose_xyz: Path,
    ignore_h: bool = True,
) -> float:
    """Calculate RMSD between two XYZ structures (no rotation/translation).

    Args:
        reference_xyz: Reference structure XYZ file.
        pose_xyz: Docked pose XYZ file.
        ignore_h: If True, skip hydrogen atoms.

    Returns:
        RMSD value in Angstrom.
    """
    ref_coords = _read_xyz_coords(reference_xyz, ignore_h=ignore_h)
    pose_coords = _read_xyz_coords(pose_xyz, ignore_h=ignore_h)

    if len(ref_coords) != len(pose_coords):
        logger.warning(
            "Atom count mismatch: ref=%d, pose=%d. Using min.",
            len(ref_coords), len(pose_coords),
        )
        n = min(len(ref_coords), len(pose_coords))
        ref_coords = ref_coords[:n]
        pose_coords = pose_coords[:n]

    if len(ref_coords) == 0:
        return 0.0

    ref = np.array(ref_coords)
    pose = np.array(pose_coords)
    diff = ref - pose
    return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))


def _read_xyz_coords(xyz_path: Path, ignore_h: bool = True) -> list[list[float]]:
    """Read coordinates from an XYZ file, optionally skipping H atoms."""
    coords = []
    with open(xyz_path) as f:
        for _ in range(2):
            next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            if ignore_h and parts[0] == "H":
                continue
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return coords


# ===================================================================
# Pose format conversion
# ===================================================================

def write_pose_pdb(
    xyz_path: Path,
    pdb_path: Path,
    graph: nx.Graph | None = None,
    atom_index_mapping: dict | None = None,
    residue_name: str = "UNK",
) -> Path:
    """Convert a pose XYZ file to PDB format.

    Args:
        xyz_path: Input XYZ file.
        pdb_path: Output PDB file.
        graph: If provided, write CONECT records from graph edges.
        atom_index_mapping: If provided, use for CONECT records.
        residue_name: Three-letter residue code.

    Returns:
        The pdb_path.
    """
    import networkx as nx

    atoms = []
    with open(xyz_path) as f:
        for _ in range(2):
            next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                atoms.append((parts[0], [float(parts[1]), float(parts[2]), float(parts[3])]))

    with open(pdb_path, "w") as f:
        for i, (el, xyz) in enumerate(atoms, 1):
            atom_type = f"{el}{i - 1}"
            f.write(
                f"HETATM{i:>5} {atom_type:>3}  {residue_name} A   1    "
                f"{xyz[0]:>8.3f}{xyz[1]:>8.3f}{xyz[2]:>8.3f}"
                f"  1.00  0.00          {el:>2}\n"
            )
        if graph is not None and atom_index_mapping is not None:
            for a1, a2 in graph.edges():
                idx1 = atom_index_mapping.get(a1, {}).get("pdbqt_index")
                idx2 = atom_index_mapping.get(a2, {}).get("pdbqt_index")
                if idx1 is not None and idx2 is not None:
                    f.write(f"CONECT {idx1 + 1:>4} {idx2 + 1:>4}\n")
        f.write("ENDMDL\n")

    return pdb_path
