"""Module 6: Docking results analysis.

Extract binding energies, ligand efficiencies, interacting residues,
and (optionally) RMSD values from AutoDock4 output.

Pure Python — no external tool dependencies (except numpy/scipy/networkx).
"""

import logging
import re
import subprocess
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np
from networkx.algorithms.isomorphism.isomorphvf2 import GraphMatcher
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from .xyz2graph import ATOMIC_RADII

logger = logging.getLogger(__name__)

# The factor xyz2graph applies to summed covalent radii when ligand prep builds
# the molecular graph. The RMSD perceives bonds by the same rule, so the
# molecule it matches is the molecule the pipeline docked.
_BOND_TOLERANCE = 1.4


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
        Dict with keys, every per-pose list in the order AutoDock ran the poses,
        which is not the order of their scores:
        - ``binding_energies``: list of floats (kcal/mol).
        - ``binding_efficiencies``: list of floats (kcal/mol per heavy atom).
        - ``interacting_residues``: list of lists of (residue_name, residue_id) tuples.
        - ``best_pose_index``: 0-based index of the lowest binding energy, the
          first of equals; pose file ``_{best_pose_index + 1}``. None if no
          analysed pose has an energy.
        - ``rmsd_values``: list of floats (only if reference_xyz is provided);
          atoms are paired by chemistry, see ``calculate_rmsd``.
        - ``rmsd_stats``: dict with mean, std, var (only if reference_xyz is provided).
    """
    if num_poses is None:
        num_poses = len(pose_xyz_paths)

    # 1. Extract binding energies
    energies = extract_binding_energies(dlg_path, num_poses)
    efficiencies = (
        [e / n_heavy_atoms for e in energies] if n_heavy_atoms > 0 else energies
    )

    # 2. Interacting residues
    residues_per_pose = []
    for pose_path in pose_xyz_paths[:num_poses]:
        res = extract_interacting_residues(pose_path, protein_pdb, cutoff=cutoff)
        residues_per_pose.append(res)

    # The energy pattern also matches the log's cluster summary, so choose only
    # among energies that belong to a pose that was analysed.
    scored = energies[: len(residues_per_pose)]
    result: dict = {
        "binding_energies": energies,
        "binding_efficiencies": efficiencies,
        "interacting_residues": residues_per_pose,
        "best_pose_index": (
            min(range(len(scored)), key=scored.__getitem__) if scored else None
        ),
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
                protein_coords.append(
                    [
                        float(parts[6]),
                        float(parts[7]),
                        float(parts[8]),
                    ]
                )

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
    """RMSD of a docked pose from a reference geometry of the same molecule.

    The pose is scored where it lies, with no rotation or translation: a docked
    pose is right only if it sits in the right place in the receptor.

    Atoms are paired by chemistry, not by their line in the file. AutoDock Run
    writes poses in the ligand PDBQT's torsion-tree order, which is not the
    order of an input XYZ, so pairing by line compares unrelated atoms. Both
    structures become bond graphs, and of the pairings that keep every element
    and every bond, the one with the lowest RMSD is used. Taking the lowest
    makes the value symmetry-corrected: a ring flipped onto itself has not moved.

    Pairing each atom with the nearest atom of its element instead, as upstream
    MetalDock does, allows pairings no bond pattern does, and reports a
    misplaced pose as closer to the reference than it is.

    Args:
        reference_xyz: Reference structure XYZ file.
        pose_xyz: Docked pose XYZ file.
        ignore_h: If True, skip hydrogen atoms.

    Returns:
        RMSD value in Angstrom.

    Raises:
        ValueError: If the two files do not hold the same molecule: other
            elements, another atom count, or the same atoms bonded otherwise.
        RuntimeError: If the molecule has so many symmetric equivalents that
            the search gives up; no RMSD is better than one that may be wrong.
    """
    ref_elements, ref_coords = _read_xyz_atoms(reference_xyz, ignore_h=ignore_h)
    pose_elements, pose_coords = _read_xyz_atoms(pose_xyz, ignore_h=ignore_h)
    names = f"{Path(pose_xyz).name} and {Path(reference_xyz).name}"

    if Counter(ref_elements) != Counter(pose_elements):
        raise ValueError(
            f"{names} do not hold the same molecule: "
            f"{_formula(pose_elements)} against {_formula(ref_elements)}."
        )
    if not ref_elements:
        raise ValueError(f"{names} hold no atoms to compare.")

    ref_graph = _bond_graph(ref_elements, ref_coords)
    pose_graph = _bond_graph(pose_elements, pose_coords)
    ref_classes, pose_classes = _symmetry_classes(ref_graph, pose_graph)
    bonded_differently = ValueError(
        f"{names} do not hold the same molecule: the same atoms, bonded differently."
    )
    if Counter(ref_classes) != Counter(pose_classes):
        raise bonded_differently

    squared = cdist(ref_coords, pose_coords, "sqeuclidean")
    matcher = _ClosestPairing(ref_graph, pose_graph, squared, ref_classes, pose_classes)
    try:
        for pairing in matcher.isomorphisms_iter():
            total = sum(squared[r, p] for r, p in pairing.items())
            matcher.best = min(matcher.best, total)
    except _SearchLimitReached:
        raise RuntimeError(
            f"Gave up pairing the atoms of {names} after {_SEARCH_LIMIT:,} "
            "candidate pairs: the molecule has too many symmetric equivalents "
            "to search. No RMSD is reported rather than one that may be too high."
        ) from None

    if not np.isfinite(matcher.best):
        raise bonded_differently
    return float(np.sqrt(matcher.best / len(ref_elements)))


# The most extensions one pairing search may weigh. The searches this package
# meets weigh a few hundred; a search that needs this many has met a symmetry
# the bound below cannot cut, and would otherwise run for hours.
_SEARCH_LIMIT = 20_000


class _SearchLimitReached(Exception):
    pass


class _ClosestPairing(GraphMatcher):
    """Search the bond-preserving pairings of two structures for the closest.

    VF2 builds a pairing one atom pair at a time and asks
    ``semantic_feasibility`` whether each extension may stand. Here that also
    rejects an extension that cannot beat the best complete pairing found so
    far. Its bound is the squared deviation of the pairs made so far plus the
    cheapest assignment of the atoms still unpaired, where an atom may go only
    to an atom of its symmetry class that lies as many bonds from each paired
    atom as it does. An isomorphism keeps every such distance, so no completion
    can undercut the bound, and an extension that leaves some atom nowhere to
    go is abandoned at once.

    Extensions are tried in order of their bound, so the first pairing found
    follows the cheapest assignment, and most symmetric alternatives are cut
    off before they are built. Nothing that could be better is ever cut, so
    the result is the exact minimum.
    """

    def __init__(
        self,
        reference: nx.Graph,
        pose: nx.Graph,
        squared: np.ndarray,
        ref_classes: list[int],
        pose_classes: list[int],
    ):
        super().__init__(reference, pose)
        self.squared = squared
        self.ref_classes = np.array(ref_classes)
        self.pose_classes = np.array(pose_classes)
        self.ref_hops = _bond_distances(reference)
        self.pose_hops = _bond_distances(pose)
        self.best = np.inf
        self.weighed = 0
        self.bounds: dict[tuple[int, int], float] = {}

    def candidate_pairs_iter(self):
        pairs = [
            (r, p)
            for r, p in super().candidate_pairs_iter()
            if self.ref_classes[r] == self.pose_classes[p]
            and self.syntactic_feasibility(r, p)
        ]
        for pair in pairs:
            self.bounds[pair] = self._bound(*pair)
        pairs.sort(key=self.bounds.__getitem__)
        return iter(pairs)

    def semantic_feasibility(self, ref_atom, pose_atom):
        return self.bounds[ref_atom, pose_atom] < self.best

    def _bound(self, ref_atom: int, pose_atom: int) -> float:
        """Least squared deviation of any pairing that extends this one."""
        self.weighed += 1
        if self.weighed > _SEARCH_LIMIT:
            raise _SearchLimitReached
        paired_ref = [*self.core_1, ref_atom]
        paired_pose = [*self.core_1.values(), pose_atom]
        total = self.squared[paired_ref, paired_pose].sum()
        free_ref = np.setdiff1d(np.arange(len(self.squared)), paired_ref)
        free_pose = np.setdiff1d(np.arange(len(self.squared)), paired_pose)
        for cls in np.unique(self.ref_classes[free_ref]):
            refs = free_ref[self.ref_classes[free_ref] == cls]
            poses = free_pose[self.pose_classes[free_pose] == cls]
            consistent = (
                self.ref_hops[np.ix_(refs, paired_ref)][:, None, :]
                == self.pose_hops[np.ix_(poses, paired_pose)][None, :, :]
            ).all(axis=2)
            cost = np.where(consistent, self.squared[np.ix_(refs, poses)], np.inf)
            try:
                rows, cols = linear_sum_assignment(cost)
            except ValueError:  # some atom has nowhere left to go
                return np.inf
            total += cost[rows, cols].sum()
        return float(total)


def _bond_distances(graph: nx.Graph) -> np.ndarray:
    """Bonds on the shortest path between each two atoms; -1 if none joins them."""
    hops = np.full((len(graph), len(graph)), -1, dtype=int)
    for source, lengths in nx.all_pairs_shortest_path_length(graph):
        for target, length in lengths.items():
            hops[source, target] = length
    return hops


def _symmetry_classes(first: nx.Graph, second: nx.Graph) -> tuple[list[int], list[int]]:
    """Colour refinement, run on both graphs with one palette.

    Each atom starts with its element as its colour; each round recolours it
    by its colour and its neighbours' colours, until no class splits further.
    An isomorphism can only pair atoms of the same final colour, so the colours
    of the two graphs must also come out in the same numbers.
    """
    graphs = (first, second)
    colours = [[g.nodes[n]["element"] for n in range(len(g))] for g in graphs]
    classes = 0
    while True:
        keys = [
            [(c[n], tuple(sorted(c[m] for m in g[n]))) for n in range(len(g))]
            for g, c in zip(graphs, colours)
        ]
        palette = {
            key: i
            for i, key in enumerate(sorted(set(keys[0]) | set(keys[1]), key=repr))
        }
        colours = [[palette[key] for key in graph_keys] for graph_keys in keys]
        if len(palette) == classes:
            return colours[0], colours[1]
        classes = len(palette)


def _bond_graph(elements: list[str], coords: np.ndarray) -> nx.Graph:
    """Bond graph of one structure, by the covalent-radius rule of ligand prep."""
    unknown = sorted({el for el in elements if el not in ATOMIC_RADII})
    if unknown:
        raise ValueError(
            f"No covalent radius for {', '.join(unknown)}, so its bonds cannot be told."
        )
    radii = np.array([ATOMIC_RADII[el] for el in elements])
    distance = cdist(coords, coords)
    limit = _BOND_TOLERANCE * (radii[:, None] + radii[None, :])
    bonded = np.triu((distance > 0.1) & (distance < limit), k=1)
    graph = nx.Graph()
    graph.add_nodes_from((i, {"element": el}) for i, el in enumerate(elements))
    graph.add_edges_from(zip(*(idx.tolist() for idx in np.nonzero(bonded))))
    return graph


def _formula(elements: list[str]) -> str:
    counts = Counter(elements)
    return " ".join(f"{el}{n if n > 1 else ''}" for el, n in sorted(counts.items()))


def _read_xyz_atoms(
    xyz_path: Path, ignore_h: bool = True
) -> tuple[list[str], np.ndarray]:
    """Read elements and coordinates from an XYZ file, optionally skipping H atoms."""
    elements, coords = [], []
    with open(xyz_path) as f:
        for _ in range(2):
            next(f)
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            element = parts[0].capitalize()
            if ignore_h and element == "H":
                continue
            elements.append(element)
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elements, np.array(coords, dtype=float).reshape(-1, 3)


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
                atoms.append(
                    (parts[0], [float(parts[1]), float(parts[2]), float(parts[3])])
                )

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
