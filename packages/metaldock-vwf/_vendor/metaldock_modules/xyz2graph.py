"""Molecular graph construction from XYZ files.

Adapted from https://github.com/zotko/xyz2graph — builds an adjacency list
from atomic coordinates and covalent radii, then converts to a NetworkX graph.
"""

from pathlib import Path

import networkx as nx
import numpy as np

# ---------------------------------------------------------------------------
# Covalent radii (Angstrom) — used to decide bond connectivity
# ---------------------------------------------------------------------------
ATOMIC_RADII: dict[str, float] = {
    "Ac": 1.88, "Ag": 1.59, "Al": 1.35, "Am": 1.51, "As": 1.21,
    "Au": 1.50, "B":  0.83, "Ba": 1.34, "Be": 0.35, "Bi": 1.54,
    "Br": 1.21, "C":  0.68, "Ca": 0.99, "Cd": 1.69, "Ce": 1.83,
    "Cl": 0.99, "Co": 1.33, "Cr": 1.35, "Cs": 1.67, "Cu": 1.52,
    "D":  0.23, "Dy": 1.75, "Er": 1.73, "Eu": 1.99, "F":  0.64,
    "Fe": 1.34, "Ga": 1.22, "Gd": 1.79, "Ge": 1.17, "H":  0.23,
    "Hf": 1.57, "Hg": 1.70, "Ho": 1.74, "I":  1.40, "In": 1.63,
    "Ir": 1.32, "K":  1.33, "La": 1.87, "Li": 0.68, "Lu": 1.72,
    "Mg": 1.10, "Mn": 1.35, "Mo": 1.47, "N":  0.68, "Na": 0.97,
    "Nb": 1.48, "Nd": 1.81, "Ni": 1.50, "Np": 1.55, "O":  0.68,
    "Os": 1.37, "P":  1.05, "Pa": 1.61, "Pb": 1.54, "Pd": 1.50,
    "Pm": 1.80, "Po": 1.68, "Pr": 1.82, "Pt": 1.50, "Pu": 1.53,
    "Ra": 1.90, "Rb": 1.47, "Re": 1.35, "Rh": 1.45, "Ru": 1.40,
    "S":  1.02, "Sb": 1.46, "Sc": 1.44, "Se": 1.22, "Si": 1.20,
    "Sm": 1.80, "Sn": 1.46, "Sr": 1.12, "Ta": 1.43, "Tb": 1.76,
    "Tc": 1.35, "Te": 1.47, "Th": 1.79, "Ti": 1.47, "Tl": 1.55,
    "Tm": 1.72, "U":  1.58, "V":  1.33, "W":  1.37, "Y":  1.78,
    "Yb": 1.94, "Zn": 1.45, "Zr": 1.56,
}


class MolGraph:
    """Lightweight molecular graph built from XYZ coordinates and covalent radii."""

    __slots__ = [
        "elements", "x", "y", "z",
        "adj_list", "atomic_radii", "bond_lengths", "adj_matrix",
    ]

    def __init__(self) -> None:
        self.elements: list[str] = []
        self.x: list[float] = []
        self.y: list[float] = []
        self.z: list[float] = []
        self.adj_list: dict[int, set[int]] = {}
        self.atomic_radii: list[float] = []
        self.bond_lengths: dict[frozenset, float] = {}
        self.adj_matrix: np.ndarray | None = None

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    def read_xyz(self, file_path: str | Path) -> None:
        """Parse an XYZ file and build the adjacency list."""
        with open(file_path) as f:
            for _ in range(2):          # skip atom-count and comment lines
                next(f)
            for line in f:
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                self.elements.append(parts[0])
                self.x.append(float(parts[1]))
                self.y.append(float(parts[2]))
                self.z.append(float(parts[3]))
        self.atomic_radii = [ATOMIC_RADII[el] for el in self.elements]
        self._generate_adjacency_list()

    def _generate_adjacency_list(self) -> None:
        """Compute adjacency from pairwise distances vs. summed covalent radii."""
        xyz = np.stack((self.x, self.y, self.z), axis=-1)
        diffs = xyz[:, np.newaxis, :] - xyz
        distances = np.sqrt(np.einsum("ijk,ijk->ij", diffs, diffs))

        radii = np.array(self.atomic_radii)
        distance_bond = (radii[:, np.newaxis] + radii) * 1.4

        adj_matrix = np.logical_and(0.1 < distances, distance_bond > distances).astype(int)

        for i, j in zip(*np.nonzero(adj_matrix)):
            self.adj_list.setdefault(i, set()).add(j)
            self.adj_list.setdefault(j, set()).add(i)
            self.bond_lengths[frozenset([i, j])] = round(distance_bond[i, j], 5)

        self.adj_matrix = adj_matrix

    # ------------------------------------------------------------------
    # Iteration helpers
    # ------------------------------------------------------------------
    def edges(self):
        """Yield unique (i, j) edges."""
        seen: set[frozenset] = set()
        for node, neighbours in self.adj_list.items():
            for nb in neighbours:
                edge = frozenset([node, nb])
                if edge not in seen:
                    seen.add(edge)
                    yield node, nb

    def __len__(self) -> int:
        return len(self.elements)

    def __getitem__(self, idx: int):
        return self.elements[idx], (self.x[idx], self.y[idx], self.z[idx])


def to_networkx_graph(mol_graph: MolGraph) -> nx.Graph:
    """Convert a MolGraph to a NetworkX graph.

    Node attributes: ``element`` (str), ``xyz`` (tuple of 3 floats).
    Edge attributes: ``length`` (float — the covalent-radius bond-length threshold).
    """
    G = nx.Graph(mol_graph.adj_list)
    node_attrs = {
        i: {"element": el, "xyz": xyz}
        for i, (el, xyz) in enumerate(mol_graph)
    }
    nx.set_node_attributes(G, node_attrs)
    edge_attrs = {
        edge: {"length": length}
        for edge, length in mol_graph.bond_lengths.items()
    }
    nx.set_edge_attributes(G, edge_attrs)
    return G


def build_molecular_graph(xyz_path: str | Path) -> nx.Graph:
    """One-shot helper: read XYZ → return NetworkX graph."""
    mg = MolGraph()
    mg.read_xyz(xyz_path)
    return to_networkx_graph(mg)
