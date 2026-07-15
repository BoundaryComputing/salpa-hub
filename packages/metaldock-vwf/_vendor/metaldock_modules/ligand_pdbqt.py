"""Module 4: Ligand PDBQT file generation.

Convert an enriched molecular graph (with charges and bond orders)
to AutoDock4 PDBQT format with ROOT/BRANCH/ENDBRANCH structure.

Pure Python — no external tool dependencies.
"""

import logging
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.optimize import differential_evolution

logger = logging.getLogger(__name__)


# ===================================================================
# Public interface
# ===================================================================

def create_ligand_pdbqt(
    graph: nx.Graph,
    metal_symbol: str,
    output_path: Path,
    vacant_site: bool = True,
    max_torsions: int = 32,
    freeze_coordination_sphere: bool = True,
) -> Path:
    """Generate a PDBQT file from an enriched molecular graph.

    Args:
        graph: NetworkX graph with node attrs 'element', 'xyz', 'charge'
               and edge attr 'bond_order'.
        metal_symbol: Symbol of the metal atom (e.g. 'Ru').
        output_path: Where to write the .pdbqt file.
        vacant_site: If True, add a dummy atom at the vacant coordination site.
        max_torsions: Maximum number of active torsions. AutoDock4 limit is 32.
            If the molecule has more rotatable bonds than this, bonds closest
            to the metal center are frozen first.
        freeze_coordination_sphere: If True, freeze bonds between atoms that
            directly coordinate the metal (within 2 bonds of metal). This is
            chemically correct for chelators like DOTA where the coordination
            cage is rigid when the metal is bound.

    Returns:
        The output_path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Find the metal atom
    metal_atom = _find_metal_atom(graph, metal_symbol)

    # Get ligand subgraphs
    ligand_graphs = _obtain_ligand_subgraphs(graph, metal_atom)

    # Find rotatable bonds
    rotatable_bonds = _find_rotatable_bonds(graph, ligand_graphs)

    # Freeze coordination sphere bonds
    if freeze_coordination_sphere:
        frozen = _find_coordination_sphere_bonds(graph, metal_atom)
        if frozen:
            logger.info(
                "Freezing %d bonds in metal coordination sphere", len(frozen)
            )
            rotatable_bonds -= frozen

    # Enforce max_torsions by freezing bonds closest to the metal
    if len(rotatable_bonds) > max_torsions:
        logger.warning(
            "Rotatable bonds (%d) exceed max_torsions (%d). "
            "Freezing %d bonds closest to the metal center.",
            len(rotatable_bonds), max_torsions,
            len(rotatable_bonds) - max_torsions,
        )
        rotatable_bonds = _limit_torsions(
            graph, rotatable_bonds, metal_atom, max_torsions
        )

    logger.info("Active torsions: %d", len(rotatable_bonds))

    # Create branches
    branches, branch_connections = _create_ligand_branches(
        graph, ligand_graphs, metal_atom, rotatable_bonds
    )

    # Build atom index mapping
    atom_index_mapping = _initialize_atom_index_mapping(graph, metal_atom, vacant_site)

    # Add dummy atom if needed
    if vacant_site:
        _add_dummy_atom(graph, metal_atom, atom_index_mapping)

    # Convert elements to AutoDock types
    _convert_elements_to_autodock(graph, atom_index_mapping)

    # Map old indices to PDBQT indices
    _assign_pdbqt_indices(branches, atom_index_mapping, metal_atom, vacant_site)

    # Write
    _write_pdbqt(graph, branches, branch_connections, atom_index_mapping,
                 len(ligand_graphs), output_path)

    logger.info("Created ligand PDBQT → %s", output_path)
    return output_path


# ===================================================================
# Metal atom identification
# ===================================================================

def _find_metal_atom(graph: nx.Graph, metal_symbol: str) -> int:
    for node, data in graph.nodes(data=True):
        if data.get("element") == metal_symbol:
            return node
    raise ValueError(f"Metal atom '{metal_symbol}' not found in graph")


# ===================================================================
# Ligand subgraph extraction
# ===================================================================

def _obtain_ligand_subgraphs(graph: nx.Graph, metal_atom: int) -> list[nx.Graph]:
    """Extract individual ligand subgraphs from the metal complex."""
    metal_neighbors = list(graph.neighbors(metal_atom))
    visited: set[int] = set()
    subgraphs = []

    for atom in metal_neighbors:
        if atom not in visited:
            sg_atoms = _bfs_ligand(graph, atom, visited, metal_atom)
            sg_atoms.add(metal_atom)
            subgraphs.append(graph.subgraph(sg_atoms).copy())

    return subgraphs


def _bfs_ligand(graph: nx.Graph, start: int, visited: set[int], metal_atom: int) -> set[int]:
    """BFS to get all atoms in a ligand (stopping at the metal)."""
    queue = [start]
    atoms: set[int] = set()
    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)
        atoms.add(current)
        for nb in graph.neighbors(current):
            if nb != metal_atom and nb not in visited:
                queue.append(nb)
    return atoms


# ===================================================================
# Coordination sphere freezing and torsion limiting
# ===================================================================

def _find_coordination_sphere_bonds(
    graph: nx.Graph, metal_atom: int,
    max_bond_distance: int = 2,
    spatial_cutoff: float = 5.0,
) -> set[tuple[int, int]]:
    """Identify rotatable bonds within the metal coordination sphere.

    Uses two complementary criteria:

    1. **Graph distance**: Bonds where BOTH endpoints are within
       *max_bond_distance* hops of the metal in the molecular graph.
    2. **Spatial distance**: Bonds where EITHER endpoint is within
       *spatial_cutoff* Angstrom of the metal atom. This catches
       weakly-coordinating groups (e.g. carboxylate O in DOTA) whose
       Mayer bond orders to the metal fell below the threshold but are
       still geometrically part of the coordination cage.

    Args:
        graph: Molecular graph.
        metal_atom: Node index of the metal atom.
        max_bond_distance: Maximum graph hops from metal (criterion 1).
        spatial_cutoff: Maximum Cartesian distance in Angstrom (criterion 2).

    Returns:
        Set of (min, max) bond tuples to freeze.
    """
    # Criterion 1: BFS for graph-distance-based coordination atoms
    coord_atoms: set[int] = {metal_atom}
    frontier = {metal_atom}
    for _ in range(max_bond_distance):
        next_frontier: set[int] = set()
        for atom in frontier:
            for nb in graph.neighbors(atom):
                if nb not in coord_atoms:
                    coord_atoms.add(nb)
                    next_frontier.add(nb)
        frontier = next_frontier

    # Criterion 2: Spatial proximity to the metal
    metal_xyz = np.asarray(graph.nodes[metal_atom].get("xyz", (0, 0, 0)))
    spatial_atoms: set[int] = set()
    for node, data in graph.nodes(data=True):
        xyz = data.get("xyz")
        if xyz is not None:
            dist = float(np.linalg.norm(np.asarray(xyz) - metal_xyz))
            if dist <= spatial_cutoff:
                spatial_atoms.add(node)

    # Union: atoms near the metal by either criterion
    all_coord = coord_atoms | spatial_atoms

    # Freeze bonds where BOTH endpoints are in the coordination sphere
    frozen: set[tuple[int, int]] = set()
    for u, v in graph.edges():
        if u in all_coord and v in all_coord:
            bond = (min(u, v), max(u, v))
            frozen.add(bond)

    if spatial_atoms - coord_atoms:
        logger.info(
            "  Spatial criterion added %d atoms beyond graph-distance criterion",
            len(spatial_atoms - coord_atoms),
        )

    return frozen


def _limit_torsions(
    graph: nx.Graph,
    rotatable_bonds: set[tuple[int, int]],
    metal_atom: int,
    max_torsions: int,
) -> set[tuple[int, int]]:
    """Reduce rotatable bonds to at most *max_torsions*.

    Keeps bonds farthest from the metal center (peripheral linker bonds
    are most important for docking; coordination sphere bonds are rigid).

    Args:
        graph: Molecular graph.
        rotatable_bonds: Current set of rotatable bonds.
        metal_atom: Node index of the metal atom.
        max_torsions: Target maximum.

    Returns:
        Reduced set of rotatable bonds.
    """
    if len(rotatable_bonds) <= max_torsions:
        return rotatable_bonds

    # Compute shortest path distance from metal for each bond
    try:
        distances = nx.single_source_shortest_path_length(graph, metal_atom)
    except nx.NetworkXError:
        distances = {}

    def bond_distance(bond: tuple[int, int]) -> float:
        d1 = distances.get(bond[0], 0)
        d2 = distances.get(bond[1], 0)
        return (d1 + d2) / 2.0

    # Sort bonds by distance from metal (descending = keep farthest first)
    ranked = sorted(rotatable_bonds, key=bond_distance, reverse=True)
    kept = set(ranked[:max_torsions])

    frozen_count = len(rotatable_bonds) - len(kept)
    logger.info(
        "Kept %d torsions (farthest from metal), froze %d (closest to metal)",
        len(kept), frozen_count,
    )
    return kept


# ===================================================================
# Rotatable bond detection
# ===================================================================

def _find_rotatable_bonds(graph: nx.Graph, ligand_graphs: list[nx.Graph]) -> set[tuple[int, int]]:
    """Find all rotatable bonds across all ligand subgraphs."""
    bonds: set[tuple[int, int]] = set()
    for lg in ligand_graphs:
        bonds.update(_generate_proper_dihedrals(graph, lg))
    return {(min(a, b), max(a, b)) for a, b in bonds}


def _generate_proper_dihedrals(graph: nx.Graph, ligand_graph: nx.Graph) -> set[tuple[int, int]]:
    """Find central bonds of proper dihedral angles."""
    proper = set()
    for atom in ligand_graph.nodes():
        dihedrals = _find_dihedrals(ligand_graph, atom, 3)
        for dih in dihedrals:
            a2, a3 = dih[1], dih[2]
            if _is_valid_bond_order(graph, a2, a3) and (a2, a3) not in proper:
                proper.add((a2, a3))
    return proper


def _find_dihedrals(graph: nx.Graph, node: int, depth: int) -> list[list[int]]:
    """Recursively find paths of given depth from a node."""
    if depth == 0:
        return [[node]]
    result = []
    for nb in graph.neighbors(node):
        for path in _find_dihedrals(graph, nb, depth - 1):
            if node not in path:
                result.append([node] + path)
    return result


def _is_valid_bond_order(graph: nx.Graph, a1: int, a2: int) -> bool:
    """Check if the bond order is consistent with a single bond (0.8 < BO < 1.2)."""
    edge_data = graph.get_edge_data(a1, a2)
    if edge_data is None:
        return False
    bo = edge_data.get("bond_order", 0)
    return 0.8 < bo < 1.2


def _is_rotatable(rotatable_bonds: set, a1: int, a2: int) -> bool:
    return (min(a1, a2), max(a1, a2)) in rotatable_bonds


def _is_not_in_ring(graph: nx.Graph, a1: int, a2: int) -> bool:
    """True if the bond is not part of any ring."""
    for cycle in nx.cycle_basis(graph):
        if a1 in cycle and a2 in cycle:
            return False
    return True


def _is_dihedral_bond(graph: nx.Graph, a1: int, a2: int) -> bool:
    """True if both atoms have ≥2 non-H neighbors (proper dihedral center)."""
    n1 = [n for n in graph.neighbors(a1) if graph.nodes[n].get("element") != "H"]
    n2 = [n for n in graph.neighbors(a2) if graph.nodes[n].get("element") != "H"]
    return len(n1) >= 2 and len(n2) >= 2


# ===================================================================
# Branch creation (DFS)
# ===================================================================

def _create_ligand_branches(
    graph: nx.Graph,
    ligand_graphs: list[nx.Graph],
    metal_atom: int,
    rotatable_bonds: set[tuple[int, int]],
) -> tuple[dict, dict]:
    """Build PDBQT branch structure via DFS over each ligand subgraph."""
    branches: dict[str, list[int]] = {"ROOT": [metal_atom]}
    branch_connections: dict[str, dict] = {}
    connection_idx = 0

    for ligand_idx, lg in enumerate(ligand_graphs):
        def dfs(atom: int, current_branch: str, visited: set[int]) -> None:
            nonlocal connection_idx
            visited.add(atom)
            for nb in lg.neighbors(atom):
                if nb in visited:
                    continue
                if (_is_rotatable(rotatable_bonds, atom, nb)
                        and _is_not_in_ring(graph, atom, nb)
                        and _is_dihedral_bond(graph, atom, nb)):
                    new_branch = f"Branch_{len(branches)}_ligand_{ligand_idx}"
                    branches[new_branch] = [nb]
                    branch_connections[f"connection_{connection_idx}"] = {
                        "from_name": current_branch,
                        "to_name": new_branch,
                        "from_atom": atom,
                        "to_atom": nb,
                    }
                    connection_idx += 1
                    dfs(nb, new_branch, visited)
                else:
                    branches[current_branch].append(nb)
                    dfs(nb, current_branch, visited)

        dfs(metal_atom, "ROOT", set())

    return branches, branch_connections


# ===================================================================
# Atom index mapping and AutoDock element conversion
# ===================================================================

def _initialize_atom_index_mapping(
    graph: nx.Graph, metal_atom: int, vacant_site: bool
) -> dict[int, dict]:
    mapping: dict[int, dict] = {}
    if vacant_site:
        mapping[-1] = {"pdbqt_index": None, "autodock_element": None}
    for node in graph.nodes():
        mapping[node] = {"pdbqt_index": None, "autodock_element": None}
    return mapping


def _convert_elements_to_autodock(graph: nx.Graph, mapping: dict) -> None:
    """Map element symbols to AutoDock atom types."""
    cycles = nx.cycle_basis(graph)

    for node in graph.nodes():
        el = graph.nodes[node].get("element", "")

        if el == "H":
            ad_el = "H"
            for nb in graph.neighbors(node):
                if graph.nodes[nb].get("element") in ("O", "N", "S"):
                    ad_el = "HD"
                    break
        elif el == "C":
            ad_el = "C"
            if any(node in cycle for cycle in cycles):
                for nb in graph.neighbors(node):
                    if graph.nodes[nb].get("element") == "H":
                        continue
                    if graph.edges[node, nb].get("bond_order", 0) > 1.10:
                        ad_el = "A"
        elif el == "O":
            ad_el = "OA"
        elif el == "N":
            ad_el = "N" if len(list(graph.neighbors(node))) == 3 else "NA"
        elif el == "S":
            ad_el = "SA"
        else:
            ad_el = el

        if node in mapping:
            mapping[node]["autodock_element"] = ad_el


def _assign_pdbqt_indices(
    branches: dict, mapping: dict, metal_atom: int, vacant_site: bool
) -> None:
    """Assign sequential PDBQT atom indices following ROOT then branches."""
    idx = 0

    if vacant_site:
        # Insert dummy atom right after the metal in ROOT
        root = branches["ROOT"]
        metal_pos = root.index(metal_atom)
        root.insert(metal_pos + 1, -1)

    for atom in branches["ROOT"]:
        mapping[atom]["pdbqt_index"] = idx
        idx += 1

    ligand_branches = {k: v for k, v in branches.items() if k != "ROOT"}
    ligand_branches = sorted(ligand_branches.items(), key=lambda x: int(x[0].split("_")[1]))

    for _, atoms in ligand_branches:
        for atom in atoms:
            mapping[atom]["pdbqt_index"] = idx
            idx += 1


# ===================================================================
# Dummy atom placement (vacant coordination site)
# ===================================================================

def _add_dummy_atom(graph: nx.Graph, metal_atom: int, mapping: dict) -> None:
    """Add a dummy atom (DD) at the vacant coordination site of the metal."""
    metal_neighbors = list(graph.neighbors(metal_atom))
    ligand_positions = []
    cycles = nx.cycle_basis(graph)

    for nb in metal_neighbors:
        in_cycle = False
        for cycle in cycles:
            if nb in cycle:
                ring_nbs = [n for n in metal_neighbors if n in cycle]
                if len(ring_nbs) > 1:
                    positions = [graph.nodes[n]["xyz"] for n in ring_nbs]
                    ligand_positions.append(np.mean(positions, axis=0))
                    in_cycle = True
                    break
        if not in_cycle:
            ligand_positions.append(graph.nodes[nb]["xyz"])

    metal_pos = np.asarray(graph.nodes[metal_atom]["xyz"])
    h_pos = _find_max_distance_position(ligand_positions, metal_pos)

    graph.add_node(-1, element="DD", xyz=tuple(h_pos), charge=0)
    mapping[-1] = {"pdbqt_index": None, "autodock_element": "DD"}


def _spherical_to_cartesian(radius: float, theta: float, phi: float, center) -> np.ndarray:
    sin_theta = np.sin(theta)
    x = radius * sin_theta * np.cos(phi) + center[0]
    y = radius * sin_theta * np.sin(phi) + center[1]
    z = radius * np.cos(theta) + center[2]
    return np.array([x, y, z])


def _find_max_distance_position(ligand_positions: list, metal_pos) -> np.ndarray:
    """Find the point on a unit sphere (centered at metal) maximizing distance to ligands."""
    metal_pos = np.asarray(metal_pos)
    shifted = np.asarray(ligand_positions) - metal_pos

    def objective(params):
        theta, phi = params
        pt = _spherical_to_cartesian(1, theta, phi, [0, 0, 0])
        return -np.sum(np.linalg.norm(pt - shifted, axis=1))

    result = differential_evolution(
        objective,
        bounds=[(0, np.pi), (0, 2 * np.pi)],
        maxiter=5000,
        tol=1e-8,
    )
    theta_opt, phi_opt = result.x
    return _spherical_to_cartesian(1, theta_opt, phi_opt, metal_pos)


# ===================================================================
# PDBQT file writing
# ===================================================================

def _write_pdbqt(
    graph: nx.Graph,
    branches: dict,
    branch_connections: dict,
    mapping: dict,
    n_ligands: int,
    output_path: Path,
) -> None:
    """Write the PDBQT file with ROOT/BRANCH/ENDBRANCH structure."""
    with open(output_path, "w") as f:
        # ROOT
        f.write("ROOT\n")
        for atom in branches["ROOT"]:
            _write_atom_line(f, graph, atom, mapping)
        f.write("ENDROOT\n")

        # BRANCHES per ligand
        for lig_idx in range(n_ligands):
            lig_branches = {
                k: v for k, v in branches.items()
                if f"ligand_{lig_idx}" in k
            }
            lig_branches = sorted(
                lig_branches.items(),
                key=lambda x: int(x[0].split("_")[1]),
            )
            endbranch_lines = []

            for branch_name, branch_atoms in lig_branches:
                conn = next(
                    (v for v in branch_connections.values() if v["to_name"] == branch_name),
                    None,
                )
                if conn:
                    from_idx = mapping[conn["from_atom"]]["pdbqt_index"] + 1
                    to_idx = mapping[conn["to_atom"]]["pdbqt_index"] + 1
                    f.write(f"BRANCH {from_idx} {to_idx}\n")

                for atom in branch_atoms:
                    _write_atom_line(f, graph, atom, mapping)

                if conn:
                    endbranch_lines.append(f"ENDBRANCH {from_idx} {to_idx}\n")

            for line in reversed(endbranch_lines):
                f.write(line)


def _write_atom_line(f, graph: nx.Graph, atom: int, mapping: dict) -> None:
    """Write a single ATOM line in PDBQT format."""
    el = graph.nodes[atom].get("element", "X")
    display_el = "H" if el == "DD" else el
    ad_el = mapping[atom]["autodock_element"] or el
    xyz = graph.nodes[atom]["xyz"]
    charge = float(graph.nodes[atom].get("charge", 0))
    pdbqt_idx = mapping[atom]["pdbqt_index"] + 1

    f.write(
        f"ATOM   {pdbqt_idx:>4} {display_el:<2}    LIG A   1    "
        f"{xyz[0]:>7.3f} {xyz[1]:>7.3f} {xyz[2]:>7.3f}"
        f"  0.00  0.00    {charge:>6.3f} {ad_el:<2}\n"
    )
