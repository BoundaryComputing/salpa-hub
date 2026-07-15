"""Shared utilities for MetalDock modules.

Graph serialization/deserialization, Lennard-Jones parameter sets,
metal validation, and common helpers.
"""

import json
import logging
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Supported metals and their optimized LJ well-depth parameters
# Format: {METAL_SYMBOL: [e_NA, e_OA, e_SA, e_HD]}
# ---------------------------------------------------------------------------
STANDARD_LJ_PARAMS: dict[str, list[float]] = {
    "V":  [4.696, 6.825, 5.658, 3.984],
    "CR": [6.371, 1.998, 0.144, 3.625],
    "CO": [5.280, 0.050, 6.673, 5.929],
    "NI": [0.630, 2.732, 4.462, 2.820],
    "CU": [4.696, 1.277, 6.791, 1.114],
    "MO": [1.330, 0.014, 0.168, 5.620],
    "RU": [6.936, 2.796, 4.295, 6.357],
    "RH": [5.559, 2.056, 0.573, 5.471],
    "PD": [4.688, 0.845, 5.574, 3.159],
    "RE": [6.738, 0.645, 3.309, 4.502],
    "OS": [5.958, 0.135, 4.102, 6.589],
    "PT": [6.532, 2.020, 6.332, 1.844],
}

# Metals that already have parameters built into AutoDock4
INTERNAL_PARAM_METALS = frozenset({"FE", "ZN", "MN"})

SUPPORTED_METALS = frozenset(STANDARD_LJ_PARAMS.keys()) | INTERNAL_PARAM_METALS

DATA_DIR = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# Metal validation
# ---------------------------------------------------------------------------
def validate_metal_symbol(metal_symbol: str) -> str:
    """Return the upper-cased metal symbol after validation.

    Raises ValueError if the symbol is not supported.
    """
    sym = metal_symbol.strip().upper()
    if sym not in SUPPORTED_METALS:
        raise ValueError(
            f"Unsupported metal symbol '{metal_symbol}'. "
            f"Supported: {sorted(SUPPORTED_METALS)}"
        )
    return sym


def get_lj_params(metal_symbol: str, custom_params: list[float] | None = None) -> list[float] | None:
    """Return [e_NA, e_OA, e_SA, e_HD] for the given metal.

    Returns None for metals with internal AutoDock4 parameters (Fe, Zn, Mn).
    If *custom_params* is provided (length-4 list), it is returned directly.
    """
    sym = validate_metal_symbol(metal_symbol)
    if sym in INTERNAL_PARAM_METALS:
        return None
    if custom_params is not None:
        if len(custom_params) != 4:
            raise ValueError("custom_params must be a length-4 list [e_NA, e_OA, e_SA, e_HD]")
        return list(custom_params)
    return STANDARD_LJ_PARAMS.get(sym, [5.0, 5.0, 5.0, 5.0])


# ---------------------------------------------------------------------------
# XYZ helpers
# ---------------------------------------------------------------------------
def find_metal_in_xyz(xyz_path: str | Path, metal_symbol: str) -> bool:
    """Return True if *metal_symbol* appears in the XYZ file."""
    sym = metal_symbol.strip()
    with open(xyz_path) as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if parts and parts[0] == sym:
                return True
    return False


def count_heavy_atoms(xyz_path: str | Path) -> int:
    """Count non-hydrogen atoms in an XYZ file."""
    count = 0
    with open(xyz_path) as f:
        for line in f.readlines()[2:]:
            parts = line.strip().split()
            if parts and parts[0] != "H":
                count += 1
    return count


# ---------------------------------------------------------------------------
# Graph serialization (NetworkX <-> JSON)
# ---------------------------------------------------------------------------
def _convert_for_json(obj: Any) -> Any:
    """Recursively convert numpy types to native Python for JSON serialization."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, dict):
        return {_convert_for_json(k): _convert_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        converted = [_convert_for_json(i) for i in obj]
        return converted if isinstance(obj, list) else tuple(converted)
    if isinstance(obj, frozenset):
        return sorted(_convert_for_json(i) for i in obj)
    return obj


def graph_to_json(G: nx.Graph) -> str:
    """Serialize a NetworkX molecular graph to a JSON string."""
    data = nx.node_link_data(G)
    data = _convert_for_json(data)
    return json.dumps(data, indent=2)


def graph_from_json(json_str: str) -> nx.Graph:
    """Deserialize a NetworkX molecular graph from a JSON string."""
    data = json.loads(json_str)
    return nx.node_link_graph(data)


def save_graph(G: nx.Graph, path: str | Path) -> Path:
    """Write a molecular graph to a JSON file. Returns the path."""
    path = Path(path)
    path.write_text(graph_to_json(G))
    logger.info("Saved graph to %s", path)
    return path


def load_graph(path: str | Path) -> nx.Graph:
    """Load a molecular graph from a JSON file."""
    path = Path(path)
    return graph_from_json(path.read_text())


# ---------------------------------------------------------------------------
# Parameter file path helper
# ---------------------------------------------------------------------------
def default_parameter_file() -> Path:
    """Return the path to the bundled metal_dock.dat AutoDock4 parameter file."""
    return DATA_DIR / "metal_dock.dat"
