"""MetalDock Ligand Prep — wraps metaldock_modules.ligand_prep.

Canonicalize a metal-complex XYZ with OpenBabel and build its molecular graph.
The graph is saved as JSON (utils.save_graph) and its path forwarded as
``graph_json`` so the QM-charges and PDBQT nodes can reload it. Also records the
metal symbol and heavy-atom count for downstream nodes.
"""

import os
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FolderParameter, StringParameter,
)
from bocoflow_core.stream_logger import stream_log

_NODE_DIR = os.path.dirname(os.path.abspath(__file__))


def _ensure_metaldock_modules():
    try:
        import metaldock_modules  # noqa: F401  installed in env / already on path
        return
    except ImportError:
        pass
    candidates = [
        os.environ.get("METALDOCK_SRC"),
        os.path.join(_NODE_DIR, "scripts"),
        os.path.join(_NODE_DIR, "..", "_vendor"),   # bundled with the package (release)
        os.path.abspath(os.path.join(_NODE_DIR, "..", "..", "..", "src")),
    ]
    for c in candidates:
        if c and os.path.isdir(os.path.join(c, "metaldock_modules")):
            if c not in sys.path:
                sys.path.insert(0, c)
            return c
    raise NodeException(
        "setup",
        "Cannot locate the metaldock_modules package. Set METALDOCK_SRC.",
    )


def _merge_predecessors(predecessor_data):
    """Flatten all upstream node `data` dicts into one carry-forward dict."""
    merged = {}
    for pred in (predecessor_data or []):
        if not isinstance(pred, dict):
            continue
        scope = pred.get("data") if isinstance(pred.get("data"), dict) else pred
        for k, v in scope.items():
            if k not in ("success", "message", "metadata", "files"):
                merged[k] = v
    return merged


class MdockLigandPrep(Node):
    """Canonicalize ligand XYZ (OpenBabel) → build molecular graph (JSON)."""

    category = "Metal Docking"
    tags = ["metaldock", "ligand-prep", "openbabel", "molecular-graph", "metal-complex"]

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name", default="complex",
            docstring="Case name. Leave default to inherit from a predecessor.",
        ),
        "xyz_file": FileParameterEdit(
            "Ligand XYZ", docstring="Metal-complex geometry (.xyz).",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Working dir for canonical XYZ + molecular-graph JSON.",
        ),
        "metal_symbol": StringParameter(
            "Metal Symbol", default="",
            docstring="Metal element symbol (e.g. Re, Ru, Cu). Forwarded downstream.",
        ),
        "obabel_path": StringParameter(
            "obabel Executable", default="obabel",
            docstring="OpenBabel CLI binary (used unless the API fallback is forced).",
        ),
        "use_openbabel_python_api": BooleanParameter(
            "Force OpenBabel Python API", default=False,
            docstring="Skip the obabel CLI and use the Python API directly. "
                      "Useful when mgltools shadows the conda-forge obabel binary.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting ligand prep...", node_id=self.node_id, progress=0)
        try:
            _ensure_metaldock_modules()
            from pathlib import Path
            from metaldock_modules import ligand_prep
            from metaldock_modules.utils import (
                save_graph, count_heavy_atoms, validate_metal_symbol,
            )

            carried = _merge_predecessors(predecessor_data)
            result = NodeResult()

            case_name = (
                flow_vars["case_name"].get_value()
                or carried.get("case_name")
                or "complex"
            )
            xyz_file = self.resolve_path(flow_vars["xyz_file"].get_value())
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            metal_symbol = (flow_vars["metal_symbol"].get_value() or "").strip()
            if not metal_symbol:
                raise NodeException("setup", "Metal Symbol is required (e.g. Re, Ru, Cu).")
            validate_metal_symbol(metal_symbol)        # rejects unsupported metals
            metal_symbol = metal_symbol.capitalize()   # element case ("RE"→"Re") to match graph

            obabel_path = flow_vars["obabel_path"].get_value() or "obabel"
            use_api = bool(flow_vars["use_openbabel_python_api"].get_value())

            stem = Path(xyz_file).stem
            canonical_xyz = Path(output_dir) / f"{stem}_c.xyz"

            stream_log("Canonicalizing XYZ + building graph...",
                       node_id=self.node_id, progress=40)
            ligand_prep.canonicalize_xyz(
                Path(xyz_file), canonical_xyz,
                obabel_path=obabel_path, use_python_api=use_api,
            )
            graph = ligand_prep.build_graph_from_xyz(canonical_xyz)

            graph_json = Path(output_dir) / "mol_graph.json"
            save_graph(graph, graph_json)
            n_heavy = count_heavy_atoms(canonical_xyz)

            result.data = {
                **carried,
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "metal_symbol": metal_symbol,
                "canonical_xyz": self.format_output_path(str(canonical_xyz)),
                "graph_json": self.format_output_path(str(graph_json)),
                "n_heavy_atoms": n_heavy,
                "n_atoms": graph.number_of_nodes(),
                "n_bonds": graph.number_of_edges(),
            }
            result.metadata["case_name"] = case_name
            result.files["input"] = {"xyz": self.format_output_path(xyz_file)}
            result.files["output"] = {
                "canonical_xyz": result.data["canonical_xyz"],
                "graph_json": result.data["graph_json"],
            }
            result.success = True
            result.message = (
                f"{metal_symbol} complex: {graph.number_of_nodes()} atoms, "
                f"{graph.number_of_edges()} bonds, {n_heavy} heavy"
            )
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("ligand prep", str(e))
