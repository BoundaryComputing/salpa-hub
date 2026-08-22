"""MetalDock Ligand PDBQT — wraps metaldock_modules.ligand_pdbqt.create_ligand_pdbqt.

Convert the charge-enriched molecular graph into an AutoDock PDBQT file with a
ROOT/BRANCH torsion tree. Coordination-sphere bonds are frozen by default, and
the torsion count is capped at AutoDock4's 32-torsion limit by freezing the
bonds closest to the metal first.
"""

import os
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FolderParameter, IntegerParameter,
    StringParameter,
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
    raise NodeException("setup", "Cannot locate metaldock_modules. Set METALDOCK_SRC.")


def _merge_predecessors(predecessor_data):
    merged = {}
    for pred in (predecessor_data or []):
        if not isinstance(pred, dict):
            continue
        scope = pred.get("data") if isinstance(pred.get("data"), dict) else pred
        for k, v in scope.items():
            if k not in ("success", "message", "metadata", "files"):
                merged[k] = v
    return merged


# The values that make this node run against its own demo_data, declared once and
# read by `salpa smoke` and the shipped 1JZI workflow template alike. A parameter's
# type gives its shape and never its value: nothing can infer that the metal here is
# Re, or that the docking box belongs at the metal's coordinates.
DEMO_CONFIG = {
    "case_name": "1jzi_re",
    "output_dir": "pdbqt",
    # Inherited from mdock_qm_charges in the pipeline; declared so this node can
    # be run alone. This is a real GFN1-xTB run's output, not a fabricated graph.
    "graph_json": "demo_data/1jzi_re_enriched_graph.json",
    "metal_symbol": "Re",
    "vacant_site": True,           # 1JZI is the vacant-coordination-sphere case
    "max_torsions": 32,
    "freeze_coordination_sphere": True,
}


class MdockLigandPdbqt(Node):
    """Enriched graph → ROOT/BRANCH PDBQT with metal-aware torsion freezing."""

    category = "Metal Docking"
    tags = ["metaldock", "pdbqt", "torsion-tree", "metal-complex", "autodock"]

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory", docstring="Working dir for the ligand PDBQT.",
        ),
        "graph_json": FileParameterEdit(
            "Enriched Graph (.json)", default="",
            docstring="Charge-enriched graph from mdock_qm_charges. Leave empty to "
                      "auto-discover `graph_json` from a predecessor.",
            optional=True,
        ),
        "metal_symbol": StringParameter(
            "Metal Symbol", default="",
            docstring="Metal element symbol. Leave empty to inherit from a predecessor.",
            optional=True,
        ),
        "vacant_site": BooleanParameter(
            "Add Vacant-Site Dummy", default=True,
            docstring="Add a dummy atom (DD) at the vacant coordination site.",
        ),
        "max_torsions": IntegerParameter(
            "Max Active Torsions", default=32,
            docstring="AutoDock4 hard limit is 32. Excess bonds nearest the metal "
                      "are frozen first.",
        ),
        "freeze_coordination_sphere": BooleanParameter(
            "Freeze Coordination Sphere", default=True,
            docstring="Freeze bonds within 2 bonds of the metal (rigid cage for "
                      "chelators like DOTA when the metal is bound).",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting ligand PDBQT...", node_id=self.node_id, progress=0)
        try:
            _ensure_metaldock_modules()
            from pathlib import Path
            from metaldock_modules import ligand_pdbqt
            from metaldock_modules.utils import load_graph, validate_metal_symbol

            carried = _merge_predecessors(predecessor_data)
            result = NodeResult()

            case_name = (
                flow_vars["case_name"].get_value() or carried.get("case_name") or "complex"
            )
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            graph_ref = flow_vars["graph_json"].get_value() or carried.get("graph_json")
            if not graph_ref:
                raise NodeException("setup", "No graph_json (param or predecessor).")
            graph_json = self.resolve_path(graph_ref)

            metal_symbol = (
                (flow_vars["metal_symbol"].get_value() or "").strip()
                or carried.get("metal_symbol")
            )
            if not metal_symbol:
                raise NodeException("setup", "Metal Symbol required (param or predecessor).")
            validate_metal_symbol(metal_symbol)        # rejects unsupported metals
            metal_symbol = metal_symbol.capitalize()   # element case to match graph

            vacant_site = bool(flow_vars["vacant_site"].get_value())
            max_torsions = int(flow_vars["max_torsions"].get_value() or 32)
            freeze = bool(flow_vars["freeze_coordination_sphere"].get_value())

            graph = load_graph(Path(graph_json))
            output_pdbqt = Path(output_dir) / f"{case_name}_ligand.pdbqt"

            stream_log("Building ROOT/BRANCH torsion tree...",
                       node_id=self.node_id, progress=40)
            ligand_pdbqt.create_ligand_pdbqt(
                graph=graph,
                metal_symbol=metal_symbol,
                output_path=output_pdbqt,
                vacant_site=vacant_site,
                max_torsions=max_torsions,
                freeze_coordination_sphere=freeze,
            )

            n_atom_lines = sum(
                1 for ln in output_pdbqt.read_text().splitlines()
                if ln.startswith("ATOM")
            )

            result.data = {
                **carried,
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "metal_symbol": metal_symbol,
                "ligand_pdbqt": self.format_output_path(str(output_pdbqt)),
                "vacant_site": vacant_site,
            }
            result.metadata["case_name"] = case_name
            result.files["output"] = {"ligand_pdbqt": result.data["ligand_pdbqt"]}
            result.success = True
            result.message = f"Ligand PDBQT: {n_atom_lines} atoms → {output_pdbqt.name}"
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("ligand pdbqt", str(e))
