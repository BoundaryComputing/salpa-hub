"""MetalDock QM Charges — wraps metaldock_modules.qm_charges.run_qm_and_enrich_graph.

Compute CM5 partial charges and bond orders on the metal complex and enrich the
molecular graph with them. Reloads the graph from the ``graph_json`` produced by
mdock_ligand_prep, re-saves the enriched graph, and forwards the new
``graph_json`` downstream.

Four engines. ``xtb`` (GFN1) is the default because it is the only one that
needs no user-supplied binary: it installs from conda-forge with the node and
finishes in seconds. ORCA is the DFT accuracy reference and is what published
numbers should come from, but it must be downloaded and registered for
separately; Gaussian and ADF are commercial. Selecting anything other than xtb
means supplying that program yourself.
"""

import os
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FloatParameter, FolderParameter,
    IntegerParameter, SelectParameter, StringParameter, TextParameter,
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
    "output_dir": "qm",
    # Normally inherited from mdock_ligand_prep; declared so this node can also
    # be run — and smoke-tested — on its own.
    "graph_json": "demo_data/1jzi_re_mol_graph.json",
    "xyz_file": "demo_data/1jzi_D_REP_c.xyz",
    # xtb, not orca: it is the only engine that needs no user-supplied binary, so
    # this is the only setting under which the demo runs unattended. Its CM5
    # charges land within 0.05 e of the ORCA reference on the metal — the stored
    # comparison is demo_data/1jzi_re_orca_reference_graph.json.
    "engine": "xtb",
    "geom_opt": False,       # single point; geometry optimization is far slower
    "charge": 1,             # the Re complex is a cation
    "spin": 0.0,             # closed shell
    "ncpu": 4,
}


class MdockQmCharges(Node):
    """DFT (ORCA/Gaussian/ADF) → CM5 charges + bond orders → enrich graph."""

    category = "Metal Docking"
    tags = ["metaldock", "qm-charges", "cm5", "orca", "dft", "metal-complex"]

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory", docstring="Working dir for QM scratch + enriched graph.",
        ),
        "graph_json": FileParameterEdit(
            "Molecular Graph (.json)", default="",
            docstring="Graph from mdock_ligand_prep. Leave empty to auto-discover "
                      "`graph_json` from a predecessor.",
            optional=True,
        ),
        "xyz_file": FileParameterEdit(
            "Canonical XYZ", default="",
            docstring="Canonicalized geometry. Leave empty to auto-discover "
                      "`canonical_xyz` from a predecessor.",
            optional=True,
        ),
        "engine": SelectParameter(
            "QM Engine", options=["xtb", "orca", "gaussian", "adf"], default="xtb",
            docstring="Charge backend. 'xtb' (GFN1, semi-empirical) ships with "
                      "this node and runs in seconds — use it to get a pipeline "
                      "working and to screen. 'orca' is DFT and the accuracy "
                      "reference for transition metals, but you must download "
                      "ORCA yourself. Gaussian/ADF are commercial.",
        ),
        "geom_opt": BooleanParameter(
            "Geometry Optimization", default=False,
            docstring="Optimize geometry before charges (expensive). Off = single-point.",
        ),
        "charge": IntegerParameter(
            "Total Charge", default=0, docstring="Net charge of the complex.",
        ),
        "spin": FloatParameter(
            "Spin (unpaired e⁻)", default=0.0,
            docstring="Number of unpaired electrons.",
        ),
        "ncpu": IntegerParameter("CPU Cores", default=4),
        # ── xTB ───────────────────────────────────────────────────────
        "xtb_path": FolderParameter(
            "xtb Directory (optional)", default="",
            docstring="xtb install dir. Leave empty to use the xtb that ships "
                      "in this node's environment.",
            optional=True,
        ),
        "xtb_solvent": StringParameter(
            "xtb Solvent (optional)", default="",
            docstring="ALPB implicit solvent, e.g. 'water'. Empty = gas phase.",
            optional=True,
        ),
        "xtb_accuracy": FloatParameter(
            "xtb Accuracy", default=1.0,
            docstring="xtb SCC accuracy; lower is tighter (default 1.0).",
        ),
        # ── ORCA ──────────────────────────────────────────────────────
        "orca_path": FolderParameter(
            "ORCA Directory (optional)", default="",
            docstring="ORCA install dir. Leave empty if orca is on PATH or "
                      "ASE_ORCA_COMMAND is set.",
            optional=True,
        ),
        "orcasimpleinput": StringParameter(
            "ORCA Simple Input", default="PBE def2-TZVP CPCM(Water)",
            docstring="ORCA `! ...` keyword line.",
        ),
        "orcablocks": TextParameter(
            "ORCA Blocks", default="",
            docstring="Extra ORCA `%...end` blocks (optional).", optional=True,
        ),
        # ── Gaussian / ADF ────────────────────────────────────────────
        "functional": StringParameter(
            "Functional (Gaussian/ADF)", default="PBE", optional=True,
        ),
        "basis_set": StringParameter(
            "Basis Set (Gaussian/ADF)", default="def2-TZVP", optional=True,
        ),
        "solvent": StringParameter(
            "Solvent (Gaussian/ADF)", default="", optional=True,
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting QM charges...", node_id=self.node_id, progress=0)
        try:
            _ensure_metaldock_modules()
            from pathlib import Path
            from metaldock_modules import qm_charges
            from metaldock_modules.utils import load_graph, save_graph

            carried = _merge_predecessors(predecessor_data)
            result = NodeResult()

            case_name = (
                flow_vars["case_name"].get_value() or carried.get("case_name") or "complex"
            )
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            graph_ref = flow_vars["graph_json"].get_value() or carried.get("graph_json")
            xyz_ref = flow_vars["xyz_file"].get_value() or carried.get("canonical_xyz")
            if not graph_ref:
                raise NodeException("setup", "No graph_json (param or predecessor).")
            if not xyz_ref:
                raise NodeException("setup", "No canonical XYZ (param or predecessor).")
            graph_json = self.resolve_path(graph_ref)
            xyz_file = self.resolve_path(xyz_ref)

            engine = (flow_vars["engine"].get_value() or "xtb").lower()
            geom_opt = bool(flow_vars["geom_opt"].get_value())
            charge = int(flow_vars["charge"].get_value() or 0)
            spin = float(flow_vars["spin"].get_value() or 0.0)
            ncpu = int(flow_vars["ncpu"].get_value() or 1)

            orca_path = flow_vars["orca_path"].get_value()
            orca_path = self.resolve_path(orca_path) if orca_path else None

            xtb_path = flow_vars["xtb_path"].get_value()
            xtb_path = self.resolve_path(xtb_path) if xtb_path else None

            graph = load_graph(Path(graph_json))

            stream_log(f"Running {engine.upper()} ({'opt' if geom_opt else 'single-point'})...",
                       node_id=self.node_id, progress=20)
            qm = qm_charges.run_qm_and_enrich_graph(
                graph=graph,
                xyz_path=Path(xyz_file),
                output_dir=Path(output_dir),
                engine=engine,
                geom_opt=geom_opt,
                charge=charge,
                spin=spin,
                ncpu=ncpu,
                orca_path=orca_path,
                orcasimpleinput=flow_vars["orcasimpleinput"].get_value()
                or "PBE def2-TZVP CPCM(Water)",
                orcablocks=flow_vars["orcablocks"].get_value() or "",
                functional=flow_vars["functional"].get_value() or "PBE",
                basis_set=flow_vars["basis_set"].get_value() or "def2-TZVP",
                solvent=flow_vars["solvent"].get_value() or "",
                xtb_path=xtb_path,
                xtb_accuracy=float(flow_vars["xtb_accuracy"].get_value() or 1.0),
                xtb_solvent=flow_vars["xtb_solvent"].get_value() or "",
            )

            enriched_json = Path(output_dir) / "enriched_graph.json"
            save_graph(qm["graph"], enriched_json)

            result.data = {
                **carried,
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "graph_json": self.format_output_path(str(enriched_json)),
                "qm_energy": str(qm.get("energy", "")),
                "qm_run_type": qm.get("run_type", ""),
                "qm_engine": engine,
            }
            if qm.get("output_xyz"):
                result.data["qm_output_xyz"] = self.format_output_path(str(qm["output_xyz"]))
            result.metadata["case_name"] = case_name
            result.files["output"] = {"graph_json": result.data["graph_json"]}
            result.success = True
            result.message = (
                f"{engine.upper()} {qm.get('run_type', '')}: enriched graph with "
                f"CM5 charges (energy {qm.get('energy', 'n/a')})"
            )
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("qm charges", str(e))
