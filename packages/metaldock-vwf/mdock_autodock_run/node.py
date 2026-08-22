"""MetalDock AutoDock Run — wraps metaldock_modules.autodock_run.run_autodock.

Generate the AutoGrid/AutoDock parameter files, run autogrid4 then autodock4,
and extract docked poses. Pulls the ligand PDBQT, receptor PDBQT, molecular
graph, and metal symbol from upstream nodes; forwards the DLG and pose paths
(plus the cleaned protein PDB and heavy-atom count) to the analysis node.
"""

import os
import shutil
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


def _mgltools_env_prefix():
    """The sibling pixi environment holding MGLTools, or None.

    MGLTools is installed into its OWN environment (see the package pixi.toml):
    its conda package replaces `bin/python` with Python 2.7, and BoCoFlow
    launches nodes as `pixi run python -m bocoflow_core.node_runner`, so sharing
    one environment stops every node in this package from starting at all.

    Node code therefore runs in `default`, and MGLTools sits beside it:

        <env root>/.pixi/envs/default    <- sys.prefix, where this code runs
        <env root>/.pixi/envs/mgltools   <- what we are looking for
    """
    sibling = os.path.join(os.path.dirname(sys.prefix), "mgltools")
    return sibling if os.path.isdir(sibling) else None


def _find_pythonsh():
    """MGLTools' launcher, which runs the Python 2.7 its scripts need.

    Not on PATH once MGLTools lives in its own environment, so look there first
    and only then fall back to PATH (covers a hand-rolled MGLTools install).
    """
    env = _mgltools_env_prefix()
    if env:
        candidate = os.path.join(env, "bin", "pythonsh")
        if os.path.isfile(candidate):
            return candidate
    return shutil.which("pythonsh") or "pythonsh"


def _find_mgltools_dir():
    """Locate the MGLTools AutoDockTools/Utilities24 directory.

    Checks the dedicated MGLTools environment first, then this node's own
    prefix — the latter still works for anyone who installed MGLTools into the
    single shared environment the older layout used.
    """
    roots = [r for r in (_mgltools_env_prefix(), sys.prefix) if r]
    for root in roots:
        candidates = [
            os.path.join(root, "MGLToolsPckgs", "AutoDockTools", "Utilities24"),
            os.path.join(root, "lib", "python2.7", "site-packages",
                         "AutoDockTools", "Utilities24"),
        ]
        for c in candidates:
            if os.path.isfile(os.path.join(c, 'prepare_gpf4.py')):
                return c
    p = shutil.which('prepare_gpf4.py')
    return os.path.dirname(p) if p else None


def _parse_vec3(text):
    """Parse 'x,y,z' (or 'x y z') into [float, float, float], or None if empty."""
    if not text:
        return None
    parts = [p for p in text.replace(",", " ").split() if p]
    if len(parts) != 3:
        raise NodeException("setup", f"Expected 3 comma-separated numbers, got: {text!r}")
    return [float(p) for p in parts]


# The values that make this node run against its own demo_data, declared once and
# read by `salpa smoke` and the shipped 1JZI workflow template alike. A parameter's
# type gives its shape and never its value: nothing can infer that the metal here is
# Re, or that the docking box belongs at the metal's coordinates.
DEMO_CONFIG = {
    "case_name": "1jzi_re",
    "output_dir": "docking",
    "num_poses": 10,
    # All three are inherited from upstream in the pipeline. Declared so the
    # docking step — the one with the longest setup chain — can be tried alone.
    "ligand_pdbqt": "demo_data/1jzi_re_ligand.pdbqt",
    "receptor_pdbqt": "demo_data/clean_1jzi.pdbqt",
    "graph_json": "demo_data/1jzi_re_enriched_graph.json",
    "metal_symbol": "Re",
    # Centred on the crystallographic Re position. Left empty the box would centre
    # itself on the metal, which is nearly the same thing but not reproducible.
    "box_center": "1.65,-7.803,27.176",
    "box_size": "20,20,20",
}


class MdockAutodockRun(Node):
    """GPF/DPF → autogrid4 → autodock4 → extract poses."""

    category = "Metal Docking"
    tags = ["metaldock", "autodock4", "autogrid4", "docking", "metal-complex"]

    OPTIONS = {
        "case_name": StringParameter("Case Name", default="complex"),
        "output_dir": FolderParameter(
            "Output Directory", docstring="Working dir for grid/docking files + poses.",
        ),
        "ligand_pdbqt": FileParameterEdit(
            "Ligand PDBQT", default="",
            docstring="Leave empty to auto-discover `ligand_pdbqt` from a predecessor.",
            optional=True,
        ),
        "receptor_pdbqt": FileParameterEdit(
            "Receptor PDBQT", default="",
            docstring="Leave empty to auto-discover `receptor_pdbqt` from a predecessor.",
            optional=True,
        ),
        "graph_json": FileParameterEdit(
            "Molecular Graph (.json)", default="",
            docstring="Used for the box center (metal position). Leave empty to "
                      "auto-discover `graph_json`.",
            optional=True,
        ),
        "metal_symbol": StringParameter(
            "Metal Symbol", default="", optional=True,
            docstring="Leave empty to inherit from a predecessor.",
        ),
        "parameter_file": FileParameterEdit(
            "AutoDock4 Parameter File", default="",
            docstring="metal_dock.dat. Leave empty to use the bundled default.",
            optional=True,
        ),
        "num_poses": IntegerParameter("Number of Poses", default=10),
        "box_center": StringParameter(
            "Box Center (x,y,z)", default="", optional=True,
            docstring="Grid box center in Å. Leave empty to center on the metal atom.",
        ),
        "box_size": StringParameter(
            "Box Size (x,y,z)", default="20,20,20",
            docstring="Grid box dimensions in Å.",
        ),
        "python_path": StringParameter("MGLTools Python", default="pythonsh"),
        "mgltools_dir": FolderParameter(
            "MGLTools Dir (optional)", default="", optional=True,
            docstring="AutoDockTools/Utilities24. Leave empty to auto-detect.",
        ),
        "autogrid4_path": StringParameter("autogrid4 Executable", default="autogrid4"),
        "autodock4_path": StringParameter("autodock4 Executable", default="autodock4"),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting AutoDock4 run...", node_id=self.node_id, progress=0)
        try:
            _ensure_metaldock_modules()
            from pathlib import Path
            from metaldock_modules import autodock_run
            from metaldock_modules.utils import (
                load_graph, validate_metal_symbol, default_parameter_file,
            )

            carried = _merge_predecessors(predecessor_data)
            result = NodeResult()

            case_name = (
                flow_vars["case_name"].get_value() or carried.get("case_name") or "complex"
            )
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            def _input(param_key, *carry_keys):
                ref = flow_vars[param_key].get_value()
                if not ref:
                    for ck in carry_keys:
                        if carried.get(ck):
                            ref = carried[ck]
                            break
                return self.resolve_path(ref) if ref else None

            ligand_pdbqt = _input("ligand_pdbqt", "ligand_pdbqt")
            receptor_pdbqt = _input("receptor_pdbqt", "receptor_pdbqt")
            graph_json = _input("graph_json", "graph_json")
            if not ligand_pdbqt:
                raise NodeException("setup", "No ligand_pdbqt (param or predecessor).")
            if not receptor_pdbqt:
                raise NodeException("setup", "No receptor_pdbqt (param or predecessor).")
            if not graph_json:
                raise NodeException("setup", "No graph_json (param or predecessor).")

            metal_symbol = (
                (flow_vars["metal_symbol"].get_value() or "").strip()
                or carried.get("metal_symbol")
            )
            if not metal_symbol:
                raise NodeException("setup", "Metal Symbol required (param or predecessor).")
            validate_metal_symbol(metal_symbol)        # rejects unsupported metals
            metal_symbol = metal_symbol.capitalize()   # element case to match graph

            param_ref = flow_vars["parameter_file"].get_value()
            parameter_file = (
                Path(self.resolve_path(param_ref)) if param_ref
                else default_parameter_file()
            )

            num_poses = int(flow_vars["num_poses"].get_value() or 10)
            box_center = _parse_vec3(flow_vars["box_center"].get_value())
            box_size = _parse_vec3(flow_vars["box_size"].get_value() or "20,20,20")

            python_path = flow_vars["python_path"].get_value() or "pythonsh"
            if python_path == "pythonsh":
                # the shown default, not a real PATH entry — resolve it
                python_path = _find_pythonsh()
            mgltools_dir = flow_vars["mgltools_dir"].get_value()
            mgltools_dir = self.resolve_path(mgltools_dir) if mgltools_dir else None
            if not mgltools_dir:
                mgltools_dir = _find_mgltools_dir()
            if not mgltools_dir:
                raise NodeException("setup", "MGLTools not found; set the MGLTools Dir param.")

            graph = load_graph(Path(graph_json))

            stream_log("Running autogrid4 + autodock4...",
                       node_id=self.node_id, progress=25)
            dock = autodock_run.run_autodock(
                ligand_pdbqt=Path(ligand_pdbqt),
                receptor_pdbqt=Path(receptor_pdbqt),
                output_dir=Path(output_dir),
                graph=graph,
                metal_symbol=metal_symbol,
                parameter_file=parameter_file,
                num_poses=num_poses,
                box_center=box_center,
                box_size=box_size,
                autogrid4_path=flow_vars["autogrid4_path"].get_value() or "autogrid4",
                autodock4_path=flow_vars["autodock4_path"].get_value() or "autodock4",
                python_path=python_path,
                mgltools_dir=mgltools_dir,
            )

            pose_xyz = [self.format_output_path(str(p)) for p in dock["pose_xyz_paths"]]
            pose_pdbqt = [self.format_output_path(str(p)) for p in dock["pose_pdbqt_paths"]]

            result.data = {
                **carried,
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "metal_symbol": metal_symbol,
                "dlg_path": self.format_output_path(str(dock["dlg_path"])),
                "pose_xyz_paths": pose_xyz,
                "pose_pdbqt_paths": pose_pdbqt,
                "num_poses": len(pose_xyz),
            }
            result.metadata["case_name"] = case_name
            result.files["output"] = {
                "dlg": result.data["dlg_path"],
                "poses": pose_xyz,
            }
            result.success = True
            result.message = f"Docked {len(pose_xyz)} pose(s) → {Path(dock['dlg_path']).name}"
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("autodock run", str(e))
