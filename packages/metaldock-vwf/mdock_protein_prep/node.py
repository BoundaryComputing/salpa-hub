"""MetalDock Protein Prep — wraps metaldock_modules.protein_prep.prepare_protein.

Clean a protein PDB, protonate it with pdb2pqr, and convert to PDBQT with
AutoDockTools' prepare_receptor4.py. This is the entry node of the metaldock
pipeline; it forwards ``receptor_pdbqt`` / ``cleaned_pdb`` down the chain so the
later docking + analysis nodes can find them without re-wiring.
"""

import os
import shutil
import sys

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter, FileParameterEdit, FloatParameter, FolderParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log

_NODE_DIR = os.path.dirname(os.path.abspath(__file__))


def _ensure_metaldock_modules():
    """Put the metaldock_modules package on sys.path.

    Resolution: METALDOCK_SRC env var > node-bundled scripts/ > repo src/.
    """
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
        "Cannot locate the metaldock_modules package. Set the METALDOCK_SRC "
        "environment variable to a directory containing metaldock_modules/.",
    )


def _find_mgltools_dir():
    """Locate the MGLTools AutoDockTools/Utilities24 directory in the env."""
    candidates = [
        os.path.join(sys.prefix, "MGLToolsPckgs", "AutoDockTools", "Utilities24"),
        os.path.join(sys.prefix, "lib", "python2.7", "site-packages",
                     "AutoDockTools", "Utilities24"),
    ]
    for c in candidates:
        if os.path.isfile(os.path.join(c, "prepare_receptor4.py")):
            return c
    p = shutil.which("prepare_receptor4.py")
    return os.path.dirname(p) if p else None


class MdockProteinPrep(Node):
    """Clean → protonate (pdb2pqr) → PDBQT (prepare_receptor4) a protein."""

    category = "Metal Docking"
    tags = ["metaldock", "protein-prep", "pdb2pqr", "pdbqt", "mgltools"]

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name", default="complex",
            docstring="Name for this docking case (forwarded downstream).",
        ),
        "pdb_file": FileParameterEdit(
            "Protein PDB", docstring="Input protein structure (.pdb).",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Working dir for protonated/cleaned PDB + receptor PDBQT.",
        ),
        "ph": FloatParameter(
            "Protonation pH", default=7.4,
            docstring="pH used by pdb2pqr for protonation.",
        ),
        "clean": BooleanParameter(
            "Remove HETATM", default=True,
            docstring="Strip ligands/cofactors/waters (HETATM lines) before PDBQT.",
        ),
        "pdb2pqr_path": StringParameter(
            "pdb2pqr Executable", default="pdb2pqr30",
            docstring="pdb2pqr binary name or path.",
        ),
        "python_path": StringParameter(
            "MGLTools Python", default="pythonsh",
            docstring="Interpreter for prepare_receptor4.py (MGLTools' pythonsh).",
        ),
        "mgltools_dir": FolderParameter(
            "MGLTools Dir (optional)",
            docstring="AutoDockTools/Utilities24 dir. Leave empty to auto-detect "
                      "from the active environment.",
            optional=True,
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting protein prep...", node_id=self.node_id, progress=0)
        try:
            _ensure_metaldock_modules()
            from pathlib import Path
            from metaldock_modules import protein_prep

            result = NodeResult()
            case_name = flow_vars["case_name"].get_value() or "complex"
            pdb_file = self.resolve_path(flow_vars["pdb_file"].get_value())
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(output_dir, exist_ok=True)

            ph = float(flow_vars["ph"].get_value() or 7.4)
            clean = bool(flow_vars["clean"].get_value())
            pdb2pqr_path = flow_vars["pdb2pqr_path"].get_value() or "pdb2pqr30"
            python_path = flow_vars["python_path"].get_value() or "pythonsh"

            mgltools_dir = flow_vars["mgltools_dir"].get_value()
            mgltools_dir = self.resolve_path(mgltools_dir) if mgltools_dir else None
            if not mgltools_dir:
                mgltools_dir = _find_mgltools_dir()
            if not mgltools_dir:
                raise NodeException(
                    "setup",
                    "MGLTools not found. Install mgltools in the environment or "
                    "set the MGLTools Dir parameter.",
                )

            stream_log("Running pdb2pqr + prepare_receptor4...",
                       node_id=self.node_id, progress=30)
            prep = protein_prep.prepare_protein(
                pdb_path=Path(pdb_file),
                output_dir=Path(output_dir),
                ph=ph,
                clean=clean,
                pdb2pqr_path=pdb2pqr_path,
                python_path=python_path,
                mgltools_dir=mgltools_dir,
            )

            receptor_pdbqt = str(prep["pdbqt"])
            cleaned_pdb = str(prep["cleaned_pdb"])
            protonated_pdb = (
                str(prep["protonated_pdb"]) if prep.get("protonated_pdb") else ""
            )

            result.data = {
                "case_name": case_name,
                "working_path": self.format_output_path(output_dir),
                "receptor_pdbqt": self.format_output_path(receptor_pdbqt),
                "cleaned_pdb": self.format_output_path(cleaned_pdb),
                "protonated_pdb": (
                    self.format_output_path(protonated_pdb) if protonated_pdb else ""
                ),
            }
            result.metadata["case_name"] = case_name
            result.files["input"] = {"pdb": self.format_output_path(pdb_file)}
            result.files["output"] = {
                "receptor_pdbqt": result.data["receptor_pdbqt"],
                "cleaned_pdb": result.data["cleaned_pdb"],
            }
            result.success = True
            result.message = f"Receptor PDBQT ready: {os.path.basename(receptor_pdbqt)}"
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            raise NodeException("protein prep", str(e))
