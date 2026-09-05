"""ligand-pdbqt-meeko — docking-ready PDBQT ligands from an SDF, with Meeko.

Salpa node wrapper. The computation itself lives in core.py; this file only
handles Salpa I/O — reading parameters, resolving paths, streaming progress,
and shaping the NodeResult.
"""

import os

from bocoflow_core.node import Node, NodeException, NodeResult, upstream
from bocoflow_core.parameters import (
    BooleanParameter,
    FileParameterEdit,
    FolderParameter,
)
from bocoflow_core.stream_logger import stream_log

# Standard 3-stage import. Do not replace this with a plain `from .core import` —
# each stage covers a real context, and stage 3 is what lets the SERVER import
# node.py (to read OPTIONS for the UI) without your heavy deps installed.
try:
    from .core import process                    # package context (pytest/direct)
except ImportError:
    try:
        from core import process                 # at runtime: node dir is on sys.path
    except ImportError:
        process = None                           # server env: heavy deps absent


#: How to run this node — the values `salpa smoke` feeds it, and the ones your own
#: node.py test uses. ONE declaration, two readers, so they cannot drift.
#:
#: Strings starting with `demo_data/` are resolved relative to this directory.
#: Anything you leave out is guessed; output folders get the run's case directory.
DEMO_CONFIG = {
    "input_sdf": "demo_data/ligands.sdf",
    "rigid_macrocycles": False,
    "hydrate": False,
    "flexible_amides": False,
}


class LigandPdbqtMeeko(Node):
    """Docking-ready PDBQT ligands from an SDF, with Meeko."""

    # Mirrors meta.toml. The registry reads meta.toml; these attributes are the
    # in-code fallback and keep node.py readable on its own.
    category = "Docking"
    tags = ["pdbqt", "meeko", "autodock", "vina", "ligand-prep", "sdf", "docking"]

    # NOTE: do NOT add force_to_run — it is inherited from Node.BASE_OPTIONS.
    OPTIONS = {
        "input_sdf": FileParameterEdit(
            "Ligands (SDF)",
            docstring="An SDF of 3D ligands, hydrogens preferably explicit — what "
                      "smiles-to-3d writes. Leave empty to take the SDF carried from a "
                      "predecessor node.",
        ),
        "rigid_macrocycles": BooleanParameter(
            "Rigid macrocycles", default=False,
            docstring="Keep macrocycle rings rigid. Meeko's default makes them flexible "
                      "with glue atoms, which AutoDock Vina handles and AutoDock4 does not.",
        ),
        "hydrate": BooleanParameter(
            "Hydrated docking", default=False,
            docstring="Add Meeko's explicit water sites for hydrated docking.",
        ),
        "flexible_amides": BooleanParameter(
            "Flexible amides", default=False,
            docstring="Let amide C–N bonds rotate. Off by default, as in AutoDock.",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory where one <name>.pdbqt per molecule is written.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        """Run the node.

        Args:
            predecessor_data: one entry per upstream node, each that node's
                `result["data"]` — already unwrapped, so read its keys directly.
                None when the node runs on its own, which it should be able to:
                every input a parameter, this only filling one in when left blank.
            flow_vars: dict of parameter name -> Parameter. Read with .get_value().

        Returns:
            str: JSON from NodeResult.to_json()
        """
        stream_log("Starting ligand PDBQT preparation...", node_id=self.node_id, progress=0)
        try:
            # Stage-3 import landed (see above) — the real deps are missing.
            if process is None:
                raise NodeException(
                    "setup",
                    "core.py could not be imported. Check that this node's pixi "
                    "environment installed meeko and rdkit (`pixi install`).",
                )

            carried = upstream(predecessor_data)
            result = NodeResult()

            sdf_val = flow_vars["input_sdf"].get_value()
            sdf_path = self.resolve_path(sdf_val) if sdf_val else ""
            # Fall back to the SDF an upstream node wrote — smiles-to-3d puts it under
            # `output_file`; any node may name one under `sdf_file`.
            if not sdf_path:
                for key in ("sdf_file", "output_file"):
                    if carried.get(key):
                        sdf_path = self.resolve_path(carried[key])
                        break
            if not sdf_path:
                raise NodeException(
                    "input",
                    "No SDF provided (via the Ligands parameter or a predecessor node).",
                )
            if not os.path.isfile(sdf_path):
                raise NodeException("input", f"SDF not found: {sdf_path}")

            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)

            stream_log("Preparing ligands with Meeko...", node_id=self.node_id, progress=30)
            try:
                out = process(
                    sdf_path, out_dir,
                    rigid_macrocycles=bool(flow_vars["rigid_macrocycles"].get_value()),
                    hydrate=bool(flow_vars["hydrate"].get_value()),
                    flexible_amides=bool(flow_vars["flexible_amides"].get_value()),
                )
            except ValueError as e:
                # A record Meeko refuses is the input's problem, named as such —
                # not skipped, and not an "execution" failure.
                raise NodeException("input", str(e))

            for lig in out["ligands"]:
                stream_log(
                    f"  {lig['name']}: {lig['n_atoms_pdbqt']} PDBQT atoms from "
                    f"{lig['n_atoms_input']}, {lig['torsions']} torsions, charge "
                    f"{lig['total_charge']:+.3f}"
                    + (" (hydrogens added)" if lig["hydrogens_added"] else ""),
                    node_id=self.node_id,
                )

            result.data = {
                **carried,                 # carry upstream keys forward
                "ligands": out["ligands"],
                "pdbqt_files": [self.format_output_path(p) for p in out["pdbqt_files"]],
                "n_in": out["n_in"],
                "n_out": out["n_out"],
                "output_dir": self.format_output_path(out_dir),
            }
            result.files["output"] = {
                lig["name"]: self.format_output_path(lig["file"]) for lig in out["ligands"]
            }
            result.success = True
            result.message = f"{out['n_out']} ligand(s) prepared → PDBQT"
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise                          # keep the original stage tag
        except Exception as e:
            raise NodeException("execution", str(e))
