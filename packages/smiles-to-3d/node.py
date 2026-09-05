"""smiles-to-3d — 3D conformers from SMILES, written as SDF (RDKit ETKDG + MMFF94).

Salpa node wrapper. The computation itself lives in core.py; this file only
handles Salpa I/O — reading parameters, resolving paths, streaming progress,
and shaping the NodeResult.
"""

import os

from bocoflow_core.node import Node, NodeException, NodeResult, upstream
from bocoflow_core.parameters import (
    FileParameterEdit,
    FolderParameter,
    IntegerParameter,
    TextParameter,
)
from bocoflow_core.stream_logger import stream_log

# Standard 3-stage import. Do not replace this with a plain `from .core import` —
# each stage covers a real context, and stage 3 is what lets the SERVER import
# node.py (to read OPTIONS for the UI) without your heavy deps installed.
try:
    from .core import process, read_smiles          # package context (pytest/direct)
except ImportError:
    try:
        from core import process, read_smiles        # at runtime: node dir is on sys.path
    except ImportError:
        process = read_smiles = None                 # server env: heavy deps absent


#: How to run this node — the values `salpa smoke` feeds it, and the ones your own
#: node.py test uses. ONE declaration, two readers, so they cannot drift.
#:
#: Strings starting with `demo_data/` are resolved relative to this directory.
#: Anything you leave out is guessed; output folders get the run's case directory.
DEMO_CONFIG = {
    "input_file": "demo_data/ligands.smi",
    "num_conformers": 1,
    "random_seed": 42,
}


class SmilesTo3d(Node):
    """3D conformers from SMILES, written as SDF (RDKit ETKDG + MMFF94)."""

    # Mirrors meta.toml. The registry reads meta.toml; these attributes are the
    # in-code fallback and keep node.py readable on its own.
    category = "Chemistry"
    tags = ["smiles", "3d", "conformer", "sdf", "rdkit", "ligand-prep", "docking"]

    # NOTE: do NOT add force_to_run — it is inherited from Node.BASE_OPTIONS.
    OPTIONS = {
        "inline_text": TextParameter(
            "SMILES (one per line)", default="",
            docstring="One molecule per line as SMILES, with an optional name after a "
                      "space, e.g. `CCO ethanol`. Leave empty to use an input file or "
                      "SMILES carried from a predecessor node.",
        ),
        "input_file": FileParameterEdit(
            "SMILES file (optional)",
            docstring="A .smi or .txt file: one SMILES per line, optional name after "
                      "whitespace, `#` comments allowed.",
        ),
        "num_conformers": IntegerParameter(
            "Conformers per molecule", default=1,
            docstring="How many 3D conformers to embed and optimise for each molecule. "
                      "All of them are written; the summary reports the lowest MMFF94 energy.",
        ),
        "random_seed": IntegerParameter(
            "Random seed", default=42,
            docstring="Seed for the ETKDG embedding. The same seed gives the same "
                      "coordinates, so a run can be repeated exactly.",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory where conformers.sdf is written.",
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
        stream_log("Starting SMILES → 3D...", node_id=self.node_id, progress=0)
        try:
            # Stage-3 import landed (see above) — the real deps are missing.
            if process is None:
                raise NodeException(
                    "setup",
                    "core.py could not be imported. Check that this node's pixi "
                    "environment installed rdkit (`pixi install`).",
                )

            carried = upstream(predecessor_data)
            result = NodeResult()

            text = flow_vars["inline_text"].get_value() or ""
            file_val = flow_vars["input_file"].get_value()
            file_path = self.resolve_path(file_val) if file_val else ""

            records = read_smiles(text, file_path)
            # Fall back to SMILES carried from an upstream node: a list of strings,
            # or the `molecules` summary another run of this node produced.
            if not records and carried.get("smiles"):
                carried_smiles = carried["smiles"]
                if isinstance(carried_smiles, str):
                    carried_smiles = [carried_smiles]
                records = read_smiles("\n".join(carried_smiles))
            if not records:
                raise NodeException(
                    "input",
                    "No SMILES provided (via text, file, or a predecessor node).",
                )

            num_conformers = int(flow_vars["num_conformers"].get_value() or 1)
            if num_conformers < 1:
                raise NodeException("input", "Conformers per molecule must be at least 1.")
            seed = int(flow_vars["random_seed"].get_value() or 0)

            out_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            os.makedirs(out_dir, exist_ok=True)
            out_file = os.path.join(out_dir, "conformers.sdf")

            stream_log(
                f"Embedding {len(records)} molecule(s), {num_conformers} conformer(s) each...",
                node_id=self.node_id, progress=30,
            )
            try:
                out = process(records, out_file, num_conformers=num_conformers, seed=seed)
            except ValueError as e:
                # A SMILES RDKit refuses is the user's input problem, named as such —
                # not skipped, and not an "execution" failure.
                raise NodeException("input", str(e))

            for m in out["molecules"]:
                stream_log(
                    f"  {m['name']}: {m['n_heavy_atoms']} heavy atoms, {m['n_atoms']} with H, "
                    f"lowest MMFF94 energy {m['lowest_energy_kcal_mol']} kcal/mol",
                    node_id=self.node_id,
                )

            result.data = {
                **carried,                 # carry upstream keys forward
                "molecules": out["molecules"],
                "smiles": [m["smiles"] for m in out["molecules"]],
                "n_in": out["n_in"],
                "n_out": out["n_out"],
                "output_file": self.format_output_path(out_file),
            }
            result.files["output"] = {"conformers_sdf": result.data["output_file"]}
            result.success = True
            result.message = (
                f"{out['n_out']} molecule(s), {num_conformers} conformer(s) each → conformers.sdf"
            )
            stream_log(result.message, node_id=self.node_id, progress=100)
            return result.to_json()

        except NodeException:
            raise                          # keep the original stage tag
        except Exception as e:
            raise NodeException("execution", str(e))
