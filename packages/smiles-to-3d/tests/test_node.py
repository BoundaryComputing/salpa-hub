"""Tests for the smiles-to-3d node.

Two levels, and the split is load-bearing:

  * core.py is pure Python and imports in any env; its RDKit-using tests skip
    cleanly where rdkit is absent, and RUN under `pixi run test`.
  * node.py imports bocoflow_core (the node-authoring API). `pixi run test`
    provides it via the `test` env's bocoflow-core-sdk, so node.py's tests RUN
    there too. They carry @needs_runtime so they skip GRACEFULLY (not error) in a
    bare env without the SDK.
"""

import json
import math
import os

import pytest

from core import conformer_coordinates, embed_conformers, process, read_smiles

try:
    from node import DEMO_CONFIG, SmilesTo3d
    _NODE_OK = True
except Exception:  # node.py imports bocoflow_core (runtime-provided)
    _NODE_OK = False

needs_runtime = pytest.mark.skipif(
    not _NODE_OK,
    reason="bocoflow_core not importable — run `pixi run test` (its test env provides bocoflow-core-sdk)",
)
needs_rdkit = pytest.mark.skipif(
    pytest.importorskip("rdkit", reason="rdkit not installed") is None, reason="rdkit not installed"
)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEMO = os.path.join(HERE, "demo_data", "ligands.smi")

# Hand-worked from the molecular formulas, not copied from the node's output:
# atoms with hydrogens, and heavy atoms.
EXPECTED = {
    "ethanol":  (9, 3),    # C2H6O
    "aspirin":  (21, 13),  # C9H8O4
    "caffeine": (24, 14),  # C8H10N4O2
}


# ---- pure parsing (runs anywhere) ------------------------------------------

def test_read_smiles_parses_names_comments_and_blank_lines():
    text = "# a comment\nCCO ethanol\n\n  c1ccccc1  \nCC(=O)O acetic acid\n"
    assert read_smiles(text) == [
        ("CCO", "ethanol"), ("c1ccccc1", "mol2"), ("CC(=O)O", "acetic acid"),
    ]


def test_read_smiles_from_the_demo_file():
    assert [name for _, name in read_smiles("", DEMO)] == ["ethanol", "aspirin", "caffeine"]


# ---- the chemistry (needs rdkit: present under `pixi run test`) -------------

@needs_rdkit
def test_atom_counts_match_the_formulas():
    for smiles, name in read_smiles("", DEMO):
        mol, energies = embed_conformers(smiles, 1, 42)
        assert (mol.GetNumAtoms(), mol.GetNumHeavyAtoms()) == EXPECTED[name], name
        assert len(energies) == 1 and math.isfinite(energies[0])


@needs_rdkit
def test_an_unparseable_smiles_is_refused_not_skipped():
    with pytest.raises(ValueError, match="could not parse"):
        embed_conformers("this is not smiles", 1, 42)


@needs_rdkit
def test_the_same_seed_gives_the_same_coordinates():
    assert conformer_coordinates("CC(=O)Oc1ccccc1C(=O)O", 42) == conformer_coordinates("CC(=O)Oc1ccccc1C(=O)O", 42)


@needs_rdkit
def test_process_writes_every_conformer_with_its_energy(tmp_path):
    out = str(tmp_path / "out.sdf")
    result = process(read_smiles("", DEMO), out, num_conformers=2, seed=7)
    assert result["n_in"] == result["n_out"] == 3
    assert [m["n_conformers"] for m in result["molecules"]] == [2, 2, 2]
    from rdkit import Chem
    records = [m for m in Chem.SDMolSupplier(out, removeHs=False) if m is not None]
    assert len(records) == 6                                   # 3 molecules × 2 conformers
    assert records[0].GetProp("_Name") == "ethanol"
    assert records[0].GetNumAtoms() == 9                        # hydrogens kept in the file
    assert math.isfinite(float(records[0].GetProp("MMFF94_energy_kcal_mol")))
    # the summary's lowest energy is one of the energies written for that molecule
    written = [float(r.GetProp("MMFF94_energy_kcal_mol")) for r in records[:2]]
    assert result["molecules"][0]["lowest_energy_kcal_mol"] == round(min(written), 4)


# ---- node wiring (needs the bocoflow_core runtime) -------------------------

def _resolve(value):
    """`demo_data/...` resolves against the node directory, as DEMO_CONFIG declares."""
    if isinstance(value, str) and value.startswith("demo_data/"):
        return os.path.join(HERE, value)
    return value


@needs_runtime
def test_node_instantiation_and_options():
    node = SmilesTo3d({"node_id": "test"})
    for key in ("inline_text", "input_file", "num_conformers", "random_seed", "output_dir"):
        assert key in node.OPTIONS


@needs_runtime
@needs_rdkit
def test_execute_end_to_end(tmp_path):
    class _P:  # minimal Parameter stub; OPTIONS objects are class-level and shared
        def __init__(self, v):
            self._v = v

        def get_value(self):
            return self._v

    node = SmilesTo3d({"node_id": "test"})
    flow_vars = {k: _P(_resolve(v)) for k, v in DEMO_CONFIG.items()}
    flow_vars["inline_text"] = _P("")
    flow_vars["output_dir"] = _P(str(tmp_path))
    result = json.loads(node.execute(None, flow_vars))
    assert result["success"] is True
    assert result["data"]["n_out"] == 3
    assert [m["n_atoms"] for m in result["data"]["molecules"]] == [9, 21, 24]
    assert os.path.exists(tmp_path / "conformers.sdf")


@needs_runtime
def test_no_input_fails_loudly():
    class _P:
        def __init__(self, v):
            self._v = v

        def get_value(self):
            return self._v

    node = SmilesTo3d({"node_id": "test"})
    flow_vars = {
        "inline_text": _P(""), "input_file": _P(""), "num_conformers": _P(1),
        "random_seed": _P(42), "output_dir": _P("."),
    }
    with pytest.raises(Exception, match="No SMILES"):
        node.execute(None, flow_vars)
