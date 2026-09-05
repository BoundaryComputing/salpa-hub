"""Tests for the ligand-pdbqt-meeko node.

Two levels, and the split is load-bearing:

  * core.py is pure Python and imports in any env; its Meeko/RDKit tests skip
    cleanly where those are absent and RUN under `pixi run test`.
  * node.py imports bocoflow_core (the node-authoring API). `pixi run test`
    provides it via the `test` env's bocoflow-core-sdk, so node.py's tests RUN
    there too. They carry @needs_runtime so they skip GRACEFULLY (not error) in a
    bare env without the SDK.
"""

import json
import os

import pytest

from core import process, read_sdf, safe_name

try:
    from node import DEMO_CONFIG, LigandPdbqtMeeko
    _NODE_OK = True
except Exception:  # node.py imports bocoflow_core (runtime-provided)
    _NODE_OK = False

needs_runtime = pytest.mark.skipif(
    not _NODE_OK,
    reason="bocoflow_core not importable — run `pixi run test` (its test env provides bocoflow-core-sdk)",
)
needs_meeko = pytest.mark.skipif(
    pytest.importorskip("meeko", reason="meeko not installed") is None, reason="meeko not installed"
)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEMO = os.path.join(HERE, "demo_data", "ligands.sdf")

# Hand-worked from AutoDock's rules, not copied from the node's output: PDBQT atoms
# (heavy + polar H, non-polar H merged) and TORSDOF. See demo_data/README.md.
EXPECTED = {
    "ethanol":  (4, 1),
    "aspirin":  (14, 4),
    "caffeine": (14, 0),
}


# ---- pure helpers (run anywhere) -------------------------------------------

def test_safe_name_keeps_titles_usable_as_file_names():
    assert safe_name("aspirin", 0) == "aspirin"
    assert safe_name("my ligand (2)", 3) == "my_ligand_2"
    assert safe_name("", 3) == "ligand4"


# ---- the chemistry (needs meeko + rdkit: present under `pixi run test`) -----

@needs_meeko
def test_every_demo_ligand_gets_the_expected_atoms_and_torsions(tmp_path):
    out = process(DEMO, str(tmp_path))
    assert out["n_in"] == out["n_out"] == 3
    for lig in out["ligands"]:
        assert (lig["n_atoms_pdbqt"], lig["torsions"]) == EXPECTED[lig["name"]], lig
        assert abs(lig["total_charge"]) < 0.01, lig          # neutral molecules
        assert lig["hydrogens_added"] is False                 # the SDF already carries them
        assert os.path.isfile(lig["file"])
    text = open(os.path.join(tmp_path, "aspirin.pdbqt"), encoding="utf-8").read()
    assert text.startswith("REMARK") or text.startswith("ROOT")
    assert "TORSDOF 4" in text


@needs_meeko
def test_hydrogens_are_added_when_the_sdf_lacks_them_and_the_result_says_so(tmp_path):
    from rdkit import Chem
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    Chem.AllChem = __import__("rdkit.Chem.AllChem", fromlist=["AllChem"])
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=1)
    bare = Chem.RemoveHs(mol)                                  # 3D coordinates, no hydrogens
    bare.SetProp("_Name", "ethanol_noH")
    sdf = tmp_path / "noh.sdf"
    w = Chem.SDWriter(str(sdf)); w.write(bare); w.close()
    out = process(str(sdf), str(tmp_path / "out"))
    lig = out["ligands"][0]
    assert lig["hydrogens_added"] is True
    assert lig["n_atoms_pdbqt"] == 4                           # same PDBQT as the explicit-H input


@needs_meeko
def test_a_record_without_3d_coordinates_is_refused_not_skipped(tmp_path):
    from rdkit import Chem
    flat = Chem.MolFromSmiles("CC(=O)O"); flat.SetProp("_Name", "acetic_acid_2d")
    Chem.rdDepictor.Compute2DCoords(flat)
    sdf = tmp_path / "flat.sdf"
    w = Chem.SDWriter(str(sdf)); w.write(flat); w.close()
    with pytest.raises(ValueError, match="no 3D coordinates"):
        process(str(sdf), str(tmp_path / "out"))


@needs_meeko
def test_an_unreadable_sdf_is_refused():
    with pytest.raises(ValueError):
        read_sdf(os.path.join(HERE, "demo_data", "README.md"))


# ---- node wiring (needs the bocoflow_core runtime) -------------------------

def _resolve(value):
    if isinstance(value, str) and value.startswith("demo_data/"):
        return os.path.join(HERE, value)
    return value


class _P:  # minimal Parameter stub; OPTIONS objects are class-level and shared
    def __init__(self, v):
        self._v = v

    def get_value(self):
        return self._v


@needs_runtime
def test_node_instantiation_and_options():
    node = LigandPdbqtMeeko({"node_id": "test"})
    for key in ("input_sdf", "rigid_macrocycles", "hydrate", "flexible_amides", "output_dir"):
        assert key in node.OPTIONS


@needs_runtime
@needs_meeko
def test_execute_end_to_end(tmp_path):
    node = LigandPdbqtMeeko({"node_id": "test"})
    flow_vars = {k: _P(_resolve(v)) for k, v in DEMO_CONFIG.items()}
    flow_vars["output_dir"] = _P(str(tmp_path))
    result = json.loads(node.execute(None, flow_vars))
    assert result["success"] is True
    assert result["data"]["n_out"] == 3
    assert [lig["torsions"] for lig in result["data"]["ligands"]] == [1, 4, 0]
    assert sorted(os.listdir(tmp_path)) == ["aspirin.pdbqt", "caffeine.pdbqt", "ethanol.pdbqt"]


@needs_runtime
@needs_meeko
def test_the_sdf_carried_from_smiles_to_3d_is_used_when_the_parameter_is_blank(tmp_path):
    node = LigandPdbqtMeeko({"node_id": "test"})
    flow_vars = {k: _P(v) for k, v in DEMO_CONFIG.items()}
    flow_vars["input_sdf"] = _P("")
    flow_vars["output_dir"] = _P(str(tmp_path))
    carried = [{"output_file": _resolve("demo_data/ligands.sdf"), "smiles": ["CCO"]}]
    result = json.loads(node.execute(carried, flow_vars))
    assert result["success"] is True and result["data"]["n_out"] == 3
    assert result["data"]["smiles"] == ["CCO"]                  # upstream keys carried forward


@needs_runtime
def test_no_input_fails_loudly():
    node = LigandPdbqtMeeko({"node_id": "test"})
    flow_vars = {k: _P(v) for k, v in DEMO_CONFIG.items()}
    flow_vars["input_sdf"] = _P("")
    flow_vars["output_dir"] = _P(".")
    with pytest.raises(Exception, match="No SDF"):
        node.execute(None, flow_vars)
