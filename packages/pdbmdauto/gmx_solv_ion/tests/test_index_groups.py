"""
The index groups Solvate & Ionize keeps from the index it was given (pdbmdauto 1.2.9).

After genion the step makes a fresh index with make_ndx's default groups, which now name the
water and the ions, and appends what the given index added to those defaults. Up to 1.2.8 that
was OriHeavy and OriBackBone and nothing else, so a group made upstream for the solute, such as
a protein-ligand complex's Protein_LIG, did not reach the MD run, and a production run could
not couple the protein with its ligand. Now every group the fresh defaults do not make again is
kept, in the order it was given. Solvent and ions are added after the solute, so a group of
solute atoms keeps its atom numbers.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parent.parent.parent
_spec = importlib.util.spec_from_file_location(
    "solv_ion_core_groups", PACKAGE / "gmx_solv_ion" / "core.py"
)
core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(core)

VACUUM_PROTEIN = """[ System ]
   1    2    3    4    5    6
[ Protein ]
   1    2    3    4    5    6
[ OriHeavy ]
   1    3    5
[ OriBackBone ]
   1    2
"""

VACUUM_COMPLEX = """[ System ]
   1    2    3    4    5    6    7    8
[ Protein ]
   1    2    3    4    5    6
[ Other ]
   7    8
[ LIG ]
   7    8
[ Protein_LIG ]
   1    2    3    4    5    6    7    8
[ OriHeavy ]
   1    3    5    7
[ OriBackBone ]
   1    2
"""

# what make_ndx writes for the solvated complex: the defaults, waters and ions included
SOLVATED_DEFAULTS = """[ System ]
   1    2    3    4    5    6    7    8    9   10   11   12
[ Protein ]
   1    2    3    4    5    6
[ Other ]
   7    8
[ LIG ]
   7    8
[ NA ]
  12
[ Water ]
   9   10   11
[ SOL ]
   9   10   11
[ non-Water ]
   1    2    3    4    5    6    7    8   12
[ Ion ]
  12
[ Water_and_ions ]
   9   10   11   12
"""


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_groups_are_read_with_their_atoms_in_order(tmp_path):
    groups = core._read_groups(write(tmp_path, "a.ndx", VACUUM_COMPLEX))
    assert [name for name, _ in groups] == [
        "System",
        "Protein",
        "Other",
        "LIG",
        "Protein_LIG",
        "OriHeavy",
        "OriBackBone",
    ]
    assert dict(groups)["OriHeavy"] == "[ OriHeavy ]\n   1    3    5    7\n"


def test_a_protein_keeps_exactly_what_it_kept_before(tmp_path):
    # the groups 1.2.8 kept, byte for byte: OriHeavy and OriBackBone
    old = core._read_groups(write(tmp_path, "a.ndx", VACUUM_PROTEIN))
    defaults = core._group_names(write(tmp_path, "b.ndx", SOLVATED_DEFAULTS))
    assert core._carried_groups(old, defaults) == (
        "[ OriHeavy ]\n   1    3    5\n[ OriBackBone ]\n   1    2\n"
    )


def test_a_complex_keeps_its_protein_ligand_group(tmp_path):
    old = core._read_groups(write(tmp_path, "a.ndx", VACUUM_COMPLEX))
    defaults = core._group_names(write(tmp_path, "b.ndx", SOLVATED_DEFAULTS))
    carried = core._carried_groups(old, defaults)
    assert carried == (
        "[ Protein_LIG ]\n   1    2    3    4    5    6    7    8\n"
        "[ OriHeavy ]\n   1    3    5    7\n"
        "[ OriBackBone ]\n   1    2\n"
    )


def test_a_default_group_is_made_again_not_carried(tmp_path):
    old = core._read_groups(write(tmp_path, "a.ndx", VACUUM_COMPLEX))
    defaults = core._group_names(write(tmp_path, "b.ndx", SOLVATED_DEFAULTS))
    carried = core._carried_groups(old, defaults)
    for made_again in ("[ System ]", "[ Protein ]", "[ Other ]", "[ LIG ]"):
        assert made_again not in carried


def test_no_index_carries_nothing(tmp_path):
    assert core._read_groups(str(tmp_path / "none.ndx")) == []
    assert core._carried_groups([], ["System"]) == ""


def test_the_step_appends_the_carried_groups_after_the_defaults(tmp_path, monkeypatch):
    """The whole step with GROMACS replaced: make_ndx writes the defaults over the given
    index, and the step appends Protein_LIG, OriHeavy and OriBackBone after them."""
    gro = write(tmp_path, "in.gro", "x\n0\n   1.0   1.0   1.0\n")
    top = write(tmp_path, "topol.top", "; topology\n")
    mdp = write(tmp_path, "ions.mdp", "; ions\n")
    ndx = write(tmp_path, "index.ndx", VACUUM_COMPLEX)

    def fake_gmx(argv, cwd=None, stdin_text=None):
        tool = argv[1]
        if tool == "editconf":
            Path(argv[argv.index("-o") + 1]).write_text("x\n0\n   5.0   5.0   5.0\n")
        elif tool in ("solvate", "genion"):
            Path(argv[argv.index("-o") + 1]).write_text("x\n0\n   5.0   5.0   5.0\n")
        elif tool == "make_ndx":
            Path(argv[argv.index("-o") + 1]).write_text(SOLVATED_DEFAULTS)
        return 0, ""

    monkeypatch.setattr(core, "_run_gmx", fake_gmx)
    result = core.process_solv_ion(gro, top, mdp, ndx, str(tmp_path), "case")
    assert result.success, result.log
    names = core._group_names(result.output_ndx)
    assert names == [
        "System",
        "Protein",
        "Other",
        "LIG",
        "NA",
        "Water",
        "SOL",
        "non-Water",
        "Ion",
        "Water_and_ions",
        "Protein_LIG",
        "OriHeavy",
        "OriBackBone",
    ]
    assert "Protein_LIG, OriHeavy, OriBackBone" in result.log
