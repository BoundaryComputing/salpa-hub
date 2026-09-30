"""
Tests for the simulation box Solvate & Ionize asks editconf for.

Up to 1.2.5 the template used a fixed 5 nm cube. That leaves 4Z8J, the reference protein
(diameter about 4.4 nm), 0.56 nm from the box along its longest axis, and 0.29 nm once it
rotates, against the 1.0 nm cutoffs of the shipped .mdp files: in a longer run the protein
meets its own periodic image, and a larger protein does not fit at all. The default is now
"auto": a cube sized from the protein's diameter plus a padding on every side, which stays
safe as the protein rotates. Three lengths still give an explicit box.
"""

import importlib.util
from pathlib import Path
from unittest.mock import patch

import pytest

PACKAGE = Path(__file__).resolve().parent.parent.parent
_spec = importlib.util.spec_from_file_location(
    "solv_ion_core", PACKAGE / "gmx_solv_ion" / "core.py"
)
core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(core)


def test_auto_is_a_cube_padded_around_the_protein():
    assert core.editconf_box_args("auto", 1.0) == ["-bt", "cubic", "-d", "1.0", "-c"]


def test_the_padding_is_the_callers():
    assert core.editconf_box_args("auto", 1.5) == ["-bt", "cubic", "-d", "1.5", "-c"]


def test_empty_and_zeros_mean_auto():
    # "0 0 0" was documented as the automatic box; it keeps meaning that.
    assert core.editconf_box_args("", 1.0) == core.editconf_box_args("auto", 1.0)
    assert core.editconf_box_args("0 0 0", 1.0) == core.editconf_box_args("auto", 1.0)


def test_three_lengths_are_an_explicit_box():
    assert core.editconf_box_args("7 7 8.5", 1.0) == ["-box", "7.0", "7.0", "8.5", "-c"]


@pytest.mark.parametrize("bad", ["5 5", "a b c", "5 -5 5", "5 5 5 5"])
def test_anything_else_is_refused(bad):
    with pytest.raises(ValueError, match="Box Size"):
        core.editconf_box_args(bad, 1.0)


@pytest.mark.parametrize("padding", [0, -1.0])
def test_a_padding_that_is_not_positive_is_refused(padding):
    with pytest.raises(ValueError, match="Box Padding"):
        core.editconf_box_args("auto", padding)


def test_the_step_runs_editconf_with_the_padded_box(tmp_path):
    gro = tmp_path / "in.gro"
    top = tmp_path / "topol.top"
    gro.write_text("x\n0\n   1.0   1.0   1.0\n")
    top.write_text("; topology\n")
    commands = []

    def fake_run(cmd, cwd=None, **kwargs):
        commands.append(cmd)
        return 1, "stop after editconf"  # the step stops at the first failure

    with patch.object(core, "_run_gmx", side_effect=fake_run):
        result = core.process_solv_ion(
            gro_file=str(gro),
            top_file=str(top),
            mdp_file="",
            ndx_file="",
            output_dir=str(tmp_path),
            case_name="case",
        )

    assert not result.success
    editconf = commands[0]
    assert editconf[:2] == ["gmx", "editconf"]
    assert editconf[-5:] == ["-bt", "cubic", "-d", "1.0", "-c"]
