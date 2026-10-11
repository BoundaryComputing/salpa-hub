"""Solvate & Ionize passes on the .mdp files an upstream step offers (pdbmdauto 1.2.9).

See gmx_md_relax/tests/test_passes_mdp_files.py: `mdp_files` maps a run label to the .mdp for
that run, and GROMACS MD Run takes the entry for its own label, so every step in between hands
the mapping on unchanged.
"""

import json
import os
import sys

import pytest

try:
    from bocoflow_core.node import NodeException  # noqa: F401
except ImportError:
    pytest.skip(
        "bocoflow_core is not installed: run in the package's test environment",
        allow_module_level=True,
    )

PACKAGE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PACKAGE)

from gmx_solv_ion import core  # noqa: E402
from gmx_solv_ion import node as node_module  # noqa: E402

MDP_FILES = {"npt": "abs:/case/gmx/npt.mdp", "md": "abs:/case/gmx/production.mdp"}


class Value:
    def __init__(self, value):
        self.value = value

    def get_value(self):
        return self.value


def run(tmp_path, monkeypatch, upstream):
    for name in ("complex.gro", "topol.top", "index.ndx"):
        (tmp_path / name).write_text("")
    done = core.SolvIonResult(
        output_gro=str(tmp_path / "ion.gro"),
        output_top=str(tmp_path / "topol.top"),
        output_ndx=str(tmp_path / "index.ndx"),
        success=True,
    )
    monkeypatch.setattr(node_module, "process_solv_ion", lambda **kw: done)
    node = node_module.GmxSolvIon({"node_id": "n1", "node_key": "GmxSolvIon"})
    node.format_output_path = lambda p: p
    flow_vars = {
        "case_name": Value("demo"),
        "run_label": Value("md"),
        "force_field": Value("amber99sb"),
        "water_model": Value("tip3p"),
        "box_size": Value("auto"),
        "box_padding": Value(1.0),
        "ion_conc": Value(0.15),
        "input_top_file": Value(""),
        "input_gro_file": Value(""),
        "input_mdp_file": Value(str(tmp_path / "ions.mdp")),
    }
    (tmp_path / "ions.mdp").write_text("")
    predecessor = {
        "working_path": str(tmp_path),
        "output_gro": str(tmp_path / "complex.gro"),
        "output_top": str(tmp_path / "topol.top"),
        "output_ndx": str(tmp_path / "index.ndx"),
        **upstream,
    }
    return json.loads(node.execute([predecessor], flow_vars))


def test_the_mapping_is_passed_on_unchanged(tmp_path, monkeypatch):
    data = run(tmp_path, monkeypatch, {"mdp_files": MDP_FILES})["data"]
    assert data["mdp_files"] == MDP_FILES


def test_nothing_is_added_when_nothing_was_offered(tmp_path, monkeypatch):
    data = run(tmp_path, monkeypatch, {})["data"]
    assert "mdp_files" not in data
