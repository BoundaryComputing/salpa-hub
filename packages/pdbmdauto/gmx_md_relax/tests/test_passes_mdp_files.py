"""Relaxation passes on the .mdp files an upstream step offers for later runs (pdbmdauto 1.2.9).

A step before the relaxations can offer `mdp_files`, a mapping from a run label to the .mdp
for that run (a step that adds a ligand can offer one for "npt" and one for "md"). GROMACS
MD Run takes the entry for its own run label, so every step between them hands the mapping on
unchanged. Relaxation reads only its own predecessor's data, so it must pass it on itself.
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

from gmx_md_relax import core  # noqa: E402
from gmx_md_relax import node as node_module  # noqa: E402

MDP_FILES = {"npt": "abs:/case/gmx/npt.mdp", "md": "abs:/case/gmx/production.mdp"}


class Value:
    def __init__(self, value):
        self.value = value

    def get_value(self):
        return self.value


def run(tmp_path, monkeypatch, upstream):
    done = core.RelaxResult()
    done.success = True
    done.output_gro = str(tmp_path / "em.gro")
    done.steps_completed = ["em"]
    done.max_force = 900.0
    done.em_safe = True
    monkeypatch.setattr(node_module, "process_single_step", lambda **kw: done)
    for name in ("ion.gro", "topol.top", "index.ndx"):
        (tmp_path / name).write_text("")
    node = node_module.GmxMdRelax({"node_id": "n1", "node_key": "GmxMdRelax"})
    node.format_output_path = lambda p: p
    flow_vars = {
        "protocol": Value("em_only"),
        "case_name": Value("demo"),
        "input_top_file": Value(""),
        "input_gro_file": Value(""),
        "input_mdp_file": Value(""),
        "input_ndx_file": Value(""),
    }
    return json.loads(
        node.execute([{"working_path": str(tmp_path), **upstream}], flow_vars)
    )


def test_the_mapping_is_passed_on_unchanged(tmp_path, monkeypatch):
    data = run(tmp_path, monkeypatch, {"mdp_files": MDP_FILES})["data"]
    assert data["mdp_files"] == MDP_FILES


def test_nothing_is_added_when_nothing_was_offered(tmp_path, monkeypatch):
    data = run(tmp_path, monkeypatch, {})["data"]
    assert "mdp_files" not in data
