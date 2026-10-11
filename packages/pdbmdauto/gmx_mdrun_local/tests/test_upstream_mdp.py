"""GROMACS MD Run takes its .mdp from upstream when its own is left at the default (1.2.9).

An upstream step can offer `mdp_files`, a mapping from a run label to the .mdp for that run.
When the node's Parameters File is left at its default, the run uses the entry for its own Run
Label; a file the person sets always wins, and a label the mapping does not name keeps the
bundled default. So a protein-ligand complex's equilibration ("npt") and production ("md") each
find the file written for them, with no path in the workflow.

The node also hands on what a second run needs, so one MD Run can follow another: its final
structure as `output_gro` (mdrun writes the velocities into it), the topology and index it used,
its folder, and the `mdp_files` mapping it was given.
"""

import json
import os
import sys
from unittest.mock import patch

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

from gmx_mdrun_local import core  # noqa: E402
from gmx_mdrun_local import node as node_module  # noqa: E402

DEFAULT_MDP = node_module.GmxMdRunLocal.OPTIONS["input_mdp_file"].default
MDP_FILES = {"npt": "abs:/case/gmx/npt.mdp", "md": "abs:/case/gmx/production.mdp"}
RELAXED = {
    "case_name": "case",
    "working_path": "abs:/case/gmx",
    "output_gro": "abs:/case/gmx/em.gro",
    "output_top": "abs:/case/gmx/topol.top",
    "output_ndx": "abs:/case/gmx/index.ndx",
}


class Value:
    def __init__(self, value):
        self.value = value

    def get_value(self):
        return self.value


def flow_vars(run_label="md", mdp=DEFAULT_MDP):
    return {
        "case_name": Value(""),
        "run_label": Value(run_label),
        "input_top_file": Value(""),
        "input_gro_file": Value(""),
        "input_mdp_file": Value(mdp),
        "input_ndx_file": Value(""),
        "output_folder": Value(""),
        "num_threads": Value(1),
        "max_warnings": Value(0),
        "verbose": Value(False),
        "force_to_run": Value(False),
    }


def run(predecessor, fv):
    """Execute with GROMACS replaced; returns (the core call's kwargs, the result data)."""
    label = fv["run_label"].get_value()
    simulated = core.SimulationResult(
        success=True,
        message="done",
        tpr_file=f"/case/gmx/{label}.tpr",
        gro_file=f"/case/gmx/{label}.gro",
        xtc_file=f"/case/gmx/{label}.xtc",
    )
    with (
        patch.object(node_module, "check_gromacs_available", return_value=True),
        patch.object(node_module, "run_md_simulation", return_value=simulated) as sim,
        patch.object(
            node_module.GmxMdRunLocal,
            "resolve_path",
            side_effect=lambda p: (
                p[4:] if isinstance(p, str) and p.startswith("abs:") else p
            ),
        ),
    ):
        node = node_module.GmxMdRunLocal({"node_id": "n1", "node_key": "GmxMdRunLocal"})
        node.format_output_path = lambda p: f"abs:{p}"
        data = json.loads(node.execute([predecessor] if predecessor else [], fv))[
            "data"
        ]
    return sim.call_args[1], data


def test_production_takes_the_upstream_file_for_its_label():
    called, _ = run({**RELAXED, "mdp_files": MDP_FILES}, flow_vars("md"))
    assert called["mdp_file"] == "/case/gmx/production.mdp"


def test_equilibration_takes_the_file_for_its_own_label():
    called, _ = run({**RELAXED, "mdp_files": MDP_FILES}, flow_vars("npt"))
    assert called["mdp_file"] == "/case/gmx/npt.mdp"


def test_a_file_the_person_sets_wins():
    called, _ = run(
        {**RELAXED, "mdp_files": MDP_FILES}, flow_vars("md", "abs:/mine/md.mdp")
    )
    assert called["mdp_file"] == "/mine/md.mdp"


def test_an_unnamed_label_keeps_the_bundled_default():
    called, _ = run({**RELAXED, "mdp_files": MDP_FILES}, flow_vars("nvt"))
    assert called["mdp_file"] == DEFAULT_MDP


def test_without_a_mapping_nothing_changes():
    called, data = run(RELAXED, flow_vars("md"))
    assert called["mdp_file"] == DEFAULT_MDP
    assert "mdp_files" not in data


def test_the_run_hands_on_what_a_second_run_needs():
    _, data = run({**RELAXED, "mdp_files": MDP_FILES}, flow_vars("npt"))
    assert data["output_gro"] == "abs:/case/gmx/npt.gro"
    assert data["output_top"] == "abs:/case/gmx/topol.top"
    assert data["output_ndx"] == "abs:/case/gmx/index.ndx"
    assert data["working_path"] == "abs:/case/gmx"
    assert data["mdp_files"] == MDP_FILES


def test_a_second_run_follows_the_first():
    _, first = run({**RELAXED, "mdp_files": MDP_FILES}, flow_vars("npt"))
    called, second = run(first, flow_vars("md"))
    assert called["gro_file"] == "/case/gmx/npt.gro"
    assert called["top_file"] == "/case/gmx/topol.top"
    assert called["ndx_file"] == "/case/gmx/index.ndx"
    assert called["mdp_file"] == "/case/gmx/production.mdp"
    assert second["output_files"]["trajectory"] == "abs:/case/gmx/md.xtc"
