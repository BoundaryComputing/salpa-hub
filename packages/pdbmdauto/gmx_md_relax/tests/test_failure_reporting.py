"""A failed relaxation step says why, in the UI and in the error it raises.

GROMACS opens its output with a banner of some 500 characters and closes with its error.
The node raised with the log's first 500 characters, so a failure showed the banner and
never the reason. Found running the template on 1AKI inside Salpa (2026-09-30):
mdrun ended with signal 6 at nvt_fixOriBackbone, and nothing on screen said why.

The same wording must still read as a blown-up step to Salpa's pipeline test, which accepts
one "mdrun(<step>) failed (rc=N)" after its own grompp succeeded; MDRUN_FAILED below is its
pattern.
"""

import os
import re
import sys

import pytest

try:
    from bocoflow_core.node import NodeException
except ImportError:
    pytest.skip(
        "bocoflow_core is not installed: run in the package's test environment",
        allow_module_level=True,
    )

PACKAGE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PACKAGE)

from gmx_md_relax import core  # noqa: E402
from gmx_md_relax import node as node_module  # noqa: E402

MDRUN_FAILED = re.compile(
    r"mdrun\(([^)]+)\) failed \(rc=-?\d+\)|mdrun\(([^)]+)\): rc=-?[1-9]\d*|\bmdrun failed:"
)
BANNER = "\n                :-) GROMACS - gmx mdrun, 2026.3-conda_forge (-:\n\n" + (
    "Executable:   /x/.pixi/envs/default/bin.ARM_NEON_ASIMD/gmx\n" * 12
)
REASON = (
    "Fatal error:\nThe domain decomposition grid has shifted too much in the Y-direction "
    "around cell 0 4 0."
)


class Value:
    """A flow variable, as node_runner hands one to execute()."""

    def __init__(self, value):
        self.value = value

    def get_value(self):
        return self.value


def _fake_run(mdrun_rc, mdrun_out):
    def run(argv, cwd=None, timeout=3600, stdin_text=None):
        if argv[1] == "grompp":
            return 0, "grompp ok"
        return mdrun_rc, mdrun_out
    return run


class TestTheStepsMessage:
    def test_it_carries_the_end_of_gromacss_output(self, tmp_path, monkeypatch):
        monkeypatch.setattr(core, "_run", _fake_run(1, BANNER + "step 500\n" + REASON))
        ok, _, log = core._run_grompp_mdrun(
            "x.mdp", "x.gro", "x.top", "x.ndx", str(tmp_path), "nvt_fixOriBackbone"
        )
        assert not ok
        assert log.startswith("mdrun(nvt_fixOriBackbone) failed (rc=1):")
        assert REASON in log

    def test_a_long_output_keeps_its_end(self, tmp_path, monkeypatch):
        monkeypatch.setattr(core, "_run", _fake_run(1, BANNER * 20 + REASON))
        _, _, log = core._run_grompp_mdrun(
            "x.mdp", "x.gro", "x.top", "x.ndx", str(tmp_path), "mm1"
        )
        assert REASON in log
        assert len(log) < 2200

    def test_a_signal_is_named(self, tmp_path, monkeypatch):
        monkeypatch.setattr(core, "_run", _fake_run(-6, BANNER))
        _, _, log = core._run_grompp_mdrun(
            "x.mdp", "x.gro", "x.top", "x.ndx", str(tmp_path), "nvt_fixOriBackbone"
        )
        assert log.startswith(
            "mdrun(nvt_fixOriBackbone) failed (rc=-6), killed by signal 6 (SIGABRT):"
        )

    @pytest.mark.parametrize("rc", [1, -6])
    def test_app_test_20_still_reads_it_as_a_blown_up_step(self, tmp_path, monkeypatch, rc):
        monkeypatch.setattr(core, "_run", _fake_run(rc, BANNER + REASON))
        _, _, log = core._run_grompp_mdrun(
            "x.mdp", "x.gro", "x.top", "x.ndx", str(tmp_path), "mm1"
        )
        match = MDRUN_FAILED.search(log)
        assert match and match.group(1) == "mm1"


def test_the_node_streams_and_raises_the_reason(tmp_path, monkeypatch):
    log = "nvt_fixOri: OK\nmdrun(nvt_fixOriBackbone) failed (rc=-6), killed by signal 6 " \
          "(SIGABRT):\n[...]\n" + "banner line\n" * 60 + REASON
    failed = core.RelaxResult()
    failed.log = log
    monkeypatch.setattr(node_module, "process_full_4step", lambda **kw: failed)
    streamed = []
    import bocoflow_core.stream_logger as stream_logger

    monkeypatch.setattr(stream_logger, "stream_log", lambda msg, **kw: streamed.append((msg, kw)))
    for name in ("ion.gro", "topol.top", "index.ndx"):
        (tmp_path / name).write_text("")

    node = node_module.GmxMdRelax({"node_id": "n1", "node_key": "GmxMdRelax"})
    flow_vars = {
        "protocol": Value("full_4step"),
        "case_name": Value("demo"),
        "input_top_file": Value(""),
        "input_gro_file": Value(""),
        "input_mdp_file": Value(""),
        "input_ndx_file": Value(""),
    }
    with pytest.raises(NodeException) as raised:
        node.execute([{"working_path": str(tmp_path)}], flow_vars)

    assert REASON in str(raised.value)
    assert "nvt_fixOri: OK" in str(raised.value)
    errors = [msg for msg, kw in streamed if kw.get("level") == "error"]
    assert any(REASON in msg for msg in errors), "the UI was never shown the failure"
