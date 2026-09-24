"""A failed step reaches the UI, and the node raises with its reason.

For a subprocess node the UI shows what stream_log() sends and the error the
node raises. log_message() reaches only the log file, and that is where the
step's log used to go. The error quoted the log's first 500 characters, while
PDB2PQR and GROMACS both print the reason last.
"""

import os
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

from pka_gmx_em import core  # noqa: E402
from pka_gmx_em import node as node_module  # noqa: E402

DEMO_PDB = os.path.join(PACKAGE, "pdb2pqr", "demo_data", "mini.pdb")
REASON = "UnicodeDecodeError: 'utf-8' codec can't decode byte 0xd6 in position 19"


class Value:
    """A flow variable, as node_runner hands one to execute()."""

    def __init__(self, value):
        self.value = value

    def get_value(self):
        return self.value


def flow_vars(tmp_path):
    return {
        "case_name": Value("demo"),
        "input_pdb": Value(DEMO_PDB),
        "output_dir": Value(str(tmp_path / "gmx")),
        "force_field": Value("amber99sb"),
        "water_model": Value("tip3p"),
        "box_distance": Value(2.0),
        "em_steps": Value(50),
        "ph": Value(7.0),
        "run_pdb2pqr": Value(True),
    }


def test_a_failed_step_is_streamed_and_raised_with_its_reason(tmp_path, monkeypatch):
    # A long log with the reason at the end, as a tool's output has it.
    log = "pdb2pqr: FAILED (rc=1)\n" + "INFO: banner line\n" * 60 + REASON
    failed = core.PkaGmxEmResult(success=False, log=log)
    monkeypatch.setattr(node_module, "process_pka_gmx_em", lambda **kw: failed)
    streamed = []
    monkeypatch.setattr(
        node_module, "stream_log", lambda msg, **kw: streamed.append((msg, kw))
    )

    node = node_module.PkaGmxEm({"node_id": "n1", "node_key": "PkaGmxEm"})
    with pytest.raises(NodeException) as raised:
        node.execute([{}], flow_vars(tmp_path))

    assert REASON in str(raised.value)
    errors = [msg for msg, kw in streamed if kw.get("level") == "error"]
    assert any(REASON in msg for msg in errors), "the UI was never shown the failure"
