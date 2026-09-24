"""A failed PDB2PQR fails the pKa + GROMACS EM step.

WHAT WENT WRONG

With "Run PDB2PQR" on, which is the default, `process_pka_gmx_em` treated a
PDB2PQR failure as non-fatal. It gave pdb2gmx the unprotonated structure, let
pdb2gmx's own hydrogen analysis choose the states, and the node reported
success with "0 protonation change(s)". The protonation was not the one asked
for. The only record was a log line that went to log_message(), which the UI
does not show, and it quoted the first 400 characters of PDB2PQR's output:
the banner, not the error, which PDB2PQR prints last.

Seen for real on 2026-09-24: one Latin-1 byte in a REMARK line makes PDB2PQR
3.6.1 exit 1 with a UnicodeDecodeError, and GROMACS accepts the same file.
"""

import os
import shutil
import sys

import pytest

PACKAGE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PACKAGE)

from pka_gmx_em import core  # noqa: E402

DEMO_PDB = os.path.join(PACKAGE, "pdb2pqr", "demo_data", "mini.pdb")

# What PDB2PQR prints when it cannot read a file: its banner first, the reason last.
BANNER = "INFO:PDB2PQR v3.6.1: biomolecular structure conversion software.\n" * 8
ERROR = "UnicodeDecodeError: 'utf-8' codec can't decode byte 0xd6 in position 19"


def run_with_stubs(tmp_path, monkeypatch, run_pdb2pqr=True):
    """Run the step with PDB2PQR failing and every gmx call stubbed.

    Returns the result and the argv of each command, in order.
    """
    calls = []

    def fake_run(argv, cwd=None, timeout=300, stdin_text=None):
        calls.append(argv)
        if argv[0] == "pdb2pqr":
            return 1, BANNER + ERROR
        return 0, ""

    monkeypatch.setattr(core, "_run", fake_run)
    result = core.process_pka_gmx_em(
        DEMO_PDB, str(tmp_path / "gmx"), "demo", run_pdb2pqr=run_pdb2pqr
    )
    return result, calls


def test_a_failed_pdb2pqr_stops_the_step_before_pdb2gmx(tmp_path, monkeypatch):
    result, calls = run_with_stubs(tmp_path, monkeypatch)
    ran = [" ".join(argv[:2]) for argv in calls]
    assert ran == ["pdb2pqr --ff"], f"went on to run {ran[1:]} after PDB2PQR failed"
    assert result.success is False


def test_the_log_carries_pdb2pqr_own_error(tmp_path, monkeypatch):
    result, _ = run_with_stubs(tmp_path, monkeypatch)
    assert "pdb2pqr: FAILED (rc=1)" in result.log
    assert ERROR in result.log
    # And it names the way to choose pdb2gmx's own states deliberately.
    assert "Run PDB2PQR" in result.log


def test_turning_pdb2pqr_off_still_hands_pdb2gmx_the_input(tmp_path, monkeypatch):
    """The opt-out is unchanged: no PDB2PQR run, and pdb2gmx gets the input."""
    _, calls = run_with_stubs(tmp_path, monkeypatch, run_pdb2pqr=False)
    assert calls[0][:4] == ["gmx", "pdb2gmx", "-f", DEMO_PDB]
    assert all(argv[0] != "pdb2pqr" for argv in calls)


@pytest.mark.skipif(
    shutil.which("gmx") is None or os.name == "nt",
    reason="needs GROMACS: run in the package's environment",
)
def test_with_gromacs_no_topology_is_written_from_an_unprotonated_input(
    tmp_path, monkeypatch
):
    """The same, with the real gmx and a pdb2pqr that fails as 3.6.1 does."""
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake = fake_bin / "pdb2pqr"
    fake.write_text(f'#!/bin/sh\necho "{ERROR}" >&2\nexit 1\n')
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")

    out = tmp_path / "gmx"
    result = core.process_pka_gmx_em(DEMO_PDB, str(out), "demo", em_steps=50)

    assert not (out / "pdb2gmx.top").exists(), "pdb2gmx ran without PDB2PQR's states"
    assert result.success is False
    assert ERROR in result.log
