"""A small system's relaxation runs as one rank; a solvated one keeps GROMACS's choice.

On a 10-core Mac GROMACS split 1AKI in vacuum (1,960 atoms) into ten domains, and such
runs failed now and then with a domain decomposition error that one rank cannot raise
(see SINGLE_RANK_BELOW_ATOMS in core.py for the cases and the timings).
"""

import os
import sys

PACKAGE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PACKAGE)

from gmx_md_relax import core  # noqa: E402


def _gro(tmp_path, natoms, name="in.gro"):
    path = tmp_path / name
    path.write_text(f"title\n{natoms:>5}\n   1.00000   1.00000   1.00000\n")
    return str(path)


def test_a_vacuum_protein_runs_as_one_rank(tmp_path):
    assert core._rank_args(_gro(tmp_path, 1960)) == ["-ntmpi", "1"]


def test_a_solvated_system_keeps_gromacss_choice(tmp_path):
    assert core._rank_args(_gro(tmp_path, 25355)) == []


def test_an_unreadable_file_keeps_gromacss_choice(tmp_path):
    assert core._rank_args(str(tmp_path / "missing.gro")) == []
    bad = tmp_path / "bad.gro"
    bad.write_text("title\nnot a number\n")
    assert core._rank_args(str(bad)) == []


def test_mdrun_is_called_with_it(tmp_path, monkeypatch):
    calls = []

    def run(argv, cwd=None, timeout=3600, stdin_text=None):
        calls.append(argv)
        if argv[1] == "mdrun":
            (tmp_path / "nvt_fixOri.gro").write_text("")
        return 0, ""

    monkeypatch.setattr(core, "_run", run)
    ok, _, _ = core._run_grompp_mdrun(
        "x.mdp", _gro(tmp_path, 1960), "x.top", "x.ndx", str(tmp_path), "nvt_fixOri"
    )
    assert ok
    assert calls[1] == ["gmx", "mdrun", "-v", "-deffnm", "nvt_fixOri", "-ntmpi", "1"]
