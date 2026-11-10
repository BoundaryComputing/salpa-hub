"""A small system's MD run is one rank; a solvated one keeps GROMACS's choice (1.2.10).

On the 10-core Mac Mini, Number of Threads left at 0 had thread-MPI start ten ranks, and the
bundled water box (1,530 atoms, a 2.5 nm cube) cannot be cut into ten domains: mdrun stopped with
"There is no domain decomposition for 10 ranks that is compatible with the given box". GMX MD
Relaxation has run a small system as one rank since 1.2.7, for the same cause (see
SINGLE_RANK_BELOW_ATOMS in core.py).
"""

import os
import subprocess
import sys

PACKAGE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PACKAGE)

from gmx_mdrun_local import core  # noqa: E402

WATER_BOX = os.path.join(PACKAGE, "gmx_mdrun_local", "demo_data", "water_box")


def _gro(tmp_path, natoms, name="in.gro"):
    path = tmp_path / name
    path.write_text(f"title\n{natoms:>5}\n   1.00000   1.00000   1.00000\n")
    return str(path)


def test_the_bundled_water_box_runs_as_one_rank():
    assert core._atom_count(os.path.join(WATER_BOX, "water_box.gro")) == 1530
    assert core.rank_args(os.path.join(WATER_BOX, "water_box.gro")) == ["-ntmpi", "1"]


def test_a_solvated_system_keeps_gromacss_choice(tmp_path):
    assert core.rank_args(_gro(tmp_path, 25355)) == []


def test_an_unreadable_file_keeps_gromacss_choice(tmp_path):
    assert core.rank_args(str(tmp_path / "missing.gro")) == []
    bad = tmp_path / "bad.gro"
    bad.write_text("title\nnot a number\n")
    assert core.rank_args(str(bad)) == []


def _run_with(
    tmp_path, monkeypatch, natoms, num_threads=0, gro_name="in.gro", absolute=True
):
    """run_md_simulation with grompp and mdrun replaced; returns the mdrun argv."""
    for name in ("x.mdp", "x.top"):
        (tmp_path / name).write_text("")
    gro = _gro(tmp_path, natoms, gro_name)
    calls = []

    def run(cmd, cwd=None, capture_output=True, text=True):
        calls.append(cmd)
        if cmd[1] == "grompp":
            (tmp_path / "md.tpr").write_text("")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(core.subprocess, "run", run)
    result = core.run_md_simulation(
        top_file=str(tmp_path / "x.top"),
        gro_file=gro if absolute else gro_name,
        mdp_file=str(tmp_path / "x.mdp"),
        working_dir=str(tmp_path),
        num_threads=num_threads,
    )
    assert result.success, result.message
    return calls[1]


def test_mdrun_is_called_with_it(tmp_path, monkeypatch):
    assert _run_with(tmp_path, monkeypatch, 1530) == [
        "gmx",
        "mdrun",
        "-deffnm",
        "md",
        "-v",
        "-ntmpi",
        "1",
    ]


def test_a_thread_count_is_kept_and_its_threads_are_one_ranks(tmp_path, monkeypatch):
    # -nt 10 alone failed as 0 did; with -ntmpi 1 the ten threads are one rank's OpenMP threads
    assert _run_with(tmp_path, monkeypatch, 1530, num_threads=10) == [
        "gmx",
        "mdrun",
        "-deffnm",
        "md",
        "-nt",
        "10",
        "-v",
        "-ntmpi",
        "1",
    ]


def test_a_relative_structure_is_read_from_the_working_folder(tmp_path, monkeypatch):
    assert _run_with(tmp_path, monkeypatch, 1530, absolute=False)[-2:] == [
        "-ntmpi",
        "1",
    ]


def test_a_solvated_system_is_run_as_before(tmp_path, monkeypatch):
    assert _run_with(tmp_path, monkeypatch, 25355) == [
        "gmx",
        "mdrun",
        "-deffnm",
        "md",
        "-v",
    ]
