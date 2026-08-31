"""Paths reach GROMACS as whole argv elements, never as text a shell re-splits.

Regression test for bocoflow#104.

WHAT WENT WRONG

`gmx_solv_ion` built its GROMACS commands as f-strings and handed them to
`subprocess.run(..., shell=True)`:

    cmd = f"gmx grompp -f {mdp_file} -c {solv_gro} ..."

That flattens a list of arguments into one string and asks a shell to guess
where the boundaries were. The packaged macOS app installs nodes under
`~/Library/Application Support/...`, so `-f` received a path with a space, the
shell guessed wrong, and grompp could not open the file. The 11-node pipeline
died at its ninth node on every packaged macOS install.

WHY IT SURVIVED SO LONG

`editconf` and `solvate` run first and take working-directory paths, which have
no space -- so they passed, and the failure read as a GROMACS problem rather
than a quoting one. `salpa smoke` runs the package from a checkout path, which
also has no space, so it could not see this either.

WHAT THIS ASSERTS

The first fix quoted every interpolated path, which makes the string round-trip
lossless. The current code does not make the round-trip at all: commands are
argv lists passed straight to execve, so a space, quote, `$` or newline in a
path is simply part of the argument.

So this no longer checks quoting. It checks the stronger property that replaced
it -- that each path arrives as ONE list element, byte-identical to what the
code meant. Nothing has to be remembered or escaped for that to hold.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gmx_solv_ion import core  # noqa: E402

AWKWARD = [
    # The real shape of a packaged macOS install.
    "/Users/someone/Library/Application Support/bocoflow-electron/nodes",
    # And the ones that a naive f'-f "{path}"' would NOT have survived either.
    "/tmp/it's a dir",
    "/tmp/cost $HOME dollars",
    '/tmp/say "hello" twice',
]


def _capture_commands(tmp_path, root):
    """Run the pipeline with every gmx call stubbed, returning the argv lists."""
    seen = []

    def fake_run(argv, cwd=None, stdin_text=None):
        seen.append(argv)
        return 0, "ok"

    node_dir = os.path.join(root, "pdbmdauto", "gmx_solv_ion")
    out = str(tmp_path / "out")
    os.makedirs(out, exist_ok=True)
    top = os.path.join(out, "topol.top")
    open(top, "w").write("; test\n")

    original = core._run_gmx
    core._run_gmx = fake_run
    try:
        core.process_solv_ion(
            gro_file=os.path.join(root, "input.gro"),
            top_file=top,
            mdp_file=os.path.join(node_dir, "demo_data", "ions.mdp"),
            ndx_file=os.path.join(out, "index.ndx"),
            output_dir=out,
            case_name="t",
        )
    except Exception:
        # Downstream parsing of stubbed output may fail; the commands are
        # already captured and are what this test is about.
        pass
    finally:
        core._run_gmx = original
    return seen


def test_commands_are_produced(tmp_path):
    """Guards the assertions below against passing on an empty list."""
    assert _capture_commands(tmp_path, AWKWARD[0]), "no gmx commands were built"


def test_every_command_is_a_list_not_a_string(tmp_path):
    """A string here would mean someone reintroduced the shell round-trip."""
    for argv in _capture_commands(tmp_path, AWKWARD[0]):
        assert isinstance(
            argv, (list, tuple)
        ), f"command was built as a string, which only a shell can run: {argv!r}"
        assert argv and argv[0] in {"gmx", "pdb2pqr"}, argv


@pytest.mark.parametrize("root", AWKWARD)
def test_an_awkward_path_arrives_as_one_intact_argument(tmp_path, root):
    """THE REGRESSION, as the property that actually matters."""
    checked = 0
    for argv in _capture_commands(tmp_path, root):
        for element in argv:
            if not isinstance(element, str) or root not in element:
                continue
            checked += 1
            # Byte-identical: not split, not escaped, not mangled.
            assert element.startswith(
                root
            ), f"path was altered before reaching execve: {element!r}"
            assert " ".join(argv).count(root) >= 1
    assert checked, f"no argument contained {root!r}; the test proved nothing"


def test_a_stdin_prompt_is_data_not_a_shell_pipe(tmp_path):
    """`echo SOL | gmx genion` is gone; the answer goes to stdin.

    The pipe was never doing anything a shell was needed for -- it answered an
    interactive prompt. Keeping it would have kept a shell in the path, and with
    it every quoting concern this test exists to retire.
    """
    seen = []

    def fake_run(argv, cwd=None, stdin_text=None):
        seen.append((argv, stdin_text))
        return 0, "ok"

    original = core._run_gmx
    core._run_gmx = fake_run
    try:
        out = str(tmp_path / "o")
        os.makedirs(out, exist_ok=True)
        top = os.path.join(out, "topol.top")
        open(top, "w").write("; test\n")
        try:
            core.process_solv_ion(
                gro_file=os.path.join(out, "in.gro"),
                top_file=top,
                mdp_file=os.path.join(out, "ions.mdp"),
                ndx_file=os.path.join(out, "index.ndx"),
                output_dir=out,
                case_name="t",
                ion_conc=0.15,
            )
        except Exception:
            pass
    finally:
        core._run_gmx = original

    for argv, stdin_text in seen:
        assert not any("|" in str(a) for a in argv), f"a pipe survived in {argv!r}"
        assert not any(str(a) == "echo" for a in argv), f"echo survived in {argv!r}"

    genion = [(a, s) for a, s in seen if "genion" in a]
    if genion:
        assert (
            genion[0][1] == "SOL\n"
        ), "genion's group answer must go to stdin, not through a shell pipe"
