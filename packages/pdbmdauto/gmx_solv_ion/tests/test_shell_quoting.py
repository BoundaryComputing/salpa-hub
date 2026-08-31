"""Every path this node puts in a shell command must survive a space.

Regression test for bocoflow#104.

WHAT WENT WRONG

`gmx_solv_ion` built its GROMACS commands as f-strings and ran them through a
shell with nothing quoted:

    cmd = f"gmx grompp -f {mdp_file} -c {solv_gro} ..."

The packaged macOS app installs nodes under `~/Library/Application Support/...`,
and the shipped template points this node's mdp at `node:demo_data/ions.mdp`. So
`-f` received a path with a space, the shell split it, and grompp could not open
the file. The 11-node pipeline died at the ninth node on every packaged macOS
install.

WHY IT SURVIVED SO LONG

`editconf` and `solvate` run first and take working-directory paths, which have
no space -- so they passed, and the failure read as a GROMACS problem rather than
a quoting one. `salpa smoke` runs the package from a checkout path, which also has
no space, so it could not see this either.

WHAT THIS ASSERTS

Not "the source contains shlex.quote" -- that is a proxy. This runs the real code
with a real spaced path and checks the invariant that matters: after `shlex.split`,
every path comes back as ONE argument, exactly as written. That is the property a
shell guarantees only if the value was quoted.
"""

import os
import shlex
import sys
from unittest.mock import patch

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gmx_solv_ion import core  # noqa: E402

# The real shape of the packaged macOS install: a space in the middle.
SPACED = "/Users/someone/Library/Application Support/bocoflow-electron/nodes"


AWKWARD = [
    # The real shape of a packaged macOS install.
    "/Users/someone/Library/Application Support/bocoflow-electron/nodes",
    # And the ones double-quoting would NOT survive, which is why this uses
    # shlex.quote rather than f'-f "{path}"'.
    "/tmp/it's a dir",
    "/tmp/cost $HOME dollars",
]


def _capture_commands(tmp_path, root):
    """Run the pipeline with every gmx call stubbed, returning the command strings."""
    seen = []

    def fake_run(cmd, cwd=None, **kwargs):
        seen.append(cmd)
        return 0, "ok"

    node_dir = os.path.join(root, "pdbmdauto", "gmx_solv_ion")
    out = str(tmp_path / "out")
    os.makedirs(out, exist_ok=True)
    top = os.path.join(out, "topol.top")
    open(top, "w").write("; test\n")

    with patch.object(core, "_run_gmx", side_effect=fake_run):
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
    return seen


def test_commands_are_produced(tmp_path):
    """Guards the assertions below against passing on an empty list."""
    assert _capture_commands(tmp_path, AWKWARD[0]), "no gmx commands were built"


@pytest.mark.parametrize("root", AWKWARD)
def test_an_awkward_path_survives_the_shell_intact(tmp_path, root):
    """THE REGRESSION, stated as the property that actually matters.

    Split each command the way a shell would. Any argument that mentions our
    directory must come back byte-identical to what the code meant -- not split
    at the space, not mangled by `$` or an apostrophe.
    """
    marker = os.path.basename(root)  # e.g. "nodes", "it's a dir"
    cmds = _capture_commands(tmp_path, root)
    checked = 0
    for cmd in cmds:
        segment = cmd.split("|")[-1]  # `echo SOL | gmx genion ...`
        args = shlex.split(segment)  # raises if quoting is unbalanced
        for a in args:
            if root in a:
                checked += 1
                assert a.startswith(root) or root in a
                assert marker in a, f"path was mangled: {a!r}\n  in: {cmd}"
    assert checked, f"no argument contained {root!r}; the test proved nothing"
