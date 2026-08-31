"""No node in this package may interpolate a path into a shell command unquoted.

Regression guard for bocoflow#104, which broke the packaged macOS app for three
releases: `~/Library/Application Support/...` contains a space, the shell split
it, and GROMACS could not open the file.

WHY THIS EXISTS ALONGSIDE THE BEHAVIOURAL TEST

`gmx_solv_ion/tests/test_shell_quoting.py` runs the real code against a real
spaced path -- that is the honest test, and it is the one that proves the fix.
But the same defect was present in FIVE nodes, and a behavioural test only ever
covers the node someone thought to write it for. This walks the whole package
and fails on any command string that interpolates a value without quoting it,
including in a node added tomorrow.

It derives its own file list by walking the tree, so a new node is covered the
moment it lands -- there is nothing to remember to update.
"""

import ast
import re
import sys
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parent.parent

# A command line starts with one of these verbs.
COMMAND_VERB = re.compile(r"^\s*(gmx|echo|pdb2pqr|obabel|python)\b")


def _starts_a_command(node: ast.JoinedStr) -> bool:
    """Read the f-string's own first literal chunk, not its unparsed source.

    Matching against `ast.unparse(node)` looked simpler and was version-dependent.
    An f-string whose interpolation contains ANOTHER f-string round-trips as
    triple-quoted on some CPythons and single-quoted on others, so a regex
    anchored on the leading quote matched on one interpreter and not the next.

    Measured on the exact line this hid -- gmx_md_relax's `gmx grompp`, whose
    `-po` argument nests `f"mdout_{run_label}.mdp"` inside the outer f-string:

        py3.9.6   ->  f'''gmx grompp ...   detector MISSED it
        py3.12.3  ->  f'gmx grompp ...     detector caught it
        py3.13.7  ->  f'''gmx grompp ...   detector MISSED it

    So 3.12 was the outlier, and a full run on 3.13 reported clean while six
    unquoted paths sat in the package. Note that both quote characters being
    present is NOT enough on its own -- the same command without the nested
    f-string unparses identically on all three. The nesting is the trigger.

    The AST itself carries no quoting, so ask it directly.
    """
    for value in node.values:
        if isinstance(value, ast.Constant) and isinstance(value.value, str):
            return bool(COMMAND_VERB.match(value.value))
        return False  # starts with an interpolation: not a literal command
    return False


# Interpolations that are numbers or enums by construction, never paths.
NON_PATH = re.compile(
    r"^(rc|returncode|stdout|stderr|ph|ion_conc|scale_fill|"
    r"box_distance|box\[\d\]|n_\w+|\w*_conc|\w*_size)$"
)

QUOTED = re.compile(r"\b(_q|quote|shlex\.quote)\s*\(")

# Directories that are not this package's source. `.venv` matters: an author who
# builds a virtualenv inside the package would otherwise have this guard scan
# pip's vendored code and report findings there -- which is exactly what happened
# on the first Linux run.
NOT_SOURCE = {
    "tests",
    ".pixi",
    ".venv",
    "venv",
    "env",
    "site-packages",
    "node_modules",
    "build",
    "dist",
    "__pycache__",
}


def _skip(path) -> bool:
    return any(part in NOT_SOURCE for part in path.parts)


def _command_fstrings():
    """Yield (path, lineno, expr) for every interpolation in a command f-string."""
    for py in sorted(PKG.rglob("*.py")):
        if _skip(py):
            continue
        try:
            tree = ast.parse(py.read_text(encoding="utf-8"))
        except SyntaxError as exc:  # a broken node is a different test's problem
            pytest.fail(f"{py.relative_to(PKG)} does not parse: {exc}")
        for node in ast.walk(tree):
            if not isinstance(node, ast.JoinedStr):
                continue
            if not _starts_a_command(node):
                continue
            for value in node.values:
                if isinstance(value, ast.FormattedValue):
                    yield py, node.lineno, ast.unparse(value.value)


def test_the_guard_actually_finds_commands():
    """Without this, an over-narrow regex would make the real test vacuous."""
    found = list(_command_fstrings())
    assert len(found) >= 10, (
        f"only found {len(found)} interpolations in command strings; "
        "COMMAND_START is probably too narrow and the guard below proves nothing"
    )


def test_no_path_reaches_a_shell_unquoted():
    """THE GUARD. Every non-numeric interpolation must go through shlex.quote."""
    offenders = [
        f"  {py.relative_to(PKG)}:{lineno}  {{{expr}}}"
        for py, lineno, expr in _command_fstrings()
        if not QUOTED.search(expr) and not NON_PATH.match(expr)
    ]
    assert not offenders, (
        "these values are interpolated into a shell command without quoting.\n"
        "A path containing a space -- which every packaged macOS install has --\n"
        "will be split by the shell. Wrap them in shlex.quote:\n" + "\n".join(offenders)
    )


def test_shell_true_callers_exist_and_are_known():
    """shell=True is what makes quoting load-bearing; keep the list visible."""
    shelly = []
    for py in sorted(PKG.rglob("*.py")):
        if _skip(py):
            continue
        tree = ast.parse(py.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if "subprocess" not in ast.unparse(node.func):
                continue
            if any(
                k.arg == "shell" and getattr(k.value, "value", None) is True
                for k in node.keywords
            ):
                shelly.append(f"{py.relative_to(PKG)}:{node.lineno}")
    # Not a cap -- a tripwire. If this grows, the new caller needs quoting too,
    # and whoever added it should have read this file.
    assert shelly, "no shell=True callers found; the quoting guard may be misaimed"
    assert len(shelly) <= 8, (
        "new shell=True caller(s) -- confirm every path is shlex.quote'd, then "
        f"raise this number:\n  " + "\n  ".join(shelly)
    )


def test_detection_survives_a_command_holding_both_quote_characters():
    """The guard's own blind spot, pinned.

    An f-string whose interpolation nests ANOTHER f-string round-trips through
    `ast.unparse` with different quoting on different CPythons -- triple-quoted
    on 3.9 and 3.13, single-quoted on 3.12 -- so the original detector, a regex
    over the unparsed source, matched on one interpreter and not the others.

    The effect was not theoretical: this guard passed on 3.13 while the same
    package failed on 3.12, and six unquoted paths in gmx_md_relax survived a
    full local run. A guard whose verdict depends on the interpreter is worse
    than no guard, because it is trusted.
    """
    # The nested f-string in the -po argument is what makes unparse's quoting
    # diverge. Without it this fixture passes on every version and proves nothing.
    src = (
        "import os\n"
        "def f(mdp_file, gro_file, output_dir, run_label):\n"
        "    cmd = (\n"
        "        f'gmx grompp -f \"{mdp_file}\" -c \"{gro_file}\" '\n"
        "        f'-po \"{os.path.join(output_dir, f\"mdout_{run_label}.mdp\")}\" -maxwarn 10'\n"
        "    )\n"
    )
    tree = ast.parse(src)
    found = [
        ast.unparse(v.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.JoinedStr) and _starts_a_command(node)
        for v in node.values
        if isinstance(v, ast.FormattedValue)
    ]
    assert "mdp_file" in found and "gro_file" in found, (
        f"the detector missed a mixed-quote command on Python "
        f"{sys.version_info.major}.{sys.version_info.minor}; found {found}"
    )
