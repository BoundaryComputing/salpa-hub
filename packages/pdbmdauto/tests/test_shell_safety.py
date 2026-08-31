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
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parent.parent

# An f-string that starts with one of these is a command line, not log text.
COMMAND_START = re.compile(r"""^f?['"]\s*(gmx|echo|pdb2pqr|obabel|python)\b""")

# Interpolations that are numbers or enums by construction, never paths.
NON_PATH = re.compile(
    r"^(rc|returncode|stdout|stderr|ph|ion_conc|scale_fill|"
    r"box_distance|box\[\d\]|n_\w+|\w*_conc|\w*_size)$"
)

QUOTED = re.compile(r"\b(_q|quote|shlex\.quote)\s*\(")


def _command_fstrings():
    """Yield (path, lineno, expr) for every interpolation in a command f-string."""
    for py in sorted(PKG.rglob("*.py")):
        if "tests" in py.parts or ".pixi" in py.parts:
            continue
        try:
            tree = ast.parse(py.read_text(encoding="utf-8"))
        except SyntaxError as exc:  # a broken node is a different test's problem
            pytest.fail(f"{py.relative_to(PKG)} does not parse: {exc}")
        for node in ast.walk(tree):
            if not isinstance(node, ast.JoinedStr):
                continue
            if not COMMAND_START.search(ast.unparse(node).strip()):
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
        if "tests" in py.parts or ".pixi" in py.parts:
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
