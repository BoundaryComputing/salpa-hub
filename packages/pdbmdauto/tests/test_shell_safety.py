"""No node in this package may hand a command to a shell.

Regression guard for bocoflow#104, which broke the packaged macOS app for three
releases.

THE INVARIANT CHANGED, AND GOT STRONGER

The first fix quoted every path interpolated into a shell command, and this file
checked that quoting. That check had two weaknesses, and both bit:

  1. It could only ever be as good as its idea of "looks like a command". The
     first version regexed `ast.unparse()` output, whose quoting style varies by
     CPython version, so it passed on 3.9 and 3.13 while failing on 3.12 -- and
     hid six unquoted paths in gmx_md_relax.
  2. Quoting is a thing an author must remember on every new line, forever.

The package no longer builds command strings at all. Commands are argv lists
passed to execve, so a space, quote, `$` or newline in a path is simply part of
the argument. That makes the invariant binary and unmissable, which is what this
file now asserts: NO `shell=True`, anywhere.

There is nothing to pattern-match and nothing to get subtly wrong. A future
author who reaches for a shell trips this immediately, with the reason attached.
"""

import ast
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parent.parent

# Directories that are not this package's source. `.venv` matters: an author who
# builds a virtualenv inside the package would otherwise have this guard scan
# pip's vendored code and report findings there -- which is what happened on the
# first Linux run of the previous version.
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


def _sources():
    for py in sorted(PKG.rglob("*.py")):
        if not any(part in NOT_SOURCE for part in py.parts):
            yield py


def _subprocess_calls():
    """Yield (path, lineno, call) for every subprocess.* call in the package."""
    for py in _sources():
        try:
            tree = ast.parse(py.read_text(encoding="utf-8"))
        except SyntaxError as exc:
            pytest.fail(f"{py.relative_to(PKG)} does not parse: {exc}")
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and "subprocess" in ast.unparse(node.func):
                yield py, node.lineno, node


def test_the_guard_actually_finds_subprocess_calls():
    """Without this, a broken walker would make the real test vacuous."""
    found = list(_subprocess_calls())
    assert len(found) >= 5, (
        f"only found {len(found)} subprocess calls; the walker is probably "
        "broken and the assertion below proves nothing"
    )


def test_nothing_runs_through_a_shell():
    """THE GUARD. A shell re-splits a flattened command; argv never does."""
    offenders = []
    for py, lineno, call in _subprocess_calls():
        for kw in call.keywords:
            if kw.arg == "shell" and getattr(kw.value, "value", None) is True:
                offenders.append(f"  {py.relative_to(PKG)}:{lineno}")
    assert not offenders, (
        "these calls run a command through a shell:\n"
        + "\n".join(offenders)
        + "\n\nBuild the command as an argv list instead and drop shell=True.\n"
        "A shell flattens the arguments and re-splits them by guessing where\n"
        "the boundaries were; every packaged macOS install has a space in its\n"
        "node path, so it guesses wrong (bocoflow#104). If you needed the shell\n"
        "for `echo X | cmd`, pass input='X\\n' instead -- that pipe only ever\n"
        "answered an interactive prompt."
    )


def test_no_command_is_built_by_string_concatenation():
    """The shape that precedes a shell: an f-string that looks like a command.

    Catches a reintroduction one step earlier than shell=True does -- someone
    assembling `f"gmx grompp -f {path}"` has already lost the argument
    boundaries, whether or not they have reached for a shell yet.
    """
    import re

    # A command line, not log text. Requiring a FLAG as well as the verb is what
    # separates `f"gmx grompp -f {mdp}"` from `f"gmx make_ndx failed (rc={rc})"`,
    # which the verb alone does not -- this check reported all four log messages
    # in the package on its first run.
    #
    # This is a heuristic, and that is precisely why it is the SECONDARY guard.
    # test_nothing_runs_through_a_shell above is binary and cannot be fooled;
    # this one exists only to catch a reintroduction one step earlier.
    verb = re.compile(r"^\s*(gmx|pdb2pqr|obabel)\s+\S+\s+-{1,2}[a-zA-Z]")
    offenders = []
    for py in _sources():
        for node in ast.walk(ast.parse(py.read_text(encoding="utf-8"))):
            if not isinstance(node, ast.JoinedStr):
                continue
            # Read the f-string's own first literal chunk from the AST. Do NOT
            # regex ast.unparse() output: CPython's quoting choice varies by
            # version, which is exactly how the previous guard passed on 3.13
            # and failed on 3.12 against the same file.
            first = next(
                (
                    v.value
                    for v in node.values
                    if isinstance(v, ast.Constant) and isinstance(v.value, str)
                ),
                None,
            )
            if first and verb.match(first):
                offenders.append(
                    f"  {py.relative_to(PKG)}:{node.lineno}  {first[:50]!r}"
                )
    assert not offenders, (
        "these build a command line as a string:\n"
        + "\n".join(offenders)
        + '\n\nUse a list: ["gmx", "grompp", "-f", mdp_file, ...]'
    )
