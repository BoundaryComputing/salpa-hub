"""Pure-Python helpers for md_solvate_packmol.

Kept separate from node.py so they can be unit-tested without
importing bocoflow_core (which is only available inside the runtime
pixi environment, not in the default test env).
"""

from __future__ import annotations

import os
import re
import shutil
import sys


def get_from_predecessors(predecessor_data, key):
    """First-match lookup across predecessor dicts."""
    for pred in (predecessor_data or []):
        if pred and key in pred:
            return pred[key]
    return None


# Note: an earlier draft of this module also exported `ensure_in_workdir`
# for copying resolved input files into the node's output_dir before
# invoking external tools. That turned out to be unnecessary —
# packmol-memgen (and tleap / parmed / ORCA / antechamber / xtb / GROMACS
# subcommands generally) all accept absolute paths for their input file
# arguments. After `Node.resolve_path()` returns a real filesystem path,
# you can pass it straight to the tool. See
# `dev-notes/node-io-helpers-duplication.md` for the discussion that led
# to removing the helper.


def detect_amberhome():
    """Locate AMBERHOME from either an env var or via the packmol-memgen
    binary's path. In the metalparm_vwf pixi env, AMBERHOME is the env
    root (two dirs above bin/packmol-memgen).
    """
    if os.environ.get("AMBERHOME"):
        return os.environ["AMBERHOME"]
    bin_path = shutil.which("packmol-memgen")
    if bin_path:
        return os.path.dirname(os.path.dirname(bin_path))
    # Last resort: sys.prefix (works for pixi/conda envs)
    return sys.prefix


def parse_solvent_counts(log_text):
    """Parse molecule counts per solvent + ion from packmol-memgen's log.

    The log emits lines like:
        Solvent: MOH    | molecules:    3036
        Solvent: WAT    | molecules:    4540
        Adding K+    : 1
    (Exact format varies by packmol-memgen version; we use regex tolerant
    matchers and fall back to None per-component if not found.)

    Returns: dict of {code: count}. Empty dict if nothing matches.
    """
    counts = {}
    # Solvent lines: tolerate "Solvent: <CODE>" with various surrounding
    # whitespace and the "molecules" word optionally followed by ":".
    # Leading whitespace allowed (logs may be indented).
    for m in re.finditer(
        r"^\s*[Ss]olvent[: \t]+(\w+)[^0-9\n]*?(\d+)\s*$",
        log_text, re.MULTILINE,
    ):
        counts[m.group(1).upper()] = int(m.group(2))
    # Ions: "Adding K+    : 1" or "Adding Cl-    : 5". The ion code must
    # start with an uppercase letter (so "Adding solvent layer..." is
    # filtered out — that doesn't end with a digit on the same line).
    # Leading whitespace allowed.
    for m in re.finditer(
        r"^\s*[Aa]dding\s+([A-Z][a-z]?[+-])\s*:?\s*(\d+)\s*$",
        log_text, re.MULTILINE,
    ):
        counts[m.group(1)] = int(m.group(2))
    return counts


def parse_rst7_box(rst7_path):
    """Pull the box vector from an AMBER .rst7 / .ncrst file.

    .rst7 ASCII layout:
        Line 1: title
        Line 2: natom (int) [optional time (float)]
        Lines 3..3+ceil(natom/2)-1: coords, 6 floats per line (2 atoms × 3)
        Line 3+ceil(natom/2): box `lx ly lz alpha beta gamma` (if periodic)

    The box line is at a deterministic offset given natom; a tail line
    with 6 numbers in a non-periodic file is a coordinate row, NOT a box.

    Returns (lx, ly, lz) in Å, or None if no box present / parsing fails.
    """
    try:
        with open(rst7_path) as f:
            lines = [ln for ln in f.read().splitlines() if ln.strip()]
        if len(lines) < 3:
            return None
        # Parse natom from line 2 (first token).
        natom = int(lines[1].split()[0])
        # Expected number of coord lines (2 atoms per line, round up).
        n_coord_lines = (natom + 1) // 2
        box_idx = 2 + n_coord_lines  # 0-indexed
        if box_idx >= len(lines):
            return None  # no box line present
        toks = lines[box_idx].split()
        if len(toks) != 6:
            return None
        return tuple(float(t) for t in toks[:3])
    except (OSError, ValueError, IndexError):
        pass
    return None
