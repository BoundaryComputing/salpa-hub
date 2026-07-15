"""Tiny shared helpers for md_solvate_gmx — predecessor lookup +
AMBERHOME detection. Identical-in-spirit to md_solvate_packmol/helpers.py;
see dev-notes/node-io-helpers-duplication.md for why we follow the
per-node pattern instead of promoting these to bocoflow_core."""
from __future__ import annotations

import os
import shutil
import sys
from typing import Any, Iterable, Optional


def get_from_predecessors(predecessor_data: Optional[Iterable], key: str) -> Any:
    """First-match lookup across predecessor dicts."""
    for pred in (predecessor_data or []):
        if pred and key in pred:
            return pred[key]
    return None


def detect_amberhome() -> str:
    """Locate AMBERHOME from env or via the tleap binary's path."""
    if os.environ.get("AMBERHOME"):
        return os.environ["AMBERHOME"]
    bin_path = shutil.which("tleap")
    if bin_path:
        # <env>/bin/tleap → <env>
        return os.path.dirname(os.path.dirname(bin_path))
    # Last resort: sys.prefix (works for pixi/conda envs)
    return sys.prefix
