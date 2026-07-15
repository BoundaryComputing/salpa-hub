"""Pure-Python helper for AMBER → GROMACS format conversion.

Wraps ParmEd's `Structure.save()` to emit the canonical pair
(`<prefix>.top`, `<prefix>.gro`) that GROMACS' `grompp` expects. ParmEd
ships with AmberTools and is the upstream-supported conversion path —
easyPARM's own ``amber_converter.py:71-83`` uses the exact same calls.

Why this is in its own helper instead of inline in node.py:
  - importable from tests without the bocoflow_core dependency
  - lets us add a tiny round-trip sanity check (atom count match)
    that would otherwise need a full Node test harness
"""
from __future__ import annotations

import os
from typing import Optional


def normalize_itp_basename(name: str) -> str:
    """Validate + normalize a user-supplied ITP filename.

    The .itp must live next to the .top, so we accept a basename only —
    no path separators. A trailing ``.itp`` is stripped before re-
    appending so ``"complex"`` and ``"complex.itp"`` both resolve to
    ``"complex.itp"``.

    Raises ValueError on a path-separator or empty input.
    """
    n = (name or "").strip()
    if not n:
        raise ValueError("ITP filename is empty")
    if "/" in n or "\\" in n:
        raise ValueError(
            f"ITP filename must be a basename (no path separators): {name!r}"
        )
    if n.lower().endswith(".itp"):
        n = n[:-4]
    if not n:
        raise ValueError(f"ITP filename has no stem: {name!r}")
    return f"{n}.itp"


def _split_top_into_itp(top_path: str, itp_path: str) -> None:
    """Move the [ moleculetype ] block out of <top_path> into <itp_path>;
    replace it in <top_path> with ``#include "<basename>.itp"``.

    Block boundaries (directive-aware, comment- and whitespace-tolerant):
      start = first line whose stripped non-comment content is exactly
              ``[ moleculetype ]`` (with any internal whitespace).
      end   = first line that is the ``[ system ]`` directive after
              start. The moleculetype block owns everything from its
              directive line up to (but not including) ``[ system ]``;
              ParmEd always emits ``[ system ]`` once after the molecule
              block.

    The ``[ atomtypes ]`` block (which appears BEFORE ``[ moleculetype ]``
    in ParmEd output) stays in the .top — a stand-alone .itp can't
    redeclare atomtypes (GROMACS requires them at top-level scope before
    any moleculetype directive), so leaving them in the master .top
    keeps the .itp portable as an ``#include``-able fragment.
    """
    with open(top_path, "r") as fh:
        lines = fh.readlines()

    def directive_of(line: str) -> Optional[str]:
        # Strip leading whitespace; ignore comment-only lines (";...").
        stripped = line.lstrip()
        if not stripped or stripped.startswith(";"):
            return None
        if not stripped.startswith("["):
            return None
        # "[ moleculetype ]  ; comment" is valid GROMACS — extract just
        # the bracketed token and lower it.
        end = stripped.find("]")
        if end == -1:
            return None
        return stripped[1:end].strip().lower()

    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        d = directive_of(line)
        if d is None:
            continue
        if start_idx is None and d == "moleculetype":
            start_idx = i
        elif start_idx is not None and d == "system":
            end_idx = i
            break

    if start_idx is None or end_idx is None or start_idx >= end_idx:
        raise RuntimeError(
            "could not locate [ moleculetype ]…[ system ] block in "
            f"{top_path}; refusing to split (ParmEd output format may "
            "have changed)"
        )

    itp_basename = os.path.basename(itp_path)
    top_basename = os.path.basename(top_path)
    block = lines[start_idx:end_idx]
    header = (
        f"; generated from {top_basename} by ep_amber_to_gromacs "
        f"(metalparm-vwf v1.11.0)\n"
    )
    with open(itp_path, "w") as fh:
        fh.write(header)
        fh.writelines(block)
        if not block or not block[-1].endswith("\n"):
            fh.write("\n")

    # Replace the moleculetype block in the .top with a single #include.
    new_lines = (
        lines[:start_idx]
        + [f'#include "{itp_basename}"\n', "\n"]
        + lines[end_idx:]
    )
    with open(top_path, "w") as fh:
        fh.writelines(new_lines)


def convert_amber_to_gromacs(
    prmtop_path: str,
    rst7_path: str,
    output_prefix: str,
    *,
    add_box_if_absent: bool = True,
    box_padding: float = 10.0,
    itp_filename: Optional[str] = None,
) -> dict:
    """Convert AMBER prmtop+rst7 → GROMACS top+gro using ParmEd.

    Args:
      prmtop_path: AMBER topology file (.prmtop / .parm7).
      rst7_path:   AMBER coordinates (.rst7 / .inpcrd / .ncrst).
      output_prefix: written as ``<prefix>.top`` and ``<prefix>.gro``.
      add_box_if_absent: if the prmtop has no periodic box, add a cubic
        box sized to the molecule's extent + 2 × box_padding (Å). This
        is what easyPARM's amber_converter does. Set False if you'll
        solvate later via ``gmx editconf`` / ``gmx solvate``.
      box_padding: half-width of cubic-box padding around the molecule
        (Å). Only applied when add_box_if_absent=True and the prmtop
        is non-periodic. Default 10 Å matches AMBER tutorial conventions.
      itp_filename: optional basename. When set, after ParmEd writes the
        monolithic .top, the [ moleculetype ] block is moved into a
        sibling ``<itp_filename>.itp`` and the .top is rewritten with
        ``#include "<itp_filename>.itp"`` — so the master .top can be
        consumed alongside other molecules' topologies. ``[ atomtypes ]``,
        ``[ defaults ]``, ``[ system ]``, and ``[ molecules ]`` stay in
        the .top. Both ``"complex"`` and ``"complex.itp"`` resolve to
        ``complex.itp``. Path separators are rejected (basename only).
        None or empty → no split (v1.10.0 behavior).

    Returns:
      Stats dict with paths and atom/box info for logging:
        {"top": str, "gro": str, "n_atoms": int, "has_box": bool,
         "box_added": bool, "itp": str | None}

    Raises:
      FileNotFoundError if prmtop or rst7 missing.
      ValueError if ``itp_filename`` is invalid (has a path separator).
      RuntimeError if the ITP split can't locate the moleculetype block.
      ImportError if ParmEd isn't available in the env (run inside
        the metalparm_vwf pixi env).
    """
    if not os.path.isfile(prmtop_path):
        raise FileNotFoundError(f"prmtop not found: {prmtop_path}")
    if not os.path.isfile(rst7_path):
        raise FileNotFoundError(f"rst7 not found: {rst7_path}")

    # Validate the ITP basename up front so we fail before invoking
    # ParmEd if the user's input is malformed.
    itp_basename: Optional[str] = None
    if itp_filename:
        itp_basename = normalize_itp_basename(itp_filename)

    try:
        import parmed as pmd
    except ImportError as e:  # pragma: no cover — env-dependent
        raise ImportError(
            "ParmEd is required for AMBER → GROMACS conversion. "
            "Run inside the metalparm_vwf pixi env."
        ) from e

    s = pmd.load_file(prmtop_path, rst7_path)
    n_atoms = len(s.atoms)

    has_box = s.box is not None
    box_added = False
    if not has_box and add_box_if_absent:
        # Cubic box sized to molecular extent + padding on each side.
        # This matches easyPARM's amber_converter pattern.
        coords = s.coordinates
        extent = float((coords.max(axis=0) - coords.min(axis=0)).max())
        edge = extent + 2 * box_padding
        s.box = [edge, edge, edge, 90.0, 90.0, 90.0]
        box_added = True

    top_path = f"{output_prefix}.top"
    gro_path = f"{output_prefix}.gro"
    # Order matters: write .gro first so a save error doesn't leave a
    # stale .top from a prior run.
    s.save(gro_path, format="gro", overwrite=True)
    s.save(top_path, format="gromacs", overwrite=True)

    itp_path: Optional[str] = None
    if itp_basename:
        itp_path = os.path.join(os.path.dirname(top_path) or ".", itp_basename)
        _split_top_into_itp(top_path, itp_path)

    return {
        "top": top_path,
        "gro": gro_path,
        "n_atoms": n_atoms,
        "has_box": has_box,
        "box_added": box_added,
        "itp": itp_path,
    }
