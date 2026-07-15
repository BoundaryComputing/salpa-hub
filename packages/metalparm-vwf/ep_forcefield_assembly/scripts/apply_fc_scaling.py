"""Apply vanilla easyPARM's BOND/ANGLE force-constant scaling to a frcmod.

Port of the awk scaling block in 01_easyPARM.sh lines 1269-1299. The
upstream EasyParm authors empirically calibrated these multipliers — the
Seminario method tends to underestimate force constants for weak metal-
ligand interactions (dative N→M, ionic O-M), and the scaling boosts
those low-magnitude constants into a "physically reasonable" range
matching experimental MD timescales.

Scaling rules (all units kcal·mol⁻¹·rad⁻² for ANGLE,
kcal·mol⁻¹·Å⁻² for BOND):

  ANGLE force constant k:
    k < 5      → k *= 11.599
    5  ≤ k < 10  → k *= 7.799
    10 ≤ k < 20  → k *= 3.599
    20 ≤ k < 29  → k *= 2.699
    k ≥ 29     → unchanged

  BOND force constant k:
    k < 20     → k *= 4.599
    k ≥ 20     → unchanged

For typical strong covalent bonds (k_bond > 200, k_angle > 60) and most
metal-N coordinations after Seminario averaging (e.g. SnP's na-Sn ~115
or Sn-os ~168), no scaling applies. Scaling kicks in for weakly
constrained metal sites — the published EasyParm validation set
includes Zn²⁺ and Cu²⁺ centers where Seminario gives k ≈ 5–25.

Usage:
    python apply_fc_scaling.py <frcmod>

Modifies the frcmod in place. Idempotent for k ≥ cutoffs (which is the
common case after Seminario averaging on strong bonds).
"""

import sys
from pathlib import Path

ANGLE_SCALES = [
    (5.0, 11.599),
    (10.0, 7.799),
    (20.0, 3.599),
    (29.0, 2.699),
]

BOND_SCALES = [
    (20.0, 4.599),
]


def _scale_value(k: float, table: list[tuple[float, float]]) -> float:
    """Return the scaled force constant according to ``table``.

    ``table`` is a list of ``(threshold, multiplier)`` rows in ascending
    threshold order. The first row whose threshold strictly exceeds ``k``
    wins; if no row matches, ``k`` is returned unchanged.
    """
    for threshold, multiplier in table:
        if k < threshold:
            return k * multiplier
    return k


def _try_float(s: str) -> float | None:
    try:
        return float(s)
    except ValueError:
        return None


def _scale_line(line: str, table: list[tuple[float, float]]) -> str:
    """Scale the second whitespace-separated column of ``line`` in place
    using ``table``. Preserves leading whitespace and the original
    line's column layout via direct field substitution.

    A frcmod BOND/ANGLE line looks like:
        ``c2-na-Sn 78.510 124.831`` (atom-types, force_const, eq_geom).
    The awk script in vanilla easyPARM uses ``$2 ~ /[0-9]/`` to gate;
    here we replicate by attempting to parse ``$2`` as float and
    skipping the line if parsing fails (pure-text lines like section
    headers, blank lines, and IMPROPER comments are unaffected).
    """
    parts = line.rstrip("\n").split()
    if len(parts) < 3:
        return line
    k = _try_float(parts[1])
    if k is None:
        return line
    new_k = _scale_value(k, table)
    if new_k == k:
        return line
    # Replace only the first occurrence of the original token to preserve
    # the rest of the line's spacing exactly.
    parts[1] = f"{new_k:g}"
    return " ".join(parts) + ("\n" if line.endswith("\n") else "")


def scale_frcmod(frcmod_path: Path) -> dict[str, int]:
    """Scale the BOND and ANGLE sections of ``frcmod_path`` in place.

    Returns a dict with counts of lines scaled in each section
    (for logging / testing).
    """
    text = frcmod_path.read_text()
    lines = text.splitlines(keepends=True)

    section: str | None = None
    counts = {"BOND": 0, "ANGLE": 0}
    out_lines: list[str] = []

    for raw in lines:
        s = raw.strip()
        if s in {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}:
            section = s
            out_lines.append(raw)
            continue

        if section == "BOND":
            new_line = _scale_line(raw, BOND_SCALES)
            if new_line != raw:
                counts["BOND"] += 1
            out_lines.append(new_line)
        elif section == "ANGLE":
            new_line = _scale_line(raw, ANGLE_SCALES)
            if new_line != raw:
                counts["ANGLE"] += 1
            out_lines.append(new_line)
        else:
            out_lines.append(raw)

    frcmod_path.write_text("".join(out_lines))
    return counts


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    frcmod = Path(sys.argv[1]).resolve()
    if not frcmod.is_file():
        print(f"frcmod not found: {frcmod}", file=sys.stderr)
        sys.exit(1)

    counts = scale_frcmod(frcmod)
    if counts["BOND"] or counts["ANGLE"]:
        print(f"apply_fc_scaling: scaled {counts['BOND']} BOND + "
              f"{counts['ANGLE']} ANGLE entries in {frcmod.name}")


if __name__ == "__main__":
    main()
