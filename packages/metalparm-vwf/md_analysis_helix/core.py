"""Core logic for md_analysis_helix — per-frame α-helix content of a
peptide along an MD trajectory, via DSSP.

Split into pure helpers (CSV / summary statistics — unit-testable with
no heavy deps) and ``run_helix_analysis`` (the MDAnalysis-backed pass,
which imports MDAnalysis lazily so this module stays importable for
server-side node introspection).

DSSP secondary structure is taken from MDAnalysis' simplified 3-state
assignment ('H' helix, 'E' strand, '-' loop); the 'H' state merges the
α/3-10/π helix types, so "helix fraction" here is total helix content.
"""
from __future__ import annotations

import csv
from statistics import mean, pstdev


CSV_HEADER = ["frame", "time_ps", "n_residues", "n_helix", "frac_helix"]


def summarize(fracs: list[float]) -> dict:
    """Summary statistics of a per-frame helix-fraction series."""
    if not fracs:
        return {"n_frames": 0}
    return {
        "n_frames": len(fracs),
        "mean_frac_helix": round(mean(fracs), 4),
        "std_frac_helix": round(pstdev(fracs) if len(fracs) > 1 else 0.0, 4),
        "min_frac_helix": round(min(fracs), 4),
        "max_frac_helix": round(max(fracs), 4),
    }


def helix_csv_rows(frames, times, counts, fracs, n_residues) -> list[list]:
    """Build CSV rows (one per analysed frame) from the parallel series."""
    rows = []
    for fr, t, c, f in zip(frames, times, counts, fracs):
        rows.append([int(fr), round(float(t), 3), int(n_residues),
                     int(c), round(float(f), 4)])
    return rows


def write_csv(path: str, header: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def run_helix_analysis(topology: str, trajectory: str, selection: str,
                       stride: int) -> dict:
    """Run DSSP over the trajectory and return per-frame helix content.

    Args:
      topology:   .tpr / .gro / .pdb — anything MDAnalysis reads for
                  topology (a .tpr gives the richest atom metadata).
      trajectory: .xtc / .trr.
      selection:  MDAnalysis selection string for the peptide
                  (e.g. ``protein``).
      stride:     analyse every Nth frame (>=1).

    Returns a dict with parallel series ``frames``, ``times`` (ps),
    ``counts`` (helix residues), ``fracs`` (helix fraction), the
    residue count ``n_residues``, and ``per_residue_propensity`` (mean
    helix occupancy per residue, aligned with ``resids``).
    """
    import MDAnalysis as mda
    from MDAnalysis.analysis.dssp import DSSP

    u = mda.Universe(topology, trajectory)
    peptide = u.select_atoms(selection)
    if peptide.n_atoms == 0:
        raise ValueError(
            f"peptide selection {selection!r} matched no atoms")

    # DSSP requires every residue to carry a complete N/CA/C/O backbone.
    # Drop caps (ACE, NME, …) and any other backbone-incomplete residue
    # the selection swept in, so the user need not hand-tune the string.
    bb = {"N", "CA", "C", "O"}
    keep = [bb.issubset(set(r.atoms.names)) for r in peptide.residues]
    if not any(keep):
        raise ValueError(
            f"selection {selection!r} has no residue with a complete "
            f"N/CA/C/O backbone — DSSP cannot run")
    peptide = peptide.residues[keep].atoms
    n_residues = peptide.n_residues

    dssp = DSSP(peptide).run(step=max(1, int(stride)))
    ss = dssp.results.dssp                      # (n_frames, n_residues) chars
    helix = (ss == "H")
    counts = helix.sum(axis=1).tolist()
    fracs = (helix.sum(axis=1) / n_residues).tolist()
    per_residue = helix.mean(axis=0).tolist()   # mean occupancy per residue

    return {
        "frames": [int(x) for x in dssp.frames],
        "times": [float(x) for x in dssp.times],
        "counts": counts,
        "fracs": fracs,
        "n_residues": int(n_residues),
        "resids": [int(r) for r in peptide.residues.resids],
        "resnames": [str(n) for n in peptide.residues.resnames],
        "per_residue_propensity": [round(float(p), 4) for p in per_residue],
    }


def per_frame_helix_matrix(topology: str, trajectory: str, selection: str,
                           stride: int = 1) -> dict:
    """Per-frame, per-residue DSSP helix assignment (for time-resolved
    colouring — e.g. the trajectory viz animates residue colour by the
    instantaneous helix state rather than the whole-trajectory average).

    Mirrors :func:`run_helix_analysis`'s residue selection + backbone keep
    filter, but returns the full boolean matrix instead of collapsing it.

    Returns ``{resids, resnames, frames, matrix}`` where ``matrix`` is a
    list of per-frame rows, each a list of 0/1 (1 = DSSP 'H') aligned with
    ``resids``.
    """
    import MDAnalysis as mda
    from MDAnalysis.analysis.dssp import DSSP

    u = mda.Universe(topology, trajectory)
    peptide = u.select_atoms(selection)
    if peptide.n_atoms == 0:
        raise ValueError(f"peptide selection {selection!r} matched no atoms")
    bb = {"N", "CA", "C", "O"}
    keep = [bb.issubset(set(r.atoms.names)) for r in peptide.residues]
    peptide = peptide.residues[keep].atoms

    dssp = DSSP(peptide).run(step=max(1, int(stride)))
    helix = (dssp.results.dssp == "H")
    return {
        "resids": [int(r) for r in peptide.residues.resids],
        "resnames": [str(n) for n in peptide.residues.resnames],
        "frames": [int(x) for x in dssp.frames],
        "matrix": [[int(v) for v in row] for row in helix],
    }
