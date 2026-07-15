"""Core logic for md_analysis_distance — per-frame distance from a set
of probe atoms (Tyr/His quenchers) to the metal centre (Sn) along an
MD trajectory.

Pure helpers (CSV / summary — unit-testable, no heavy deps) are split
from ``run_distance_analysis`` (the MDAnalysis-backed pass, which
imports MDAnalysis lazily). Distances are minimum-image (PBC-aware)
and reported in Ångström.
"""
from __future__ import annotations

import csv
from statistics import mean, pstdev


def probe_label(resname: str, resid: int, atom: str) -> str:
    """Stable per-probe-atom column label, e.g. ``TYR2_OH``."""
    return f"{resname}{resid}_{atom}"


def distance_csv_header(labels: list[str]) -> list[str]:
    return ["frame", "time_ps"] + [f"{lbl}_A" for lbl in labels] + ["min_A"]


def distance_csv_rows(frames, times, dist_matrix) -> list[list]:
    """Rows from parallel series. ``dist_matrix`` is a list of per-frame
    lists of probe distances (Å); the per-frame minimum is appended."""
    rows = []
    for fr, t, dists in zip(frames, times, dist_matrix):
        d = [round(float(x), 3) for x in dists]
        rows.append([int(fr), round(float(t), 3), *d,
                     round(min(d), 3) if d else None])
    return rows


def write_csv(path: str, header: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def summarize(frames, dist_matrix, labels) -> dict:
    """Summary of a distance run: the closest approach over the whole
    trajectory (value / probe / frame) and the per-frame-minimum series
    statistics."""
    if not dist_matrix or not labels:
        return {"n_frames": 0}
    per_frame_min = [min(d) for d in dist_matrix]
    # global closest approach
    best_val, best_probe, best_frame = None, None, None
    for fr, dists in zip(frames, dist_matrix):
        for lbl, dv in zip(labels, dists):
            if best_val is None or dv < best_val:
                best_val, best_probe, best_frame = dv, lbl, fr
    return {
        "n_frames": len(dist_matrix),
        "n_probes": len(labels),
        "closest_approach_A": round(float(best_val), 3),
        "closest_probe": best_probe,
        "closest_frame": int(best_frame),
        "mean_min_dist_A": round(mean(per_frame_min), 3),
        "std_min_dist_A": round(
            pstdev(per_frame_min) if len(per_frame_min) > 1 else 0.0, 3),
        "max_min_dist_A": round(max(per_frame_min), 3),
    }


def run_distance_analysis(topology: str, trajectory: str, metal_sel: str,
                          probe_sel: str, stride: int) -> dict:
    """Per-frame minimum-image distances from each probe atom to the
    metal centre.

    Args:
      topology:   .tpr / .gro / .pdb.
      trajectory: .xtc / .trr.
      metal_sel:  MDAnalysis selection for the metal centre — must
                  resolve to exactly one atom.
      probe_sel:  MDAnalysis selection for the probe atoms (e.g. the
                  Tyr hydroxyl O and His ring N).
      stride:     analyse every Nth frame (>=1).

    Returns ``frames``, ``times`` (ps), ``labels`` (per-probe column
    labels), and ``dist_matrix`` (list of per-frame distance lists, Å).
    """
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    u = mda.Universe(topology, trajectory)
    metal = u.select_atoms(metal_sel)
    if metal.n_atoms != 1:
        raise ValueError(
            f"metal selection {metal_sel!r} matched {metal.n_atoms} atoms "
            f"— it must resolve to exactly one (the metal centre)")
    probes = u.select_atoms(probe_sel)
    if probes.n_atoms == 0:
        raise ValueError(
            f"probe selection {probe_sel!r} matched no atoms")

    labels = [probe_label(str(a.resname), int(a.resid), str(a.name))
              for a in probes]

    frames, times, dist_matrix = [], [], []
    step = max(1, int(stride))
    for ts in u.trajectory[::step]:
        d = distance_array(probes.positions, metal.positions,
                           box=ts.dimensions)
        dist_matrix.append([float(x) for x in d[:, 0]])
        frames.append(int(ts.frame))
        times.append(float(ts.time))

    return {
        "frames": frames,
        "times": times,
        "labels": labels,
        "dist_matrix": dist_matrix,
    }
