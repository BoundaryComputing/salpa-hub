"""Core logic for md_traj_center — PBC-correct an MD trajectory so the
solute (the metallopeptide) is whole and centred before analysis.

Why this matters: GROMACS writes coordinates wrapped into the periodic
box, so a molecule that drifts across a box face is split — half its
atoms on one side, half on the other. Geometry-sensitive analyses
break on a split molecule: DSSP mis-assigns secondary structure
because backbone H-bond distances are wrong, a radius of gyration
explodes, a centre-of-mass jumps. (A minimum-image *distance* — like
md_analysis_distance — is immune, but helix content is not.)

This module makes the solute whole again (unwrap across PBC using the
bond graph), centres it in the box, and wraps the solvent/membrane
back in — the standard ``gmx trjconv -pbc whole`` + ``-center`` +
``-pbc mol`` recipe, done with MDAnalysis transformations so no
GROMACS binary is needed and the result is identical to what the
analysis nodes (also MDAnalysis) would see.

Pure helpers (selection / summary — unit-testable, no heavy deps) are
split from ``run_traj_center`` (the MDAnalysis-backed pass, which
imports MDAnalysis lazily so this module stays importable for
server-side node introspection).
"""
from __future__ import annotations

import os


def summarize(n_in: int, n_out: int, n_solute: int, n_other: int,
              stride: int) -> dict:
    """Summary of a centring run."""
    return {
        "n_frames_in": int(n_in),
        "n_frames_out": int(n_out),
        "stride": int(stride),
        "n_solute_atoms": int(n_solute),
        "n_other_atoms": int(n_other),
    }


def run_traj_center(topology: str, trajectory: str, selection: str,
                    out_traj: str, out_gro: str | None = None,
                    stride: int = 1, extract_first: bool = False) -> dict:
    """Make the solute whole, centre it in the box, wrap everything else.

    Args:
      topology:   .tpr / .gro / .pdb. A ``.tpr`` is strongly preferred —
                  unwrapping needs the bond graph, and a .tpr carries
                  bonds; a bare .gro does not.
      trajectory: .xtc / .trr to process.
      selection:  MDAnalysis selection for the solute to keep whole and
                  centre (e.g. ``protein`` or ``not resname WAT``).
      out_traj:   path for the PBC-corrected .xtc.
      out_gro:    optional path for a .gro of the first corrected frame
                  (a convenient topology for the centred trajectory).
      stride:     write every Nth frame (>=1).
      extract_first: two-pass mode. When True (and the solute is smaller
                  than the full system) stream the input trajectory once
                  to write a SOLUTE-ONLY temp xtc + PDB (with CONECT
                  bonds), then run the unwrap + centre on that small
                  system. Output is solute-only — perfect for downstream
                  analyses (DSSP, distances) that only need the solute,
                  and orders of magnitude faster for membrane systems
                  where the bulk of every frame is lipid/water I/O.

    Returns the ``summarize`` dict plus ``out_traj`` / ``out_gro``.
    """
    import MDAnalysis as mda
    from MDAnalysis import transformations as trans

    u = mda.Universe(topology, trajectory)
    seed = u.select_atoms(selection)
    if seed.n_atoms == 0:
        raise ValueError(f"solute selection {selection!r} matched no atoms")
    if not hasattr(u, "bonds") or len(u.bonds) == 0:
        raise ValueError(
            "the topology carries no bonds — unwrapping cannot run. "
            "Use a .tpr (it stores the bond graph); a bare .gro does not.")

    # Expand the selection to whole connected molecules. The SnP
    # fragment (porphyrin + Sn) is bonded to the peptide GLU but is not
    # 'protein' — selecting only 'protein' would leave the fragment in
    # `others` and wrap() would move it independently, splitting the
    # metallopeptide. `same fragment as` pulls in the whole molecule.
    solute = u.select_atoms(f"same fragment as ({selection})")

    # --- two-pass extract-first optimisation ----------------------------
    # For large systems (membrane + water) the per-frame xtc I/O is
    # dominated by environment atoms that the analysis never touches.
    # Stream the trajectory once, write a tiny solute-only xtc + a PDB
    # (with CONECT for bonds), then recurse on that small system to do
    # the actual unwrap/centre. ~90× smaller xtc for a DPPC system; the
    # output is solute-only (helix/distance analyses use only the solute).
    if extract_first and solute.n_atoms < u.atoms.n_atoms:
        out_dir = os.path.dirname(os.path.abspath(out_traj)) or "."
        extract_pdb = os.path.join(out_dir, "_subset.pdb")
        extract_xtc = os.path.join(out_dir, "_subset.xtc")
        solute.write(extract_pdb, bonds="conect")           # topology w/ bonds
        step = max(1, int(stride))
        with mda.Writer(extract_xtc, solute.n_atoms) as w:
            for _ in u.trajectory[::step]:
                w.write(solute)
        # Recurse on the small system. selection='all' since the subset
        # universe contains only the solute; stride is already applied.
        out = run_traj_center(extract_pdb, extract_xtc, "all", out_traj,
                              out_gro=out_gro, stride=1, extract_first=False)
        out["extracted_first"] = True
        out["extract_pdb"] = extract_pdb
        out["extract_xtc"] = extract_xtc
        return out

    others = u.atoms - solute
    workflow = [trans.unwrap(u.atoms),
                trans.center_in_box(solute, center="mass", wrap=False)]
    if others.n_atoms:
        workflow.append(trans.wrap(others))
    u.trajectory.add_transformations(*workflow)

    step = max(1, int(stride))
    n_in = u.trajectory.n_frames
    n_out = 0
    first_written = False
    with mda.Writer(out_traj, u.atoms.n_atoms) as w:
        for ts in u.trajectory[::step]:
            w.write(u.atoms)
            n_out += 1
            if out_gro and not first_written:
                u.atoms.write(out_gro)
                first_written = True

    summary = summarize(n_in, n_out, solute.n_atoms, others.n_atoms, step)
    summary["out_traj"] = os.path.abspath(out_traj)
    summary["out_gro"] = (os.path.abspath(out_gro)
                          if out_gro and first_written else None)
    return summary
