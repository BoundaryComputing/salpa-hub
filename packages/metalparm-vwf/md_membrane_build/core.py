"""Core logic for md_membrane_build — embed a dry metallopeptide
GROMACS topology in a DPPC bilayer and produce a solvated GROMACS
topology for membrane MD, **without ever re-deriving the solute**.

Why the indirection (same lesson as md_solvate_gmx): packmol-memgen's
``--parametrize`` step runs tleap, which re-types the solute against
the standard AMBER library — for a fragment-fused residue (the SnP
GLU that has lost OE2 and gained a bond to the porphyrin N) tleap
auto-completes the "missing" OE2 on top of the fragment, exploding MD
step 0. So this node:

  1. runs packmol-memgen *without* ``--parametrize`` — geometry only:
     orient the peptide, build the DPPC bilayer, add water + ions,
     write a packed PDB;
  2. splits the packed PDB into the solute part and the
     membrane part (DPPC + water + ions — all *standard* residues);
  3. tleap-parametrises only the membrane part (leaprc.lipid21 +
     leaprc.water.opc) — standard residues, so no re-derivation hazard;
  4. ParmEd-concatenates the preserved solute topology (from the
     GROMACS .top) with the membrane topology, coordinates from the
     same packmol-memgen frame;
  5. exports GROMACS .top + .gro (+ optional .itp split).

The solute topology is the GROMACS .top from ep_amber_to_gromacs and
is never round-tripped through tleap.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Residue names that belong to the *membrane* (everything packmol-memgen
# adds): Lipid21 DPPC is a split-residue model (PC headgroup + two PA
# tails), water is WAT, ions are monatomic. Anything else in the packed
# PDB is the solute.
DPPC_RESNAMES = {"PC", "PA", "DPPC"}
WATER_RESNAMES = {"WAT", "HOH", "SOL"}
ION_RESNAMES = {"NA", "CL", "K", "NA+", "CL-", "K+", "Na+", "Cl-", "K+"}
MEMBRANE_RESNAMES = DPPC_RESNAMES | WATER_RESNAMES | ION_RESNAMES


def load_solute_structure(top_path: str, gro_path: str) -> object:
    """Load the dry metallopeptide GROMACS top+gro into a ParmEd
    Structure — preserves every moleculetype verbatim (no re-typing)."""
    import parmed as pmd
    return pmd.load_file(top_path, xyz=gro_path)


def write_solute_pdb(structure: object, path: str) -> str:
    """Write the solute Structure to a PDB for packmol-memgen input."""
    structure.save(path, format="pdb", overwrite=True)
    return path


def orient_peptide_along_z(structure: object) -> object:
    """Rotate the solute in place so the peptide helix axis aligns with
    z — the membrane normal — for transmembrane embedding.

    The helix axis is the first principal component of the peptide Cα
    atoms. This is done geometrically because MEMEMBED (packmol-memgen's
    default orienter) cannot handle a metallopeptide — the non-standard
    SnP fragment residue makes it fail to emit ``*_EMBED.pdb``. So the
    node always pre-orients here and runs packmol-memgen ``--preoriented``.

    No-op (returns unchanged) if fewer than 3 Cα atoms are present.
    """
    import numpy as np

    ca = [i for i, a in enumerate(structure.atoms) if a.name == "CA"]
    if len(ca) < 3:
        return structure
    coords = np.asarray(structure.coordinates, dtype=float)
    centre = coords[ca].mean(axis=0)
    centred_ca = coords[ca] - centre
    # principal axis = eigenvector of the largest eigenvalue
    _, evecs = np.linalg.eigh(centred_ca.T @ centred_ca)
    axis = evecs[:, -1]
    axis = axis / np.linalg.norm(axis)

    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(axis, z)
    s = float(np.linalg.norm(v))
    c = float(np.dot(axis, z))
    if s < 1e-8:                       # already (anti)parallel to z
        rot = np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    else:                              # Rodrigues rotation axis→z
        vx = np.array([[0.0, -v[2], v[1]],
                       [v[2], 0.0, -v[0]],
                       [-v[1], v[0], 0.0]])
        rot = np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))

    structure.coordinates = (coords - centre) @ rot.T
    return structure


def run_packmol_memgen(
    solute_pdb: str,
    work_dir: str,
    *,
    lipid: str = "DPPC",
    dist: float = 17.5,
    dist_wat: float = 17.5,
    saltcon: float = 0.15,
    xy_box_A: Optional[float] = None,
    nloop_all: Optional[int] = None,
    overwrite: bool = True,
) -> str:
    """Run packmol-memgen for **geometry only** (no ``--parametrize``).

    Builds a ``lipid`` bilayer around ``solute_pdb``, adds a water layer
    of ``dist_wat`` Å each side and NaCl at ``saltcon`` M (0 ⇒ counter-
    ions only). Returns the path to the packed ``bilayer_*.pdb``.

    Always passes ``--preoriented`` — MEMEMBED cannot orient a
    metallopeptide, so the solute is pre-oriented geometrically (see
    ``orient_peptide_along_z``) before this call. Also ``--notprotonate
    --nottrim``: the solute already carries a complete hydrogen set from
    the metallopeptide topology; re-running ``reduce`` would corrupt
    the fragment.

    ``xy_box_A`` (optional) passes ``--distxy_fix`` to force the XY box
    to that exact width in Å. Use when packmol's auto-sized XY (driven by
    ``dist`` + solute extent) is too small for the lipid count and the
    all-together packing loop won't converge.

    ``nloop_all`` (optional) passes ``--nloop_all`` to bump the
    all-together packing iteration cap above packmol-memgen's default.
    """
    os.makedirs(work_dir, exist_ok=True)
    local_pdb = os.path.join(work_dir, os.path.basename(solute_pdb))
    if os.path.abspath(local_pdb) != os.path.abspath(solute_pdb):
        shutil.copy2(solute_pdb, local_pdb)

    cmd = [
        "packmol-memgen",
        "--pdb", os.path.basename(local_pdb),
        "--lipids", lipid,
        "--dist", f"{dist:.2f}",
        "--dist_wat", f"{dist_wat:.2f}",
        "--preoriented",
        "--notprotonate", "--nottrim",
        "--noprogress",
    ]
    if xy_box_A is not None and xy_box_A > 0:
        cmd += ["--distxy_fix", f"{xy_box_A:.2f}"]
    if nloop_all is not None and nloop_all > 0:
        cmd += ["--nloop_all", str(int(nloop_all))]
    if saltcon and saltcon > 0:
        cmd += ["--salt", "--saltcon", f"{saltcon:.3f}"]
    if overwrite:
        cmd.append("--overwrite")

    proc = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
    packed = _find_packed_pdb(work_dir, os.path.basename(local_pdb))
    if packed is None:
        tail = (proc.stderr or proc.stdout or "")[-2500:]
        raise RuntimeError(
            "packmol-memgen produced no bilayer_*.pdb. Tail of output:\n"
            + tail)
    return packed


def _find_packed_pdb(work_dir: str, solute_basename: str) -> Optional[str]:
    """Locate packmol-memgen's packed output (``bilayer_<solute>.pdb``)."""
    stem = Path(solute_basename).stem
    candidates = [
        os.path.join(work_dir, f"bilayer_{stem}.pdb"),
        os.path.join(work_dir, f"bilayer_{solute_basename}"),
    ]
    for c in candidates:
        if os.path.isfile(c) and os.path.getsize(c) > 0:
            return c
    # Fallback: any bilayer_*.pdb in the dir.
    hits = sorted(Path(work_dir).glob("bilayer_*.pdb"),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    return str(hits[0]) if hits else None


def _resname(line: str) -> str:
    return line[17:20].strip()


def split_packed_pdb(packed_pdb: str, work_dir: str,
                      ) -> Tuple[str, str, Dict[str, int]]:
    """Split the packed PDB into a solute PDB and a membrane PDB.

    Membrane = residues in ``MEMBRANE_RESNAMES`` (DPPC split-residues,
    water, ions). Solute = everything else, in packed order. Returns
    ``(solute_pdb, membrane_pdb, counts)`` where counts reports the
    number of solute / lipid / water / ion atoms.
    """
    solute_lines: List[str] = []
    membrane_lines: List[str] = []
    counts = {"solute": 0, "lipid": 0, "water": 0, "ion": 0}
    cryst = None
    last = None   # the list the previous ATOM went to — TER follows it
    for raw in Path(packed_pdb).read_text().splitlines():
        if raw.startswith("CRYST1"):
            cryst = raw
            continue
        if raw.startswith("TER"):
            # Preserve molecule separators — tleap connects consecutive
            # non-TER-separated residues, so dropping TER would bond
            # adjacent lipids into one chain (split-residue chirality
            # error). Route each TER to the side its preceding atom went.
            if last is not None:
                last.append("TER")
            continue
        if not raw.startswith(("ATOM  ", "HETATM")):
            continue
        rn = _resname(raw)
        if rn in MEMBRANE_RESNAMES:
            membrane_lines.append(raw)
            last = membrane_lines
            if rn in DPPC_RESNAMES:
                counts["lipid"] += 1
            elif rn in WATER_RESNAMES:
                counts["water"] += 1
            else:
                counts["ion"] += 1
        else:
            solute_lines.append(raw)
            last = solute_lines
            counts["solute"] += 1
    if not solute_lines:
        raise RuntimeError(
            "no solute atoms found in the packed PDB — every residue "
            "matched a membrane/solvent name; check the input solute")
    if not membrane_lines:
        raise RuntimeError(
            "no membrane atoms found in the packed PDB — packmol-memgen "
            "did not add lipids/water")

    solute_pdb = os.path.join(work_dir, "packed_solute.pdb")
    membrane_pdb = os.path.join(work_dir, "packed_membrane.pdb")
    head = [cryst] if cryst else []
    Path(solute_pdb).write_text("\n".join(head + solute_lines + ["END"]) + "\n")
    Path(membrane_pdb).write_text(
        "\n".join(head + membrane_lines + ["END"]) + "\n")
    return solute_pdb, membrane_pdb, counts


def read_packmol_tolerance(packmol_inp: str, default: float = 2.0) -> float:
    """Read the ``tolerance`` (minimum inter-atom distance, Å) packmol
    was told to pack with. Used as the periodic-box margin so periodic
    images clear each other by the same distance — see
    ``assemble_membrane_system``. Falls back to ``default`` (2.0 Å,
    packmol-memgen's own default) if the keyword is absent.
    """
    try:
        for raw in Path(packmol_inp).read_text().splitlines():
            s = raw.split()
            if len(s) >= 2 and s[0] == "tolerance":
                return float(s[1])
    except (OSError, ValueError):
        pass
    return default


def count_membrane_residues(system) -> Dict[str, int]:
    """Count membrane *molecules* in an assembled ParmEd system by
    residue name.

    Lipid21 DPPC is a split-residue model (one PC headgroup + two PA
    tails per lipid), so the PC residue count is the lipid molecule
    count; water is one residue per molecule. ``split_packed_pdb``'s
    ``counts`` are *atom* counts — use this for molecule counts.

    Returns ``lipid`` / ``water`` / ``ion`` molecule counts plus the
    raw per-resname tally under ``by_resname``.
    """
    from collections import Counter
    tally = Counter(r.name.strip() for r in system.residues)
    n_lipid = tally.get("PC", 0) or tally.get("DPPC", 0)
    n_water = sum(tally.get(w, 0) for w in (WATER_RESNAMES | {"OPC"}))
    n_ion = sum(c for n, c in tally.items() if n in ION_RESNAMES)
    return {"lipid": n_lipid, "water": n_water, "ion": n_ion,
            "by_resname": dict(tally)}


def parametrize_membrane(membrane_pdb: str, work_dir: str) -> Tuple[str, str]:
    """tleap-parametrise the membrane part (DPPC + water + ions).

    Every residue here is *standard* — Lipid21 DPPC, OPC water,
    monatomic ions — so tleap parametrises it cleanly; the
    re-derivation hazard only ever applied to the non-standard solute,
    which this node never sends through tleap.

    Returns ``(prmtop, rst7)``.
    """
    script = "\n".join([
        "source leaprc.lipid21",
        "source leaprc.water.opc",
        "loadamberparams frcmod.ionslm_126_opc",
        f"mem = loadpdb {os.path.basename(membrane_pdb)}",
        "saveamberparm mem membrane.prmtop membrane.rst7",
        "quit",
        "",
    ])
    tleap_in = os.path.join(work_dir, "membrane.tleap")
    Path(tleap_in).write_text(script)
    proc = subprocess.run(
        ["tleap", "-f", "membrane.tleap"],
        cwd=work_dir, capture_output=True, text=True)
    prmtop = os.path.join(work_dir, "membrane.prmtop")
    rst7 = os.path.join(work_dir, "membrane.rst7")
    if not (os.path.isfile(prmtop) and os.path.getsize(prmtop) > 0):
        tail = (proc.stdout or "")[-2500:]
        raise RuntimeError(
            "tleap did not parametrise the membrane. Tail of leap output:\n"
            + tail)
    return prmtop, rst7


def assemble_membrane_system(
    solute_structure: object,
    packed_solute_pdb: str,
    membrane_prmtop: str,
    membrane_rst7: str,
    margin: float = 2.0,
) -> object:
    """Concatenate the preserved solute topology with the tleap-built
    membrane topology, coordinates taken from the common packmol-memgen
    frame, and set a clash-free periodic box.

    The solute's coordinates are overlaid from ``packed_solute_pdb``
    (packmol-memgen placed the rigid solute in the bilayer; with
    ``--notprotonate --nottrim`` its atom count + order are unchanged
    from the input, so the overlay is a direct index assignment). The
    membrane Structure already carries the packed coordinates via the
    tleap rst7.

    **The periodic box.** packmol packs *non-periodically*: it enforces
    the pairwise ``tolerance`` only between explicitly-packed pairs,
    never between an atom and a periodic image, and its ``inside box``
    is a *soft* (penalty) constraint atoms leak ~1–2 Å past. So the
    box cannot be taken from the ``inside box`` constraint — sized to
    that, atoms near opposite faces become periodic neighbours ~0 Å
    apart and the r⁻¹² wall makes the first minimisation step diverge.

    Instead the box is the *actual packed-coordinate bounding box* plus
    ``margin`` per dimension, with the coordinates shifted to sit
    centred in it. Because every atom then lies within an extent of
    ``box - margin``, no atom can come within ``margin`` of any periodic
    image — set ``margin`` to the packmol tolerance and the periodic
    gap is no tighter than the closest in-box contact packmol itself
    allowed. The box is a few % larger than packmol-memgen intended
    (the loose pack + margin); the thin vacuum seam closes within the
    first ps of NPT. packmol-memgen emits no CRYST1, so this
    coordinate-derived box is the cell.
    """
    import numpy as np
    import parmed as pmd

    packed_solute = pmd.load_file(packed_solute_pdb)
    if len(packed_solute.atoms) != len(solute_structure.atoms):
        raise RuntimeError(
            f"solute atom-count mismatch: GROMACS topology has "
            f"{len(solute_structure.atoms)}, packed solute PDB has "
            f"{len(packed_solute.atoms)} — packmol-memgen must not "
            f"reorder/reprotonate the solute (check --notprotonate)")
    solute_structure.coordinates = packed_solute.coordinates

    membrane = pmd.load_file(membrane_prmtop, xyz=membrane_rst7)
    system = solute_structure + membrane

    m = float(margin)
    coords = np.asarray(system.coordinates, dtype=float)
    cmin = coords.min(axis=0)
    cmax = coords.max(axis=0)
    # shift so the packed system sits in [m/2, extent + m/2]: a margin/2
    # vacuum strip on every face ⇒ periodic images clear by >= margin.
    system.coordinates = coords - cmin + m / 2.0
    ext = cmax - cmin
    system.box = [ext[0] + m, ext[1] + m, ext[2] + m, 90.0, 90.0, 90.0]
    return system


def save_gromacs_outputs(
    structure: object,
    output_dir: str,
    *,
    prefix: str = "complex",
    itp_filename: Optional[str] = "metallopeptide_mem",
) -> Dict[str, Optional[str]]:
    """Write the assembled Structure as GROMACS .top + .gro, optionally
    splitting the metallopeptide moleculetype into a sibling .itp.
    Reuses ep_amber_to_gromacs's split helper (same package)."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    top_path = out_dir / f"{prefix}.top"
    gro_path = out_dir / f"{prefix}.gro"

    structure.save(str(gro_path), format="gro", overwrite=True)
    structure.save(str(top_path), format="gromacs", overwrite=True)

    itp_path: Optional[Path] = None
    if itp_filename:
        try:
            from ..ep_amber_to_gromacs.core import (
                _split_top_into_itp, normalize_itp_basename,
            )
        except (ImportError, ValueError):
            # Robust fallback: load the sibling node's core.py by file path
            # under a UNIQUE module name (a bare `from core import ...` would
            # resolve to THIS node's own already-imported `core`).
            import importlib.util as _ilu
            _sib = Path(__file__).resolve().parent.parent / "ep_amber_to_gromacs" / "core.py"
            _spec = _ilu.spec_from_file_location("ep_amber_to_gromacs_core", str(_sib))
            _mod = _ilu.module_from_spec(_spec)
            _spec.loader.exec_module(_mod)
            _split_top_into_itp = _mod._split_top_into_itp
            normalize_itp_basename = _mod.normalize_itp_basename
        itp_path = out_dir / normalize_itp_basename(itp_filename)
        _split_top_into_itp(str(top_path), str(itp_path))

    return {
        "top": str(top_path),
        "gro": str(gro_path),
        "itp": str(itp_path) if itp_path else None,
    }


def compute_solute_charge(structure: object) -> float:
    """Net solute charge (sum of atom partial charges)."""
    return float(sum(a.charge for a in structure.atoms))
