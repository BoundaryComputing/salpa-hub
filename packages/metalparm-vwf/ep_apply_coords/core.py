"""Coordinate-applier: write final complex.rst7 + complex.pdb from a
topology-only complex.prmtop plus aligned source PDBs.

The sister node to ``ep_fragment_fuse_topology``. Given a prmtop whose
atom order is ``peptide_residues + fragment_residues`` (the canonical
output of ``saveamberparm`` after ``combine { pep frag }`` + tleap
``remove`` operations), this node:

  1. loads the prmtop via ParmEd → ordered list of (residue_idx,
     atom_name) tuples;
  2. parses the source peptide PDB (residues 1..N) and the source
     fragment PDB (residue 1 in the fragment's own numbering, mapped
     to N+1 in the combined frame);
  3. for each prmtop atom, looks up the matching atom in the
     appropriate source PDB by ``(residue_idx_in_source, atom_name)``;
  4. assigns the source's XYZ to the prmtop atom;
  5. writes the final ``complex.rst7`` and ``complex.pdb``.

Atoms that tleap removed (caps OE2/HE2 on peptide GLU, CM/HM1-3/CAP/OAP
on fragment) are simply not iterated — they aren't in the prmtop.

Why split off this functionality:
  - Decouples re-runs: tweak interface_bonds → re-fuse topology only;
    re-run alignment → re-run coords only.
  - Reusable topology: one prmtop → many starting coordinate sets
    (docking poses, replicas, decoys for free-energy work).
  - DAG honesty: the ``fragment_align → ep_apply_coords`` edge is now
    a coordinate dependency, not pseudo-topology.

Pure-Python helpers live here so they can be unit-tested without ParmEd
where possible. The ParmEd-dependent rst7/pdb writer is kept thin.
"""
from __future__ import annotations

import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple


def parse_pdb_atoms(pdb_path: str) -> List[Dict]:
    """Parse ATOM/HETATM records from a PDB. Returns ordered list of dicts:
    ``{resseq, resname, atom_name, x, y, z, element}``.

    No external deps — works without ParmEd or BioPython. Coordinates are
    floats in Å; resseq is the integer residue sequence number from cols
    23–26 (1-based, as written by tleap savepdb / fragment_align's tleap
    transform output).
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB not found: {pdb_path}")
    out: List[Dict] = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            elem = line[76:78].strip() if len(line) >= 78 else ""
            if not elem:
                # Fallback: prefix of atom name (strip digits, keep first
                # 1-2 letters).
                stripped = "".join(c for c in atom_name if c.isalpha())
                elem = stripped[:2] if len(stripped) >= 2 and stripped[1].islower() else stripped[:1]
            out.append({
                "resseq": resseq, "resname": resname, "atom_name": atom_name,
                "x": x, "y": y, "z": z, "element": elem,
            })
    return out


def build_pdb_atom_index(pdb_atoms: List[Dict]) -> Dict[Tuple[int, str], Dict]:
    """Index PDB atoms by ``(resseq, atom_name)`` for O(1) lookup.

    Raises ValueError on duplicate keys (a malformed PDB with two atoms
    named 'CA' on the same residue would silently lose one in a dict).
    """
    idx: Dict[Tuple[int, str], Dict] = {}
    dupes = []
    for a in pdb_atoms:
        key = (a["resseq"], a["atom_name"])
        if key in idx:
            dupes.append(key)
        else:
            idx[key] = a
    if dupes:
        raise ValueError(
            f"PDB has duplicate (resseq, atom_name) keys: {dupes[:5]}"
            + (" ..." if len(dupes) > 5 else "")
        )
    return idx


def map_prmtop_to_source_pdbs(
    prmtop_atoms: List[Tuple[int, str]],
    peptide_pdb_atoms: List[Dict],
    fragment_pdb_atoms: List[Dict],
    *,
    pep_residues: int,
) -> List[Tuple[float, float, float]]:
    """For each prmtop atom (1-based resid, atom_name), look up coords in
    the right source PDB.

    Args:
      prmtop_atoms: ordered list of ``(residue_idx, atom_name)`` tuples
        in tleap's combined-unit numbering — peptide residues are
        1..pep_residues; fragment residues are pep_residues+1.. .
      peptide_pdb_atoms: pre-parsed peptide PDB atoms (from
        ``parse_pdb_atoms``).
      fragment_pdb_atoms: pre-parsed aligned-fragment PDB atoms.
      pep_residues: how many residues belong to the peptide. Used to
        switch source PDBs at the boundary.

    Returns:
      List of (x, y, z) tuples, one per prmtop atom, in prmtop order.

    Raises:
      ValueError if any prmtop atom has no matching atom in the
      appropriate source PDB. The error names the missing atom so the
      caller can debug atom-name-mismatch issues
      (e.g., antechamber rename of NH2 → N1).
    """
    pep_idx = build_pdb_atom_index(peptide_pdb_atoms)
    frag_idx = build_pdb_atom_index(fragment_pdb_atoms)

    coords: List[Tuple[float, float, float]] = []
    missing: List[Tuple[int, str]] = []
    for resid, atom_name in prmtop_atoms:
        if resid <= pep_residues:
            src = pep_idx
            lookup_key = (resid, atom_name)
            src_label = "peptide.pdb"
        else:
            src = frag_idx
            # The fragment PDB has its own residue numbering (typically
            # starting at 1); convert from combined-unit resid back.
            frag_resid = resid - pep_residues
            lookup_key = (frag_resid, atom_name)
            src_label = "fragment.pdb"
        match = src.get(lookup_key)
        if match is None:
            missing.append((resid, atom_name))
            continue
        coords.append((match["x"], match["y"], match["z"]))

    if missing:
        # Group by source PDB for a more actionable error
        by_src: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
        for resid, name in missing:
            src_label = "peptide.pdb" if resid <= pep_residues else "fragment.pdb"
            by_src[src_label].append((resid, name))
        msg_parts = [
            f"{len(items)} atom(s) missing in {src}: {items[:8]}"
            + (" ..." if len(items) > 8 else "")
            for src, items in by_src.items()
        ]
        raise ValueError(
            "apply_coords: prmtop atoms not found in source PDBs. "
            + " ; ".join(msg_parts)
        )
    return coords


def apply_coords(
    prmtop_path: str,
    peptide_pdb_path: str,
    fragment_pdb_path: str,
    output_prefix: str,
    *,
    pep_residues: int,
) -> dict:
    """Load prmtop + source PDBs, transfer coords, write rst7 + pdb.

    Args:
      prmtop_path: AMBER topology from ep_fragment_fuse_topology.
      peptide_pdb_path: aligned peptide PDB (typically the one written by
        peptide_builder; fragment_align doesn't move it).
      fragment_pdb_path: aligned fragment PDB (from fragment_align).
      output_prefix: produces ``<prefix>.rst7`` and ``<prefix>.pdb``.
      pep_residues: residue count of the peptide unit (forwarded by
        upstream nodes).

    Returns:
      Stats dict for logging::

        {"rst7": str, "pdb": str, "n_atoms": int,
         "pep_atoms": int, "frag_atoms": int}

    Raises:
      FileNotFoundError if prmtop / peptide.pdb / fragment.pdb missing.
      ValueError if any prmtop atom has no matching source-PDB atom.
      ImportError if ParmEd isn't available in the env.
    """
    if not os.path.isfile(prmtop_path):
        raise FileNotFoundError(f"prmtop not found: {prmtop_path}")
    if not os.path.isfile(peptide_pdb_path):
        raise FileNotFoundError(f"peptide PDB not found: {peptide_pdb_path}")
    if not os.path.isfile(fragment_pdb_path):
        raise FileNotFoundError(f"fragment PDB not found: {fragment_pdb_path}")
    if pep_residues <= 0:
        raise ValueError(f"pep_residues must be positive (got {pep_residues!r})")

    try:
        import parmed as pmd
    except ImportError as e:  # pragma: no cover — env-dependent
        raise ImportError(
            "ParmEd is required for ep_apply_coords. "
            "Run inside the metalparm_vwf pixi env."
        ) from e

    s = pmd.load_file(prmtop_path)
    # Build (residue_idx_1based, atom_name) tuples in prmtop order.
    prmtop_atoms: List[Tuple[int, str]] = [
        (atom.residue.idx + 1, atom.name) for atom in s.atoms
    ]

    peptide_atoms = parse_pdb_atoms(peptide_pdb_path)
    fragment_atoms = parse_pdb_atoms(fragment_pdb_path)

    coords = map_prmtop_to_source_pdbs(
        prmtop_atoms, peptide_atoms, fragment_atoms,
        pep_residues=pep_residues,
    )

    n_atoms = len(s.atoms)
    pep_atoms_n = sum(1 for r, _ in prmtop_atoms if r <= pep_residues)
    frag_atoms_n = n_atoms - pep_atoms_n

    # Apply coords. ParmEd expects a flat (n_atoms, 3) array.
    import numpy as np
    s.coordinates = np.array(coords, dtype=float)

    rst7_path = f"{output_prefix}.rst7"
    pdb_path = f"{output_prefix}.pdb"
    s.save(rst7_path, format="rst7", overwrite=True)
    s.save(pdb_path, format="pdb", overwrite=True)

    return {
        "rst7": rst7_path,
        "pdb": pdb_path,
        "n_atoms": n_atoms,
        "pep_atoms": pep_atoms_n,
        "frag_atoms": frag_atoms_n,
    }
