"""Pure-Python geometry primitives for fragment_align.

What this does (Avogadro / GaussView style — deterministic, no FF, no QM):

    Given a peptide PDB with an "anchor" atom (e.g. GLU.CD) and a fragment
    PDB with a corresponding anchor (e.g. SnP.NH2 → mol.N1 after antechamber
    renaming), compute a rigid-body transformation that places the fragment
    so its anchor sits at the right amide-bond distance + angle from the
    peptide anchor.

Algorithm:

1. Compute the "outward" direction at the peptide anchor — the direction
   into which the new bond should extend, based on the anchor's
   hybridization and its existing bonds (after the user-specified
   ``pep_remove`` atoms are dropped). For sp2 carbonyl GLU.CD with CG
   and OE1 still bonded: outward = -unit(unit(CD->CG) + unit(CD->OE1))
   (the bisector pointing away from the existing bonds).

2. Compute the fragment anchor's current bond axis — the direction it's
   currently pointing AWAY from its existing bonds (after ``frag_remove``
   atoms are dropped). Same hybridization formula, applied to the fragment.

3. Rotation = the rotation that aligns the fragment axis ANTI-parallel to
   the peptide outward direction (so the bond points correctly). Rodrigues'
   formula on the cross-product axis, signed angle.

4. Translation = move the fragment so its anchor sits at
   pep_anchor + L * pep_outward, where L is a standard bond length lookup.

5. (Optional) secondary rotation around the new bond axis to minimize
   atom-pair clashes between peptide and fragment. Cheap O(n²) scan over a
   handful of candidate angles.

This is the same math GUI molecular editors use when you "join" two fragments
through a bond — it gives a "decent initial structure" without needing any
energy minimization. Subsequent MD prep (solvate, ions, equilibration) cleans
up any small residual clashes.

The module is pure numpy + biopython; no bocoflow_core, no tleap, so it can
be unit-tested standalone. The node.py wrapper drives tleap separately to
apply the computed transformation to the fragment OFF library.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np
from Bio.PDB import PDBParser


# ---------------------------------------------------------------------------
# Standard parameters
# ---------------------------------------------------------------------------

# Common interface bond lengths (Å). Keyed by an opaque label that callers
# pick via interface_bond["bond_kind"], or auto-inferred from the anchor
# atom elements if the field is missing (see _infer_bond_kind).
STANDARD_BOND_LENGTHS: dict[str, float] = {
    # Organic / peptide single bonds
    "C-N_amide":      1.33,
    "C-N_amine":      1.47,
    "C-N_imine":      1.28,
    "C-O_ester":      1.34,
    "C-S_thioether":  1.81,
    "S-S_disulfide":  2.05,
    "C-C_single":     1.54,
    "Calpha-N_pep":   1.47,
    # Metal-ligand defaults — only used when the anchor is the metal itself.
    # Values are typical AMBER/SDB-MM equilibrium distances at standard
    # protonation; for high-spin / low-spin variants override per-bond via
    # interface_bond["bond_length"].
    "Zn-N":           2.05,    # zinc-imidazole / zinc-amine
    "Zn-S":           2.32,    # zinc-thiolate (Cys)
    "Zn-O":           1.95,    # zinc-aspartate
    "Sn-N":           2.10,    # tin-pyrrole / tin-amide
    "Sn-O":           2.00,    # tin-methoxy axial
    "Fe-N":           2.00,    # heme Fe-NHis (axial), porphyrin Fe-Npy
    "Fe-S":           2.30,    # rubredoxin / iron-sulfur Fe-Scys
    "Fe-O":           2.00,    # ferritin Fe-O / Fe-OH
    "Cu-N":           2.00,    # cupredoxin / Cu-amine
    "Cu-S":           2.30,    # plastocyanin Cu-Scys
    "Cu-O":           1.95,
    "Ni-N":           2.05,
    "Ni-S":           2.20,
    "Mg-O":           2.10,    # rubisco / kinase Mg-Asp
    "Mg-N":           2.20,
    "Ca-O":           2.40,    # calcium-binding loops
    "Mn-O":           2.10,    # photosystem II
    "Mn-N":           2.20,
    "Co-N":           1.95,    # cobalamin Co-N
    "Co-O":           1.95,
    "Pt-N":           2.00,    # cisplatin-style Pt-NHis / Pt-NMe
    "Pt-S":           2.30,    # Pt-thioether
    "Pt-Cl":          2.30,
    "Ru-N":           2.05,    # Ru-bpy / Ru-imidazole
    "Ru-O":           2.05,
}

# Hybridization at standard residue/atom anchors. Used to compute the
# "outward" direction. Extend this for new anchor types as needed.
HYBRIDIZATION_TABLE: dict[tuple[str, str], str] = {
    # GLU side-chain carbonyl (after OE2/HE2 removed it's an open sp2 site)
    ("GLU", "CD"):     "sp2_open",
    ("ASP", "CG"):     "sp2_open",
    # Cys thiol — sp3, becomes 1-bonded after HG removal
    ("CYS", "SG"):     "sp3_open",
    # His N-epsilon — sp2 imidazole, becomes 2-bonded after H removal
    ("HIS", "NE2"):    "sp2_imidazole",
    ("HIS", "ND1"):    "sp2_imidazole",
    # Lys epsilon-amine — sp3 (will form imine/amide with cofactor)
    ("LYS", "NZ"):     "sp3_open",
    # Tyr hydroxyl
    ("TYR", "OH"):     "sp3_open",
    # Ser/Thr hydroxyl
    ("SER", "OG"):     "sp3_open",
    ("THR", "OG1"):    "sp3_open",
}

# Default fragment-side hybridizations keyed by (resname or "*", atom_name).
# Fragment atoms after antechamber are typically renamed (NH2 → N1), so we
# also support keying by element-only or by metal symbol. The interface_bond
# can always override via `frag_hybrid` if the anchor's role differs from
# the default for that name.
FRAGMENT_HYBRIDIZATION_TABLE: dict[tuple[str, str], str] = {
    # SnP-builder original names (before antechamber renaming)
    ("*", "NH2"):      "sp2_open",          # aniline N pointing toward the peptide
    ("*", "N1"):       "sp2_open",          # antechamber-renamed equivalent (1st N in lib)
    # Common metal anchors — caller is expected to pass an appropriate
    # `bonded_positions` set. For metals coordinating the peptide, the
    # outward direction is the next vacant site of the metal coordination
    # sphere (varies with geometry — see compute_outward_direction notes).
    ("*", "ZN"):       "metal_tetrahedral", # tetrahedral by default (Zn-finger 2C2H)
    ("*", "ZN1"):      "metal_tetrahedral",
    ("*", "FE"):       "metal_axial",       # heme Fe (planar porphyrin + axials)
    ("*", "FE1"):      "metal_axial",
    ("*", "CU"):       "metal_axial",       # cupredoxins
    ("*", "CU1"):      "metal_axial",
    ("*", "NI"):       "metal_tetrahedral",
    ("*", "MG"):       "metal_octahedral",
    ("*", "CA"):       "metal_octahedral",
    ("*", "MN"):       "metal_octahedral",
    ("*", "CO"):       "metal_octahedral",
    ("*", "PT"):       "metal_square_planar",  # cisplatin-like
    ("*", "PT1"):      "metal_square_planar",
    ("*", "RU"):       "metal_octahedral",
    ("*", "SN"):       "metal_axial",       # SnP porphyrin axial site
    ("*", "SN1"):      "metal_axial",
}


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AtomCoord:
    name: str
    element: str
    resname: str
    resseq: int
    coord: np.ndarray   # shape (3,)


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

def parse_pdb_atoms(pdb_path: str) -> list[AtomCoord]:
    """Read a PDB file and return a list of AtomCoord objects.

    Lightweight parser — reads ATOM/HETATM lines only. Preserves atom names
    (including PDB whitespace stripping) and resnames so callers can look
    atoms up by (residue, atom).
    """
    if not pdb_path or not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"pdb_path not found: {pdb_path}")
    out: list[AtomCoord] = []
    with open(pdb_path) as f:
        for raw in f:
            if not raw.startswith(("ATOM  ", "HETATM")):
                continue
            name = raw[12:16].strip()
            resname = raw[17:20].strip()
            try:
                resseq = int(raw[22:26].strip())
            except ValueError:
                resseq = 0
            x = float(raw[30:38])
            y = float(raw[38:46])
            z = float(raw[46:54])
            element = raw[76:78].strip() or name[0]
            out.append(AtomCoord(
                name=name, element=element, resname=resname,
                resseq=resseq, coord=np.array([x, y, z], dtype=float),
            ))
    return out


def find_residue_atoms(atoms: list[AtomCoord], resseq: int,
                       resname: Optional[str] = None) -> list[AtomCoord]:
    return [a for a in atoms
            if a.resseq == resseq
            and (resname is None or a.resname == resname)]


def find_atom(atoms: list[AtomCoord], resseq: int, atom_name: str,
              resname: Optional[str] = None) -> AtomCoord:
    for a in atoms:
        if a.resseq != resseq:
            continue
        if resname is not None and a.resname != resname:
            continue
        if a.name == atom_name:
            return a
    raise ValueError(
        f"atom not found: resseq={resseq} resname={resname} atom={atom_name}"
    )


# ---------------------------------------------------------------------------
# Outward-direction computation
# ---------------------------------------------------------------------------

def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n < 1e-9:
        raise ValueError("zero-length vector — cannot unit-normalize")
    return v / n


def compute_outward_direction(
    anchor_pos: np.ndarray,
    bonded_positions: list[np.ndarray],
    hybridization: str,
) -> np.ndarray:
    """Return a unit vector pointing AWAY from the anchor's existing bonds,
    into which the new bond should extend.

    hybridization values:
      - "sp2_open"        : 1 vacant sp2 site, 2 existing bonds (e.g. GLU-CD
                            after OE2/HE2 removed; aniline N after caps removed).
                            outward = bisector pointing away from the two
                            existing bond directions.
      - "sp3_open"        : 1 vacant sp3 site, 1-3 existing bonds. Tetrahedral
                            completion.
      - "sp2_imidazole"   : Same as sp2_open but for ring N (e.g. HIS-NE2).
      - "metal_axial"     : axial site on a planar/octahedral metal — caller
                            should pass the ring/plane normal as
                            bonded_positions[0]; outward = +normal direction.
    """
    if not bonded_positions:
        raise ValueError("need at least one bonded position to compute outward")

    if hybridization in ("sp2_open", "sp2_imidazole"):
        # outward = -unit(sum of unit vectors from anchor to each bonded atom)
        #         = away from the bisector of existing bonds
        if len(bonded_positions) < 1:
            raise ValueError("sp2_open needs >= 1 bonded position")
        s = np.zeros(3)
        for b in bonded_positions:
            s += _unit(b - anchor_pos)
        return _unit(-s)

    if hybridization == "sp3_open":
        # Tetrahedral completion — outward fills the missing tetrahedron vertex.
        # For 1 bonded atom: outward = -unit(b->anchor) (just the back-axis;
        # the other 2 H's that complete the tetrahedron aren't placed by us).
        # For 3 bonded atoms: outward = -unit(sum of unit vectors), then
        # tetrahedral tilt off-axis.
        s = np.zeros(3)
        for b in bonded_positions:
            s += _unit(b - anchor_pos)
        if np.linalg.norm(s) < 1e-6:
            # bonded atoms are symmetric around anchor — pick an arbitrary axis
            return np.array([0.0, 0.0, 1.0])
        return _unit(-s)

    if hybridization == "metal_axial":
        # Heme / planar-coord: caller passes the plane normal as
        # bonded_positions[0] directly. The "outward" is +normal (the next
        # axial coordination site). For the OPPOSITE axial site, negate.
        return _unit(bonded_positions[0])

    if hybridization == "metal_tetrahedral":
        # 4-coordinate metal (e.g. Zn-finger 2C2H). Given existing ligand
        # positions, the next ligand goes at -unit(sum_of_existing_unit_vectors).
        # If 3 ligands already occupy 3 of the 4 vertices, this points to
        # the 4th. If fewer ligands are present, returns the bisector of
        # the existing ones (still a sensible "open" direction).
        if not bonded_positions:
            raise ValueError("metal_tetrahedral needs >= 1 bonded position")
        s = np.zeros(3)
        for b in bonded_positions:
            s += _unit(b - anchor_pos)
        if np.linalg.norm(s) < 1e-6:
            return np.array([0.0, 0.0, 1.0])
        return _unit(-s)

    if hybridization == "metal_square_planar":
        # 4 ligands in a plane (e.g. Pt-N4). Same outward formula as
        # tetrahedral when 3 ligands are present — the open site lies in
        # the plane of the other 3. With 4 already filled, there is no
        # in-plane outward direction; return the plane normal as a
        # fallback (axial approach).
        if not bonded_positions:
            raise ValueError("metal_square_planar needs >= 1 bonded position")
        if len(bonded_positions) >= 3:
            # Estimate plane normal from 3 ligands
            v1 = bonded_positions[1] - bonded_positions[0]
            v2 = bonded_positions[2] - bonded_positions[0]
            normal = np.cross(v1, v2)
            if np.linalg.norm(normal) > 1e-6:
                # If 4 ligands present, return normal (axial); otherwise
                # in-plane outward (-bisector projected into plane)
                if len(bonded_positions) == 4:
                    return _unit(normal)
        s = np.zeros(3)
        for b in bonded_positions:
            s += _unit(b - anchor_pos)
        if np.linalg.norm(s) < 1e-6:
            return np.array([0.0, 0.0, 1.0])
        return _unit(-s)

    if hybridization == "metal_octahedral":
        # 6-coordinate metal (e.g. Ca, Mg, Mn, Co, Ru). Same algorithmic
        # form as tetrahedral — outward = -unit(sum). With 5 ligands
        # present, this points to the 6th vertex.
        if not bonded_positions:
            raise ValueError("metal_octahedral needs >= 1 bonded position")
        s = np.zeros(3)
        for b in bonded_positions:
            s += _unit(b - anchor_pos)
        if np.linalg.norm(s) < 1e-6:
            return np.array([0.0, 0.0, 1.0])
        return _unit(-s)

    raise ValueError(f"unknown hybridization: {hybridization!r}")


# ---------------------------------------------------------------------------
# Rigid-body transformation (Rodrigues / cross-product alignment)
# ---------------------------------------------------------------------------

def rotation_aligning(v_from: np.ndarray, v_to: np.ndarray) -> np.ndarray:
    """Return a 3×3 rotation matrix that takes v_from → v_to (both unit).

    Uses Rodrigues' formula on the cross-product axis. Handles the
    parallel/antiparallel edge cases.
    """
    a = _unit(v_from)
    b = _unit(v_to)
    cos_t = float(np.dot(a, b))
    if cos_t > 1.0 - 1e-9:
        return np.eye(3)
    if cos_t < -1.0 + 1e-9:
        # 180° rotation: pick any axis perpendicular to a
        helper = np.eye(3)[int(np.argmin(np.abs(a)))]
        axis = _unit(np.cross(a, helper))
        return _rodrigues(axis, math.pi)
    axis = _unit(np.cross(a, b))
    angle = math.acos(cos_t)
    return _rodrigues(axis, angle)


def _rodrigues(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotation matrix for `angle` radians around unit `axis`."""
    a = _unit(axis)
    K = np.array([
        [0.0,    -a[2],  a[1]],
        [a[2],    0.0,  -a[0]],
        [-a[1],   a[0],  0.0],
    ])
    return np.eye(3) + math.sin(angle) * K + (1.0 - math.cos(angle)) * (K @ K)


# ---------------------------------------------------------------------------
# Bond-length lookup
# ---------------------------------------------------------------------------

def _infer_bond_kind(pep_element: str, frag_element: str) -> str:
    """Guess a STANDARD_BOND_LENGTHS key from anchor element pair."""
    pair = sorted([pep_element.upper(), frag_element.upper()])
    if pair == ["C", "N"]:
        return "C-N_amide"  # most common biological case
    if pair == ["C", "S"]:
        return "C-S_thioether"
    if pair == ["S", "S"]:
        return "S-S_disulfide"
    if pair == ["C", "O"]:
        return "C-O_ester"
    if pair == ["C", "C"]:
        return "C-C_single"
    if pair == ["N", "S"]:
        return "C-N_amine"  # rough default
    return "C-N_amide"  # final fallback


def lookup_bond_length(interface_bond: dict,
                       pep_element: str, frag_element: str,
                       overrides: Optional[dict[str, float]] = None) -> float:
    """Resolve the target bond length for an interface bond.

    Priority: explicit `interface_bond["bond_length"]` > overrides table
    keyed by `bond_kind` > STANDARD_BOND_LENGTHS keyed by `bond_kind` >
    auto-inferred from the anchor elements.
    """
    if "bond_length" in interface_bond:
        try:
            return float(interface_bond["bond_length"])
        except (TypeError, ValueError):
            pass
    kind = interface_bond.get("bond_kind") or _infer_bond_kind(pep_element, frag_element)
    if overrides and kind in overrides:
        return float(overrides[kind])
    return float(STANDARD_BOND_LENGTHS.get(kind, 1.50))


# ---------------------------------------------------------------------------
# Main transformation
# ---------------------------------------------------------------------------

def compute_rigid_transformation(
    peptide_pdb_path: str,
    fragment_pdb_path: str,
    interface_bond: dict,
    bond_length: Optional[float] = None,
    clash_optimize: bool = True,
    bond_length_overrides: Optional[dict[str, float]] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (translation, rotation) such that applying them to the fragment
    coordinates places ``frag_anchor`` at the right amide-bond geometry
    relative to ``pep_anchor`` in the peptide.

    Order of operations on the fragment: first rotate (pivoting around its
    own anchor), then translate.

    interface_bond keys (same schema as ep_fragment_fuse_topology, plus optional
    geometry hints):
      - pep_resid, pep_atom         : peptide anchor
      - frag_resid, frag_atom       : fragment anchor
      - pep_remove, frag_remove     : atoms to ignore when computing the
                                      outward direction (these are removed
                                      at fuse time anyway)
      - bond_kind (optional)        : key into STANDARD_BOND_LENGTHS
      - bond_length (optional)      : explicit override in Å
      - pep_hybrid (optional)       : override HYBRIDIZATION_TABLE
      - frag_hybrid (optional)      : override FRAGMENT_HYBRIDIZATION_TABLE
    """
    pep_atoms = parse_pdb_atoms(peptide_pdb_path)
    frag_atoms = parse_pdb_atoms(fragment_pdb_path)

    pep_resid = int(interface_bond["pep_resid"])
    pep_atom_name = interface_bond["pep_atom"]
    frag_resid = int(interface_bond["frag_resid"])
    frag_atom_name = interface_bond["frag_atom"]
    pep_remove = set(interface_bond.get("pep_remove", []))
    frag_remove = set(interface_bond.get("frag_remove", []))

    pep_anchor = find_atom(pep_atoms, pep_resid, pep_atom_name)
    frag_anchor = find_atom(frag_atoms, frag_resid, frag_atom_name)

    # Hybridization lookup
    pep_hybrid = interface_bond.get("pep_hybrid") or HYBRIDIZATION_TABLE.get(
        (pep_anchor.resname, pep_anchor.name), "sp2_open"
    )
    frag_hybrid = interface_bond.get("frag_hybrid") or (
        FRAGMENT_HYBRIDIZATION_TABLE.get((frag_anchor.resname, frag_anchor.name))
        or FRAGMENT_HYBRIDIZATION_TABLE.get(("*", frag_anchor.name))
        or "sp2_open"
    )

    # Existing-bond positions (drop the to-be-removed atoms, drop the anchor itself,
    # restrict to atoms within ~2.0 Å of the anchor — i.e., directly bonded)
    def _bonded_to(anchor: AtomCoord, all_atoms: list[AtomCoord],
                   skip_names: set[str]) -> list[np.ndarray]:
        bonded: list[np.ndarray] = []
        for a in all_atoms:
            if a.resseq != anchor.resseq:
                continue
            if a.name == anchor.name or a.name in skip_names:
                continue
            d = float(np.linalg.norm(a.coord - anchor.coord))
            if d < 2.0:  # generous covalent-bond cutoff
                bonded.append(a.coord)
        return bonded

    pep_bonded = _bonded_to(pep_anchor, pep_atoms, pep_remove)
    frag_bonded = _bonded_to(frag_anchor, frag_atoms, frag_remove)

    if not pep_bonded:
        raise ValueError(
            f"No remaining bonds at peptide anchor {pep_anchor.resname}{pep_anchor.resseq}.{pep_anchor.name}"
            f" — pep_remove may be too aggressive."
        )
    if not frag_bonded:
        raise ValueError(
            f"No remaining bonds at fragment anchor {frag_anchor.resname}{frag_anchor.resseq}.{frag_anchor.name}"
            f" — frag_remove may be too aggressive."
        )

    pep_outward = compute_outward_direction(
        pep_anchor.coord, pep_bonded, pep_hybrid,
    )
    frag_outward = compute_outward_direction(
        frag_anchor.coord, frag_bonded, frag_hybrid,
    )

    # Rotation: align frag_outward ANTI-parallel to pep_outward (so the bond
    # extends from frag_anchor BACK toward pep_anchor + L*pep_outward).
    rotation = rotation_aligning(frag_outward, -pep_outward)

    # Apply rotation around the fragment anchor
    rotated_anchor = rotation @ frag_anchor.coord

    # Translation: place rotated frag_anchor at pep_anchor + L * pep_outward
    L = bond_length if bond_length is not None else lookup_bond_length(
        interface_bond, pep_anchor.element, frag_anchor.element,
        overrides=bond_length_overrides,
    )
    target = pep_anchor.coord + L * pep_outward
    translation = target - rotated_anchor

    if clash_optimize:
        # Step 1: pick the better of the two sp-equivalent rotations.
        #
        # `rotation_aligning(frag_outward, -pep_outward)` returns the unique
        # minimum-angle rotation, but for sp2/sp3 anchors there are *two*
        # geometrically valid placements — the second is a 180° flip
        # around the bond axis. The sign of `compute_outward_direction`
        # for sp2_open is arbitrary (cross-product order), so the primary
        # rotation can land the bulk of the fragment on the wrong side
        # of the bond. The YASARA-designed SnP-peptide demo hits this:
        # the input PDB has SnP correctly positioned, but the primary
        # rotation flips the porphyrin into the helix face.
        #
        # Metric: count atom pairs within 2.0 Å, EXCLUDING the new
        # interface bond itself (the bond is at ~1.3 Å in both candidates,
        # so it would dominate any min-distance metric). The candidate
        # with fewer close contacts wins.
        bond_axis = pep_outward
        pivot = target
        flip_R = _rodrigues(bond_axis, math.pi)
        CLASH_R = 2.0  # Å — VMD's "unusual bond" cutoff is ~0.6×Σr_cov

        # Indices of atoms within 2 Å of the peptide anchor (i.e. directly
        # bonded to it including the new bond) — exclude these from the
        # clash count since they're chemically expected.
        pep_arr_all = np.asarray([a.coord for a in pep_atoms])
        pep_anchor_dists = np.linalg.norm(
            pep_arr_all - pep_anchor.coord, axis=1)
        pep_far_mask = pep_anchor_dists > CLASH_R
        pep_far = pep_arr_all[pep_far_mask]

        def _placed_with(rot: np.ndarray) -> list[np.ndarray]:
            r_anchor = rot @ frag_anchor.coord
            t = target - r_anchor
            return [rot @ a.coord + t for a in frag_atoms]

        def _clash_count(placed: list[np.ndarray]) -> int:
            """Count peptide-fragment atom pairs at < CLASH_R, excluding
            atoms directly bonded to the peptide anchor (which include
            the new amide bond and existing GLU side-chain atoms)."""
            if pep_far.size == 0 or not placed:
                return 0
            frag_arr = np.asarray(placed)
            d = np.linalg.norm(
                pep_far[:, None, :] - frag_arr[None, :, :], axis=-1)
            return int((d < CLASH_R).sum())

        original_clashes = _clash_count(_placed_with(rotation))
        flipped_clashes = _clash_count(_placed_with(flip_R @ rotation))
        if flipped_clashes < original_clashes:
            rotation = flip_R @ rotation
            rotated_anchor = rotation @ frag_anchor.coord
            translation = target - rotated_anchor

        # Step 2: secondary fine-grained scan around the bond axis (every
        # 30°) to relax any residual clashes.
        placed = [rotation @ a.coord + translation for a in frag_atoms]
        best_angle = find_clash_free_rotation(
            peptide_atoms=[a.coord for a in pep_atoms],
            fragment_atoms=placed,
            pivot=pivot,
            axis=bond_axis,
            n_steps=12,
        )
        if abs(best_angle) > 1e-6:
            secondary = _rodrigues(bond_axis, best_angle)
            rotation = secondary @ rotation
            rotated_anchor = rotation @ frag_anchor.coord
            translation = target - rotated_anchor

    return translation, rotation


def find_clash_free_rotation(
    peptide_atoms: list[np.ndarray],
    fragment_atoms: list[np.ndarray],
    pivot: np.ndarray,
    axis: np.ndarray,
    n_steps: int = 12,
) -> float:
    """Scan rotations around ``axis`` (passing through ``pivot``) and return
    the angle (radians) that maximizes the minimum peptide-fragment atom
    distance — i.e., picks the orientation with fewest steric clashes.

    Cheap O(n_steps × |pep| × |frag|) brute-force scan. n_steps=12 means
    every 30°. Returns 0.0 if no rotation improves the baseline.
    """
    if not peptide_atoms or not fragment_atoms:
        return 0.0
    pep_arr = np.asarray(peptide_atoms)
    best_angle = 0.0
    best_min_dist = -1.0
    for k in range(n_steps):
        angle = 2.0 * math.pi * k / n_steps
        if k == 0:
            rotated = np.asarray(fragment_atoms)
        else:
            R = _rodrigues(axis, angle)
            rotated = np.asarray([R @ (f - pivot) + pivot for f in fragment_atoms])
        # Min distance between any peptide and any fragment atom
        diffs = pep_arr[:, None, :] - rotated[None, :, :]
        d = np.linalg.norm(diffs, axis=-1)
        min_d = float(d.min())
        if min_d > best_min_dist:
            best_min_dist = min_d
            best_angle = angle
    return best_angle


# ---------------------------------------------------------------------------
# tleap script emission
# ---------------------------------------------------------------------------

def build_align_tleap_script(
    fragment_lib_basename: str,
    fragment_resname: str,
    translation: np.ndarray,
    rotation: np.ndarray,
    output_lib_basename: str,
    output_pdb_basename: str,
    fragment_frcmod_basename: Optional[str] = None,
) -> str:
    """Emit a tleap script that loads the fragment lib, applies translate +
    transform to the unit, and saves the aligned lib + PDB.

    `fragment_resname` is the unit name inside the lib (e.g. "mol" or "SNP").

    tleap's transform order is row-major; we emit the rotation matrix as
    {{ r11 r12 r13 } { r21 r22 r23 } { r31 r32 r33 }}. Translation comes
    BEFORE transform — tleap applies them in script order, so we transform
    (rotate around origin) and then translate to the target.

    NOTE: tleap's `transform` rotates around the origin (no pivot option).
    To rotate around the fragment anchor, the caller must include the
    pivot offset in the translation: pre-translate by -pivot, transform,
    translate by (target_pos). For our use the rotation matrix already
    incorporates this (compute_rigid_transformation returns the combined
    transformation), so we emit them as separate tleap commands.
    """
    r = rotation
    dx, dy, dz = float(translation[0]), float(translation[1]), float(translation[2])
    lines = [
        "logFile align.leap.log",
    ]
    if fragment_frcmod_basename:
        lines.append(f"loadamberparams {fragment_frcmod_basename}")
    lines.extend([
        f"loadoff {fragment_lib_basename}",
        f"transform {fragment_resname} {{"
        f" {{ {r[0,0]:.10f} {r[0,1]:.10f} {r[0,2]:.10f} }}"
        f" {{ {r[1,0]:.10f} {r[1,1]:.10f} {r[1,2]:.10f} }}"
        f" {{ {r[2,0]:.10f} {r[2,1]:.10f} {r[2,2]:.10f} }} }}",
        f"translate {fragment_resname} {{ {dx:.6f} {dy:.6f} {dz:.6f} }}",
        f"saveoff {fragment_resname} {output_lib_basename}",
        f"savepdb {fragment_resname} {output_pdb_basename}",
        "quit",
    ])
    return "\n".join(lines) + "\n"
