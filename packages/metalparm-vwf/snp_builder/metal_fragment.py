"""Shared geometry + I/O helpers for metal-fragment builders.

Extracted from `snp_builder/core.py` (originally porphyrin-specific) so a
future ``zn_finger_builder``, ``heme_builder``, etc. can reuse the same
primitives without forking the file. The split is:

  - **This module (``metal_fragment.py``)** — anything that does *not*
    depend on porphyrin-specific knowledge: data records, geometry
    primitives (planar SVD, perpendicular unit vectors, tetrahedral H
    placement), residue extraction, the metal-atom swap, axial ligand
    placement, peptide-side cleanup (YASARA atom-name rename, HIS
    tautomer inference), and file writers.

  - **``core.py``** — the SnP-specific orchestrator that knows about the
    4 pyrrole nitrogens, the aniline-NH cap, and the porphyrin ring
    plane.

When a second cofactor builder lands (Zn-finger, heme, …), it should:
  1. Either import from this module (when packaged together) or copy it
     into its own node dir (when standalone-installed); the public
     surface is small enough that copy-paste is acceptable until the
     marketplace install process is settled.
  2. Implement its own scaffold extraction (analogous to porphyrin's
     "find 4 pyrrole N → SVD → centroid + normal") and its own cap.
  3. Reuse ``swap_metal_at_centroid``, ``place_axial_ligand`` (when
     applicable), ``extract_peptide``, ``write_outputs``, and the
     ``BuildResult`` + ``AtomRec`` data classes verbatim.

Pure numpy + Bio.PDB (no RDKit). Sn-N / Zn-N dative bonds confuse RDKit
valence; we do geometric edits directly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom

# ---------------------------------------------------------------------------
# Generic chemistry constants (Angstroms / radians)
# ---------------------------------------------------------------------------
O_C_METHOXY = 1.42       # O-CH3
C_H_METHYL = 1.09        # C-H
N_C_AMIDE = 1.33         # N-C(=O) amide
C_O_CARBONYL = 1.23      # C=O
C_C_METHYL = 1.51        # C(=O)-CH3
N_C_METHYL = 1.45        # N-CH3 (NME cap methyl)
N_H_AMIDE = 1.01         # N-H amide
TETRAHEDRAL_ANGLE = math.radians(109.471)
TRIGONAL_ANGLE = math.radians(120.0)


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class AtomRec:
    """One PDB/XYZ atom: name, element symbol, coords."""
    name: str
    element: str
    coord: np.ndarray   # shape (3,)

    @classmethod
    def from_bio(cls, a: Atom) -> "AtomRec":
        return cls(name=a.get_name().strip(),
                   element=(a.element or a.get_name().strip()[0]).upper(),
                   coord=np.asarray(a.coord, dtype=float))


@dataclass
class BuildResult:
    """Output of a metal-fragment build: atom list + the cap atoms that
    should be removed at fuse time (so the fuse JSON's ``frag_remove``
    field can reference them by name)."""
    atoms: list[AtomRec]
    cap_atoms: list[str] = field(default_factory=list)

    def to_xyz(self, comment: str = "metal fragment") -> str:
        lines = [str(len(self.atoms)), comment]
        for a in self.atoms:
            lines.append(f"{a.element:<3s} {a.coord[0]:15.8f} "
                         f"{a.coord[1]:15.8f} {a.coord[2]:15.8f}")
        return "\n".join(lines) + "\n"

    def to_pdb(self, resname: str = "FRG") -> str:
        out = []
        for i, a in enumerate(self.atoms, start=1):
            name = pdb_atom_name(a.name, a.element)
            out.append(
                f"HETATM{i:5d} {name} {resname:>3s} A   1    "
                f"{a.coord[0]:8.3f}{a.coord[1]:8.3f}{a.coord[2]:8.3f}"
                f"  1.00  0.00          {a.element:>2s}"
            )
        out.append("END")
        return "\n".join(out) + "\n"


def pdb_atom_name(name: str, element: str) -> str:
    """Right/left-pad an atom name into the 4-column PDB atom-name field.

    PDB convention: 1-letter elements get a leading-space pad
    (``" N  "``); 2-letter elements take the full field (``"Sn  "``).
    Names already 4 chars are kept verbatim.
    """
    name = name.strip()
    if len(element) == 1 and len(name) < 4:
        return f" {name:<3s}"
    return f"{name:<4s}"


# ---------------------------------------------------------------------------
# Geometry primitives
# ---------------------------------------------------------------------------

def ring_plane(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(centroid, unit_normal)`` for a set of ~coplanar points.

    Normal is the singular vector with the smallest singular value of the
    centered point cloud. Sign is arbitrary but deterministic given the
    input order. Used for porphyrin (4 pyrrole N) but generalizes to any
    planar set of ≥3 ligands.
    """
    centroid = points.mean(axis=0)
    _, _, vh = np.linalg.svd(points - centroid, full_matrices=False)
    normal = vh[-1] / np.linalg.norm(vh[-1])
    return centroid, normal


def perpendicular_unit(v: np.ndarray) -> np.ndarray:
    """Return an arbitrary unit vector perpendicular to v."""
    v = v / np.linalg.norm(v)
    helper = np.eye(3)[int(np.argmin(np.abs(v)))]
    perp = np.cross(v, helper)
    return perp / np.linalg.norm(perp)


def tetrahedral_h(anchor: np.ndarray, bonded_to: np.ndarray,
                  bond_len: float, n: int) -> list[np.ndarray]:
    """Place ``n`` hydrogens around ``anchor`` at tetrahedral-ish positions,
    with the already-existing bond going to ``bonded_to`` occupying one
    vertex.

    n in {1, 3}. For n=3 we emit 3 H's spaced 120° around the axis from
    ``bonded_to`` to ``anchor``, tilted by ``180 - tetrahedral``.
    """
    axis = anchor - bonded_to
    axis /= np.linalg.norm(axis)
    tilt = math.pi - TETRAHEDRAL_ANGLE
    perp = perpendicular_unit(axis)
    result = []
    for k in range(n):
        phi = 2 * math.pi * k / n
        perp_rot = (perp * math.cos(phi)
                    + np.cross(axis, perp) * math.sin(phi)
                    + axis * np.dot(axis, perp) * (1 - math.cos(phi)))
        direction = axis * math.cos(tilt) + perp_rot * math.sin(tilt)
        result.append(anchor + bond_len * direction)
    return result


# ---------------------------------------------------------------------------
# PDB residue extraction + metal swap
# ---------------------------------------------------------------------------

def extract_residue_atoms(input_pdb: str | Path,
                          resname: str = "UNK") -> list[AtomRec]:
    """Read a PDB and return all atoms of the named residue as AtomRec.

    Used as the "find the cofactor" step. ``resname`` is the 3-char
    residue name in the input PDB (typically ``UNK`` for YASARA-extracted
    metal sites).
    """
    structure = PDBParser(QUIET=True).get_structure("in", str(input_pdb))
    out: list[AtomRec] = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname().strip() == resname:
                    for atom in res:
                        out.append(AtomRec.from_bio(atom))
    if not out:
        raise ValueError(f"No {resname} residue found in {input_pdb}")
    return out


def swap_metal_at_centroid(atoms: list[AtomRec],
                           centroid: np.ndarray,
                           metal_in: str,
                           metal_out: str) -> list[AtomRec]:
    """Drop the placeholder metal (atom name ``metal_in``) and append a
    new metal (atom name ``metal_out``, element ``metal_out.capitalize()``)
    at ``centroid``.

    Atom NAME stays in PDB convention (4-char field, all-caps for metals
    like "SN", "ZN"). Atom ELEMENT must be standard symbol case ("Sn",
    "Zn") so element-symbol-keyed lookups downstream (e.g. covalent radii
    in 02_get_bond_angle.py) work correctly.

    Returns a NEW list (does not mutate the input). Raises if the
    placeholder metal is not present.
    """
    out = [a for a in atoms if a.name != metal_in]
    if len(out) == len(atoms):
        raise ValueError(f"Metal atom {metal_in} not found")
    out.append(AtomRec(name=metal_out,
                       element=metal_out.capitalize(),
                       coord=centroid.copy()))
    return out


# ---------------------------------------------------------------------------
# Axial ligand placement
# ---------------------------------------------------------------------------

def place_axial_ligand(centroid: np.ndarray,
                       normal: np.ndarray,
                       ligand: str,
                       bond_len: float,
                       *,
                       sign: int = +1,
                       tag: str = "1") -> list[AtomRec]:
    """Place one axial ligand at ``centroid + sign * bond_len * normal``.

    Supported ligands:
      - ``"OMe"``: methoxy (O + C + 3H)
      - ``"OH"``:  hydroxyl (O + H)
      - ``"Cl"``:  chloride (Cl, replaces O — Sn-Cl ≈ 2.4 Å)
      - ``"none"``: no atoms emitted

    ``tag`` is appended to atom names to disambiguate the two axials
    (e.g. ``"1"`` and ``"2"`` for the upper and lower axial in a
    porphyrin).
    """
    if ligand == "none":
        return []
    if ligand == "Cl":
        # Replace the would-be O with Cl directly. 2.4 Å is the Sn-Cl
        # equilibrium; for other metals adjust the literal at the call site.
        return [AtomRec(name=f"CL{tag}", element="CL",
                        coord=centroid + sign * 2.40 * normal)]

    out: list[AtomRec] = []
    o_pos = centroid + sign * bond_len * normal
    out.append(AtomRec(name=f"O{tag}", element="O", coord=o_pos))

    if ligand == "OMe":
        c_pos = o_pos + sign * O_C_METHOXY * normal
        out.append(AtomRec(name=f"CM{tag}", element="C", coord=c_pos))
        for i, h_pos in enumerate(
                tetrahedral_h(c_pos, o_pos, C_H_METHYL, 3), start=1):
            out.append(AtomRec(name=f"H{tag}{i}", element="H", coord=h_pos))
    elif ligand == "OH":
        h_pos = o_pos + sign * C_H_METHYL * normal
        out.append(AtomRec(name=f"HO{tag}", element="H", coord=h_pos))
    else:
        raise ValueError(f"Unknown axial ligand={ligand!r}")
    return out


# ---------------------------------------------------------------------------
# Peptide-side cleanups
# ---------------------------------------------------------------------------

# YASARA-style PDBs use non-standard atom names for protonated acid /
# hydroxyl Hs. Map them to ff19SB names so tleap's residue templates
# apply cleanly when loadpdb-ing the extracted peptide.
YASARA_ATOM_RENAMES: dict[tuple[str, str], str] = {
    ("GLU", "COOH"): "HE2",
    ("ASP", "COOH"): "HD2",
    ("SER", "HO"):   "HG",
    ("THR", "HO"):   "HG1",
    ("TYR", "HO"):   "HH",
    ("CYS", "HS"):   "HG",
}


def _rename_atom_in_pdb_line(line: str) -> str:
    """If ``(resname, atom)`` on this PDB line is in YASARA_ATOM_RENAMES,
    rewrite the atom-name field (cols 13-16) to the ff19SB equivalent.
    """
    if not line.startswith(("ATOM  ", "HETATM")):
        return line
    resname = line[17:20].strip()
    atom = line[12:16].strip()
    new_atom = YASARA_ATOM_RENAMES.get((resname, atom))
    if new_atom is None:
        return line
    if len(new_atom) <= 3:
        formatted = f" {new_atom:<3s}"
    else:
        formatted = f"{new_atom:<4s}"
    return line[:12] + formatted + line[16:]


def infer_his_tautomer_renames(input_pdb: str | Path) -> dict[tuple[str, str], str]:
    """Inspect a PDB and decide per-HIS-residue whether to rewrite its
    name to HID, HIE, or HIP based on which protonation hydrogens are
    present.

    Default ff19SB behaviour: a residue named ``HIS`` is loaded as NHIE
    (epsilon-protonated, has HE2 on NE2). If the input PDB carries an
    HD1 atom on a HIS residue (delta-protonated tautomer), tleap fails
    with *"Atom .R<NHIE 1>.A<HD1 20> does not have a type"*.

    Decision matrix:
      - HD1 present, HE2 absent  → rename HIS → HID  (delta tautomer)
      - HD1 absent,  HE2 present → leave as HIS (ff19SB default = HIE)
      - both present              → rename HIS → HIP  (doubly protonated)
      - neither                   → leave as HIS (tleap rebuilds protons)
    """
    src = Path(input_pdb)
    if not src.is_file():
        raise FileNotFoundError(f"input_pdb not found: {src}")

    his_atoms: dict[tuple[str, str], set[str]] = {}
    for raw in src.read_text().splitlines():
        if not raw.startswith(("ATOM  ", "HETATM")):
            continue
        if raw[17:20].strip() != "HIS":
            continue
        chain = raw[21]
        resseq = raw[22:26].strip()
        atom = raw[12:16].strip()
        his_atoms.setdefault((chain, resseq), set()).add(atom)

    decisions: dict[tuple[str, str], str] = {}
    for key, atoms in his_atoms.items():
        has_hd1 = "HD1" in atoms
        has_he2 = "HE2" in atoms
        if has_hd1 and has_he2:
            decisions[key] = "HIP"
        elif has_hd1:
            decisions[key] = "HID"
    return decisions


# ---------------------------------------------------------------------------
# Terminal capping (ACE / NME) — geometric placement via internal coordinates
# ---------------------------------------------------------------------------

def _place_from_internal(a: np.ndarray, b: np.ndarray, c: np.ndarray,
                         dist: float, angle_deg: float,
                         dihedral_deg: float) -> np.ndarray:
    """NeRF-style placement: return the position of atom D bonded to C,
    with ``|C-D| = dist``, bond angle ``B-C-D = angle_deg``, and dihedral
    ``A-B-C-D = dihedral_deg`` (all anchored on the three given points).
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    theta = math.radians(angle_deg)
    phi = math.radians(dihedral_deg)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d_local = (
        (-dist * math.cos(theta)) * bc
        + (dist * math.sin(theta) * math.cos(phi)) * m
        + (dist * math.sin(theta) * math.sin(phi)) * n
    )
    return c + d_local


def _fmt_atom_line(serial: int, name: str, element: str, resname: str,
                   chain: str, resseq: int, coord: np.ndarray) -> str:
    """Format one ``ATOM`` record in fixed PDB columns."""
    name4 = pdb_atom_name(name, element)
    return (
        f"ATOM  {serial:>5d} {name4} {resname:>3s} {chain:1s}{resseq:>4d}    "
        f"{coord[0]:>8.3f}{coord[1]:>8.3f}{coord[2]:>8.3f}"
        f"  1.00  0.00          {element:>2s}"
    )


def _backbone_coords(res_lines: list[str]) -> dict[str, np.ndarray]:
    """Map atom name → coord for one residue's PDB lines."""
    out: dict[str, np.ndarray] = {}
    for ln in res_lines:
        nm = ln[12:16].strip()
        out[nm] = np.array(
            [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return out


def _build_ace_cap(bb: dict[str, np.ndarray]) -> list[tuple[str, str, np.ndarray]]:
    """Build the 6 atoms of an ACE N-cap (acetyl) for a residue whose
    backbone N/CA/C coords are in ``bb``. Returns (name, element, coord).

    The cap's carbonyl C bonds to the residue N; methyl + carbonyl O are
    placed in the amide plane; methyl Hs via the tetrahedral helper.
    Geometry continues a helical-ish backbone (φ ≈ -60°); exact values
    are not critical — the structure is energy-minimized downstream.
    """
    for need in ("N", "CA", "C"):
        if need not in bb:
            raise ValueError(
                f"cannot ACE-cap: first residue missing backbone atom {need!r}")
    n, ca, c = bb["N"], bb["CA"], bb["C"]
    # carbonyl C: bonded to N, angle CA-N-C', dihedral C-CA-N-C' = φ
    cac = _place_from_internal(c, ca, n, N_C_AMIDE, 121.9, -60.0)
    # carbonyl O: in the amide plane, trans to CA across the C'-N bond
    o = _place_from_internal(ca, n, cac, C_O_CARBONYL, 122.9, 180.0)
    # methyl C: the third sp2 substituent of C', cis to CA (opposite O)
    ch3 = _place_from_internal(ca, n, cac, C_C_METHYL, 116.6, 0.0)
    # ff19SB ACE template names: carbonyl C/O, methyl CH3, methyl H1/H2/H3.
    atoms = [("C", "C", cac), ("O", "O", o), ("CH3", "C", ch3)]
    for i, h in enumerate(tetrahedral_h(ch3, cac, C_H_METHYL, 3), start=1):
        atoms.append((f"H{i}", "H", h))
    return atoms


def _build_nme_cap(bb: dict[str, np.ndarray]) -> list[tuple[str, str, np.ndarray]]:
    """Build the 6 atoms of an NME C-cap (N-methylamide) for a residue
    whose backbone N/CA/C coords are in ``bb``. Returns (name, element,
    coord). The cap N bonds to the residue C; amide H + methyl + methyl
    Hs follow. Geometry continues a helical-ish backbone (ψ ≈ -47°).
    """
    for need in ("N", "CA", "C"):
        if need not in bb:
            raise ValueError(
                f"cannot NME-cap: last residue missing backbone atom {need!r}")
    n, ca, c = bb["N"], bb["CA"], bb["C"]
    # cap N: bonded to C, angle CA-C-N, dihedral N-CA-C-N' = ψ
    nme_n = _place_from_internal(n, ca, c, N_C_AMIDE, 116.6, -47.0)
    # methyl C: trans to CA across the C-N' bond
    ch3 = _place_from_internal(ca, c, nme_n, N_C_METHYL, 121.5, 180.0)
    # amide H: third substituent of the sp2 N, cis to CA (opposite methyl)
    h = _place_from_internal(ca, c, nme_n, N_H_AMIDE, 119.0, 0.0)
    # ff19SB NME template names: amide N/H, methyl carbon `C`, methyl H1/H2/H3.
    atoms = [("N", "N", nme_n), ("H", "H", h), ("C", "C", ch3)]
    for i, hh in enumerate(tetrahedral_h(ch3, nme_n, C_H_METHYL, 3), start=1):
        atoms.append((f"H{i}", "H", hh))
    return atoms


# N-terminal protons of a free (charged) N-terminus — dropped when an
# ACE cap is added so the capped backbone N keeps only its amide H
# (tleap rebuilds that single H from the mid-chain template).
_NTERM_PROTONS = {"H1", "H2", "H3", "HT1", "HT2", "HT3"}
# C-terminal extra oxygens — dropped when an NME cap is added.
_CTERM_OXYGENS = {"OXT", "OT1", "OT2"}


def _apply_terminal_caps(atom_lines: list[str]) -> list[str]:
    """Given the kept ATOM/HETATM lines of a single peptide chain (in
    file order), prepend an ACE cap before the first residue and append
    an NME cap after the last, and drop the now-redundant
    charged-terminus atoms.

    Every residue is **renumbered sequentially from 1** (ACE = 1,
    peptide = 2…n+1, NME = n+2) and every serial is renumbered too. The
    sequential residue numbering matters: tleap addresses residues by
    their stored sequence number (``unit.N``), so a non-contiguous
    numbering (e.g. ACE at resSeq 0) would make ``cpx.7`` resolve to the
    wrong residue at fuse time.

    Returns the new line list (ATOM records only; no TER/END).
    """
    if not atom_lines:
        raise ValueError("cannot cap an empty peptide")

    def _resseq(ln: str) -> str:
        return ln[22:26].strip()

    order: list[str] = []
    for ln in atom_lines:
        rs = _resseq(ln)
        if rs not in order:
            order.append(rs)
    first_res, last_res = order[0], order[-1]
    first_lines = [ln for ln in atom_lines if _resseq(ln) == first_res]
    last_lines = [ln for ln in atom_lines if _resseq(ln) == last_res]

    chain = atom_lines[0][21]
    ace = _build_ace_cap(_backbone_coords(first_lines))
    nme = _build_nme_cap(_backbone_coords(last_lines))

    # Sequential residue numbers: ACE=1, peptide=2..n+1, NME=n+2.
    renum = {rs: i + 2 for i, rs in enumerate(order)}
    nme_num = len(order) + 2

    out: list[str] = []
    serial = 0

    def _emit_cap(cap, resname, resseq):
        nonlocal serial
        for name, element, coord in cap:
            serial += 1
            out.append(_fmt_atom_line(serial, name, element, resname,
                                      chain, resseq, coord))

    _emit_cap(ace, "ACE", 1)
    for ln in atom_lines:
        nm = ln[12:16].strip()
        rs = _resseq(ln)
        if rs == first_res and nm in _NTERM_PROTONS:
            continue  # ACE cap replaces the charged N-terminus
        if rs == last_res and nm in _CTERM_OXYGENS:
            continue  # NME cap replaces the charged C-terminus
        serial += 1
        out.append(ln[:6] + f"{serial:>5d}" + ln[11:22]
                   + f"{renum[rs]:>4d}" + ln[26:])
    _emit_cap(nme, "NME", nme_num)
    return out


def extract_peptide(input_pdb: str | Path,
                    output_path: str | Path,
                    *,
                    unk_resname: str = "UNK",
                    rename_yasara_atoms: bool = True,
                    rename_his_tautomers: bool = True,
                    residue_range: tuple[int, int] | None = None,
                    cap_termini: bool = False) -> int:
    """Write everything in ``input_pdb`` except the ``unk_resname``
    residue to ``output_path``. Returns the number of distinct residues
    written.

    Used by metal-fragment builder nodes to produce a peptide-only PDB
    whose coordinates are in the same frame as the metal fragment — so
    the fused complex starts with chemically sensible interface geometry
    without any minimization step.

    Preserves ATOM/HETATM lines verbatim (modulo the YASARA atom-name
    and HIS-tautomer fixes); keeps TER records; drops everything else
    (HEADER, REMARK, LINK, CONECT, etc.) so tleap loadpdb sees only
    standard residue records.

    Args:
      rename_yasara_atoms: rename non-standard YASARA atom names (e.g.
        ``GLU.COOH → GLU.HE2``) so the output is ff19SB-template-clean.
      rename_his_tautomers: rewrite HIS residue names to HID/HIP based
        on which protonation hydrogens are present (HD1, HE2) so tleap
        loads the correct ff19SB template.
      residue_range: optional ``(lo, hi)`` inclusive span of PDB
        residue-sequence numbers. When given, only residues whose
        ``resSeq`` parses to an int in ``[lo, hi]`` are written — used
        to carve a sub-peptide out of a longer chain (e.g. ``(1, 7)``
        for the Case 1 heptapeptide = residues 1–7 of ``snpp.pdb``).
        Original ``TER`` records are dropped and a single synthetic
        ``TER`` is emitted after the kept atoms, so the truncated chain
        terminates cleanly for ``tleap loadpdb``. ``None`` = whole chain.
      cap_termini: when True, geometrically place an **ACE** N-cap
        before the first written residue and an **NME** C-cap after the
        last (`Ac-…-NH-CH₃`), and drop the now-redundant charged-terminus
        atoms (N-terminal H1/H2/H3, C-terminal OXT). Used for the Case 1
        heptapeptide so the folding MD starts from a properly capped
        peptide rather than a zwitterion. Assumes a single chain. The
        returned residue count includes the 2 caps.
    """
    src = Path(input_pdb)
    if not src.is_file():
        raise FileNotFoundError(f"input_pdb not found: {src}")
    if residue_range is not None and residue_range[0] > residue_range[1]:
        raise ValueError(
            f"residue_range lo > hi: {residue_range}")

    his_renames: dict[tuple[str, str], str] = {}
    if rename_his_tautomers:
        his_renames = infer_his_tautomer_renames(src)

    kept_lines: list[str] = []
    seen_residues: set[tuple[str, str, str]] = set()
    for raw in src.read_text().splitlines():
        if raw.startswith(("ATOM  ", "HETATM")):
            resname = raw[17:20].strip()
            if resname == unk_resname:
                continue
            chain = raw[21]
            resseq = raw[22:26].strip()
            if residue_range is not None:
                try:
                    resnum = int(resseq)
                except ValueError:
                    continue  # non-numeric resSeq can't be range-matched
                if not (residue_range[0] <= resnum <= residue_range[1]):
                    continue
            new_resname = his_renames.get((chain, resseq))
            if new_resname is not None and resname == "HIS":
                raw = raw[:17] + f"{new_resname:<3s}" + raw[20:]
                resname = new_resname
            seen_residues.add((chain, resseq, resname))
            if rename_yasara_atoms:
                raw = _rename_atom_in_pdb_line(raw)
            kept_lines.append(raw)
        elif raw.startswith("TER"):
            # In range mode the original TER (at the chain's true end) is
            # outside the kept span — emit one synthetic TER after the loop.
            if residue_range is None:
                kept_lines.append(raw)
    if residue_range is not None:
        if not seen_residues:
            raise ValueError(
                f"residue_range {residue_range} matched no residues in "
                f"{src.name} — check the residue-number span")
        kept_lines.append("TER")
    kept_lines.append("END")

    n_residues = len(seen_residues)
    if cap_termini:
        # Re-derive from the kept ATOM records, then prepend ACE / append
        # NME. The capped chain gets a single trailing TER.
        atom_lines = [l for l in kept_lines
                      if l.startswith(("ATOM  ", "HETATM"))]
        kept_lines = _apply_terminal_caps(atom_lines) + ["TER", "END"]
        n_residues += 2  # ACE + NME

    dst = Path(output_path)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(kept_lines) + "\n")
    return n_residues


# ---------------------------------------------------------------------------
# File writers
# ---------------------------------------------------------------------------

def write_outputs(result: BuildResult,
                  out_dir: str | Path,
                  basename: str,
                  resname: str) -> dict[str, str]:
    """Write ``<basename>.xyz`` + ``<basename>.pdb`` into ``out_dir``.

    Returns a dict with the resolved paths and the cap-atom list (for
    forwarding into BoCoFlow result.data).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    xyz_path = out_dir / f"{basename}.xyz"
    pdb_path = out_dir / f"{basename}.pdb"
    xyz_path.write_text(result.to_xyz(
        comment=f"{resname} fragment, {len(result.atoms)} atoms"))
    pdb_path.write_text(result.to_pdb(resname=resname))
    return {"xyz": str(xyz_path), "pdb": str(pdb_path),
            "cap_atoms": list(result.cap_atoms)}
