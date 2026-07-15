"""SnP fragment builder — convert a YASARA-style Zn-porphyrin UNK residue into
a capped Sn(IV)(OMe)2-porphyrin fragment ready for EasyParm.

This file is the **porphyrin-specific** orchestrator. The generic
geometry / I/O / peptide-cleanup helpers live in ``metal_fragment.py``;
both files together compose ``snp_builder``. When a second cofactor
builder lands (Zn-finger, heme, …), ``metal_fragment.py`` is the
copyable / re-importable shared core.

Algorithm
---------
1. Extract the UNK residue atoms via ``extract_residue_atoms``.
2. Identify the 4 pyrrole nitrogens (``PYRROLE_N_NAMES``), compute
   centroid + plane normal via SVD (``ring_plane``).
3. Drop the placeholder metal (default ``ZN``) and insert the new
   metal (default ``SN``) at the centroid (``swap_metal_at_centroid``).
4. Add 2 axial ligands along ±n_hat (``place_axial_ligand``; default
   methoxy at 2.00 Å).
5. Apply the porphyrin-aniline cap so the NH₂ nitrogen sees amide
   electronics during QM (``_apply_cap``, the porphyrin-specific bit).

Outputs:
    snp_frag.xyz  — ORCA/antechamber-ready XYZ
    snp_frag.pdb  — same atoms with residue name SNP (for visualization)
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

# Import the shared geometry / IO layer. Three-tier fallback mirrors the
# node-runner pattern in node.py: package-relative when imported as a
# package, direct-path when node_runner adds the dir to sys.path, and
# importlib.util when the test suite spec_from_file_location-loads
# core.py standalone (test_snp_builder.py:22).
try:
    from . import metal_fragment as _mf
except ImportError:
    try:
        import metal_fragment as _mf  # type: ignore[no-redef]
    except ImportError:
        import importlib.util as _ilu
        import sys as _sys
        from pathlib import Path as _Path
        _spec = _ilu.spec_from_file_location(
            "snp_builder_metal_fragment",
            _Path(__file__).with_name("metal_fragment.py"),
        )
        _mf = _ilu.module_from_spec(_spec)  # type: ignore[assignment]
        # Must register in sys.modules BEFORE exec_module so dataclass can
        # resolve the module via cls.__module__ during class construction.
        _sys.modules["snp_builder_metal_fragment"] = _mf
        _spec.loader.exec_module(_mf)  # type: ignore[union-attr]

# Re-export shared symbols for back-compat with callers that referenced
# them off snp_builder.core in the past.
AtomRec = _mf.AtomRec
BuildResult = _mf.BuildResult
C_C_METHYL = _mf.C_C_METHYL
C_H_METHYL = _mf.C_H_METHYL
C_O_CARBONYL = _mf.C_O_CARBONYL
N_C_AMIDE = _mf.N_C_AMIDE
O_C_METHOXY = _mf.O_C_METHOXY
TETRAHEDRAL_ANGLE = _mf.TETRAHEDRAL_ANGLE
TRIGONAL_ANGLE = _mf.TRIGONAL_ANGLE
YASARA_ATOM_RENAMES = _mf.YASARA_ATOM_RENAMES
extract_peptide = _mf.extract_peptide
extract_residue_atoms = _mf.extract_residue_atoms
infer_his_tautomer_renames = _mf.infer_his_tautomer_renames
_pdb_atom_name = _mf.pdb_atom_name
place_axial_ligand = _mf.place_axial_ligand
ring_plane = _mf.ring_plane
swap_metal_at_centroid = _mf.swap_metal_at_centroid
tetrahedral_h = _mf.tetrahedral_h
_generic_write_outputs = _mf.write_outputs


# ---------------------------------------------------------------------------
# Sn(IV)-porphyrin specific constants
# ---------------------------------------------------------------------------
SN_O_AXIAL = 2.00        # Sn(IV)-O (methoxy)

# 4 pyrrole N atom names in YASARA-style PDBs (the metal sits at their
# centroid; SVD of these gives the porphyrin plane normal).
PYRROLE_N_NAMES = ("N", "N1", "N2", "N3")


# ---------------------------------------------------------------------------
# Main entry — porphyrin-specific orchestrator
# ---------------------------------------------------------------------------

def build_snp_fragment(
    input_pdb: str | Path,
    *,
    unk_resname: str = "UNK",
    metal_in: str = "ZN",
    metal_out: str = "SN",
    axial_ligand: str = "OMe",      # "OMe" | "OH" | "Cl" | "none"
    axial_bond_len: float = SN_O_AXIAL,
    cap_style: str = "ACE",         # "ACE" | "H" | "NHMe"
) -> BuildResult:
    """Read ``input_pdb``, extract the ``unk_resname`` residue atoms,
    swap the placeholder metal for ``metal_out``, add 2 axial ligands
    along the porphyrin normal, apply ``cap_style`` to the aniline
    nitrogen, and return a BuildResult.
    """
    unk_atoms = extract_residue_atoms(input_pdb, resname=unk_resname)

    # Identify pyrrole nitrogens (porphyrin-specific: exactly 4 expected)
    ring_n = [a for a in unk_atoms if a.name in PYRROLE_N_NAMES]
    if len(ring_n) != 4:
        raise ValueError(
            f"Expected 4 pyrrole N ({PYRROLE_N_NAMES}), found {len(ring_n)}: "
            f"{[a.name for a in ring_n]}")
    centroid, normal = ring_plane(np.vstack([a.coord for a in ring_n]))

    # Generic: drop placeholder metal, insert new metal at centroid
    atoms = swap_metal_at_centroid(unk_atoms, centroid, metal_in, metal_out)

    # Generic: place 2 axial ligands along ±normal
    for sign, tag in [(+1, "1"), (-1, "2")]:
        atoms.extend(place_axial_ligand(centroid, normal, axial_ligand,
                                        axial_bond_len, sign=sign, tag=tag))

    # Porphyrin-specific: cap the aniline NH (atom name NH2)
    cap_names = _apply_cap(atoms, cap_style=cap_style)

    return BuildResult(atoms=atoms, cap_atoms=cap_names)


# ---------------------------------------------------------------------------
# Porphyrin-aniline cap (the only truly SnP-specific bit)
# ---------------------------------------------------------------------------

def _apply_cap(atoms: list[AtomRec], cap_style: str) -> list[str]:
    """Mutate ``atoms`` in place by adding cap atoms on the aniline N
    (atom name ``NH2``).

    The aniline N is sp² (amide-like). Its three substituents — the
    bonded phenyl C, any existing H, and the new cap — must sit at
    ~120° from each other in the amide plane. We use the N–C(phenyl)
    bond as the backbone direction and place the cap opposite any
    existing H so the geometry is truly trigonal planar.

    Returns the list of cap atom names — these are the atoms the fuse
    node will remove at runtime.

    Porphyrin-specific because:
      - depends on the ``NH2`` atom name in the input PDB (YASARA's
        convention for the aniline NH₂ on the SnP linker);
      - uses the porphyrin ring normal as the side-axis fallback when
        the NH₂ has no H to point opposite to.
    """
    nh2 = next((a for a in atoms if a.name == "NH2"), None)
    if nh2 is None:
        raise ValueError("Aniline nitrogen (atom name NH2) not found")

    # Backbone: bonded phenyl C (nearest C within ~1.6 Å)
    c_phe = min(
        (a for a in atoms if a.element == "C"),
        key=lambda a: float(np.linalg.norm(a.coord - nh2.coord)),
    )
    n_to_c = c_phe.coord - nh2.coord
    n_to_c /= np.linalg.norm(n_to_c)
    extension_dir = -n_to_c  # N away from phenyl C

    # Porphyrin normal — fallback for side axis when NH2 has no H
    ring = np.vstack([a.coord for a in atoms if a.name in PYRROLE_N_NAMES])
    _, ring_normal = ring_plane(ring)
    ring_normal /= np.linalg.norm(ring_normal)

    # Side axis in the amide plane: place new cap opposite any existing H on N
    h_on_n = [a for a in atoms if a.element == "H"
              and float(np.linalg.norm(a.coord - nh2.coord)) < 1.15]

    def _side_axis_opposite_h() -> np.ndarray:
        if h_on_n:
            h_rel = h_on_n[0].coord - nh2.coord
            h_perp = h_rel - np.dot(h_rel, n_to_c) * n_to_c
            if np.linalg.norm(h_perp) > 1e-6:
                return -h_perp / np.linalg.norm(h_perp)
        s = np.cross(n_to_c, ring_normal)
        return s / np.linalg.norm(s)

    half = math.pi - TRIGONAL_ANGLE  # 60° — rotation from N→(away from phenyl)
    side = _side_axis_opposite_h()

    if cap_style == "H":
        # Keep existing Hs; if only one, add a second at the proper sp² position.
        if len(h_on_n) < 2:
            new_h_dir = math.cos(half) * extension_dir + math.sin(half) * side
            new_h_dir /= np.linalg.norm(new_h_dir)
            atoms.append(AtomRec(name="HH2B", element="H",
                                 coord=nh2.coord + 1.01 * new_h_dir))
        return []

    # ACE and NHMe: the new C substituent sits at 120° from C_phe, opposite H.
    c_from_n = math.cos(half) * extension_dir + math.sin(half) * side
    c_from_n /= np.linalg.norm(c_from_n)

    if cap_style == "NHMe":
        cm_pos = nh2.coord + N_C_AMIDE * c_from_n
        atoms.append(AtomRec(name="CM", element="C", coord=cm_pos))
        for i, h_pos in enumerate(
                tetrahedral_h(cm_pos, nh2.coord, C_H_METHYL, 3), start=1):
            atoms.append(AtomRec(name=f"HM{i}", element="H", coord=h_pos))
        return ["CM", "HM1", "HM2", "HM3"]

    if cap_style == "ACE":
        cap_pos = nh2.coord + N_C_AMIDE * c_from_n
        atoms.append(AtomRec(name="CAP", element="C", coord=cap_pos))

        # At CAP, the amide plane is spanned by (extension_dir, side). OAP and
        # CM radiate at 120° from the CAP→N bond (= -c_from_n), i.e. 60° from
        # +c_from_n in the amide plane.
        n_plane = np.cross(extension_dir, side)
        if np.linalg.norm(n_plane) < 1e-6:
            n_plane = ring_normal
        n_plane /= np.linalg.norm(n_plane)
        perp = np.cross(n_plane, c_from_n)
        perp /= np.linalg.norm(perp)

        oap_dir = math.cos(half) * c_from_n + math.sin(half) * perp
        oap_dir /= np.linalg.norm(oap_dir)
        oap_pos = cap_pos + C_O_CARBONYL * oap_dir
        atoms.append(AtomRec(name="OAP", element="O", coord=oap_pos))

        cm_dir = math.cos(half) * c_from_n - math.sin(half) * perp
        cm_dir /= np.linalg.norm(cm_dir)
        cm_pos = cap_pos + C_C_METHYL * cm_dir
        atoms.append(AtomRec(name="CM", element="C", coord=cm_pos))

        for i, h_pos in enumerate(
                tetrahedral_h(cm_pos, cap_pos, C_H_METHYL, 3), start=1):
            atoms.append(AtomRec(name=f"HM{i}", element="H", coord=h_pos))
        return ["CM", "HM1", "HM2", "HM3", "CAP", "OAP"]

    raise ValueError(f"Unknown cap_style={cap_style!r}")


# ---------------------------------------------------------------------------
# Back-compat re-exports — older callers may still reference these
# directly off snp_builder.core. They live in metal_fragment.py now,
# but we keep the names accessible here so test imports don't break.
# ---------------------------------------------------------------------------
def write_outputs(result: BuildResult, out_dir: str | Path,
                  basename: str = "snp_frag",
                  resname: str = "SNP") -> dict[str, str]:
    """SnP-flavoured wrapper around the generic ``metal_fragment.write_outputs``.

    Default basename ``snp_frag`` and resname ``SNP`` are the SnP demo
    conventions; override either for non-SnP cases.
    """
    return _generic_write_outputs(result, out_dir, basename=basename, resname=resname)


# Alias for back-compat: older tests/imports reference _ring_plane directly.
_ring_plane = ring_plane
