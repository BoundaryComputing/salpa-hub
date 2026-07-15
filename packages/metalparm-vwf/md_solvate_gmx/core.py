"""Pure-Python core for md_solvate_gmx.

Builds a solvated GROMACS topology (`complex.top` + `complex.gro` +
`metallopeptide_solv.itp`) from a dry GROMACS topology produced by
`ep_amber_to_gromacs`, by:

  1. Building single-molecule AMBER prmtops for each solvent species
     via `tleap` (loading the AmberTools `solvents.lib`,
     `atomic_ions.lib`, `frcmod.opc`, `frcmod.ionslm_126_opc`).
  2. Loading those single-molecule prmtops into ParmEd `Structure`s.
  3. Computing solvent counts from a box-volume target + a user-
     supplied molar ratio (defaults: MeOH 14.7 mol/nm³, H₂O 33.3
     mol/nm³).
  4. Running raw `packmol` (NOT `packmol-memgen`) to place 1 solute
     + N solvents + ions in a cubic box. Only coordinates are
     produced; the topology never re-routes through tleap, so the
     fused metallopeptide residues (e.g. SnP-modified GLU6) are
     preserved bytewise.
  5. Assembling the final ParmEd `Structure` via `Structure * N`
     (replicate) + `Structure + Structure` (concatenate). The solute
     comes from the dry GROMACS topology; solvent moleculetypes are
     appended.
  6. Exporting GROMACS `.top` + `.gro`, then splitting the
     `[ moleculetype ]` block out into a sibling `.itp` (reusing
     `ep_amber_to_gromacs.core._split_top_into_itp`).

The pure-Python functions in this module are unit-testable without
the bocoflow_core dependency.
"""
from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ─── Constants ────────────────────────────────────────────────────────

# Bulk densities (molecules per nm³) for mixed-solvent count solving.
# Sources: H₂O ≈ 55.5 mol/L = 33.4 molec/nm³ at 25 °C; MeOH ≈ 24.6 mol/L
# = 14.8 molec/nm³ at 25 °C.
DEFAULT_SOLVENT_DENSITIES_PER_NM3 = {
    "WAT": 33.4,
    "MOH": 14.8,
}

# Map water model option → leaprc.water.* file name (tleap source line).
WATER_MODEL_LEAPRC = {
    "opc":     "leaprc.water.opc",
    "opc3":    "leaprc.water.opc3",
    "tip3p":   "leaprc.water.tip3p",
    "tip4pew": "leaprc.water.tip4pew",
    "spce":    "leaprc.water.spce",
}

# Map water model → unit name to copy in tleap. Most leaprcs alias WAT,
# HOH, and the model-specific name to the same Unit; we use the model
# name to be explicit (tleap variable assignment doesn't auto-resolve
# aliases consistently across all AmberTools versions).
WATER_MODEL_UNIT = {
    "opc":     "OPC",
    "opc3":    "OP3",   # tleap names: OP3 (not OPC3)
    "tip3p":   "TP3",
    "tip4pew": "T4E",
    "spce":    "SPC",   # SPC/E uses same unit as SPC; FF differs in frcmod
}


@dataclass
class SolventCounts:
    """Resolved per-species molecule counts ready for packmol."""

    wat: int = 0
    moh: int = 0
    cation: int = 0
    anion: int = 0

    def as_dict(self, cation_code: str, anion_code: str) -> Dict[str, int]:
        out = {"WAT": self.wat}
        if self.moh:
            out["MOH"] = self.moh
        if self.cation:
            out[cation_code] = self.cation
        if self.anion:
            out[anion_code] = self.anion
        return out


# ─── tleap-based solvent unit building ────────────────────────────────


def build_solvent_unit_structures(
    work_dir: str,
    *,
    water_model: str = "opc",
    want_moh: bool = True,
    cation: str = "K+",
    anion: str = "Cl-",
    tleap_bin: str = "tleap",
    env: Optional[Dict[str, str]] = None,
    timeout: int = 60,
) -> Dict[str, "object"]:
    """Build single-molecule ParmEd Structures for each solvent species.

    Args:
      work_dir: directory to write tleap scratch files (prmtop/rst7).
      water_model: one of WATER_MODEL_LEAPRC keys (e.g. "opc").
      want_moh: include MeOH (extracted from MEOHBOX, first residue).
      cation/anion: ion unit names (must exist in atomic_ions.lib).
      tleap_bin: tleap executable (default 'tleap' on PATH).
      env: subprocess env (must have AMBERHOME set if non-default).
      timeout: seconds.

    Returns:
      dict of code → parmed.Structure. Keys: "WAT" always; "MOH" if
      want_moh; cation/anion codes if non-empty.

    Raises:
      KeyError if water_model unknown.
      RuntimeError if tleap fails.
      ImportError if ParmEd not available.
    """
    if water_model not in WATER_MODEL_LEAPRC:
        raise KeyError(
            f"unknown water_model {water_model!r}; valid: "
            f"{sorted(WATER_MODEL_LEAPRC.keys())}"
        )

    try:
        import parmed as pmd
    except ImportError as e:  # pragma: no cover — env-dependent
        raise ImportError(
            "ParmEd is required (provided by ambertools in the "
            "metalparm_vwf pixi env)."
        ) from e

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    water_leaprc = WATER_MODEL_LEAPRC[water_model]
    water_unit = WATER_MODEL_UNIT[water_model]

    # Sanitize ion names for tleap variable use (K+ → Kp, Cl- → Clm).
    def var_name(code: str) -> str:
        return code.replace("+", "p").replace("-", "m")

    script_lines: List[str] = [
        # MEOHBOX (Cieplak 2001 methanol) uses Cornell-family AMBER
        # protein atom types (CT/OH/H1/HO), NOT GAFF types. Need
        # leaprc.protein.ff19SB to load the parm10 vdW + bonded params
        # for those types. Source it before everything else so later
        # leaprcs can override specific params if needed.
        "source leaprc.protein.ff19SB",
        f"source {water_leaprc}",
        "source leaprc.gaff2",
        # leaprc.water.* already loadOff solvents.lib + atomic_ions.lib,
        # but explicitly loading frcmod.opc ensures the LJ params for the
        # OPC virtual site (MW) survive even if a custom leaprc chain
        # drops them. Harmless for non-OPC models.
        "loadAmberParams frcmod.opc",
        "loadAmberParams frcmod.ionslm_126_opc",
        "",
        f"wat = copy {water_unit}",
        "saveAmberParm wat wat.prmtop wat.rst7",
    ]

    if want_moh:
        # MEOHBOX in solvents.lib is a pre-equilibrated box of ~512 MeOH
        # residues. We save it whole and let ParmEd extract residue 1.
        # (tleap has no clean way to delete N-1 residues from a Unit.)
        script_lines += [
            "",
            "moh = copy MEOHBOX",
            "saveAmberParm moh meohbox.prmtop meohbox.rst7",
        ]

    for ion_code in (cation, anion):
        if not ion_code:
            continue
        v = var_name(ion_code)
        script_lines += [
            "",
            f"ion_{v} = copy {ion_code}",
            f"saveAmberParm ion_{v} ion_{v}.prmtop ion_{v}.rst7",
        ]

    script_lines += ["", "quit", ""]

    script_path = work / "build_solvents.leap"
    script_path.write_text("\n".join(script_lines))

    log_path = work / "build_solvents.log"
    with open(log_path, "w") as lf:
        proc = subprocess.run(
            [tleap_bin, "-f", str(script_path)],
            cwd=str(work),
            env=env or os.environ.copy(),
            stdout=lf, stderr=subprocess.STDOUT,
            timeout=timeout,
        )

    if proc.returncode != 0:
        tail = log_path.read_text().splitlines()[-40:]
        raise RuntimeError(
            f"tleap failed building solvent units (rc={proc.returncode}); "
            f"see {log_path}. Tail:\n" + "\n".join(tail)
        )

    structures: Dict[str, object] = {}

    structures["WAT"] = pmd.load_file(
        str(work / "wat.prmtop"), str(work / "wat.rst7"),
    )

    if want_moh:
        meohbox_struct = pmd.load_file(
            str(work / "meohbox.prmtop"), str(work / "meohbox.rst7"),
        )
        # Extract the FIRST residue (`:1`). ParmEd's bracket selection
        # returns a Structure containing only the matched atoms (and
        # the residue + bonded structure they're attached to).
        moh_one = meohbox_struct[":1"]
        # Move the single MeOH to the origin so packmol can place it
        # freely. parmed's coordinates is an (N,3) array we overwrite.
        coords = moh_one.coordinates
        moh_one.coordinates = coords - coords.mean(axis=0)
        structures["MOH"] = moh_one

    for ion_code in (cation, anion):
        if not ion_code:
            continue
        v = var_name(ion_code)
        structures[ion_code] = pmd.load_file(
            str(work / f"ion_{v}.prmtop"),
            str(work / f"ion_{v}.rst7"),
        )

    return structures


# ─── Box and count math ───────────────────────────────────────────────


def compute_box_dimensions_nm(
    solute_coords_angstrom,
    padding_A: float = 12.0,
) -> Tuple[float, float, float]:
    """Compute cubic-box edge length covering solute extent + 2 × padding.

    The box is cubic (single edge length used for x/y/z) — packmol works
    best with right-rectangular boxes, and a cubic box wastes a bit of
    volume but keeps the per-side margin uniform.

    Args:
      solute_coords_angstrom: numpy (N,3) array in Å.
      padding_A: half-width margin on each face (Å).

    Returns:
      (lx, ly, lz) in nanometres, all equal.
    """
    import numpy as np

    coords = np.asarray(solute_coords_angstrom, dtype=float)
    extent_A = (coords.max(axis=0) - coords.min(axis=0)).max()
    edge_A = float(extent_A + 2.0 * padding_A)
    edge_nm = edge_A * 0.1
    return (edge_nm, edge_nm, edge_nm)


def parse_solvent_ratio(ratio_str: str, codes: List[str]) -> Dict[str, float]:
    """Parse '2:1' (or '1') into a code → ratio map.

    Args:
      ratio_str: colon-separated numerical values, same length as codes.
      codes: list of solvent codes, e.g. ['MOH', 'WAT'].

    Returns:
      dict of code → float ratio.

    Raises:
      ValueError if the value counts don't match codes.
    """
    parts = [p.strip() for p in (ratio_str or "").split(":") if p.strip()]
    if len(parts) != len(codes):
        raise ValueError(
            f"solvent_ratio {ratio_str!r} has {len(parts)} parts but "
            f"solvents list has {len(codes)} entries ({codes})"
        )
    try:
        values = [float(p) for p in parts]
    except ValueError as e:
        raise ValueError(
            f"solvent_ratio {ratio_str!r} has non-numeric parts"
        ) from e
    return dict(zip(codes, values))


def compute_solvent_counts(
    box_nm: Tuple[float, float, float],
    solute_volume_nm3: float,
    ratio: Dict[str, float],
    densities_per_nm3: Optional[Dict[str, float]] = None,
) -> Dict[str, int]:
    """Solve for integer molecule counts given a mixed-solvent target.

    For a single solvent: n = density × free_volume.

    For two solvents (e.g. MOH + WAT) with molar ratio r_A : r_B and
    bulk densities ρ_A, ρ_B, volume conservation gives:

        n_A / ρ_A + n_B / ρ_B = V_free
        n_A / n_B = r_A / r_B

    Solving: n_B = V_free / ((r_A/r_B) / ρ_A + 1 / ρ_B)
             n_A = (r_A/r_B) × n_B

    Args:
      box_nm: cubic box edge tuple in nm.
      solute_volume_nm3: approximate solute excluded volume.
      ratio: dict of code → numerical ratio (e.g. {"MOH": 2, "WAT": 1}).
      densities_per_nm3: optional override; defaults to
        DEFAULT_SOLVENT_DENSITIES_PER_NM3.

    Returns:
      dict of code → int molecule count. Codes with ratio 0 are omitted.
    """
    densities = dict(DEFAULT_SOLVENT_DENSITIES_PER_NM3)
    if densities_per_nm3:
        densities.update(densities_per_nm3)

    total_volume = box_nm[0] * box_nm[1] * box_nm[2]
    free_volume = max(0.0, total_volume - solute_volume_nm3)

    active = [c for c in ratio if ratio[c] > 0]
    if not active:
        return {}

    if len(active) == 1:
        c = active[0]
        if c not in densities:
            raise KeyError(
                f"no bulk density for solvent {c!r}; provide "
                f"densities_per_nm3 explicitly"
            )
        return {c: int(round(free_volume * densities[c]))}

    if len(active) == 2:
        a, b = active
        if a not in densities or b not in densities:
            missing = [c for c in (a, b) if c not in densities]
            raise KeyError(
                f"no bulk density for {missing}; provide "
                f"densities_per_nm3 explicitly"
            )
        ra, rb = ratio[a], ratio[b]
        da, db = densities[a], densities[b]
        # n_b = V_free / ((ra/rb)/da + 1/db) -- but expressed as fraction
        # of molecules, not volume. We work with "molecules", so:
        #   V_free = n_a / da + n_b / db
        #   n_a = (ra/rb) * n_b
        #   V_free = (ra/rb) * n_b / da + n_b / db = n_b * ((ra/rb)/da + 1/db)
        n_b = free_volume / ((ra / rb) / da + 1.0 / db)
        n_a = (ra / rb) * n_b
        return {
            a: int(round(n_a)),
            b: int(round(n_b)),
        }

    raise NotImplementedError(
        f"more than 2 solvents not supported (got {active})"
    )


def compute_ion_counts(
    solute_charge: float,
    total_volume_nm3: float,
    saltcon_M: float,
    cation: str,
    anion: str,
) -> Dict[str, int]:
    """Neutralise the solute charge and add salt at saltcon_M (mol/L).

    Conversion: saltcon (mol/L) × volume (L) × N_A (per mol)
              = saltcon × total_volume_nm3 × 1e-24 (L/nm³) × 6.022e23
              = saltcon × total_volume_nm3 × 0.6022

    Args:
      solute_charge: net charge of the solute (signed).
      total_volume_nm3: full box volume.
      saltcon_M: salt concentration above neutralisation; 0 = neutralise only.
      cation/anion: code strings (e.g. "K+", "Cl-").

    Returns:
      dict of code → int (only non-zero counts are included).
    """
    # Round solute charge to nearest integer for ion accounting.
    q = int(round(solute_charge))
    if q > 0:
        n_cation_neut = 0
        n_anion_neut = q
    elif q < 0:
        n_cation_neut = -q
        n_anion_neut = 0
    else:
        n_cation_neut = n_anion_neut = 0

    n_pairs = int(round(saltcon_M * total_volume_nm3 * 0.6022))

    counts: Dict[str, int] = {}
    n_c = n_cation_neut + n_pairs
    n_a = n_anion_neut + n_pairs
    if n_c > 0:
        counts[cation] = n_c
    if n_a > 0:
        counts[anion] = n_a
    return counts


# ─── PDB writing for packmol templates ────────────────────────────────


def write_single_molecule_pdbs(
    structures: Dict[str, object],
    work_dir: str,
) -> Dict[str, str]:
    """Write one-molecule PDB templates for packmol's `structure` blocks.

    For each code, calls `Structure.save(...pdb...)` then returns the
    relative basename (packmol reads file names relative to its CWD).

    Args:
      structures: dict of code → parmed.Structure (one molecule each).
      work_dir: directory to write PDBs into.

    Returns:
      dict of code → basename (e.g. {"WAT": "wat_one.pdb", ...}).
    """
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    out: Dict[str, str] = {}
    for code, s in structures.items():
        # Sanitize code for filename: K+ → Kp, Cl- → Clm.
        safe = code.replace("+", "p").replace("-", "m")
        name = f"{safe.lower()}_one.pdb"
        path = work / name
        s.save(str(path), format="pdb", overwrite=True)
        out[code] = name
    return out


def write_solute_pdb(
    solute_structure: object,
    work_dir: str,
    *,
    basename: str = "solute.pdb",
    center: bool = True,
) -> str:
    """Write the solute PDB for packmol, optionally centred at origin.

    Args:
      solute_structure: parmed.Structure for the dry solute.
      work_dir: write directory.
      basename: output basename (default 'solute.pdb').
      center: if True, translate coords so the geometric centre is at
        the origin (matches the box layout used by `write_packmol_input`).

    Returns:
      basename (relative to work_dir).
    """
    import numpy as np

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    if center:
        coords = solute_structure.coordinates.copy()
        coords -= coords.mean(axis=0)
        solute_structure.coordinates = coords

    path = work / basename
    solute_structure.save(str(path), format="pdb", overwrite=True)
    return basename


# ─── packmol input + run ──────────────────────────────────────────────


def write_packmol_input(
    work_dir: str,
    *,
    solute_pdb: str,
    solvent_pdbs: Dict[str, str],
    solvent_counts: Dict[str, int],
    box_nm: Tuple[float, float, float],
    seed: int = -1,
    tolerance: float = 2.0,
    nloop_solvent: int = 20,
    solute_radius: float = 1.5,
    output_pdb: str = "system.pdb",
    input_name: str = "pack.inp",
) -> str:
    """Write a packmol input file packing solute + solvents in a box.

    Box convention: centred at origin, ranging from -L/2 to +L/2 on
    each axis (consistent with packmol-memgen's output and with
    `write_solute_pdb(center=True)`).

    Args:
      work_dir: write directory.
      solute_pdb: basename of the solute PDB (must already exist in work_dir).
      solvent_pdbs: dict of code → basename for each solvent PDB.
      solvent_counts: dict of code → integer count.
      box_nm: cubic box edges in nm.
      seed: packmol RNG seed; -1 = packmol picks a random seed.
      tolerance: minimum atom-atom distance (Å). 2.0 is packmol-memgen's default.
      nloop_solvent: packmol's per-structure attempt count.
      solute_radius: packmol exclusion radius around solute (Å). 1.5 = packmol-memgen.
      output_pdb: filename for packmol's packed output.
      input_name: filename for the input script.

    Returns:
      str: the full input file text (also written to disk).
    """
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    lx_A = box_nm[0] * 10.0
    ly_A = box_nm[1] * 10.0
    lz_A = box_nm[2] * 10.0
    hx, hy, hz = lx_A / 2.0, ly_A / 2.0, lz_A / 2.0
    box_xyz = (
        f"{-hx:.3f} {-hy:.3f} {-hz:.3f}  "
        f"{hx:.3f} {hy:.3f} {hz:.3f}"
    )

    lines: List[str] = [
        f"tolerance {tolerance}",
        "filetype pdb",
        f"output {output_pdb}",
    ]
    if seed >= 0:
        lines.append(f"seed {seed}")
    lines += [
        "",
        "add_amber_ter",
        "amber_ter_preserve",
        "nloop 100",
        "",
        f"structure {solute_pdb}",
        "  number 1",
        "  fixed 0. 0. 0. 0. 0. 0.",
        f"  radius {solute_radius}",
        "end structure",
        "",
    ]

    # Iterate by solvent_counts insertion order — this MUST match the
    # order the caller will use to assemble the ParmEd Structure
    # (solute first, then each code in `counts.items()` order). The
    # packmol output PDB lays atoms down block-by-block, so the block
    # order here defines the atom layout that ParmEd will copy coords
    # from. Iterating `solvent_pdbs` instead would scramble the
    # mapping whenever the two dicts have different insertion order
    # (e.g. user says "MOH:WAT" but solvent_units built WAT first).
    for code, n in solvent_counts.items():
        if n <= 0:
            continue
        basename = solvent_pdbs.get(code)
        if basename is None:
            raise KeyError(
                f"solvent_counts has {code!r} but solvent_pdbs has no "
                f"matching PDB template (keys: {list(solvent_pdbs)})"
            )
        lines += [
            f"structure {basename}",
            f"  nloop {nloop_solvent}",
            f"  number {n}",
            f"  inside box {box_xyz}",
            "end structure",
            "",
        ]

    text = "\n".join(lines)
    if not text.endswith("\n"):
        text += "\n"

    (work / input_name).write_text(text)
    return text


def run_packmol(
    work_dir: str,
    *,
    input_name: str = "pack.inp",
    output_name: str = "system.pdb",
    packmol_bin: str = "packmol",
    env: Optional[Dict[str, str]] = None,
    timeout: int = 600,
) -> str:
    """Run packmol with stdin redirected from the input file.

    packmol takes its directives on stdin (no -i flag), so we open the
    input file and pipe it in.

    Args:
      work_dir: where input + output live (subprocess cwd).
      input_name: input file basename.
      output_name: expected packmol output basename.
      packmol_bin: executable name.
      env: subprocess env.
      timeout: seconds.

    Returns:
      Absolute path to the produced output PDB.

    Raises:
      RuntimeError if packmol exits non-zero or output is missing.
    """
    work = Path(work_dir)
    log_path = work / "packmol.log"
    err_path = work / "packmol.err"

    with open(work / input_name, "r") as stdin_fh:
        with open(log_path, "w") as out_fh, open(err_path, "w") as err_fh:
            proc = subprocess.run(
                [packmol_bin],
                cwd=str(work),
                env=env or os.environ.copy(),
                stdin=stdin_fh,
                stdout=out_fh,
                stderr=err_fh,
                timeout=timeout,
            )

    out_path = work / output_name
    if proc.returncode != 0 or not out_path.is_file():
        tail = "\n".join(log_path.read_text().splitlines()[-40:])
        raise RuntimeError(
            f"packmol failed (rc={proc.returncode}). Tail of "
            f"{log_path}:\n{tail}"
        )
    return str(out_path)


# ─── ParmEd assembly ──────────────────────────────────────────────────


def assemble_solvated_structure(
    solute_structure: object,
    solvent_units: Dict[str, object],
    counts: Dict[str, int],
    packed_pdb_path: str,
    box_nm: Tuple[float, float, float],
) -> object:
    """Build the final solvated ParmEd Structure.

    Atom-order convention (matches `write_packmol_input` block order):
      1. solute atoms (from solute_structure, unchanged topology)
      2. for each code in `counts.keys()` in insertion order:
           solvent_units[code] × counts[code] (replicated topology)
      3. coordinates pulled from `packed_pdb_path` by index

    The `+` operator on `parmed.Structure` concatenates atoms and
    bonded params; `*` operator replicates a Structure.

    Args:
      solute_structure: dry solute as parmed.Structure (loaded from .top + .gro).
      solvent_units: dict of code → single-molecule Structure.
      counts: dict of code → int (must match write_packmol_input's input).
      packed_pdb_path: path to packmol's output PDB (atom order matches block order).
      box_nm: cubic box edge tuple in nm. Used to set the final box.

    Returns:
      parmed.Structure with all atoms + bonds; coords from packmol.
    """
    import parmed as pmd

    # Start with the solute as the base structure. `+` builds new
    # Structures iteratively.
    system = solute_structure

    for code, n in counts.items():
        if n <= 0:
            continue
        unit = solvent_units.get(code)
        if unit is None:
            raise KeyError(
                f"counts requested {n!r} of {code} but no Structure provided"
            )
        # `Structure * N` replicates the topology N times. Each copy
        # carries a copy of the coordinates from the original; those
        # will be overwritten from the packmol PDB below.
        replicated = unit * n
        system = system + replicated

    # Pull coordinates from packmol's packed PDB. The PDB's atoms are
    # in the same order as the assembled system because we built both
    # from the same code-ordered blocks.
    packed = pmd.load_file(packed_pdb_path)
    if len(packed.atoms) != len(system.atoms):
        raise RuntimeError(
            f"atom count mismatch: assembled system has "
            f"{len(system.atoms)} atoms, packmol produced "
            f"{len(packed.atoms)}. Block order or counts mismatched."
        )
    system.coordinates = packed.coordinates

    # Re-centre coords to the box interior (translate so min corner is
    # at origin), and set the box vectors.
    coords = system.coordinates.copy()
    coords -= coords.min(axis=0)
    system.coordinates = coords

    lx_A = box_nm[0] * 10.0
    ly_A = box_nm[1] * 10.0
    lz_A = box_nm[2] * 10.0
    system.box = [lx_A, ly_A, lz_A, 90.0, 90.0, 90.0]

    return system


# ─── GROMACS output (reuses ep_amber_to_gromacs split) ────────────────


def save_gromacs_outputs(
    structure: object,
    output_dir: str,
    *,
    prefix: str = "complex",
    itp_filename: Optional[str] = "metallopeptide_solv",
) -> Dict[str, str]:
    """Write the assembled Structure as GROMACS top + gro, optionally
    splitting the moleculetype block into a sibling .itp.

    Args:
      structure: parmed.Structure to serialise.
      output_dir: write directory.
      prefix: master file prefix (default 'complex' → complex.top / complex.gro).
      itp_filename: if set, split the [ moleculetype ] block out into
        this basename's .itp; otherwise leave the .top monolithic.

    Returns:
      dict with keys 'top', 'gro', 'itp' (itp = None when not split).
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    top_path = out_dir / f"{prefix}.top"
    gro_path = out_dir / f"{prefix}.gro"

    # Write .gro first so a save error doesn't leave a stale .top.
    structure.save(str(gro_path), format="gro", overwrite=True)
    structure.save(str(top_path), format="gromacs", overwrite=True)

    itp_path: Optional[Path] = None
    if itp_filename:
        # Reuse ep_amber_to_gromacs's split helper. Both modules ship in
        # the same package; relative import works whether we're under
        # the package namespace or invoked as scripts.
        try:
            from ..ep_amber_to_gromacs.core import (
                _split_top_into_itp, normalize_itp_basename,
            )
        except (ImportError, ValueError):
            # Robust fallback: load the sibling node's core.py by file path
            # under a UNIQUE module name. A bare `from core import ...` would
            # resolve to THIS node's already-imported `core` (which lacks the
            # helper) — the cause of the silent ITP-split / solvation failure.
            import importlib.util as _ilu
            _sib = Path(__file__).resolve().parent.parent / "ep_amber_to_gromacs" / "core.py"
            _spec = _ilu.spec_from_file_location("ep_amber_to_gromacs_core", str(_sib))
            _mod = _ilu.module_from_spec(_spec)
            _spec.loader.exec_module(_mod)
            _split_top_into_itp = _mod._split_top_into_itp
            normalize_itp_basename = _mod.normalize_itp_basename

        itp_basename = normalize_itp_basename(itp_filename)
        itp_path = out_dir / itp_basename
        _split_top_into_itp(str(top_path), str(itp_path))

    return {
        "top": str(top_path),
        "gro": str(gro_path),
        "itp": str(itp_path) if itp_path else None,
    }


# ─── Solute helpers ───────────────────────────────────────────────────


def load_solute_structure(top_path: str, gro_path: str) -> object:
    """Load a GROMACS top+gro pair into a ParmEd Structure.

    `parmed.load_file(top, xyz=gro)` consumes the GROMACS topology
    (atom types, charges, bonds) and overlays the coordinates from the
    .gro. The Structure preserves all moleculetypes verbatim — no
    re-typing, no auto-completion — which is the whole point of the
    md_solvate_gmx node.
    """
    import parmed as pmd

    return pmd.load_file(top_path, xyz=gro_path)


def estimate_solute_volume_nm3(structure: object, vdw_factor: float = 1.0) -> float:
    """Estimate the solute's excluded volume.

    Uses the molecular bounding-box volume as a rough upper bound (good
    enough for solvent-count math, which is order-of-magnitude).

    Args:
      structure: parmed.Structure for the solute (Å coordinates).
      vdw_factor: scale factor; 1.0 returns bounding-box volume. Lower
        values (e.g. 0.7) approximate the molecular volume better but
        risk over-packing.

    Returns:
      Volume in nm³.
    """
    import numpy as np

    coords = np.asarray(structure.coordinates, dtype=float)
    extents_A = coords.max(axis=0) - coords.min(axis=0)
    bbox_A3 = float(extents_A[0] * extents_A[1] * extents_A[2])
    return bbox_A3 * 1e-3 * vdw_factor  # Å³ → nm³


def compute_solute_charge(structure: object) -> float:
    """Sum atom partial charges → net solute charge."""
    return float(sum(a.charge for a in structure.atoms))
