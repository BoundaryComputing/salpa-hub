"""Append UFF-derived NONBON + MASS entries for metal element types.

`parmchk2` doesn't emit MASS or NONBON entries for metal elements (Ru,
Sn, Zn, Fe, Pt, …) because GAFF doesn't define them. EasyParm step 11
only *patches* existing NONBON entries via UFF data — it never *adds*
new ones, and never touches MASS. So a fresh tleap session loading
just the produced `COMPLEX.lib + COMPLEX.frcmod` always fails with:

    Error! For atom (.R<mol 1>.A<Sn1 78>) could not find vdW (or other)
    parameters for type (Sn)

This script closes that gap. For each atom type in the mol2 whose
corresponding atom-name (e.g., `Sn1` → `Sn`) is a metal element with
atomic number > 10 and *not* already in the frcmod, append:

  - a **MASS** entry (`<type>  <atomic_mass>`) so tleap knows the
    element identity for the new atom type
  - a **NONBON** entry (`<type>  R_min/2  ε`) so tleap can compute
    Lennard-Jones interactions

Both are required: tleap rejects a type with NONBON but no MASS, and
vice versa.

Reads:
    NEW_COMPLEX.mol2   — atom types + names (any element)
    metal_number.dat   — IDs of the metal atoms
    distance.dat       — bond list, used to compute coordination number
                         for picking the right UFF entry
    uff_data.txt       — UFF library (sibling of this script)
    <frcmod>           — passed as argv[2], else COMPLEX.frcmod
                         (in-place update)

Writes:
    <frcmod>           — appended NONBON lines for missing metal types
"""

import os
import re
import sys
from pathlib import Path


ATOMIC_NUMBERS = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16,
    "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23,
    "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37,
    "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51,
    "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Hf": 72,
    "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79,
    "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
}

# Standard atomic masses (amu) — IUPAC 2021 conventional values.
ATOMIC_MASSES = {
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "K": 39.098, "Ca": 40.078, "Sc": 44.956,
    "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845,
    "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723,
    "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Rb": 85.468,
    "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95,
    "Ru": 101.07, "Rh": 102.906, "Pd": 106.42, "Ag": 107.868, "Cd": 112.414,
    "In": 114.818, "Sn": 118.71, "Sb": 121.760, "Te": 127.60, "I": 126.904,
    "Cs": 132.905, "Ba": 137.327, "La": 138.905, "Hf": 178.49, "Ta": 180.948,
    "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
    "Au": 196.967, "Hg": 200.592, "Tl": 204.38, "Pb": 207.2, "Bi": 208.980,
}


def _element_from_atom_name(name: str) -> str | None:
    """Strip trailing digits, then match the longest element prefix.

    e.g. ``Sn1`` → ``Sn``; ``CA`` → ``Ca`` (atom name is uppercase but
    ``Ca`` is the element). Falls back to single-letter element.
    """
    base = re.sub(r"\d+$", "", name)
    if not base:
        return None
    two = base[:2].capitalize() if len(base) >= 2 else None
    if two and two in ATOMIC_NUMBERS:
        return two
    one = base[0].upper()
    if one in ATOMIC_NUMBERS:
        return one
    return None


def parse_mol2_atoms(mol2_path: Path) -> list[dict]:
    """Return [{id, name, type, element}, ...] for the @<TRIPOS>ATOM section."""
    atoms = []
    in_atom = False
    for line in mol2_path.read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            continue
        if in_atom and line.startswith("@<TRIPOS>"):
            break
        if not in_atom:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            aid = int(parts[0])
        except ValueError:
            continue
        name = parts[1]
        atype = parts[5].split(".")[0]  # strip MOL2 sub-typing like Sn.4
        atoms.append({
            "id": aid,
            "name": name,
            "type": atype,
            "element": _element_from_atom_name(name),
        })
    return atoms


def read_metal_ids(metal_file: Path) -> list[int]:
    if not metal_file.exists():
        return []
    return [int(s) for s in metal_file.read_text().split() if s.strip()]


def count_metal_bonds(metal_ids: set[int], distance_file: Path) -> dict[int, int]:
    bond_counts = {m: 0 for m in metal_ids}
    if not distance_file.exists():
        return bond_counts
    for line in distance_file.read_text().splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            a, b = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        if a in bond_counts:
            bond_counts[a] += 1
        if b in bond_counts:
            bond_counts[b] += 1
    return bond_counts


def read_uff_data(uff_file: Path) -> dict[str, tuple[float, float]]:
    """Map UFF atom-type → (R_min/2, ε) for AMBER NONBON.

    UFF cols: type, col1, col2, σ, ε, col5, col6
    AMBER NONBON: type  R_min/2  ε  →  (σ/2, ε)
    """
    out: dict[str, tuple[float, float]] = {}
    for line in uff_file.read_text().splitlines():
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            r_half = float(parts[3]) / 2.0
            eps = float(parts[4])
        except ValueError:
            continue
        out[parts[0]] = (r_half, eps)
    return out


def pick_uff_entry(uff: dict[str, tuple[float, float]],
                   element: str,
                   coord_number: int | None) -> tuple[float, float] | None:
    """Choose the UFF row for a metal element + coordination number.

    Strategy: prefer ``<Element><coord_number>`` (e.g. ``Sn4``, ``Fe6+2``,
    ``Ru6+2``) — UFF encodes coordination geometry in the type suffix.
    Fall back to *any* entry whose key starts with the element symbol.
    """
    candidates = [k for k in uff.keys() if k.startswith(element)]
    if not candidates:
        return None
    if coord_number is not None:
        target = str(coord_number)
        for k in candidates:
            tail = k[len(element):]
            # match e.g. "4" or "6+2" or "3+2" — first char of tail is the
            # coordination number per UFF naming convention
            if tail.startswith(target):
                return uff[k]
    return uff[candidates[0]]


def existing_nonbon_types(frcmod_path: Path) -> set[str]:
    """Return the atom types already present in the frcmod's NONBON section."""
    return _section_types(frcmod_path, "NONBON", min_parts=3)


def existing_mass_types(frcmod_path: Path) -> set[str]:
    """Return the atom types already present in the frcmod's MASS section."""
    return _section_types(frcmod_path, "MASS", min_parts=2)


def _section_types(frcmod_path: Path, section: str, min_parts: int) -> set[str]:
    types: set[str] = set()
    in_section = False
    for line in frcmod_path.read_text().splitlines():
        s = line.strip()
        if s == section:
            in_section = True
            continue
        if in_section and s in {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}:
            in_section = False
            continue
        if not in_section:
            continue
        parts = s.split()
        if len(parts) >= min_parts:
            types.add(parts[0])
    return types


def append_metal_nonbon(frcmod_path: Path, lines_to_add: list[str]) -> None:
    """Insert lines into the NONBON section, immediately after the header.

    Critical: tleap treats a blank line as a NONBON-section terminator, so
    appending at the end of the file (after the trailing blanks the
    pipeline produces) puts the entries OUTSIDE the section. Insert after
    the `NONBON` header so they're contiguous with any existing entries.
    """
    if not lines_to_add:
        return

    text = frcmod_path.read_text()
    lines = text.splitlines(keepends=True)

    for i, raw in enumerate(lines):
        if raw.strip() == "NONBON":
            new_lines = lines[: i + 1] + lines_to_add + lines[i + 1 :]
            frcmod_path.write_text("".join(new_lines))
            return

    # No NONBON section — add one at the end with our entries inline
    if not text.endswith("\n"):
        text += "\n"
    text += "\nNONBON\n" + "".join(lines_to_add)
    frcmod_path.write_text(text)


def insert_into_mass_section(frcmod_path: Path, lines_to_add: list[str]) -> None:
    """Insert lines into the MASS section (right after the `MASS` header).

    Creates the section at the top of the file if absent.
    """
    if not lines_to_add:
        return

    text = frcmod_path.read_text()
    lines = text.splitlines(keepends=True)

    for i, raw in enumerate(lines):
        if raw.strip() == "MASS":
            # Insert immediately after the MASS header.
            new_lines = lines[: i + 1] + lines_to_add + lines[i + 1 :]
            frcmod_path.write_text("".join(new_lines))
            return

    # No MASS section found — prepend one at the top, after the remark line
    # (frcmod files always start with a one-line remark).
    if not lines:
        lines = ["Remark line goes here\n"]
    head, tail = lines[:1], lines[1:]
    new_block = ["MASS\n"] + lines_to_add + ["\n"]
    frcmod_path.write_text("".join(head + new_block + tail))


def fill_metal_nonbon(work_dir: Path, frcmod_name: str = "COMPLEX.frcmod") -> dict[str, list[str]]:
    """Add MASS + NONBON entries for metal element atom types.

    Returns ``{'mass': [<lines>], 'nonbon': [<lines>]}`` — the lines added
    to each section. Both lists empty if nothing to do.
    """
    frcmod_path = work_dir / frcmod_name
    if not frcmod_path.exists():
        return {"mass": [], "nonbon": []}

    mol2_path = work_dir / "NEW_COMPLEX.mol2"
    if not mol2_path.exists():
        mol2_path = work_dir / "COMPLEX.mol2"
    if not mol2_path.exists():
        return {"mass": [], "nonbon": []}

    atoms = parse_mol2_atoms(mol2_path)
    metal_ids = set(read_metal_ids(work_dir / "metal_number.dat"))
    bond_counts = count_metal_bonds(metal_ids, work_dir / "distance.dat")

    script_dir = Path(__file__).parent
    uff_file = script_dir / "uff_data.txt"
    if not uff_file.exists():
        return {"mass": [], "nonbon": []}
    uff = read_uff_data(uff_file)

    nonbon_already = existing_nonbon_types(frcmod_path)
    mass_already = existing_mass_types(frcmod_path)

    seen_types: set[str] = set()
    nonbon_lines: list[str] = []
    mass_lines: list[str] = []
    for atom in atoms:
        atype = atom["type"]
        element = atom["element"]
        if atype in seen_types:
            continue
        if not element or ATOMIC_NUMBERS.get(element, 0) <= 10:
            continue
        # Only emit for atoms that are themselves metals (i.e., listed in
        # metal_number.dat) — the same definition step 11 uses.
        if metal_ids and atom["id"] not in metal_ids:
            continue
        coord = bond_counts.get(atom["id"])
        params = pick_uff_entry(uff, element, coord)
        if params is None:
            continue
        r_half, eps = params
        seen_types.add(atype)
        # NB format: 2 leading spaces + type (4 chars) + R_min/2 + ε
        # (matches gaff.dat / amber frcmod convention; tleap rejects
        # entries without the leading whitespace as "unknown keyword").
        if atype not in nonbon_already:
            nonbon_lines.append(f"  {atype:<8}{r_half:>10.4f} {eps:>7.4f}\n")
        if atype not in mass_already:
            mass = ATOMIC_MASSES.get(element)
            if mass is not None:
                mass_lines.append(f"{atype:<6}{mass:>8.3f}\n")

    if mass_lines:
        insert_into_mass_section(frcmod_path, mass_lines)
    if nonbon_lines:
        append_metal_nonbon(frcmod_path, nonbon_lines)
    return {"mass": mass_lines, "nonbon": nonbon_lines}


def main() -> None:
    work_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()
    frcmod_name = sys.argv[2] if len(sys.argv) > 2 else "COMPLEX.frcmod"
    added = fill_metal_nonbon(work_dir, frcmod_name)
    types_added = {l.split()[0] for l in added["nonbon"]} | {l.split()[0] for l in added["mass"]}
    if types_added:
        print(f"metal_nonbon_fill: added MASS + NONBON for {sorted(types_added)} "
              f"({len(added['mass'])} MASS, {len(added['nonbon'])} NONBON)")


if __name__ == "__main__":
    main()
