"""Step 05 — pass-through (no-ULS) variant.

The original 05_prepare_mol2_frcmod.py renames metal-coordinating atom
types to unique labels (na -> n1/n2/n3/n4, os -> o1/o2, ...) so that the
metal-bonded force constants don't pollute generic GAFF types. That isolation
is unused for single-metal-cofactor cases (the entire metalparm-vwf roadmap)
and creates a vocabulary divide that has to be reconciled by every
downstream consumer.

This pass-through variant skips the rename. NEW_COMPLEX.mol2 is identical
to COMPLEX_modified.mol2; the metadata files are written with identity
mappings so steps 06-13 continue to work unchanged. The resulting
.lib + .frcmod use GAFF type names throughout (`na-Sn`, `Sn-os`, ...) and
match by construction — no Hybridization_Info.dat sourcing needed.

Reads:
    COMPLEX_modified.mol2  — antechamber output (GAFF default types)
    metal_number.dat       — atom IDs of metal atoms (1-indexed)
    distance.dat           — bond list with distances; used to find which
                             atoms are bonded to metals

Writes:
    NEW_COMPLEX.mol2              — copy of COMPLEX_modified.mol2 (no rename)
    new_atomtype.dat              — one line per metal-bonded atom type (identity)
    metalloprotein_atomtype.dat   — `<type> <type>` per metal-bonded atom (identity)
    Hybridization_Info.dat        — empty addAtomTypes block (no new types)
"""

import shutil
import sys
from pathlib import Path


def read_metal_numbers(path: Path) -> list[int]:
    return [int(line.strip()) for line in path.read_text().splitlines() if line.strip()]


def read_distances(path: Path) -> dict[int, list[int]]:
    """Return {metal_atom_id: [bonded_atom_id, ...]}.

    distance.dat has rows: atom1 atom2 distance — symmetric. We collect
    every neighbour for every atom; the caller filters by metal IDs.
    """
    bonds: dict[int, list[int]] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            a, b = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        bonds.setdefault(a, []).append(b)
        bonds.setdefault(b, []).append(a)
    return bonds


def find_metal_bonded_atoms(metal_ids: list[int], bonds: dict[int, list[int]]) -> set[int]:
    bonded: set[int] = set()
    for m in metal_ids:
        for partner in bonds.get(m, []):
            if partner not in metal_ids:
                bonded.add(partner)
    return bonded


def parse_mol2_atoms(mol2_path: Path) -> dict[int, dict[str, str]]:
    """Return {atom_id: {'name': str, 'type': str}}."""
    atoms: dict[int, dict[str, str]] = {}
    in_atom_section = False
    for line in mol2_path.read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom_section = True
            continue
        if line.startswith("@<TRIPOS>") and in_atom_section:
            break
        if not in_atom_section:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            atom_id = int(parts[0])
        except ValueError:
            continue
        atoms[atom_id] = {"name": parts[1], "type": parts[5]}
    return atoms


def main() -> None:
    work_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()

    mol2_in = work_dir / "COMPLEX_modified.mol2"
    mol2_out = work_dir / "NEW_COMPLEX.mol2"
    metal_file = work_dir / "metal_number.dat"
    distance_file = work_dir / "distance.dat"

    new_atomtype_out = work_dir / "new_atomtype.dat"
    metalloprotein_out = work_dir / "metalloprotein_atomtype.dat"
    hybrid_out = work_dir / "Hybridization_Info.dat"

    if not mol2_in.exists():
        raise FileNotFoundError(f"missing {mol2_in}")

    shutil.copy2(mol2_in, mol2_out)

    atoms = parse_mol2_atoms(mol2_in)
    metal_ids = read_metal_numbers(metal_file) if metal_file.exists() else []
    bonds = read_distances(distance_file) if distance_file.exists() else {}
    bonded_ids = find_metal_bonded_atoms(metal_ids, bonds)

    bonded_types_in_order: list[str] = []
    seen_atom_ids: set[int] = set()
    for atom_id in sorted(bonded_ids):
        if atom_id in seen_atom_ids:
            continue
        seen_atom_ids.add(atom_id)
        type_name = atoms.get(atom_id, {}).get("type", "")
        if type_name:
            bonded_types_in_order.append(type_name)

    new_atomtype_out.write_text("".join(f"{t} \n" for t in bonded_types_in_order))
    metalloprotein_out.write_text("".join(f"{t} {t}\n" for t in bonded_types_in_order))
    hybrid_out.write_text("addAtomTypes {\n}\n")


if __name__ == "__main__":
    main()
