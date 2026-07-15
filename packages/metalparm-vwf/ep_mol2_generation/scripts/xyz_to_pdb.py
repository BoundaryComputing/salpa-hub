"""Convert an XYZ file to PDB.

Two modes:

- **Legacy** (no ``--template``): atom names are generated as
  ``<element><1-based-counter>`` (e.g. the 1st nitrogen becomes ``N1``,
  the 47th carbon becomes ``C47``). This is the original easyPARM
  behavior and matches the file layout that downstream scripts
  (``atomtype_helper.py``, ``atomtype_detector.py``) expect.

- **Template overlay** (``--template <pdb>``): atom names are read from
  the template PDB and overlaid onto the XYZ coordinates atom-by-atom
  (in order). Used to preserve meaningful labels like ``NH2``, ``CAP``,
  ``OAP``, ``CM``, ``HM1-3`` from snp_builder through QM optimization →
  antechamber → final lib. On any validation failure (atom-count
  mismatch or per-position element mismatch) the script emits a
  warning to stderr and falls back to the legacy mode — the run does
  not fail.

  **Constraint**: antechamber re-normalizes atom names whose leading
  character isn't a valid element symbol (e.g., ``Z..`` becomes ``C..``
  for a carbon atom). Names that start with the correct element symbol
  — e.g. ``NH2`` (nitrogen), ``CAP`` (carbon), ``HM1`` (hydrogen) —
  are preserved verbatim by antechamber, which is the real-world case
  for snp_builder and similar fragment builders.

Original easyPARM xyz_to_pdb.py header preserved below.
"""
###################################################################################################################
#   Automated Force Fields for Metals     /$$$$$$$   /$$$$$$  /$$$$$$$  /$$      /$$                              #
#                                        | $$__  $$ /$$__  $$| $$__  $$| $$$    /$$$                              #
# Original: Abdelazim M. A. Abdelgawwad — Universitat de València — LGPL-2.1
###################################################################################################################
import argparse
import re
import sys
from collections import defaultdict


def parse_xyz(xyz_path):
    """Yield ``(element, x, y, z)`` tuples from an XYZ file.

    Skips the first two header lines per the XMOL convention.
    """
    with open(xyz_path) as fh:
        next(fh, None)
        next(fh, None)
        for line in fh:
            parts = line.split()
            if len(parts) < 4:
                continue
            element = parts[0]
            yield element, float(parts[1]), float(parts[2]), float(parts[3])


def parse_template_pdb(pdb_path):
    """Return list of ``(atom_name, element)`` from ATOM/HETATM rows.

    Element is taken from cols 77-78 if present (PDB v3 convention),
    otherwise inferred from the leading alphabetic prefix of the atom
    name (cols 13-16). All values normalized to capitalized element
    symbols (e.g. ``"N"``, ``"Sn"``).
    """
    rows = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            element_field = line[76:78].strip() if len(line) > 76 else ""
            if element_field:
                element = element_field.capitalize()
            else:
                m = re.match(r"([A-Za-z]+)", atom_name)
                element = m.group(1).capitalize() if m else ""
            rows.append((atom_name, element))
    return rows


def _format_atom_name(name):
    """Format the atom-name field (cols 13-16) per PDB convention.

    - 4-char names occupy the full field.
    - ≤3-char names are left-padded with one space and right-padded
      to fill the 4-column field. This is what most pdb writers do
      and keeps gawk-style downstream parsers happy.
    """
    if len(name) >= 4:
        return name[:4]
    return f" {name:<3}"


def _emit_pdb_row(fh, atom_number, atom_name, x, y, z, element=""):
    """Write one ATOM record. PDB v3 columns:

    ``ATOM  ``    cols 1-6     record name (with trailing 2 spaces)
    ``%5d``       cols 7-11    atom serial number
    ``  ``        col 12       blank (alternative location)
    name field    cols 13-16   atom name (4 chars; see _format_atom_name)
    ``  ``        col 17       blank (alternate location)
    ``mol``       cols 18-20   residue name
    `` ``         col 21       blank
    ``A``/``  ``  col 22       chain ID (we use blank for legacy parity)
    `` ``         col 23-26    residue sequence number (4 cols)
    ``    ``      cols 27-30   blank + insertion code
    coords        cols 31-54   x, y, z each %8.3f
    occupancy     cols 55-60   %6.2f
    tempfactor    cols 61-66   %6.2f
    blank         cols 67-76
    element       cols 77-78   right-justified element symbol
    """
    name_field = _format_atom_name(atom_name)
    elt = (element or "").capitalize()[:2]
    fh.write(
        f"ATOM  {atom_number:5d} {name_field} mol     1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elt:>2}\n"
    )


def _legacy_xyz_to_pdb(xyz_path, pdb_path):
    """Original element+counter naming."""
    counts = defaultdict(int)
    atom_number = 1
    with open(pdb_path, "w") as fh:
        for element, x, y, z in parse_xyz(xyz_path):
            counts[element] += 1
            name = f"{element}{counts[element]}"
            _emit_pdb_row(fh, atom_number, name, x, y, z, element)
            atom_number += 1
        fh.write("END\n")


def _template_xyz_to_pdb(xyz_path, pdb_path, template_path):
    """Overlay template names onto XYZ coords. Returns True on success;
    False on validation failure (caller falls back to legacy)."""
    xyz_atoms = list(parse_xyz(xyz_path))
    template = parse_template_pdb(template_path)

    if len(xyz_atoms) != len(template):
        print(
            f"WARNING: template/XYZ atom count mismatch "
            f"(XYZ={len(xyz_atoms)}, template={len(template)}) — "
            f"falling back to element+counter naming",
            file=sys.stderr,
        )
        return False

    for i, ((xyz_elt, *_), (_, tmpl_elt)) in enumerate(zip(xyz_atoms, template), start=1):
        if xyz_elt.capitalize() != tmpl_elt.capitalize():
            print(
                f"WARNING: element mismatch at position {i} "
                f"(XYZ={xyz_elt}, template={tmpl_elt}) — "
                f"falling back to element+counter naming",
                file=sys.stderr,
            )
            return False

    with open(pdb_path, "w") as fh:
        for atom_number, ((element, x, y, z), (name, _)) in enumerate(
            zip(xyz_atoms, template), start=1
        ):
            _emit_pdb_row(fh, atom_number, name, x, y, z, element)
        fh.write("END\n")
    print(
        f"Reconstructed atom names from template: {len(template)} atoms preserved"
    )
    return True


def xyz_to_pdb(xyz_file, pdb_file, template=None):
    """Public entry point. Used as a library function from tests."""
    if template:
        if _template_xyz_to_pdb(xyz_file, pdb_file, template):
            return
    _legacy_xyz_to_pdb(xyz_file, pdb_file)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Convert XYZ to PDB.")
    p.add_argument("xyz_file")
    p.add_argument("pdb_file")
    p.add_argument(
        "--template",
        default=None,
        help="Optional PDB whose atom names overlay the XYZ coords. "
        "When omitted, names are generated as <element><counter>.",
    )
    args = p.parse_args()
    xyz_to_pdb(args.xyz_file, args.pdb_file, args.template)
