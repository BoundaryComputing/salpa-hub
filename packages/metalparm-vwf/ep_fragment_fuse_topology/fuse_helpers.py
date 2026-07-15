"""Pure-Python helpers for ep_fragment_fuse — no bocoflow_core dep so they
can be unit-tested standalone.

After the peptide_builder split, this node only fuses: it loads a peptide
PDB written by an upstream peptide_builder, plus a fragment library + frcmod
written by ep_library_generation, and stitches them with one or more
interface bonds. Sequence-handling and peptide construction live in
peptide_builder/core.py.
"""

from __future__ import annotations

import json
import os


DEFAULT_INTERFACE_BONDS = [
    {
        "pep_resid": 6,
        "pep_atom": "CD",
        "frag_resid": 1,
        "frag_atom": "NH2",
        "pep_remove": ["OE2", "HE2"],
        "frag_remove": ["CM", "HM1", "HM2", "HM3", "CAP", "OAP"],
    }
]

FF_MAP = {
    "ff19SB": "leaprc.protein.ff19SB",
    "ff14SB": "leaprc.protein.ff14SB",
}


def parse_interface_bonds(raw: str):
    """Parse the JSON-list-of-dicts from the UI; fall back to default.

    Raises ValueError (not NodeException — caller translates) on malformed input.
    """
    if not raw or not raw.strip():
        return DEFAULT_INTERFACE_BONDS
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as e:
        raise ValueError(f"interface_bonds is not valid JSON: {e}")
    if not isinstance(value, list) or not value:
        raise ValueError("interface_bonds must be a non-empty JSON list of dicts")
    required = {"pep_resid", "pep_atom", "frag_resid", "frag_atom"}
    for i, entry in enumerate(value):
        if not isinstance(entry, dict) or not required.issubset(entry):
            raise ValueError(
                f"interface_bonds[{i}] missing required keys {required - entry.keys()}")
    return value


def read_lib_unit_name(lib_path: str) -> str:
    """Return the unit name embedded in an AMBER OFF library.

    Format: line 1 is `!!index array str`, line 2 is ` "NAME"`. We strip
    quotes and whitespace to recover NAME. Used by ep_fragment_fuse to
    issue a `copy <NAME>` even when the unit was renamed by ep_library_generation
    (typically "mol") or peptide_builder (the auto-generated unit name).
    """
    if not lib_path or not os.path.isfile(lib_path):
        raise FileNotFoundError(f"lib not found: {lib_path}")
    with open(lib_path) as f:
        for raw in f:
            stripped = raw.strip()
            if not stripped or stripped.startswith("!!"):
                continue
            # First non-header, non-empty line is the unit name in quotes
            return stripped.strip('"').strip("'").strip()
    raise ValueError(f"no unit name found in {lib_path}")


def count_pdb_residues(pdb_path: str) -> int:
    """Count distinct (chain, resseq, resname) triples in a PDB file."""
    if not pdb_path or not os.path.isfile(pdb_path):
        return 0
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                chain = line[21]
                resseq = line[22:26].strip()
                resname = line[17:20].strip()
                seen.add((chain, resseq, resname))
    return len(seen)


def build_tleap_script(*, forcefield, fragment_lib, fragment_frcmod,
                        peptide_lib, fragment_resname, peptide_resname,
                        interface_bonds, pep_unit_size,
                        linkage_frcmods=None):
    """Emit the fuse tleap script as a string.

    Pure topology fuser: loads two parameterized OFF libraries (peptide +
    fragment), combines them, removes cap atoms, creates the interface
    bond(s), saves the combined AMBER topology. No structural input
    (peptide.pdb / fragment.pdb) — atom positions come from the libs.

    `pep_unit_size` is the residue count of the peptide unit, used to
    compute the fragment residue's 1-based index inside the combined unit.

    `linkage_frcmods` is an optional list of additional frcmod paths
    that fill cross-FF parameter gaps at the interface bond. Most
    GAFF2-ff19SB amide attachments need ~10 cross-types
    (BOND CO-ns, ANGLES 2C-CO-ns / O2-CO-ns / CO-ns-hn / CO-ns-ca,
    7 zero-amplitude torsions, and the canonical X-X-C-O improper)
    that neither standalone FF defines. Bundled examples live under
    ``ep_fragment_fuse_topology/demo_data/linkages/``.
    """
    lines = []
    lines.append(f"source {FF_MAP.get(forcefield, FF_MAP['ff19SB'])}")
    lines.append("source leaprc.gaff2")
    lines.append("source leaprc.water.opc")
    lines.append("")
    lines.append(f"loadamberparams {os.path.basename(fragment_frcmod)}")
    for lk in (linkage_frcmods or []):
        lines.append(f"loadamberparams {os.path.basename(lk)}")
    lines.append(f"loadoff {os.path.basename(fragment_lib)}")
    lines.append(f"loadoff {os.path.basename(peptide_lib)}")
    lines.append("")
    lines.append(f"pep = copy {peptide_resname}")
    lines.append(f"frag = copy {fragment_resname}")
    lines.append("cpx = combine { pep frag }")
    lines.append("")
    frag_first_resid = pep_unit_size + 1
    any_removes = False
    for ifb in interface_bonds:
        pep_rid = ifb["pep_resid"]
        frag_rid_in_cpx = frag_first_resid + int(ifb.get("frag_resid", 1)) - 1
        for atom in ifb.get("pep_remove", []):
            # tleap `remove A B` removes B from container A. A must be a
            # unit (cpx), B is the residue or atom path within that unit.
            lines.append(f"remove cpx cpx.{pep_rid}.{atom}")
            any_removes = True
        for atom in ifb.get("frag_remove", []):
            lines.append(f"remove cpx cpx.{frag_rid_in_cpx}.{atom}")
            any_removes = True
    if any_removes:
        lines.append("")
    for ifb in interface_bonds:
        pep_rid = ifb["pep_resid"]
        frag_rid_in_cpx = frag_first_resid + int(ifb.get("frag_resid", 1)) - 1
        lines.append(
            f"bond cpx.{pep_rid}.{ifb['pep_atom']} "
            f"cpx.{frag_rid_in_cpx}.{ifb['frag_atom']}"
        )
    lines.append("")
    lines.append("check cpx")
    lines.append("charge cpx")
    lines.append("saveamberparm cpx complex.prmtop complex.rst7")
    lines.append("savepdb cpx complex.pdb")
    lines.append("quit")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Interface charge rebalancing
# ---------------------------------------------------------------------------

def redistribute_to_integer(charges, tol=0.4):
    """Adjust a list of per-atom charges so their sum is the nearest
    integer, spreading the remainder equally over every atom.

    Fragment fusion deletes interface cap atoms (an interface bond's
    ``pep_remove`` / ``frag_remove``) so the linkage bond can form. Each
    deleted atom carried a partial charge, so the residue it was removed
    from is left with a non-integer net charge — e.g. deleting GLU's
    ``OE2`` carboxylate oxygen (-0.8188 e) from a -1 glutamate template
    leaves the residue at -0.1812 e. This repairs that residue.

    ``tol``: max accepted distance from the nearest integer. An interface
    deletion remainder is small (the SnP case is 0.18 e); a residue more
    than ``tol`` off is a genuine parameterisation error (wrong
    protonation / charge state) or an ambiguous near-0.5 rounding — raise
    ``ValueError`` rather than silently snap it to the wrong integer.

    Returns a new list (input not mutated). A list already summing to an
    integer is returned unchanged.
    """
    charges = list(charges)
    if not charges:
        return charges
    total = sum(charges)
    target = round(total)
    delta = target - total
    if abs(delta) <= 1e-9:
        return charges
    if abs(delta) > tol:
        raise ValueError(
            f"net charge {total:+.5f} is {abs(delta):.3f} e from the "
            f"nearest integer ({target}) — beyond the {tol} e tolerance; "
            f"this is a parameterisation error (e.g. wrong protonation "
            f"state), not an interface-deletion remainder")
    per_atom = delta / len(charges)
    return [q + per_atom for q in charges]


def rebalance_residue_charges(prmtop_path, tol=0.4, report_threshold=1e-4):
    """Rewrite ``prmtop_path`` in place so every residue's net charge is
    an integer.

    tleap's ``saveamberparm`` does not repair the non-integer net charge
    left by the fuse's interface-atom deletion (see
    ``redistribute_to_integer``); for a fragment-fused metallopeptide the
    whole complex inherits the remainder (the SnP demo lands at
    -0.1812 e, all of it on the linkage GLU). This loads the prmtop,
    redistributes every residue's remainder over its own atoms so the
    complex total is an exact integer, and saves.

    ``report_threshold``: every residue is normalised, but only those
    whose original imbalance exceeded this (a genuine deletion remainder,
    not float noise from prmtop charge precision) appear in the returned
    ``adjusted`` list.

    Returns a summary dict: ``total_before``, ``total_after``, and
    ``adjusted`` — a list of ``(resname, resid, before, after)`` for each
    meaningfully-adjusted residue.
    """
    import parmed

    parm = parmed.load_file(str(prmtop_path))
    total_before = sum(a.charge for a in parm.atoms)
    adjusted = []
    changed = False
    for res in parm.residues:
        before = sum(a.charge for a in res.atoms)
        new = redistribute_to_integer([a.charge for a in res.atoms], tol=tol)
        if abs(sum(new) - before) <= 1e-12:
            continue
        for atom, q in zip(res.atoms, new):
            atom.charge = q
        changed = True
        if abs(before - round(before)) > report_threshold:
            adjusted.append((res.name, res.idx + 1, before, sum(new)))
    if changed:
        parm.save(str(prmtop_path), overwrite=True)
    total_after = sum(a.charge for a in parm.atoms)
    return {
        "total_before": total_before,
        "total_after": total_after,
        "adjusted": adjusted,
    }
