"""Pure-Python helpers for peptide_builder — no bocoflow_core dep so they're
unit-testable.

Three responsibilities:
  1. Emit a tleap script that builds a standalone peptide from either a
     sequence string or a user-supplied PDB, and saves it as a canonical
     ``peptide.pdb`` (+ optional ``peptide.prmtop`` / ``peptide.rst7``).
  2. Compute the residue count of the built peptide so downstream nodes
     (notably ep_fragment_fuse_topology) can place fragment residues correctly.
  3. **PDB-mode preprocessing** (v1.8.0) — clean a user-supplied peptide
     PDB by selecting MODEL 1, filtering chains + residue range, dropping
     heteroatoms, applying YASARA → ff19SB atom renames, and inferring
     HIS tautomers (HD1/HE2 → HID/HIE/HIP). So an AlphaFold or
     ProteinMPNN peptide can be loaded directly without preprocessing.
"""

from __future__ import annotations

import os
import re


FF_MAP = {
    "ff19SB": "leaprc.protein.ff19SB",
    "ff14SB": "leaprc.protein.ff14SB",
}


def peptide_residue_count(sequence: str, n_term: str, c_term: str) -> int:
    """Return number of residues tleap will build, given the user sequence
    and requested terminal caps. ACE/NME are added if not already at the
    respective end.
    """
    tokens = [t for t in (sequence or "").split() if t]
    n = len(tokens)
    if n_term == "ACE" and tokens and tokens[0] != "ACE":
        n += 1
    if c_term == "NME" and tokens and tokens[-1] != "NME":
        n += 1
    return n


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


def normalize_sequence(sequence: str, n_term: str, c_term: str) -> list[str]:
    """Return the token list tleap will see, with caps inserted as needed."""
    tokens = [t for t in (sequence or "").split() if t]
    if n_term == "ACE" and tokens and tokens[0] != "ACE":
        tokens = ["ACE"] + tokens
    if c_term == "NME" and tokens and tokens[-1] != "NME":
        tokens = tokens + ["NME"]
    return tokens


PEPTIDE_FRCMOD_PLACEHOLDER = """\
# peptide_builder placeholder — pure ff19SB/ff14SB has no parameter overrides.
# Kept as an empty frcmod so the fuse boundary is symmetric with the metal
# side (which always ships a paired .lib + .frcmod).

MASS

BOND

ANGLE

DIHE

IMPROPER

NONBON
"""


def build_peptide_tleap_script(
    *,
    forcefield: str,
    peptide_sequence: str,
    n_term: str,
    c_term: str,
    save_topology: bool = False,
) -> str:
    """Emit a tleap script that builds a standalone peptide from a sequence
    and saves both `peptide.pdb` and `peptide.lib` (the parameterized AMBER
    OFF unit). Optionally also writes peptide.prmtop / peptide.rst7.
    """
    tokens = normalize_sequence(peptide_sequence, n_term, c_term)
    if not tokens:
        raise ValueError("peptide_sequence is empty")

    leaprc = FF_MAP.get(forcefield, FF_MAP["ff19SB"])
    lines = [
        f"source {leaprc}",
        "",
        f"pep = sequence {{ {' '.join(tokens)} }}",
        "check pep",
        "savepdb pep peptide.pdb",
        "saveoff pep peptide.lib",
    ]
    if save_topology:
        lines.append("saveamberparm pep peptide.prmtop peptide.rst7")
    lines.append("quit")
    return "\n".join(lines) + "\n"


def build_peptide_from_pdb_tleap_script(
    *,
    forcefield: str,
    peptide_pdb_basename: str,
    save_topology: bool = False,
) -> str:
    """Load a user-supplied peptide PDB into tleap with the chosen ff and
    save it as a parameterized AMBER OFF unit. The PDB must use standard
    residue/atom names recognized by ff19SB/ff14SB.
    """
    leaprc = FF_MAP.get(forcefield, FF_MAP["ff19SB"])
    lines = [
        f"source {leaprc}",
        "",
        f"pep = loadpdb {peptide_pdb_basename}",
        "check pep",
        "savepdb pep peptide.pdb",
        "saveoff pep peptide.lib",
    ]
    if save_topology:
        lines.append("saveamberparm pep peptide.prmtop peptide.rst7")
    lines.append("quit")
    return "\n".join(lines) + "\n"


def validate_user_pdb(pdb_path: str) -> int:
    """Validate a user-supplied peptide PDB. Returns residue count.

    Raises ValueError with a clear message on the common failure modes
    (missing file, no ATOM/HETATM records, zero residues).
    """
    if not pdb_path:
        raise ValueError("peptide_pdb path is empty")
    if not os.path.isfile(pdb_path):
        raise ValueError(f"peptide_pdb not found: {pdb_path}")
    n = count_pdb_residues(pdb_path)
    if n <= 0:
        raise ValueError(
            f"peptide_pdb has zero ATOM/HETATM residues: {pdb_path}"
        )
    return n


# ---------------------------------------------------------------------------
# PDB-mode preprocessing helpers (v1.8.0)
#
# Mirror of the YASARA-rename + HIS-tautomer-infer logic in
# snp_builder/metal_fragment.py. The two copies are intentional per the
# v1.7.0 cross-node-import policy: each node bundles what it needs so
# the marketplace install process doesn't have to handle shared libs.
# Keep these in sync if the underlying ff19SB conventions change.
# ---------------------------------------------------------------------------

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

    PDB atom-name field: 3-char names get a leading-space pad; 4-char
    names occupy the field fully.
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


def infer_his_tautomer_renames(pdb_text: str) -> dict[tuple[str, str], str]:
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

    Returns a dict keyed on ``(chain, resseq)`` for every HIS residue
    that needs renaming. ``HID``/``HIE``/``HIP``/``NHIE`` in the input
    are left alone (already explicit).

    Operates on PDB text rather than a file path so the same helper
    can be plugged into the streaming preprocessor below.
    """
    his_atoms: dict[tuple[str, str], set[str]] = {}
    for raw in pdb_text.splitlines():
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


# Standard ff19SB-recognized residue names (the 20 AAs + protonation /
# tautomer variants + caps). Used by ``drop_heteroatoms`` to distinguish
# a non-standard ATOM-line residue from a standard one. HETATM is always
# dropped when ``drop_heteroatoms=True``.
STANDARD_RESIDUES: frozenset[str] = frozenset({
    # Caps
    "ACE", "NME", "NHE",
    # 20 standard AAs
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # ff19SB protonation / tautomer variants
    "ASH", "GLH",                    # protonated ASP / GLU
    "LYN",                            # neutral LYS
    "HID", "HIE", "HIP",              # HIS tautomers / both-protonated
    "CYM", "CYX",                     # deprotonated / disulfide-bonded CYS
    "MSE",                            # selenomethionine — standard PDB
    # N-terminal / C-terminal forms (ff19SB also recognizes these prefixes)
    "NALA", "NARG", "NASN", "NASP", "NCYS", "NGLN", "NGLU", "NGLY", "NHIS",
    "NILE", "NLEU", "NLYS", "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP",
    "NTYR", "NVAL",
    "CALA", "CARG", "CASN", "CASP", "CCYS", "CGLN", "CGLU", "CGLY", "CHIS",
    "CILE", "CLEU", "CLYS", "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP",
    "CTYR", "CVAL",
})


def parse_residue_range(spec: str) -> tuple[int, int] | None:
    """Parse a residue-range spec like ``"5-30"`` into ``(5, 30)``.

    Returns ``None`` for empty / "all" / unparseable specs (caller treats
    as "keep everything"). Raises ``ValueError`` for clearly malformed
    specs (e.g. ``"abc"``, ``"5-3"`` reversed).
    """
    if not spec or not spec.strip() or spec.strip().lower() == "all":
        return None
    m = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", spec)
    if not m:
        raise ValueError(
            f"Invalid residue range {spec!r}; expected 'N-M' (e.g. '5-30')")
    lo, hi = int(m.group(1)), int(m.group(2))
    if lo > hi:
        raise ValueError(f"Invalid residue range {spec!r}: lo > hi")
    return (lo, hi)


def parse_chain_filter(spec: str) -> set[str] | None:
    """Parse a chain-filter spec like ``"A"`` or ``"A,B"`` into a set.

    Returns ``None`` for empty / "all" specs (caller treats as "keep
    everything").
    """
    if not spec or not spec.strip() or spec.strip().lower() == "all":
        return None
    return {c.strip() for c in spec.split(",") if c.strip()}


def peptide_pdb_preprocess(
    input_pdb: str,
    output_pdb: str,
    *,
    chain_filter: str = "",
    residue_range: str = "",
    drop_heteroatoms: bool = True,
    rename_yasara_atoms: bool = True,
    rename_his_tautomers: bool = True,
) -> dict[str, int]:
    """Clean a user-supplied peptide PDB so tleap's loadpdb accepts it.

    Pipeline (in order):
      1. Take only ``MODEL 1`` if the file has multiple models. NMR
         ensembles + AlphaFold-multimer outputs commonly have multiple
         frames; tleap doesn't handle them and we always want the first.
      2. Filter by chain (``chain_filter="A"`` or ``"A,B"``); empty /
         ``"all"`` keeps everything.
      3. Filter by residue sequence range (``residue_range="5-30"``);
         empty / ``"all"`` keeps everything.
      4. Drop ``HETATM`` records and any ``ATOM`` records with a
         residue name not in ``STANDARD_RESIDUES`` (when
         ``drop_heteroatoms=True``).
      5. Rename non-standard YASARA atom names (``GLU.COOH → HE2`` etc.)
         when ``rename_yasara_atoms=True``.
      6. Rewrite ``HIS`` residue names to ``HID`` / ``HIE`` / ``HIP``
         based on which protonation hydrogens are present, so tleap's
         ff19SB template loader picks the correct tautomer.

    Returns a stats dict ``{"residues": int, "atoms": int, "models_dropped":
    int, "chains_dropped": int, "residues_out_of_range": int,
    "het_dropped": int, "his_renamed": int}`` for logging / testing.
    Writes the cleaned PDB to ``output_pdb``.
    """
    if not os.path.isfile(input_pdb):
        raise ValueError(f"input_pdb not found: {input_pdb}")

    text = open(input_pdb).read()

    # 1. Take only MODEL 1
    model_blocks: list[str] = []
    if re.search(r"^MODEL\s+\d+", text, flags=re.MULTILINE):
        # Multi-model: split into MODEL ... ENDMDL blocks; keep first
        in_model = False
        cur: list[str] = []
        for raw in text.splitlines():
            if raw.startswith("MODEL"):
                in_model = True
                continue
            if raw.startswith("ENDMDL"):
                in_model = False
                if cur:
                    model_blocks.append("\n".join(cur))
                    break
                continue
            if in_model:
                cur.append(raw)
        if not model_blocks and cur:
            model_blocks.append("\n".join(cur))
        body = model_blocks[0] if model_blocks else ""
    else:
        body = text

    # Build set/range filters
    chain_set = parse_chain_filter(chain_filter)
    res_range = parse_residue_range(residue_range)

    # Pre-pass: collect HIS atoms for tautomer inference (uses the same
    # body the rest of the pipeline operates on, so range/chain filters
    # affect which HIS residues are considered).
    his_renames: dict[tuple[str, str], str] = {}

    # 2-5. Filter + rename atoms in a single pass
    out_lines: list[str] = []
    stats = {
        "residues": 0, "atoms": 0,
        "models_dropped": 0, "chains_dropped": 0,
        "residues_out_of_range": 0, "het_dropped": 0, "his_renamed": 0,
    }
    if re.search(r"^MODEL\s+\d+", text, flags=re.MULTILINE):
        # crudely count extra models for the stats dict
        n_models = len(re.findall(r"^MODEL\s+\d+", text, flags=re.MULTILINE))
        stats["models_dropped"] = max(n_models - 1, 0)

    seen_residues: set[tuple[str, str]] = set()

    for raw in body.splitlines():
        if not raw.startswith(("ATOM  ", "HETATM", "TER")):
            continue

        if raw.startswith("TER"):
            out_lines.append(raw)
            continue

        chain = raw[21]
        resseq_str = raw[22:26].strip()
        try:
            resseq_int = int(resseq_str)
        except ValueError:
            resseq_int = None
        resname = raw[17:20].strip()

        # Chain filter
        if chain_set is not None and chain not in chain_set:
            stats["chains_dropped"] += 1
            continue

        # Residue range filter
        if res_range is not None and resseq_int is not None:
            lo, hi = res_range
            if not (lo <= resseq_int <= hi):
                stats["residues_out_of_range"] += 1
                continue

        # Heteroatom / non-standard-residue filter.
        # YASARA-style PDBs use HETATM for *standard* residues too (e.g.
        # `HETATM ... GLU A 6 ...`), so the gate must be on the residue
        # name, not the record type. Drop a record only when its resname
        # isn't an ff19SB-recognized template name (waters, ligands,
        # post-translationally modified residues).
        if drop_heteroatoms and resname not in STANDARD_RESIDUES:
            stats["het_dropped"] += 1
            continue

        # Apply YASARA atom rename
        if rename_yasara_atoms:
            raw = _rename_atom_in_pdb_line(raw)

        out_lines.append(raw)
        seen_residues.add((chain, resseq_str))
        stats["atoms"] += 1

    # 6. HIS tautomer rename (operates on the filtered output)
    if rename_his_tautomers:
        joined_pre = "\n".join(out_lines) + "\n"
        his_renames = infer_his_tautomer_renames(joined_pre)
        if his_renames:
            stats["his_renamed"] = len(his_renames)
            renamed_lines: list[str] = []
            for raw in out_lines:
                if not raw.startswith(("ATOM  ", "HETATM")):
                    renamed_lines.append(raw)
                    continue
                if raw[17:20].strip() != "HIS":
                    renamed_lines.append(raw)
                    continue
                chain = raw[21]
                resseq = raw[22:26].strip()
                new_resname = his_renames.get((chain, resseq))
                if new_resname is None:
                    renamed_lines.append(raw)
                    continue
                renamed_lines.append(raw[:17] + f"{new_resname:<3s}" + raw[20:])
            out_lines = renamed_lines

    out_lines.append("END")
    stats["residues"] = len(seen_residues)

    os.makedirs(os.path.dirname(os.path.abspath(output_pdb)) or ".", exist_ok=True)
    with open(output_pdb, "w") as f:
        f.write("\n".join(out_lines) + "\n")
    return stats


def rename_lib_unit(lib_path: str, new_name: str) -> str:
    """Rewrite an AMBER OFF library so its unit name is `new_name`.

    tleap's `saveoff varname filename` writes the unit under tleap's
    *internal* name (often "mol" — the antechamber/sequence default), NOT
    the variable name. When two libraries with unit name "mol" are both
    `loadoff`'d in fuse, the second clobbers the first. Renaming peptide's
    unit from "mol" to "pep" (or whatever) fixes this.

    Returns the actual unit name found in the file (or `new_name` after
    rename). If the unit is already named `new_name`, leaves the file
    untouched.
    """
    with open(lib_path) as f:
        text = f.read()
    # First non-header line: ` "NAME"`
    cur = None
    for raw in text.splitlines():
        s = raw.strip()
        if not s or s.startswith("!!"):
            continue
        cur = s.strip('"').strip("'").strip()
        break
    if not cur:
        raise ValueError(f"could not find unit name in {lib_path}")
    if cur == new_name:
        return cur
    # Replace `"<cur>"` (index header) and `entry.<cur>.` (every section)
    # Be conservative: don't substring-match across boundaries.
    rewritten = text.replace(f'"{cur}"', f'"{new_name}"')
    rewritten = rewritten.replace(f'entry.{cur}.', f'entry.{new_name}.')
    with open(lib_path, "w") as f:
        f.write(rewritten)
    return new_name
