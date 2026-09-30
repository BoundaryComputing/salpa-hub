"""
gen-gmx-ndx core — generate custom GROMACS index groups.

Creates two groups that distinguish original (crystallographic) residues from
homology-modeled residues:
- OriHeavy: all heavy atoms from original residues (for full restraints)
- OriBackBone: backbone atoms from original residues (for backbone restraints)

These groups are used during vacuum MD relaxation (gmx_md_relax full_4step)
with GROMACS freezegrps to keep original atoms frozen while modeled loops relax.

Runs BEFORE solvation — reads the vacuum pdb2gmx.gro structure.

WHICH RESIDUES WERE REBUILT comes from Fix Missing Residues' record,
`rebuilt_residues.json` beside `fixed.pdb`. ProMod3 numbers every chain from 1,
while missing_residues_chain_*.csv keeps the deposited numbers, and a .gro has no
chain IDs. Up to pdbmdauto 1.2.5 the deposited numbers were matched against the
model's by number alone, so on 4Z8J the groups left out crystal residues A33-A37
and held all six rebuilt residues. The CSVs are still used, per
chain, for a PDB that keeps its deposited numbering and has no record.

gen_gmx_ndx/core.py is a copy of this file; keep the two identical.
"""

import csv
import json
import os
import subprocess
from dataclasses import dataclass

# Backbone atom names for position restraints
BACKBONE_ATOMS = {
    "CA",
    "C",
    "N",
    "O",
    "O1P",
    "P",
    "O2P",
    "O5'",
    "C5'",
    "C4'",
    "C3'",
    "O3'",
}

# Solvent and ions: never part of the protein groups, never counted as residues.
SOLVENT_AND_IONS = ("SOL", "NA", "CL", "HOH", "WAT")

# Written by Fix Missing Residues into its output folder (Merge/), beside fixed.pdb.
REBUILT_MAP_NAME = "rebuilt_residues.json"


@dataclass
class NdxResult:
    """Result of index group generation."""

    output_ndx: str = ""
    n_ori_heavy: int = 0
    n_ori_backbone: int = 0
    success: bool = False
    log: str = ""


def read_structure_atoms(path: str) -> list:
    """Parse PDB or GRO file and extract atom records.

    Uses sequential 1-based indices (GROMACS NDX convention), NOT PDB serial
    numbers — PDB serials can have duplicates across chains.

    Returns:
        List of dicts with keys: index (1-based sequential), name, resname, chain, resid, element
    """
    atoms = []
    is_gro = path.endswith(".gro")

    with open(path) as f:
        if is_gro:
            lines = f.readlines()
            # GRO format: line 0 = title, line 1 = natoms, lines 2..n+1 = atoms
            for i, line in enumerate(lines[2:]):
                if len(line.strip()) < 20:  # box line at end
                    break
                try:
                    resid = int(line[0:5].strip())
                    resname = line[5:10].strip()
                    name = line[10:15].strip()
                    atoms.append(
                        {
                            "index": i + 1,
                            "name": name,
                            "resname": resname,
                            "chain": "",
                            "resid": resid,
                            "element": "",
                        }
                    )
                except (ValueError, IndexError):
                    continue
        else:
            seq_idx = 0
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                seq_idx += 1
                try:
                    atoms.append(
                        {
                            "index": seq_idx,
                            "name": line[12:16].strip(),
                            "resname": line[17:20].strip(),
                            "chain": line[21].strip(),
                            "resid": int(line[22:26].strip()),
                            "element": line[76:78].strip() if len(line) > 76 else "",
                        }
                    )
                except (ValueError, IndexError):
                    continue

    return atoms


def read_missing_residues_csv(csv_path: str) -> set:
    """Read missing residues CSV and return set of (chain, resid) tuples."""
    missing = set()
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            missing.add((row["chain"], int(row["ssseq"])))
    return missing


def locate_rebuilt_map(
    explicit: str, case_dir: str, merge_folder: str = "Merge"
) -> str:
    """The record of rebuilt residues to use, or "" when there is none.

    An explicit path wins. Otherwise the record Fix Missing Residues writes into
    the case folder's Merge/ subfolder, which is where the pipeline puts it.
    """
    if explicit:
        return explicit if os.path.isfile(explicit) else ""
    candidate = (
        os.path.join(case_dir, merge_folder, REBUILT_MAP_NAME) if case_dir else ""
    )
    return candidate if candidate and os.path.isfile(candidate) else ""


def load_rebuilt_map(path: str) -> list:
    """Read Fix Missing Residues' record: one entry per chain, in fixed.pdb's order.

    Each entry is {"chain": id, "residues": [model numbers in order],
    "rebuilt": [the numbers ProMod3 built]}. Raises ValueError on a malformed record.
    """
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    chains = data.get("chains")
    if not isinstance(chains, list) or not chains:
        raise ValueError(f"{os.path.basename(path)} lists no chains")
    for c in chains:
        residues, rebuilt = c.get("residues"), c.get("rebuilt")
        if not isinstance(residues, list) or not isinstance(rebuilt, list):
            raise ValueError(
                f"{os.path.basename(path)}: chain {c.get('chain')!r} is incomplete"
            )
        if not set(rebuilt) <= set(residues):
            raise ValueError(
                f"{os.path.basename(path)}: chain {c.get('chain')!r} lists rebuilt residues "
                f"that are not in its model"
            )
    return chains


def residue_runs(atoms: list) -> list:
    """Group the protein atoms into residues, in file order: [((chain, resid), [atoms])].

    A new residue starts wherever the chain or the residue number changes. Solvent
    and ions are left out.
    """
    runs = []
    for atom in atoms:
        if atom["resname"] in SOLVENT_AND_IONS:
            continue
        key = (atom["chain"], atom["resid"])
        if runs and runs[-1][0] == key:
            runs[-1][1].append(atom)
        else:
            runs.append((key, [atom]))
    return runs


def _span(numbers: list) -> str:
    return f"{numbers[0]}-{numbers[-1]} ({len(numbers)})" if numbers else "none"


def rebuilt_atom_indices(atoms: list, chains: list) -> set:
    """Indices of the atoms that belong to rebuilt residues, located chain by chain.

    With chain IDs (a PDB), each chain is found by its ID. A .gro has none, so its
    residues are taken in order and split by the record's residue counts: pdb2gmx
    writes the chains in fixed.pdb's order and keeps ProMod3's numbers. Each chain's
    residue numbers must equal the record's exactly; anything else (renumbered, a
    residue missing or added, chains reordered) raises ValueError rather than guess.
    """
    runs = residue_runs(atoms)
    rebuilt = set()

    if runs and all(key[0] for key, _ in runs):
        by_chain = {}
        for (chain, resid), members in runs:
            by_chain.setdefault(chain, []).append((resid, members))
        recorded = {c["chain"] for c in chains}
        extra = sorted(set(by_chain) - recorded)
        if extra:
            raise ValueError(
                f"chains {', '.join(extra)} are in the structure but not in the record "
                f"Fix Missing Residues wrote"
            )
        segments = [(c, by_chain.get(c["chain"], [])) for c in chains]
    else:
        segments, position = [], 0
        for c in chains:
            n = len(c["residues"])
            segment = [
                (key[1], members) for key, members in runs[position : position + n]
            ]
            segments.append((c, segment))
            position += n
        if position != len(runs):
            raise ValueError(
                f"the structure has {len(runs)} residues; Fix Missing Residues built {position}"
            )

    for c, segment in segments:
        numbers = [resid for resid, _ in segment]
        if numbers != c["residues"]:
            raise ValueError(
                f"chain {c['chain']}: the structure's residues {_span(numbers)} are not the "
                f"{_span(c['residues'])} Fix Missing Residues built, so the rebuilt ones "
                f"cannot be located"
            )
        rebuilt_numbers = set(c["rebuilt"])
        for resid, members in segment:
            if resid in rebuilt_numbers:
                rebuilt.update(a["index"] for a in members)
    return rebuilt


def csv_missing_atom_indices(
    atoms: list, missing_csv_dir: str, chain_ids: list = None
) -> set:
    """Indices of atoms in residues the CSVs list as missing, matched per chain.

    Only right for a structure that keeps the DEPOSITED numbering and its chain IDs.
    A .gro has no chain IDs and carries ProMod3's numbers, so it raises ValueError.
    """
    if not any(a["chain"] for a in atoms):
        raise ValueError(
            f"no {REBUILT_MAP_NAME} from Fix Missing Residues. A .gro has no chain IDs and "
            f"carries ProMod3's residue numbers, not the deposited ones in "
            f"missing_residues_chain_*.csv, so the rebuilt residues cannot be told apart. "
            f"Re-run Fix Missing Residues (pdbmdauto 1.2.6 or later), or give the record's path"
        )
    if not chain_ids:
        chain_ids = sorted({a["chain"] for a in atoms if a["chain"]})
    missing = set()
    for chain_id in chain_ids:
        csv_path = os.path.join(
            missing_csv_dir, f"missing_residues_chain_{chain_id}.csv"
        )
        if os.path.exists(csv_path):
            missing.update(read_missing_residues_csv(csv_path))
    return {a["index"] for a in atoms if (a["chain"], a["resid"]) in missing}


def ori_group_indices(atoms: list, skip: set) -> tuple:
    """OriHeavy and OriBackBone: protein atoms outside `skip`."""
    heavy, backbone = [], []
    for atom in atoms:
        if atom["index"] in skip or atom["resname"] in SOLVENT_AND_IONS:
            continue
        element = atom["element"].upper()
        if (element and element != "H") or (
            not element and not atom["name"].startswith("H")
        ):
            heavy.append(atom["index"])
        if atom["name"] in BACKBONE_ATOMS:
            backbone.append(atom["index"])
    return heavy, backbone


def generate_ori_ndx(
    structure_path: str,
    ndx_path: str,
    missing_csv_dir: str,
    chain_ids: list = None,
    rebuilt_map_path: str = "",
) -> NdxResult:
    """Generate custom GROMACS index groups for original residues.

    Reads the vacuum structure (pdb2gmx.gro or PDB) and identifies which atoms
    belong to original vs rebuilt residues, from Fix Missing Residues' record when
    one is given, else from the missing_residues CSVs (a PDB in the deposited
    numbering only). Nothing is written unless the rebuilt residues are located.

    Args:
        structure_path: Path to vacuum structure (.gro or .pdb).
        ndx_path: Path to output NDX file.
        missing_csv_dir: Directory with missing_residues_chain_*.csv files.
        chain_ids: Chain IDs to process (default: auto-detect).
        rebuilt_map_path: Fix Missing Residues' rebuilt_residues.json ("" for none).

    Returns:
        NdxResult with group sizes and file path.
    """
    result = NdxResult()
    log_lines = []

    # Step 1: Read atoms
    atoms = read_structure_atoms(structure_path)
    if not atoms:
        result.log = "No atoms found in structure"
        return result
    log_lines.append(f"Read {len(atoms)} atoms from {os.path.basename(structure_path)}")

    # Step 2: Locate the rebuilt residues, before anything is written
    try:
        if rebuilt_map_path:
            chains = load_rebuilt_map(rebuilt_map_path)
            skip = rebuilt_atom_indices(atoms, chains)
            for c in chains:
                log_lines.append(
                    f"Chain {c['chain']}: rebuilt residues {c['rebuilt'] or 'none'}"
                )
        else:
            skip = csv_missing_atom_indices(atoms, missing_csv_dir, chain_ids)
            log_lines.append(
                f"No {REBUILT_MAP_NAME}: missing residues from the CSVs, per chain"
            )
    except (OSError, ValueError) as e:
        result.log = "\n".join(log_lines + [str(e)])
        return result

    # Step 3: Build OriHeavy and OriBackBone groups
    ori_heavy_indices, ori_backbone_indices = ori_group_indices(atoms, skip)
    log_lines.append(f"OriHeavy: {len(ori_heavy_indices)} atoms")
    log_lines.append(f"OriBackBone: {len(ori_backbone_indices)} atoms")

    # Step 4: Generate base NDX
    # "q" quits make_ndx's interactive group editor, accepting the defaults.
    # Fed to stdin rather than piped through a shell, so a path containing a
    # space (every packaged macOS install) stays one argument.
    cmd = ["gmx", "make_ndx", "-f", structure_path, "-o", ndx_path]
    rc = subprocess.run(
        cmd, input="q\n", capture_output=True, text=True, timeout=60
    ).returncode
    if rc != 0:
        result.log = "\n".join(log_lines + [f"gmx make_ndx failed (rc={rc})"])
        return result

    # Step 5: Append groups to NDX file (strip any previous ones first)
    existing_lines = []
    with open(ndx_path) as f:
        skip_group = False
        for line in f:
            if line.strip() in ("[ OriHeavy ]", "[ OriBackBone ]"):
                skip_group = True
                continue
            if line.startswith("[") and skip_group:
                skip_group = False
            if not skip_group:
                existing_lines.append(line)

    with open(ndx_path, "w") as f:
        f.writelines(existing_lines)
        f.write("\n[ OriHeavy ]\n")
        for i, idx in enumerate(ori_heavy_indices):
            f.write(f"{idx:>8d}")
            if (i + 1) % 15 == 0:
                f.write("\n")
        f.write("\n\n[ OriBackBone ]\n")
        for i, idx in enumerate(ori_backbone_indices):
            f.write(f"{idx:>8d}")
            if (i + 1) % 15 == 0:
                f.write("\n")
        f.write("\n")

    result.output_ndx = ndx_path
    result.n_ori_heavy = len(ori_heavy_indices)
    result.n_ori_backbone = len(ori_backbone_indices)
    result.success = True
    result.log = "\n".join(log_lines)

    return result
