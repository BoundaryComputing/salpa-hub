"""
gen-gmx-ndx core — generate custom GROMACS index groups.

Creates two groups that distinguish original (crystallographic) residues from
homology-modeled residues:
- OriHeavy: all heavy atoms from original residues (for full restraints)
- OriBackBone: backbone atoms from original residues (for backbone restraints)

These groups are used during vacuum MD relaxation (gmx_md_relax full_4step)
with GROMACS freezegrps to keep original atoms frozen while modeled loops relax.

Runs BEFORE solvation — reads the vacuum pdb2gmx.gro structure.
"""

import csv
import glob
import os
import subprocess
from dataclasses import dataclass
from shlex import quote as _q  # every path in a shell command goes through this


# Backbone atom names for position restraints
BACKBONE_ATOMS = {"CA", "C", "N", "O", "O1P", "P", "O2P", "O5'", "C5'", "C4'", "C3'", "O3'"}


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
                    atoms.append({
                        "index": i + 1,
                        "name": name,
                        "resname": resname,
                        "chain": "",
                        "resid": resid,
                        "element": "",
                    })
                except (ValueError, IndexError):
                    continue
        else:
            seq_idx = 0
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                seq_idx += 1
                try:
                    atoms.append({
                        "index": seq_idx,
                        "name": line[12:16].strip(),
                        "resname": line[17:20].strip(),
                        "chain": line[21].strip(),
                        "resid": int(line[22:26].strip()),
                        "element": line[76:78].strip() if len(line) > 76 else "",
                    })
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


def generate_ori_ndx(
    structure_path: str,
    ndx_path: str,
    missing_csv_dir: str,
    chain_ids: list = None,
) -> NdxResult:
    """Generate custom GROMACS index groups for original residues.

    Reads the vacuum structure (pdb2gmx.gro or PDB) and identifies which
    atoms belong to original vs modeled residues using missing_residues CSVs.

    Args:
        structure_path: Path to vacuum structure (.gro or .pdb).
        ndx_path: Path to output NDX file.
        missing_csv_dir: Directory with missing_residues_chain_*.csv files.
        chain_ids: Chain IDs to process (default: auto-detect).

    Returns:
        NdxResult with group sizes and file path.
    """
    result = NdxResult()
    log_lines = []

    # Step 1: Generate base NDX
    cmd = f'echo q | gmx make_ndx -f {_q(structure_path)} -o {_q(ndx_path)}'
    rc = subprocess.run(cmd, shell=True, capture_output=True, timeout=60).returncode
    if rc != 0:
        result.log = f"gmx make_ndx failed (rc={rc})"
        return result

    # Step 2: Read atoms
    atoms = read_structure_atoms(structure_path)
    if not atoms:
        result.log = "No atoms found in structure"
        return result
    log_lines.append(f"Read {len(atoms)} atoms from {os.path.basename(structure_path)}")

    # Step 3: Collect missing residue positions
    all_missing = set()        # (chain, resid)
    all_missing_resids = set() # resid only (for GRO files without chain IDs)

    if chain_ids is None:
        chain_ids = sorted(set(a["chain"] for a in atoms if a["chain"]))

    # For GRO files (no chain IDs), scan CSV directory
    if not chain_ids:
        csv_files = glob.glob(os.path.join(missing_csv_dir, "missing_residues_chain_*.csv"))
        chain_ids = [
            os.path.basename(f).replace("missing_residues_chain_", "").replace(".csv", "")
            for f in csv_files
        ]
        log_lines.append(f"GRO file — found chain CSVs: {chain_ids}")

    for chain_id in chain_ids:
        csv_path = os.path.join(missing_csv_dir, f"missing_residues_chain_{chain_id}.csv")
        if os.path.exists(csv_path):
            chain_missing = read_missing_residues_csv(csv_path)
            all_missing.update(chain_missing)
            all_missing_resids.update(resid for _, resid in chain_missing)
            log_lines.append(f"Chain {chain_id}: {len(chain_missing)} missing residues")

    # Step 4: Build OriHeavy and OriBackBone groups
    ori_heavy_indices = []
    ori_backbone_indices = []

    for atom in atoms:
        chain = atom["chain"]
        resid = atom["resid"]

        # Skip modeled (missing) residues
        if (chain, resid) in all_missing:
            continue
        if not chain and resid in all_missing_resids:
            continue

        # Skip water/ions (only protein atoms for OriHeavy)
        if atom["resname"] in ("SOL", "NA", "CL", "HOH", "WAT"):
            continue

        # OriHeavy: all non-hydrogen atoms
        element = atom["element"].upper()
        is_heavy = False
        if element and element != "H":
            is_heavy = True
        elif not element and not atom["name"].startswith("H"):
            is_heavy = True

        if is_heavy:
            ori_heavy_indices.append(atom["index"])

        # OriBackBone: backbone atoms
        if atom["name"] in BACKBONE_ATOMS:
            ori_backbone_indices.append(atom["index"])

    log_lines.append(f"OriHeavy: {len(ori_heavy_indices)} atoms")
    log_lines.append(f"OriBackBone: {len(ori_backbone_indices)} atoms")

    # Step 5: Append groups to NDX file (strip any previous ones first)
    existing_lines = []
    with open(ndx_path) as f:
        skip = False
        for line in f:
            if line.strip() in ("[ OriHeavy ]", "[ OriBackBone ]"):
                skip = True
                continue
            if line.startswith("[") and skip:
                skip = False
            if not skip:
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
