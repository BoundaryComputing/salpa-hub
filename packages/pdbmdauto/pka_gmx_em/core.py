"""
pka_gmx_em core — protonation bridge + GROMACS topology + energy minimization.

Replaces legacy fix_pka_gmx_em.py (Docker + gromacs_py) with direct CLI calls.

Pipeline:
1. PDB2PQR + PROPKA  — predict protonation states at target pH
2. Protonation bridge — parse PQR residue names, patch original PDB
3. gmx pdb2gmx -ignh  — generate GROMACS topology from patched PDB
4. gmx editconf        — create simulation box (triclinic)
5. Two-step EM         — steepest descent (no constraints → h-bonds)

Dependencies: pdb2pqr CLI, gromacs CLI, Python stdlib only.
No gromacs_py, no MDAnalysis, no Docker, no plotting libraries.

Protonation bridge rationale (replicating gromacs_py.prepare_top):
  PDB2PQR with --ffout AMBER renames titratable residues to GROMACS-compatible
  AMBER names (HID/HIE/HIP, ASH, GLH, CYX). These are directly recognised by
  GROMACS amber99sb-ildn without any translation.  When the patched PDB is fed
  to pdb2gmx with -ignh, pdb2gmx reads the residue names, skips its own
  interactive histidine analysis, and assigns correct force-field parameters.
"""

import os
import re
import subprocess
from dataclasses import dataclass, field

# ── Constants ──────────────────────────────────────────────────────────────

# Protonation-state residue names that PDB2PQR may write (--ffout AMBER)
_AMBER_PROTONATION = {"HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM"}

# CHARMM equivalents (for --ffout CHARMM)
_CHARMM_PROTONATION = {"HSD", "HSE", "HSP", "ASPP", "GLUP", "CYX", "CYM"}


# ── Result dataclass ──────────────────────────────────────────────────────

@dataclass
class PkaGmxEmResult:
    """Result of the protonation + topology + EM pipeline."""

    em_gro: str = ""
    em_top: str = ""
    pdb2gmx_gro: str = ""
    patched_pdb: str = ""
    pqr_file: str = ""
    em_max_force: float = 0.0
    protonation_changes: dict = field(default_factory=dict)
    success: bool = False
    log: str = ""


# ── Helpers ───────────────────────────────────────────────────────────────

def _run(argv, cwd=None, timeout=300, stdin_text=None):
    """Run a command as an argv list. Returns (returncode, combined output).

    NO SHELL. Arguments are passed to execve as-is, so a space, quote or `$` in
    a path is simply part of the argument. Building one string and letting a
    shell re-split it is what broke the pipeline on packaged macOS, where every
    node lives under `~/Library/Application Support/...` (bocoflow#104).

    `stdin_text` replaces the `echo q | gmx make_ndx` idiom: the pipe only ever
    answered an interactive prompt, and stdin does that without a shell.
    """
    r = subprocess.run(
        argv, input=stdin_text, capture_output=True, text=True,
        cwd=cwd, timeout=timeout,
    )
    return r.returncode, (r.stdout or "") + "\n" + (r.stderr or "")


def _extract_max_force(log_path):
    """Extract final 'Maximum force' value from a GROMACS log file."""
    max_force = 0.0
    try:
        with open(log_path) as fh:
            for line in fh:
                m = re.search(r"Maximum force\s*=\s*([\d.eE+\-]+)", line)
                if m:
                    max_force = float(m.group(1))
    except (OSError, ValueError):
        pass
    return max_force


def _purge_backup_files(directory):
    """Remove GROMACS backup files (#*) and step*.pdb temporaries."""
    for fname in os.listdir(directory):
        if re.match(r"^#.*|^step.*\.pdb$", fname):
            try:
                os.remove(os.path.join(directory, fname))
            except OSError:
                pass


def _write_em_mdp(path, constraints="none", nsteps=1000):
    """Write a minimal steepest-descent energy-minimisation MDP file."""
    with open(path, "w") as fh:
        fh.write(
            f"; Energy minimisation — constraints={constraints}\n"
            f"integrator  = steep\n"
            f"nsteps      = {nsteps}\n"
            f"emtol       = 1000.0\n"
            f"emstep      = 0.01\n"
            f"nstxout     = 0\n"
            f"nstvout     = 0\n"
            f"nstfout     = 0\n"
            f"nstlog      = 100\n"
            f"nstenergy   = 100\n"
            f"nstlist     = 10\n"
            f"cutoff-scheme = Verlet\n"
            f"coulombtype = PME\n"
            f"rcoulomb    = 1.0\n"
            f"rvdw        = 1.0\n"
            f"pbc         = xyz\n"
            f"constraints = {constraints}\n"
        )


# ── Protonation bridge ───────────────────────────────────────────────────
#
# This is the key piece that gromacs_py.prepare_top() does internally:
#   1. Run PDB2PQR  →  PQR file with force-field-named residues
#   2. Parse PQR    →  extract (chain, resSeq) → new residue name
#   3. Patch PDB    →  apply renamed residues to the original PDB
#   4. Feed to pdb2gmx -ignh  →  topology with correct protonation
#

def _parse_pqr_residues(pqr_path, ff_type="AMBER"):
    """Parse PQR file and return residues that PDB2PQR renamed for protonation.

    Returns:
        dict of (chain_id, res_seq_str) → new_residue_name
        Only entries where the name is a known protonation-state name.
    """
    known = _AMBER_PROTONATION if ff_type == "AMBER" else _CHARMM_PROTONATION
    residues = {}

    with open(pqr_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            # PQR follows PDB column layout: resName=17-19, chainID=21, resSeq=22-25
            try:
                res_name = line[17:21].strip()
                chain_id = line[21:22].strip()
                res_seq = line[22:26].strip()
            except IndexError:
                continue

            if res_name in known:
                key = (chain_id, res_seq)
                if key not in residues:
                    residues[key] = res_name

    return residues


def _detect_unrenamed_his(pqr_path):
    """Detect HIS residues that PDB2PQR left generic and assign from H atoms.

    PDB2PQR should rename all HIS to HID/HIE/HIP with --ffout AMBER, but
    as a safety net we check for any remaining "HIS" and assign based on
    which nitrogen hydrogens are present (HD1 → delta, HE2 → epsilon).

    Returns:
        dict of (chain_id, res_seq_str) → assigned AMBER name (HID/HIE/HIP)
    """
    his_atoms = {}  # (chain, resSeq) → set of atom names

    with open(pqr_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            res_name = line[17:21].strip()
            if res_name != "HIS":
                continue
            chain_id = line[21:22].strip()
            res_seq = line[22:26].strip()
            atom_name = line[12:16].strip()

            key = (chain_id, res_seq)
            his_atoms.setdefault(key, set()).add(atom_name)

    assignments = {}
    for key, atoms in his_atoms.items():
        has_hd1 = "HD1" in atoms
        has_he2 = "HE2" in atoms
        if has_hd1 and has_he2:
            assignments[key] = "HIP"   # doubly protonated (+1 charge)
        elif has_hd1:
            assignments[key] = "HID"   # delta-protonated (neutral)
        else:
            assignments[key] = "HIE"   # epsilon-protonated (neutral, most common)

    return assignments


def patch_pdb_with_protonation(original_pdb, pqr_path, output_pdb, ff_type="AMBER"):
    """Apply PDB2PQR protonation-state residue names to original PDB.

    This replicates the residue-renaming bridge that gromacs_py does between
    PDB2PQR and pdb2gmx.  The patched PDB can then be fed to pdb2gmx -ignh.

    Steps:
        1. Parse PQR for renamed protonation residues (HID/HIE/HIP, ASH, GLH, CYX)
        2. Detect any remaining generic "HIS" and assign from hydrogen atoms
        3. Replace residue names in original PDB at matching (chain, resSeq)

    Args:
        original_pdb: Path to the original PDB file.
        pqr_path:     Path to PDB2PQR output PQR file.
        output_pdb:   Path for the patched PDB output.
        ff_type:      "AMBER" or "CHARMM" — determines expected names.

    Returns:
        dict of changes: (chain, resSeq) → (old_name, new_name)
    """
    # Step 1: get renamed residues from PQR
    protonation_map = _parse_pqr_residues(pqr_path, ff_type)

    # Step 2: handle any remaining generic HIS
    his_fallback = _detect_unrenamed_his(pqr_path)
    for key, name in his_fallback.items():
        if key not in protonation_map:
            protonation_map[key] = name

    # Step 3: patch PDB
    changes = {}
    patched_lines = []

    with open(original_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                old_name = line[17:21].strip()
                chain_id = line[21:22].strip()
                res_seq = line[22:26].strip()
                key = (chain_id, res_seq)

                if key in protonation_map and old_name != protonation_map[key]:
                    new_name = protonation_map[key]

                    # PDB columns 17-19 (0-indexed): residue name field
                    # 3-char names (AMBER: HID, HIE, HIP, ASH, GLH, CYX)
                    # 4-char names (CHARMM: ASPP, GLUP) extend into column 20
                    if len(new_name) <= 3:
                        res_field = f" {new_name:<3s}"    # " HID"
                    else:
                        res_field = f"{new_name:<4s}"     # "ASPP"
                    line = line[:17] + res_field + line[21:]

                    if key not in changes:
                        changes[key] = (old_name, new_name)

            patched_lines.append(line)

    with open(output_pdb, "w") as fh:
        fh.writelines(patched_lines)

    return changes


# ── Main pipeline ─────────────────────────────────────────────────────────

def process_pka_gmx_em(
    input_pdb,
    output_dir,
    case_name,
    force_field="amber99sb",
    water_model="tip3p",
    box_distance=2.0,
    em_steps=1000,
    ph=7.0,
    run_pdb2pqr=True,
):
    """Run protonation prediction + GROMACS topology + energy minimisation.

    Args:
        input_pdb:     Input PDB (from fix_residues or merge step).
        output_dir:    Working directory for all outputs.
        case_name:     Case identifier (used in file names).
        force_field:   GROMACS force field name (amber99sb, amber99sb-ildn, …).
        water_model:   Water model for pdb2gmx (tip3p, spc, spce, tip4p).
        box_distance:  Padding distance around solute (nm) for editconf.
        em_steps:      Maximum steps per EM stage.
        ph:            Target pH for PROPKA protonation prediction.
        run_pdb2pqr:   Whether to run PDB2PQR (if False, feeds raw PDB to pdb2gmx).

    Returns:
        PkaGmxEmResult with output file paths and status.
    """
    result = PkaGmxEmResult()
    log_lines = []
    os.makedirs(output_dir, exist_ok=True)

    # Determine PDB2PQR force-field flag
    if force_field.startswith("amber"):
        pdb2pqr_ff = "AMBER"
    elif force_field.startswith("charmm"):
        pdb2pqr_ff = "CHARMM"
    elif force_field.startswith("opls"):
        pdb2pqr_ff = "AMBER"   # OPLS ≈ AMBER naming for titratable residues
    else:
        pdb2pqr_ff = "AMBER"

    # Start with original PDB; may be replaced by patched version
    gmx_input_pdb = input_pdb

    # ── Step 1: PDB2PQR + PROPKA ──────────────────────────────────────────
    if run_pdb2pqr:
        pqr_file = os.path.join(output_dir, "propka.pqr")

        cmd = [
            "pdb2pqr", "--ff", pdb2pqr_ff, "--ffout", pdb2pqr_ff,
            "--keep-chain", "--titration-state-method=propka",
            f"--with-ph={ph:.2f}", "--log-level=INFO", "--include-header",
            input_pdb, pqr_file,
        ]

        rc, out = _run(cmd, cwd=output_dir, timeout=120)

        if rc == 0 and os.path.exists(pqr_file):
            result.pqr_file = pqr_file
            log_lines.append("pdb2pqr: OK")

            # ── Step 2: Protonation bridge ────────────────────────────────
            patched_pdb = os.path.join(output_dir, "protonated.pdb")
            changes = patch_pdb_with_protonation(
                input_pdb, pqr_file, patched_pdb, ff_type=pdb2pqr_ff,
            )
            result.patched_pdb = patched_pdb
            result.protonation_changes = changes
            gmx_input_pdb = patched_pdb

            if changes:
                log_lines.append(
                    f"protonation bridge: {len(changes)} residue(s) renamed"
                )
                for (ch, seq), (old, new) in sorted(changes.items()):
                    log_lines.append(
                        f"  chain {ch or '-'} res {seq}: {old} -> {new}"
                    )
            else:
                log_lines.append(
                    "protonation bridge: no changes (standard states at this pH)"
                )
        else:
            log_lines.append(f"pdb2pqr: FAILED (rc={rc}), falling back to raw PDB")
            log_lines.append(f"  {out[:400]}")
            # Not fatal — pdb2gmx will use its own hydrogen-bond analysis

    # ── Step 3: gmx pdb2gmx ──────────────────────────────────────────────
    pdb2gmx_gro = os.path.join(output_dir, "pdb2gmx.gro")
    pdb2gmx_top = os.path.join(output_dir, "pdb2gmx.top")

    cmd = ["gmx", "pdb2gmx", "-f", gmx_input_pdb, "-o", pdb2gmx_gro,
           "-p", pdb2gmx_top, "-ff", force_field, "-water", water_model,
           "-ignh"]

    rc, out = _run(cmd, cwd=output_dir)
    log_lines.append(f"pdb2gmx: rc={rc}")
    if rc != 0:
        result.log = "\n".join(log_lines) + "\n" + out
        return result
    result.pdb2gmx_gro = pdb2gmx_gro

    # ── Step 4: gmx editconf — box ───────────────────────────────────────
    box_gro = os.path.join(output_dir, "box.gro")

    cmd = ["gmx", "editconf", "-f", pdb2gmx_gro, "-o", box_gro,
           "-bt", "triclinic", "-d", str(box_distance)]

    rc, out = _run(cmd, cwd=output_dir)
    log_lines.append(f"editconf: rc={rc}")
    if rc != 0:
        result.log = "\n".join(log_lines) + "\n" + out
        return result

    # ── Step 5a: EM stage 1 — no constraints ─────────────────────────────
    em1_mdp = os.path.join(output_dir, "em_noconstr.mdp")
    _write_em_mdp(em1_mdp, constraints="none", nsteps=em_steps)

    em1_tpr = os.path.join(output_dir, "em_noconstr.tpr")
    cmd = ["gmx", "grompp", "-f", em1_mdp, "-c", box_gro,
           "-p", pdb2gmx_top, "-o", em1_tpr, "-maxwarn", "10"]
    rc, out = _run(cmd, cwd=output_dir)
    log_lines.append(f"grompp(em1): rc={rc}")
    if rc != 0:
        result.log = "\n".join(log_lines) + "\n" + out
        return result

    cmd = ["gmx", "mdrun", "-v", "-deffnm", "em_noconstr"]
    rc, out = _run(cmd, cwd=output_dir, timeout=600)
    log_lines.append(f"mdrun(em1): rc={rc}")
    if rc != 0:
        result.log = "\n".join(log_lines) + "\n" + out
        return result

    em1_gro = os.path.join(output_dir, "em_noconstr.gro")

    # ── Step 5b: EM stage 2 — h-bonds constraints ────────────────────────
    em2_mdp = os.path.join(output_dir, "em_hbonds.mdp")
    _write_em_mdp(em2_mdp, constraints="h-bonds", nsteps=em_steps)

    em2_tpr = os.path.join(output_dir, "em_hbonds.tpr")
    cmd = ["gmx", "grompp", "-f", em2_mdp, "-c", em1_gro,
           "-p", pdb2gmx_top, "-o", em2_tpr, "-maxwarn", "10"]
    rc, out = _run(cmd, cwd=output_dir)
    log_lines.append(f"grompp(em2): rc={rc}")
    if rc != 0:
        result.log = "\n".join(log_lines) + "\n" + out
        return result

    cmd = ["gmx", "mdrun", "-v", "-deffnm", "em_hbonds"]
    rc, out = _run(cmd, cwd=output_dir, timeout=600)
    log_lines.append(f"mdrun(em2): rc={rc}")
    if rc != 0:
        result.log = "\n".join(log_lines) + "\n" + out
        return result

    em2_gro = os.path.join(output_dir, "em_hbonds.gro")
    em2_log = os.path.join(output_dir, "em_hbonds.log")

    # ── Finalise ──────────────────────────────────────────────────────────
    result.em_max_force = _extract_max_force(em2_log)
    _purge_backup_files(output_dir)

    result.em_gro = em2_gro
    result.em_top = pdb2gmx_top
    result.success = os.path.exists(em2_gro)
    result.log = "\n".join(log_lines)

    return result
