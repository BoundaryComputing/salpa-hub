"""Module 1: Protein preparation.

Clean PDB → protonate with pdb2pqr → convert to PDBQT with prepare_receptor4.

All functions take explicit paths and tool locations as arguments.
No os.chdir(), no global env vars, no god-object.
"""

import logging
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def clean_pdb(input_pdb: Path, output_pdb: Path) -> Path:
    """Remove HETATM lines (ligands, cofactors, waters) from a PDB file.

    Args:
        input_pdb: Path to the input PDB (e.g. protonated).
        output_pdb: Path for the cleaned output PDB.

    Returns:
        The output_pdb path.
    """
    with open(input_pdb) as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if "HETATM" not in line:
                fout.write(line)
    logger.info("Cleaned PDB: removed HETATM lines → %s", output_pdb)
    return output_pdb


def protonate_pdb(
    input_pdb: Path,
    output_pdb: Path,
    ph: float = 7.4,
    drop_water: bool = True,
    pdb2pqr_path: str = "pdb2pqr30",
) -> Path:
    """Protonate a protein PDB at a given pH using pdb2pqr.

    Args:
        input_pdb: Source PDB file.
        output_pdb: Destination for the protonated PDB.
        ph: Target pH for protonation.
        drop_water: If True, remove water molecules.
        pdb2pqr_path: Path or command name for the pdb2pqr executable.

    Returns:
        The output_pdb path.
    """
    if output_pdb.exists():
        logger.info("Protonated PDB already exists: %s", output_pdb)
        return output_pdb

    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        pdb2pqr_path,
        "--noopt",
        "--pdb-output", str(output_pdb),
        "--with-ph", str(ph),
    ]
    if drop_water:
        cmd.append("--drop-water")
    cmd += [str(input_pdb), str(output_pdb)]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        logger.warning("pdb2pqr stderr: %s", result.stderr)
    logger.info("Protonated PDB at pH %.1f → %s", ph, output_pdb)
    return output_pdb


def create_receptor_pdbqt(
    pdb_path: Path,
    pdbqt_path: Path,
    prepare_receptor_script: str | Path | None = None,
    python_path: str = "pythonsh",
    mgltools_dir: str | Path | None = None,
) -> Path:
    """Convert a clean PDB to PDBQT format using AutoDockTools' prepare_receptor4.py.

    Args:
        pdb_path: Input PDB file (cleaned, protonated).
        pdbqt_path: Output PDBQT file.
        prepare_receptor_script: Full path to prepare_receptor4.py.
            If None, derived from *mgltools_dir*.
        python_path: Python interpreter to run the script. Defaults to
            ``pythonsh`` (the MGLTools bundled interpreter). Use ``python3``
            only if the script is compatible.
        mgltools_dir: Path to the AutoDockTools directory containing
            prepare_receptor4.py (used only if *prepare_receptor_script* is None).

    Returns:
        The pdbqt_path.
    """
    if pdbqt_path.exists():
        logger.info("Receptor PDBQT already exists: %s", pdbqt_path)
        return pdbqt_path

    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)

    if prepare_receptor_script is None:
        if mgltools_dir is None:
            raise ValueError("Provide either prepare_receptor_script or mgltools_dir")
        prepare_receptor_script = Path(mgltools_dir) / "prepare_receptor4.py"

    cmd = [
        python_path,
        str(prepare_receptor_script),
        "-U", "nphs",
        "-A", "None",
        "-r", str(pdb_path),
        "-o", str(pdbqt_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(
            "prepare_receptor4 failed (rc=%d):\nstdout: %s\nstderr: %s",
            result.returncode, result.stdout, result.stderr,
        )
        raise RuntimeError(
            f"prepare_receptor4 failed. Check PDB format. stderr: {result.stderr}"
        )
    logger.info("Created receptor PDBQT → %s", pdbqt_path)
    return pdbqt_path


def prepare_protein(
    pdb_path: Path,
    output_dir: Path,
    ph: float = 7.4,
    clean: bool = True,
    pdb2pqr_path: str = "pdb2pqr30",
    python_path: str = "python3",
    mgltools_dir: str | Path | None = None,
) -> dict:
    """Full protein preparation pipeline: protonate → clean → PDBQT.

    Args:
        pdb_path: Input protein PDB.
        output_dir: Working directory for intermediate and output files.
        ph: Protonation pH.
        clean: If True, remove HETATM lines.
        pdb2pqr_path: pdb2pqr executable.
        python_path: Python for running MGLTools scripts.
        mgltools_dir: AutoDockTools directory.

    Returns:
        Dict with keys: ``protonated_pdb``, ``cleaned_pdb``, ``pdbqt``.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = pdb_path.stem

    # Step 1: Protonate
    protonated_pdb = output_dir / f"{stem}_protonated.pdb"
    protonate_pdb(pdb_path, protonated_pdb, ph=ph, drop_water=clean,
                  pdb2pqr_path=pdb2pqr_path)

    # Step 2: Clean (remove HETATM)
    if clean:
        cleaned_pdb = output_dir / f"clean_{stem}.pdb"
        clean_pdb(protonated_pdb, cleaned_pdb)
    else:
        cleaned_pdb = protonated_pdb

    # Step 3: PDBQT
    pdbqt = output_dir / f"clean_{stem}.pdbqt"
    create_receptor_pdbqt(
        cleaned_pdb, pdbqt,
        python_path=python_path,
        mgltools_dir=mgltools_dir,
    )

    return {
        "protonated_pdb": protonated_pdb,
        "cleaned_pdb": cleaned_pdb,
        "pdbqt": pdbqt,
    }
