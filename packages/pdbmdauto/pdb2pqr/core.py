"""
pdb2pqr core — pure Python logic, no BoCoFlow dependencies.

Provides:
- pdb2pqr executable discovery (env or PATH)
- CLI command building with all supported flags
- Subprocess execution
- PQR file statistics extraction
- PQR to PDB conversion via MDAnalysis
"""

import os
import subprocess
import sys

import MDAnalysis as mda


# ---------------------------------------------------------------------------
# Executable discovery
# ---------------------------------------------------------------------------


def find_pdb2pqr_executable(custom_path=""):
    """Find the pdb2pqr executable.

    Resolution order:
    1. Custom path (if provided and exists)
    2. Same directory as current Python interpreter
    3. Fallback to bare command name (relies on PATH)

    Args:
        custom_path: Optional user-specified path to pdb2pqr binary.

    Returns:
        Path string to the pdb2pqr executable.

    Raises:
        FileNotFoundError: If custom_path is given but does not exist.
    """
    if custom_path:
        if not os.path.exists(custom_path):
            raise FileNotFoundError(
                f"Custom pdb2pqr path not found: {custom_path}"
            )
        return custom_path

    # Auto-detect from the same environment as Python
    python_bin_dir = os.path.dirname(sys.executable)
    potential_names = (
        ["pdb2pqr.exe", "pdb2pqr"] if os.name == "nt" else ["pdb2pqr"]
    )

    for name in potential_names:
        candidate = os.path.join(python_bin_dir, name)
        if os.path.exists(candidate):
            return candidate

    # Fallback to PATH
    return "pdb2pqr.exe" if os.name == "nt" else "pdb2pqr"


# ---------------------------------------------------------------------------
# Command building
# ---------------------------------------------------------------------------


def build_pdb2pqr_command(
    pdb2pqr_cmd,
    input_pdb,
    output_pqr,
    force_field="AMBER",
    ph=7.0,
    keep_chain=True,
    optimize_hydrogens=True,
    include_header=True,
    use_propka=True,
    log_level="INFO",
):
    """Build the pdb2pqr CLI argument list.

    Args:
        pdb2pqr_cmd: Path to pdb2pqr executable.
        input_pdb: Input PDB file path.
        output_pqr: Output PQR file path.
        force_field: Force field name (AMBER, CHARMM, etc.).
        ph: pH value for protonation state assignment.
        keep_chain: Preserve chain identifiers.
        optimize_hydrogens: Optimize hydrogen positions.
        include_header: Include header in PQR output.
        use_propka: Use PROPKA for pKa predictions.
        log_level: Logging verbosity (DEBUG, INFO, WARNING, ERROR).

    Returns:
        List of command-line arguments.
    """
    cmd = [
        pdb2pqr_cmd,
        "--ff", force_field,
        "--ffout", force_field,
        "--log-level", log_level,
    ]

    if use_propka:
        cmd.extend(["--titration-state-method", "propka"])
        cmd.extend(["--with-ph", str(ph)])

    if keep_chain:
        cmd.append("--keep-chain")

    if not optimize_hydrogens:
        cmd.append("--no-optimize")

    if include_header:
        cmd.append("--include-header")

    cmd.extend([input_pdb, output_pqr])

    return cmd


# ---------------------------------------------------------------------------
# Subprocess execution
# ---------------------------------------------------------------------------


def run_pdb2pqr(cmd, working_dir):
    """Execute pdb2pqr as a subprocess.

    Args:
        cmd: Command-line argument list.
        working_dir: Working directory for the process.

    Returns:
        Tuple of (stdout_text, return_code).
    """
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=working_dir,
    )
    return result.stdout or "", result.returncode


# ---------------------------------------------------------------------------
# PQR statistics
# ---------------------------------------------------------------------------


def extract_pqr_statistics(pqr_path):
    """Parse PQR file and count atoms/hydrogens.

    Args:
        pqr_path: Path to the PQR file.

    Returns:
        Dict with keys: total_atoms, hydrogen_atoms, non_hydrogen_atoms.
    """
    atom_count = 0
    hydrogen_count = 0

    with open(pqr_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_count += 1
                atom_name = line[12:16].strip()
                if atom_name.startswith("H"):
                    hydrogen_count += 1

    return {
        "total_atoms": atom_count,
        "hydrogen_atoms": hydrogen_count,
        "non_hydrogen_atoms": atom_count - hydrogen_count,
    }


# ---------------------------------------------------------------------------
# PQR to PDB conversion
# ---------------------------------------------------------------------------


def convert_pqr_to_pdb(pqr_path, pdb_path):
    """Convert PQR to PDB using MDAnalysis.

    Args:
        pqr_path: Input PQR file path.
        pdb_path: Output PDB file path.
    """
    system = mda.Universe(pqr_path, format="PQR", in_memory=True)

    # Remove the segment IDs to avoid 'SYST' being printed
    for atom in system.atoms:
        atom.segment.segid = ""

    system.atoms.write(pdb_path)
