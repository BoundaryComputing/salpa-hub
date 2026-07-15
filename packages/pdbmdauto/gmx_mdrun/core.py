"""
GROMACS MD Run - Core Functions (Level 1)

Pure Python functions for running GROMACS molecular dynamics simulations.
These functions are independent of BoCoFlow and can be:
- Tested with standard pytest
- Used from command line
- Called from Jupyter notebooks
- Imported into other projects

For BoCoFlow integration, see node.py which wraps these functions.

Architecture:
    This module follows the Node Wrapper Mechanism pattern:
    - Level 1 (this file): Pure Python functions with no BoCoFlow dependencies
    - Level 2 (node.py): BoCoFlow wrapper that calls these functions

See: dev-notes/node-wrapper-mechanism-design.md
"""

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class GromppConfig:
    """Configuration for gmx grompp command.

    Attributes:
        mdp_file: Path to MDP parameters file
        gro_file: Path to structure/coordinate file
        top_file: Path to topology file
        output_tpr: Path for output TPR file
        ndx_file: Optional path to index file
        maxwarn: Maximum warnings to allow (default: 10)
        restraint_file: Optional path to restraint reference structure
    """

    mdp_file: str
    gro_file: str
    top_file: str
    output_tpr: str
    ndx_file: Optional[str] = None
    maxwarn: int = 10
    restraint_file: Optional[str] = None  # -r flag, defaults to gro_file


@dataclass
class MdrunConfig:
    """Configuration for gmx mdrun command.

    Attributes:
        deffnm: Default filename prefix for all output files
        num_threads: Number of OpenMP threads (0 = auto-detect)
        verbose: Show progress during simulation
        gpu_ids: GPU device IDs to use (e.g., "0" or "0,1")
        extra_args: Additional command-line arguments
    """

    deffnm: str
    num_threads: int = 0  # 0 = auto-detect
    verbose: bool = True
    gpu_ids: Optional[str] = None
    extra_args: List[str] = field(default_factory=list)


@dataclass
class SimulationResult:
    """Result of a GROMACS MD simulation.

    Attributes:
        success: Whether the simulation completed successfully
        message: Human-readable status message
        tpr_file: Path to generated TPR file (if created)
        gro_file: Path to output structure file (if created)
        xtc_file: Path to trajectory file (if created)
        edr_file: Path to energy file (if created)
        log_file: Path to log file (if created)
        grompp_returncode: Return code from grompp command
        mdrun_returncode: Return code from mdrun command
        grompp_stderr: stderr output from grompp
        mdrun_stderr: stderr output from mdrun
    """

    success: bool
    message: str
    tpr_file: Optional[str] = None
    gro_file: Optional[str] = None
    xtc_file: Optional[str] = None
    edr_file: Optional[str] = None
    log_file: Optional[str] = None
    grompp_returncode: int = 0
    mdrun_returncode: int = 0
    grompp_stderr: str = ""
    mdrun_stderr: str = ""


def build_grompp_command(config: GromppConfig) -> List[str]:
    """
    Build the gmx grompp command from configuration.

    Args:
        config: GromppConfig with all grompp settings

    Returns:
        List of command arguments for subprocess

    Example:
        >>> config = GromppConfig(
        ...     mdp_file="md.mdp",
        ...     gro_file="conf.gro",
        ...     top_file="topol.top",
        ...     output_tpr="md.tpr"
        ... )
        >>> cmd = build_grompp_command(config)
        >>> cmd[0:2]
        ['gmx', 'grompp']
    """
    cmd = [
        "gmx",
        "grompp",
        "-f",
        config.mdp_file,
        "-c",
        config.gro_file,
        "-r",
        config.restraint_file or config.gro_file,
        "-p",
        config.top_file,
        "-o",
        config.output_tpr,
        "-maxwarn",
        str(config.maxwarn),
    ]

    if config.ndx_file:
        cmd.extend(["-n", config.ndx_file])

    return cmd


def build_mdrun_command(config: MdrunConfig) -> List[str]:
    """
    Build the gmx mdrun command from configuration.

    Args:
        config: MdrunConfig with all mdrun settings

    Returns:
        List of command arguments for subprocess

    Example:
        >>> config = MdrunConfig(deffnm="md", num_threads=4)
        >>> cmd = build_mdrun_command(config)
        >>> "-nt" in cmd and "4" in cmd
        True
    """
    cmd = ["gmx", "mdrun", "-deffnm", config.deffnm]

    if config.num_threads > 0:
        cmd.extend(["-nt", str(config.num_threads)])

    if config.verbose:
        cmd.append("-v")

    if config.gpu_ids:
        cmd.extend(["-gpu_id", config.gpu_ids])

    cmd.extend(config.extra_args)

    return cmd


def run_grompp(
    config: GromppConfig,
    working_dir: str,
    capture_output: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run gmx grompp to generate TPR file.

    Args:
        config: GromppConfig with grompp settings
        working_dir: Directory to run command in
        capture_output: Whether to capture stdout/stderr

    Returns:
        CompletedProcess with return code and output

    Raises:
        FileNotFoundError: If input files don't exist
    """
    # Validate input files exist
    for filepath, name in [
        (config.mdp_file, "MDP"),
        (config.gro_file, "GRO"),
        (config.top_file, "TOP"),
    ]:
        if not os.path.isabs(filepath):
            filepath = os.path.join(working_dir, filepath)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"{name} file not found: {filepath}")

    if config.ndx_file:
        ndx_path = config.ndx_file
        if not os.path.isabs(ndx_path):
            ndx_path = os.path.join(working_dir, ndx_path)
        if not os.path.exists(ndx_path):
            raise FileNotFoundError(f"NDX file not found: {ndx_path}")

    cmd = build_grompp_command(config)

    return subprocess.run(
        cmd,
        cwd=working_dir,
        capture_output=capture_output,
        text=True,
    )


def run_mdrun(
    config: MdrunConfig,
    working_dir: str,
    capture_output: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run gmx mdrun to execute MD simulation.

    Args:
        config: MdrunConfig with mdrun settings
        working_dir: Directory to run command in
        capture_output: Whether to capture stdout/stderr

    Returns:
        CompletedProcess with return code and output

    Raises:
        FileNotFoundError: If TPR file doesn't exist
    """
    # Validate TPR file exists
    tpr_file = os.path.join(working_dir, f"{config.deffnm}.tpr")
    if not os.path.exists(tpr_file):
        raise FileNotFoundError(f"TPR file not found: {tpr_file}")

    cmd = build_mdrun_command(config)

    return subprocess.run(
        cmd,
        cwd=working_dir,
        capture_output=capture_output,
        text=True,
    )


def run_md_simulation(
    top_file: str,
    gro_file: str,
    mdp_file: str,
    working_dir: str,
    run_label: str = "md",
    ndx_file: Optional[str] = None,
    num_threads: int = 0,
    max_warnings: int = 10,
    verbose: bool = True,
) -> SimulationResult:
    """
    Run a complete GROMACS MD simulation (grompp + mdrun).

    This is the main entry point for running MD simulations.
    It's a pure Python function with no BoCoFlow dependencies.

    Args:
        top_file: Path to topology file (.top)
        gro_file: Path to structure file (.gro)
        mdp_file: Path to parameters file (.mdp)
        working_dir: Working directory for simulation
        run_label: Label for output files (default: "md")
        ndx_file: Optional path to index file (.ndx)
        num_threads: Number of threads (0 = auto)
        max_warnings: Maximum warnings for grompp
        verbose: Show verbose mdrun output

    Returns:
        SimulationResult with success status and output file paths

    Example:
        >>> result = run_md_simulation(
        ...     top_file="topol.top",
        ...     gro_file="conf.gro",
        ...     mdp_file="md.mdp",
        ...     working_dir="/path/to/simulation",
        ...     run_label="nvt"
        ... )
        >>> result.success
        True
    """
    # Output file paths
    tpr_file = os.path.join(working_dir, f"{run_label}.tpr")
    out_gro = os.path.join(working_dir, f"{run_label}.gro")
    out_xtc = os.path.join(working_dir, f"{run_label}.xtc")
    out_edr = os.path.join(working_dir, f"{run_label}.edr")
    out_log = os.path.join(working_dir, f"{run_label}.log")

    # Step 1: Run grompp
    grompp_config = GromppConfig(
        mdp_file=mdp_file,
        gro_file=gro_file,
        top_file=top_file,
        output_tpr=tpr_file,
        ndx_file=ndx_file,
        maxwarn=max_warnings,
    )

    try:
        grompp_result = run_grompp(grompp_config, working_dir)
    except FileNotFoundError as e:
        return SimulationResult(
            success=False,
            message=str(e),
            grompp_returncode=-1,
        )

    if grompp_result.returncode != 0:
        return SimulationResult(
            success=False,
            message=f"grompp failed: {grompp_result.stderr or grompp_result.stdout}",
            grompp_returncode=grompp_result.returncode,
            grompp_stderr=grompp_result.stderr,
        )

    if not os.path.exists(tpr_file):
        return SimulationResult(
            success=False,
            message="TPR file was not created",
            grompp_returncode=grompp_result.returncode,
        )

    # Step 2: Run mdrun
    mdrun_config = MdrunConfig(
        deffnm=run_label,
        num_threads=num_threads,
        verbose=verbose,
    )

    try:
        mdrun_result = run_mdrun(mdrun_config, working_dir)
    except FileNotFoundError as e:
        return SimulationResult(
            success=False,
            message=str(e),
            tpr_file=tpr_file,
            grompp_returncode=grompp_result.returncode,
            mdrun_returncode=-1,
        )

    if mdrun_result.returncode != 0:
        return SimulationResult(
            success=False,
            message=f"mdrun failed: {mdrun_result.stderr or mdrun_result.stdout}",
            tpr_file=tpr_file,
            grompp_returncode=grompp_result.returncode,
            mdrun_returncode=mdrun_result.returncode,
            mdrun_stderr=mdrun_result.stderr,
        )

    # Collect output files
    return SimulationResult(
        success=True,
        message=f"MD simulation '{run_label}' completed successfully",
        tpr_file=tpr_file if os.path.exists(tpr_file) else None,
        gro_file=out_gro if os.path.exists(out_gro) else None,
        xtc_file=out_xtc if os.path.exists(out_xtc) else None,
        edr_file=out_edr if os.path.exists(out_edr) else None,
        log_file=out_log if os.path.exists(out_log) else None,
        grompp_returncode=grompp_result.returncode,
        mdrun_returncode=mdrun_result.returncode,
    )


def check_gromacs_available() -> bool:
    """
    Check if GROMACS is installed and available in PATH.

    Returns:
        True if 'gmx' command is available, False otherwise
    """
    try:
        result = subprocess.run(
            ["gmx", "--version"],
            capture_output=True,
            text=True,
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def get_gromacs_version() -> Optional[str]:
    """
    Get the installed GROMACS version.

    Returns:
        Version string (e.g., "2023.3") or None if not available
    """
    try:
        result = subprocess.run(
            ["gmx", "--version"],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            # Parse version from output
            for line in result.stdout.split("\n"):
                if "GROMACS version" in line:
                    return line.split()[-1]
        return None
    except FileNotFoundError:
        return None
