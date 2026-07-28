"""
GROMACS MD Run Node (HPC-enabled)

This node executes GROMACS molecular dynamics simulations with support for both
local and remote (HPC/SLURM) execution using the HPCNodeBase infrastructure.

Architecture:
    - Level 1 (core.py): Pure Python functions for running GROMACS
    - Level 2 (this file): BoCoFlow wrapper with HPCNodeBase for HPC support

The node supports:
1. Local execution using core.py functions (gmx grompp + gmx mdrun)
2. Remote HPC execution via user-provided SLURM scripts

For remote execution, users provide their own SLURM script with template
variables ({{VAR_NAME}}) that get replaced at runtime.

Template Variables Available:
- {{JOB_NAME}}: Auto-generated job name
- {{REMOTE_WORK_DIR}}: Remote working directory
- {{OUTPUT_DIR}}: Output directory (same as REMOTE_WORK_DIR)
- {{WORKING_DIR}}: Working directory (same as REMOTE_WORK_DIR)
- {{NODE_ID}}: BoCoFlow node ID
- {{RUN_LABEL}}: User-specified run label (e.g., md, nvt, npt)
- {{INPUT_TOP_FILE}}: Topology file basename
- {{INPUT_GRO_FILE}}: Structure file basename
- {{INPUT_MDP_FILE}}: Parameters file basename
- {{INPUT_NDX_FILE}}: Index file basename (empty if not provided)

Example SLURM Script:
    #!/bin/bash -l
    #SBATCH --job-name={{JOB_NAME}}
    #SBATCH --output={{OUTPUT_DIR}}/slurm-%j.out
    #SBATCH --partition=gpu
    #SBATCH --nodes=2
    #SBATCH --ntasks-per-node=16
    #SBATCH --gres=gpu:2
    #SBATCH --time=48:00:00

    module load GROMACS/2023.3

    cd {{WORKING_DIR}}

    gmx grompp -f {{INPUT_MDP_FILE}} -c {{INPUT_GRO_FILE}} \\
        -r {{INPUT_GRO_FILE}} -p {{INPUT_TOP_FILE}} \\
        -n {{INPUT_NDX_FILE}} -o {{RUN_LABEL}}.tpr -maxwarn 10

    srun gmx_mpi mdrun -deffnm {{RUN_LABEL}} -v

See: dev-notes/slurm-user-script-and-base-class-analysis.md for design details.
See: dev-notes/node-wrapper-mechanism-design.md for architecture patterns.
"""

import glob
import json
import os
from datetime import datetime
from typing import Dict, List

from bocoflow_core.hpc_node import HPCNodeBase
from bocoflow_core.logger import log_message
from bocoflow_core.node import NodeException, NodeResult
from bocoflow_core.parameters import BooleanParameter, FileParameterEdit, StringParameter

try:
    from .core import run_md_simulation
except ImportError:
    # Stage 2. node_runner puts the node's directory on sys.path and imports
    # node.py as a TOP-LEVEL module, so there is no package for `.core` to be
    # relative to. Without this the next stage ran instead and every symbol
    # below was None by the time execute() called it.
    try:
        from core import run_md_simulation
    except ImportError:
        # Server environment: core deps not available.
        # OPTIONS still work — functions only called at execution time.
        run_md_simulation = None

class GmxMdRun(HPCNodeBase):
    """
    GROMACS MD Run node with local and HPC execution support.

    Executes GROMACS molecular dynamics simulations:
    1. Preprocessing with 'gmx grompp' to generate TPR file
    2. Running the simulation with 'gmx mdrun'
    3. For HPC: Manages file transfer and job monitoring via SLURM

    Input Files:
    - TOP file (.top): Topology file with molecular structure and force field
    - GRO file (.gro): Structure file with atomic coordinates
    - MDP file (.mdp): MD parameters file
    - NDX file (.ndx): Index file defining atom groups (optional)
    - ITP files (*.itp): Included topology files (auto-detected)

    Output Files:
    - TPR file (.tpr): Portable binary run input file
    - GRO file (.gro): Final structure coordinates
    - XTC file (.xtc): Compressed trajectory
    - EDR file (.edr): Energy data
    - LOG file (.log): Simulation log
    """

    # Metadata
    name = "GROMACS MD Run"
    node_key = "GmxMdRun"
    category = "simulation"
    tags = ["molecular-dynamics", "md", "gromacs", "hpc", "slurm"]

    # Connection ports (default: 1 in, 1 out from base class)
    num_in = 1
    num_out = 1


    # Node options: Combine HPC options with node-specific options
    OPTIONS = {
        # Include standard HPC options (execution_mode, hpc_profile, slurm_script, force_resubmit)
        **HPCNodeBase.HPC_OPTIONS,
        # Node-specific parameters
        "case_name": StringParameter(
            label="Case Name",
            default="",
            docstring="Name identifier for this simulation case (optional, uses predecessor if empty)",
            optional=True,
        ),
        "run_label": StringParameter(
            label="Run Label",
            default="md",
            docstring="Label for output files (e.g., md, nvt, npt, production)",
        ),
        "input_top_file": FileParameterEdit(
            label="Topology File (.top)",
            docstring="GROMACS topology file containing molecular structure and force field parameters",
        ),
        "input_gro_file": FileParameterEdit(
            label="Structure File (.gro)",
            docstring="GROMACS structure file with initial atomic coordinates",
        ),
        "input_mdp_file": FileParameterEdit(
            label="Parameters File (.mdp)",
            docstring="GROMACS MD parameters file specifying simulation settings",
        ),
        "input_ndx_file": FileParameterEdit(
            label="Index File (.ndx)",
            docstring="GROMACS index file defining atom groups",
            optional=True,
        ),
        "num_threads": StringParameter(
            label="Number of Threads",
            default="0",
            docstring="Number of OpenMP threads (0 = auto-detect)",
            optional=True,
        ),
        "max_warnings": StringParameter(
            label="Max Warnings",
            default="10",
            docstring="Maximum warnings to allow in grompp",
            optional=True,
        ),
        # force_to_run is now inherited from Node.BASE_OPTIONS (bocoflow-core/node.py)
        # "force_to_run": BooleanParameter(
        #     label="Force to Run",
        #     default=False,
        #     docstring="Execute regardless of previous execution record",
        #     optional=True,
        # ),
    }

    # ==========================================================================
    # HPCNodeBase abstract method implementations
    # ==========================================================================

    def get_input_files(self, flow_vars: dict) -> List[str]:
        """
        Return list of input files to transfer to remote HPC.

        Includes:
        - Topology file (.top)
        - Structure file (.gro)
        - Parameters file (.mdp)
        - Index file (.ndx) if provided
        - All .itp files from the working directory
        """
        files = []

        # Required input files
        top_file = flow_vars["input_top_file"].get_value()
        gro_file = flow_vars["input_gro_file"].get_value()
        mdp_file = flow_vars["input_mdp_file"].get_value()

        if top_file:
            files.append(self.resolve_path(top_file))
        if gro_file:
            files.append(self.resolve_path(gro_file))
        if mdp_file:
            files.append(self.resolve_path(mdp_file))

        # Optional index file
        ndx_file = flow_vars["input_ndx_file"].get_value()
        if ndx_file:
            files.append(self.resolve_path(ndx_file))

        # Include all .itp files from working directory (topology includes)
        if gro_file:
            working_dir = os.path.dirname(self.resolve_path(gro_file))
            itp_files = glob.glob(os.path.join(working_dir, "*.itp"))
            files.extend(itp_files)

        return files

    def get_output_files(self, flow_vars: dict) -> List[str]:
        """
        Return list of output file patterns to retrieve from HPC.
        """
        run_label = flow_vars["run_label"].get_value()
        return [
            f"{run_label}.tpr",
            f"{run_label}.gro",
            f"{run_label}.xtc",
            f"{run_label}.edr",
            f"{run_label}.log",
        ]

    def get_output_files_by_category(self, flow_vars: dict) -> Dict[str, List[str]]:
        """
        Return output files grouped by size category for selective download.

        Categories for GROMACS MD output:
        - "essential": Log and final structure (< 10 MB typically)
        - "standard": Run input and energy data (< 100 MB typically)
        - "large": Trajectory files (can be 1-100+ GB for long simulations)

        Users can enable "Download Large Files" option to also download
        trajectory files (.xtc, .trr) which can be very large.
        """
        run_label = flow_vars["run_label"].get_value()
        return {
            "essential": [
                f"{run_label}.log",  # ~1-10 MB - simulation log
                f"{run_label}.gro",  # ~1-5 MB - final structure
            ],
            "standard": [
                f"{run_label}.tpr",  # ~10-50 MB - run input (needed for analysis)
                f"{run_label}.edr",  # ~10-100 MB - energy data
            ],
            "large": [
                f"{run_label}.xtc",  # 1-100+ GB - compressed trajectory
                f"{run_label}.trr",  # 1-100+ GB - full precision trajectory
                f"{run_label}.cpt",  # ~10-100 MB - checkpoint file
            ],
        }

    def get_template_variables(self, flow_vars: dict) -> Dict[str, str]:
        """
        Return template variables for SLURM script processing.

        Variables provided:
        - RUN_LABEL: User-specified run label
        - INPUT_TOP_FILE: Basename of topology file
        - INPUT_GRO_FILE: Basename of structure file
        - INPUT_MDP_FILE: Basename of parameters file
        - INPUT_NDX_FILE: Basename of index file (empty if not provided)
        """
        run_label = flow_vars["run_label"].get_value()

        top_file = flow_vars["input_top_file"].get_value() or ""
        gro_file = flow_vars["input_gro_file"].get_value() or ""
        mdp_file = flow_vars["input_mdp_file"].get_value() or ""
        ndx_file = flow_vars["input_ndx_file"].get_value() or ""

        return {
            "RUN_LABEL": run_label,
            "INPUT_TOP_FILE": os.path.basename(top_file),
            "INPUT_GRO_FILE": os.path.basename(gro_file),
            "INPUT_MDP_FILE": os.path.basename(mdp_file),
            "INPUT_NDX_FILE": os.path.basename(ndx_file) if ndx_file else "",
        }

    def run_local(self, predecessor_data: list, flow_vars: dict) -> dict:
        """
        Execute GROMACS MD simulation locally using core.py functions.

        This method delegates to the Level 1 core functions for actual execution,
        following the Node Wrapper Mechanism pattern.
        """
        log_message("Starting local GROMACS MD execution")

        try:
            result = NodeResult()
            result.metadata["execution_time"] = datetime.now().isoformat()

            # Get parameters from predecessor and flow_vars
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}
            case_name = flow_vars["case_name"].get_value() or input_data.get(
                "case_name", "protein"
            )
            run_label = flow_vars["run_label"].get_value()

            # Get and resolve file paths
            top_file = self.resolve_path(flow_vars["input_top_file"].get_value())
            gro_file = self.resolve_path(flow_vars["input_gro_file"].get_value())
            mdp_file = self.resolve_path(flow_vars["input_mdp_file"].get_value())

            ndx_file = flow_vars["input_ndx_file"].get_value()
            if ndx_file:
                ndx_file = self.resolve_path(ndx_file)

            # Parse optional parameters
            num_threads_str = flow_vars.get("num_threads")
            num_threads = 0
            if num_threads_str:
                try:
                    num_threads = int(num_threads_str.get_value())
                except (ValueError, AttributeError):
                    num_threads = 0

            max_warnings_str = flow_vars.get("max_warnings")
            max_warnings = 10
            if max_warnings_str:
                try:
                    max_warnings = int(max_warnings_str.get_value())
                except (ValueError, AttributeError):
                    max_warnings = 10

            working_dir = os.path.dirname(gro_file)

            log_message(f"Case: {case_name}, Run label: {run_label}")
            log_message(f"Working directory: {working_dir}")
            log_message(f"Input TOP: {top_file}")
            log_message(f"Input GRO: {gro_file}")
            log_message(f"Input MDP: {mdp_file}")
            log_message(f"Input NDX: {ndx_file or 'None'}")

            # Record input files
            result.files["input"].update(
                {
                    "input_top_file": self.format_output_path(top_file),
                    "input_gro_file": self.format_output_path(gro_file),
                    "input_mdp_file": self.format_output_path(mdp_file),
                }
            )
            if ndx_file:
                result.files["input"]["input_ndx_file"] = self.format_output_path(
                    ndx_file
                )

            # Execute using Level 1 core functions
            sim_result = run_md_simulation(
                top_file=top_file,
                gro_file=gro_file,
                mdp_file=mdp_file,
                working_dir=working_dir,
                run_label=run_label,
                ndx_file=ndx_file,
                num_threads=num_threads,
                max_warnings=max_warnings,
                verbose=True,
            )

            # Check for simulation errors
            if not sim_result.success:
                raise NodeException("gromacs", sim_result.message)

            # Collect output files
            output_files = {}
            if sim_result.tpr_file:
                output_files["tpr_file"] = self.format_output_path(sim_result.tpr_file)
            if sim_result.gro_file:
                output_files["gro_file"] = self.format_output_path(sim_result.gro_file)
            if sim_result.xtc_file:
                output_files["xtc_file"] = self.format_output_path(sim_result.xtc_file)
            if sim_result.edr_file:
                output_files["edr_file"] = self.format_output_path(sim_result.edr_file)
            if sim_result.log_file:
                output_files["log_file"] = self.format_output_path(sim_result.log_file)

            # Update result
            result.files["output"].update(output_files)
            result.data.update(
                {
                    "case_name": case_name,
                    "run_label": run_label,
                    "working_path": self.format_output_path(working_dir),
                    "output_files": output_files,
                    "grompp_returncode": sim_result.grompp_returncode,
                    "mdrun_returncode": sim_result.mdrun_returncode,
                }
            )

            result.success = True
            result.message = sim_result.message

            log_message(
                f"Successfully completed GROMACS MD for {case_name}/{run_label}"
            )
            return json.loads(result.to_json())

        except NodeException:
            raise
        except Exception as e:
            log_message(f"Error in GmxMdRun: {str(e)}", level="error")
            raise NodeException("gromacs_mdrun", str(e))

    def build_hpc_result_data(self, flow_vars: dict, local_dir: str) -> dict:
        """
        Build result data after HPC job completion.

        This method is called by HPCNodeBase after files have been retrieved
        from the HPC cluster to populate the result data field.

        Args:
            flow_vars: Node flow variables
            local_dir: Local directory where output files were downloaded

        Returns:
            Dictionary with node-specific result data
        """
        log_message(f"[GmxMdRun] Building HPC result data from {local_dir}")

        run_label = flow_vars["run_label"].get_value()
        case_name = flow_vars["case_name"].get_value() or "protein"

        # Build output files dict from retrieved files
        output_files = {}
        file_mappings = [
            (".tpr", "tpr_file"),
            (".gro", "gro_file"),
            (".xtc", "xtc_file"),
            (".edr", "edr_file"),
            (".log", "log_file"),
        ]

        for ext, key in file_mappings:
            file_path = os.path.join(local_dir, f"{run_label}{ext}")
            if os.path.exists(file_path):
                output_files[key] = self.format_output_path(file_path)
                log_message(f"[GmxMdRun] Found HPC output: {key} = {file_path}")
            else:
                log_message(
                    f"[GmxMdRun] HPC output not found: {run_label}{ext}",
                    level="warning",
                )

        result_data = {
            "case_name": case_name,
            "run_label": run_label,
            "working_path": self.format_output_path(local_dir),
            "output_files": output_files,
            "execution_mode": "remote_hpc",
        }

        log_message(f"[GmxMdRun] HPC result data keys: {list(result_data.keys())}")
        return result_data
