"""
GROMACS MD Run Node (Local Only) - Level 2 Wrapper

A simplified node for running GROMACS molecular dynamics simulations locally.
Designed for demonstrations, tutorials, and testing purposes.

This node wraps the core functions from core.py, handling:
- Parameter extraction from UI (flow_vars)
- Path resolution (abs:/rel: prefixes)
- Result formatting for BoCoFlow

For the core algorithm (Level 1), see core.py which can be tested independently.
For HPC/SLURM cluster support, use the full 'gmx-mdrun' node instead.

Architecture:
    This module follows the Node Wrapper Mechanism pattern:
    - Level 1 (core.py): Pure Python functions with no BoCoFlow dependencies
    - Level 2 (this file): BoCoFlow wrapper that calls core functions

See: dev-notes/node-wrapper-mechanism-design.md
"""

import os
from datetime import datetime
from pathlib import Path

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter,
    FileParameterEdit,
    FolderParameter,
    IntegerParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import SimulationResult, check_gromacs_available, run_md_simulation
except ImportError:
    # Stage 2. node_runner puts the node's directory on sys.path and imports
    # node.py as a TOP-LEVEL module, so there is no package for `.core` to be
    # relative to. Without this the next stage ran instead and every symbol
    # below was None by the time execute() called it.
    try:
        from core import SimulationResult, check_gromacs_available, run_md_simulation
    except ImportError:
        # Server environment: core deps not available.
        # OPTIONS still work — functions only called at execution time.
        SimulationResult = check_gromacs_available = run_md_simulation = None

class GmxMdRunLocal(Node):
    """
    Simple GROMACS MD simulation node for local execution.

    This is a Level 2 wrapper that:
    1. Extracts parameters from BoCoFlow UI
    2. Resolves file paths
    3. Calls core functions (Level 1)
    4. Formats results for BoCoFlow

    The actual GROMACS execution logic is in core.py and can be
    tested independently without BoCoFlow dependencies.
    """

    # Metadata (loaded from meta.toml, these are fallbacks)
    name = "GROMACS MD Run (Local)"
    node_key = "GmxMdRunLocal"
    category = "simulation"
    tags = ["molecular-dynamics", "md", "gromacs", "local", "tutorial"]

    # Connection ports
    num_in = 1
    num_out = 1

    OPTIONS = {
        "case_name": StringParameter(
            label="Case Name",
            default="",
            docstring="Name for this simulation (optional)",
            optional=True,
        ),
        "run_label": StringParameter(
            label="Run Label",
            default="md",
            docstring="Label for output files (e.g., md, nvt, npt, em)",
        ),
        "input_top_file": FileParameterEdit(
            label="Topology File (.top)",
            default="",
            docstring="Topology. Leave empty: auto-discovers from predecessor.",
        ),
        "input_gro_file": FileParameterEdit(
            label="Structure File (.gro)",
            default="",
            docstring="Relaxed structure. Leave empty: auto-discovers from predecessor.",
        ),
        "input_mdp_file": FileParameterEdit(
            label="Parameters File (.mdp)",
            default="node:demo_data/md.mdp",
            docstring="Production MD parameters (bundled: 1000 steps, 2 ps).",
        ),
        "input_ndx_file": FileParameterEdit(
            label="Index File (.ndx)",
            default="",
            docstring="Index file. Leave empty: auto-discovers from predecessor.",
            optional=True,
        ),
        "output_folder": FolderParameter(
            label="Output Folder",
            default="",
            docstring="Output subfolder. Leave empty: creates gmx_md/ in predecessor working_path.",
            optional=True,
        ),
        "num_threads": IntegerParameter(
            label="Number of Threads",
            default=0,
            docstring="Number of threads for mdrun (0 = auto-detect)",
            optional=True,
        ),
        "max_warnings": IntegerParameter(
            label="Max Warnings",
            default=10,
            docstring="Maximum warnings allowed by grompp",
            optional=True,
        ),
        "verbose": BooleanParameter(
            label="Verbose Output",
            default=True,
            docstring="Show verbose output during simulation",
            optional=True,
        ),
        # force_to_run is now inherited from Node.BASE_OPTIONS (bocoflow-core/node.py)
        # "force_to_run": BooleanParameter(
        #     label="Force to Run",
        #     default=False,
        #     docstring="Execute regardless of previous results",
        #     optional=True,
        # ),
    }

    def execute(self, predecessor_data, flow_vars):
        """
        Execute GROMACS MD simulation locally.

        This method bridges BoCoFlow data flow with core functions:
        1. Extract parameters from flow_vars
        2. Resolve file paths
        3. Call run_md_simulation() from core.py
        4. Format result for BoCoFlow
        """
        stream_log(
            "Starting GROMACS MD Run (Local)",
            node_id=self.node_id,
            component="GmxMdRunLocal",
            progress=0,
        )

        try:
            # Check GROMACS availability
            if not check_gromacs_available():
                raise NodeException(
                    "dependency",
                    "GROMACS not found. Please install GROMACS and ensure "
                    "'gmx' is available in PATH.",
                )

            # Initialize result
            result = NodeResult()
            result.metadata["execution_time"] = datetime.now().isoformat()
            result.metadata["execution_mode"] = "local"

            # ================================================================
            # Extract parameters from flow_vars (BoCoFlow UI)
            # ================================================================
            input_data = predecessor_data[0] if predecessor_data and predecessor_data[0] else {}
            case_name = flow_vars["case_name"].get_value() or input_data.get(
                "case_name", "simulation"
            )
            run_label = flow_vars["run_label"].get_value()
            num_threads = flow_vars["num_threads"].get_value()
            max_warnings = flow_vars["max_warnings"].get_value()
            verbose = flow_vars["verbose"].get_value()

            # ================================================================
            # Resolve file paths (handle abs:/rel: prefixes)
            # Auto-discover from predecessor data if not explicitly set.
            # ================================================================
            top_file = self.resolve_path(flow_vars["input_top_file"].get_value()) or ""
            gro_file = self.resolve_path(flow_vars["input_gro_file"].get_value()) or ""
            mdp_file = self.resolve_path(flow_vars["input_mdp_file"].get_value()) or ""

            ndx_file = flow_vars["input_ndx_file"].get_value()
            ndx_file = self.resolve_path(ndx_file) if ndx_file else ""

            # Auto-discover from predecessor (e.g. gmx_md_relax output)
            if not gro_file and input_data.get("output_gro"):
                gro_file = self.resolve_path(input_data["output_gro"])
            if not top_file and input_data.get("output_top"):
                top_file = self.resolve_path(input_data["output_top"])
            if not ndx_file and input_data.get("output_ndx"):
                ndx_file = self.resolve_path(input_data["output_ndx"])

            # Auto-discover from predecessor working_path
            working_path = input_data.get("working_path", "")
            if working_path:
                search_dir = self.resolve_path(working_path)
                if os.path.isdir(search_dir):
                    for f in sorted(os.listdir(search_dir)):
                        fp = os.path.join(search_dir, f)
                        if not gro_file and (f.endswith(".gro") and "em" not in f and "box" not in f and "solv" not in f):
                            gro_file = fp
                        elif not gro_file and f.endswith("_ion.gro"):
                            gro_file = fp
                        elif not top_file and f.endswith(".top"):
                            top_file = fp
                        elif not ndx_file and f == "index.ndx":
                            ndx_file = fp

            # Use bundled default MDP if none provided
            if not mdp_file:
                node_dir = getattr(self, "_node_dir", None) or os.path.dirname(__file__)
                default_mdp = os.path.join(node_dir, "demo_data", "md.mdp")
                if os.path.exists(default_mdp):
                    mdp_file = default_mdp
                    stream_log("Using bundled default MDP: md.mdp", node_id=self.node_id)

            # Validate required files
            if not gro_file or not top_file:
                missing = [n for n, v in [("GRO", gro_file), ("TOP", top_file), ("MDP", mdp_file)] if not v]
                raise NodeException(
                    "gromacs",
                    f"Missing required files: {', '.join(missing)}. "
                    f"Provide explicit paths or connect to a predecessor node.",
                )

            # Determine output/working directory
            output_folder = flow_vars["output_folder"].get_value()
            if output_folder:
                working_dir = self.resolve_path(output_folder)
                os.makedirs(working_dir, exist_ok=True)
            else:
                working_dir = os.path.dirname(gro_file)

            # Log configuration (streamed to UI in real-time)
            stream_log(f"Case: {case_name}", node_id=self.node_id)
            stream_log(f"Run label: {run_label}", node_id=self.node_id)
            stream_log(f"Working directory: {working_dir}", node_id=self.node_id)
            stream_log(f"Topology: {Path(top_file).name}", node_id=self.node_id)
            stream_log(f"Structure: {Path(gro_file).name}", node_id=self.node_id)
            stream_log(f"Parameters: {Path(mdp_file).name}", node_id=self.node_id)
            if ndx_file:
                stream_log(f"Index: {Path(ndx_file).name}", node_id=self.node_id)

            stream_log(
                "Running grompp (preprocessing)...",
                node_id=self.node_id,
                progress=20,
            )

            # Record input files
            result.files["input"] = {
                "topology": self.format_output_path(top_file),
                "structure": self.format_output_path(gro_file),
                "parameters": self.format_output_path(mdp_file),
            }
            if ndx_file:
                result.files["input"]["index"] = self.format_output_path(ndx_file)

            # ================================================================
            # CALL CORE FUNCTION (Level 1 - Pure Python)
            # ================================================================
            stream_log(
                "Running mdrun (simulation)...",
                node_id=self.node_id,
                progress=40,
            )
            sim_result: SimulationResult = run_md_simulation(
                top_file=top_file,
                gro_file=gro_file,
                mdp_file=mdp_file,
                working_dir=working_dir,
                run_label=run_label,
                ndx_file=ndx_file,
                num_threads=num_threads,
                max_warnings=max_warnings,
                verbose=verbose,
            )
            # ================================================================

            # Check for errors
            if not sim_result.success:
                raise NodeException("gromacs", sim_result.message)

            # ================================================================
            # Format result for BoCoFlow
            # ================================================================
            # All paths are formatted relative to workflow's working directory
            stream_log(
                "Simulation complete, processing output files...",
                node_id=self.node_id,
                progress=80,
            )
            output_files = {}
            if sim_result.tpr_file:
                output_files["tpr"] = self.format_output_path(sim_result.tpr_file)
                stream_log(
                    f"Output: {Path(sim_result.tpr_file).name}", node_id=self.node_id
                )
            if sim_result.gro_file:
                output_files["structure"] = self.format_output_path(sim_result.gro_file)
                stream_log(
                    f"Output: {Path(sim_result.gro_file).name}", node_id=self.node_id
                )
            if sim_result.xtc_file:
                output_files["trajectory"] = self.format_output_path(
                    sim_result.xtc_file
                )
                stream_log(
                    f"Output: {Path(sim_result.xtc_file).name}", node_id=self.node_id
                )
            if sim_result.edr_file:
                output_files["energy"] = self.format_output_path(sim_result.edr_file)
                stream_log(
                    f"Output: {Path(sim_result.edr_file).name}", node_id=self.node_id
                )
            if sim_result.log_file:
                output_files["log"] = self.format_output_path(sim_result.log_file)
                stream_log(
                    f"Output: {Path(sim_result.log_file).name}", node_id=self.node_id
                )

            result.files["output"] = output_files

            # Result data
            result.data = {
                "case_name": case_name,
                "run_label": run_label,
                "output_files": output_files,
                "threads_used": num_threads if num_threads > 0 else "auto",
            }

            result.success = True
            result.message = sim_result.message

            stream_log(
                f"Success: {result.message}",
                node_id=self.node_id,
                progress=100,
            )
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            stream_log(f"Error: {str(e)}", level="error", node_id=self.node_id)
            raise NodeException("gromacs", str(e))
