"""
pdb2pqr — BoCoFlow node wrapper.

Converts PDB files to PQR format with protonation states assigned using PDB2PQR.
Optionally generates a protonated PDB file from the PQR output via MDAnalysis.

Features: multiple force fields, PROPKA pKa predictions, hydrogen optimization.
"""

import os
from datetime import datetime

from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter,
    FileParameterEdit,
    FloatParameter,
    FolderParameter,
    SelectParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log

try:
    from .core import (
        build_pdb2pqr_command,
        convert_pqr_to_pdb,
        extract_pqr_statistics,
        find_pdb2pqr_executable,
        run_pdb2pqr,
    )
except ImportError:
    # Stage 2. node_runner puts the node's directory on sys.path and imports
    # node.py as a TOP-LEVEL module, so there is no package for `.core` to be
    # relative to. Without this the next stage ran instead and every symbol
    # below was None by the time execute() called it.
    try:
        from core import (
            build_pdb2pqr_command,
            convert_pqr_to_pdb,
            extract_pqr_statistics,
            find_pdb2pqr_executable,
            run_pdb2pqr,
        )
    except ImportError:
        # Server environment: pdb2pqr not installed.
        # OPTIONS still work — functions only called in PIXI_SUBPROCESS.
        build_pdb2pqr_command = convert_pqr_to_pdb = None
        extract_pqr_statistics = find_pdb2pqr_executable = run_pdb2pqr = None

class PDB2PQR(Node):
    """
    Converts PDB files to PQR format with protonation states assigned.

    PDB2PQR adds missing hydrogen atoms, assigns protonation states based on
    pH, and calculates partial charges. Useful for preparing structures from
    ESMFold or other prediction tools for molecular dynamics simulations.

    Features:
    - Adds missing hydrogen atoms
    - Assigns protonation states based on pH via PROPKA
    - Supports AMBER, CHARMM, and other force fields
    - Optionally generates a protonated PDB file from PQR output

    Input: PDB file (typically from structure prediction)
    Output: PQR file with hydrogens and charges, protonated PDB file (optional)
    """

    # NOTE: Metadata (name, hashtags, num_in, num_out) comes from meta.toml.
    # NOTE: force_to_run is inherited from Node.BASE_OPTIONS — do NOT add it here.

    category = "io"
    tags = ["protein", "pdb", "processing", "protonation", "force-field"]

    name = "PDB to PQR Converter"
    node_key = "PDB2PQR"

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name",
            docstring="Name of the case/system for PQR conversion",
        ),
        "input_pdb": FileParameterEdit(
            "Input PDB File",
            docstring="PDB file to convert (e.g., from ESMFold prediction)",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Output directory for PQR and log files",
        ),
        "force_field": SelectParameter(
            "Force Field",
            default="AMBER",
            options=["AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON"],
            docstring="Force field for atom typing and charge assignment",
        ),
        "ph": FloatParameter(
            "pH Value",
            default=7.0,
            docstring="pH value for protonation state assignment",
        ),
        "keep_chain": BooleanParameter(
            "Keep Chain IDs",
            default=True,
            docstring="Preserve chain identifiers from input PDB",
        ),
        "optimize_hydrogens": BooleanParameter(
            "Optimize Hydrogens",
            default=True,
            docstring="Optimize hydrogen positions after addition",
        ),
        "include_header": BooleanParameter(
            "Include Header",
            default=True,
            docstring="Include header information in output PQR",
        ),
        "use_propka": BooleanParameter(
            "Use PROPKA",
            default=True,
            docstring="Use PROPKA for pKa predictions and titration states",
        ),
        "log_level": SelectParameter(
            "Log Level",
            default="INFO",
            options=["DEBUG", "INFO", "WARNING", "ERROR"],
            docstring="Logging verbosity for pdb2pqr",
        ),
        "custom_pdb2pqr_path": StringParameter(
            "Custom PDB2PQR Path",
            default="",
            docstring="Optional custom path to pdb2pqr executable (overrides auto-detection)",
        ),
        "generate_pdb": BooleanParameter(
            "Generate PDB from PQR",
            default=True,
            docstring="Generate a protonated PDB file from the PQR output using MDAnalysis",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute the PDB2PQR conversion."""
        stream_log("Starting PDB2PQR conversion...", node_id=self.node_id, progress=0)

        try:
            result = NodeResult()
            result.metadata.update(
                {
                    "case_name": flow_vars["case_name"].get_value(),
                    "execution_time": datetime.now().isoformat(),
                }
            )

            # --- Read parameters ---
            case_name = flow_vars["case_name"].get_value()
            input_pdb = self.resolve_path(flow_vars["input_pdb"].get_value())
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            force_field = flow_vars["force_field"].get_value()
            ph = flow_vars["ph"].get_value()
            keep_chain = flow_vars["keep_chain"].get_value()
            optimize_hydrogens = flow_vars["optimize_hydrogens"].get_value()
            include_header = flow_vars["include_header"].get_value()
            use_propka = flow_vars["use_propka"].get_value()
            log_level = flow_vars["log_level"].get_value()
            custom_pdb2pqr_path = flow_vars["custom_pdb2pqr_path"].get_value().strip()
            generate_pdb = flow_vars["generate_pdb"].get_value()

            # --- Resolve case_name from predecessor ---
            if not case_name and predecessor_data and predecessor_data[0]:
                input_data = predecessor_data[0]
                case_name = input_data.get("case_name", "protein")

            log_message(f"Case: {case_name}, Input: {input_pdb}")
            log_message(f"Output: {output_dir}")

            result.metadata.update({"output_dir": self.format_output_path(output_dir)})

            # --- Validate input file ---
            if not os.path.exists(input_pdb):
                raise NodeException(
                    "execution", f"Input PDB file not found: {input_pdb}"
                )

            os.makedirs(output_dir, exist_ok=True)

            # --- Build output file paths ---
            output_pqr = os.path.join(output_dir, f"{case_name}_structure.pqr")
            output_pdb = (
                os.path.join(output_dir, f"{case_name}_protonated.pdb")
                if generate_pdb
                else None
            )
            propka_output = (
                os.path.join(output_dir, f"{case_name}_propka.out")
                if use_propka
                else None
            )

            # --- Find pdb2pqr executable ---
            stream_log(
                "Locating pdb2pqr executable...",
                node_id=self.node_id,
                progress=10,
            )
            pdb2pqr_cmd = find_pdb2pqr_executable(custom_pdb2pqr_path)
            log_message(f"Using pdb2pqr: {pdb2pqr_cmd}")

            # --- Build and run command ---
            cmd = build_pdb2pqr_command(
                pdb2pqr_cmd,
                input_pdb,
                output_pqr,
                force_field=force_field,
                ph=ph,
                keep_chain=keep_chain,
                optimize_hydrogens=optimize_hydrogens,
                include_header=include_header,
                use_propka=use_propka,
                log_level=log_level,
            )

            stream_log(
                "Running pdb2pqr...", node_id=self.node_id, progress=20
            )
            log_message(f"Command: {' '.join(cmd)}")

            stdout, returncode = run_pdb2pqr(cmd, output_dir)

            if stdout:
                log_message(f"pdb2pqr output: {stdout}")

            if returncode != 0:
                raise NodeException(
                    "execution",
                    f"pdb2pqr failed with return code {returncode}. Output: {stdout}",
                )

            if not os.path.exists(output_pqr):
                raise NodeException(
                    "execution",
                    "PQR output file was not created. Check workflow logs for errors.",
                )

            stream_log(
                "Extracting PQR statistics...",
                node_id=self.node_id,
                progress=60,
            )

            # --- Extract statistics ---
            stats = extract_pqr_statistics(output_pqr)

            # --- Convert PQR to PDB (optional) ---
            if generate_pdb and output_pdb:
                stream_log(
                    "Converting PQR to protonated PDB...",
                    node_id=self.node_id,
                    progress=80,
                )
                log_message(f"Converting PQR to PDB: {output_pqr} -> {output_pdb}")
                convert_pqr_to_pdb(output_pqr, output_pdb)
                log_message(f"Protonated PDB file created: {output_pdb}")

            # --- Record input files ---
            result.files["input"].update(
                {"input_pdb": self.format_output_path(input_pdb)}
            )

            # --- Store processing results ---
            result.data = {
                "case_name": case_name,
                "conversion_parameters": {
                    "force_field": force_field,
                    "ph": ph,
                    "keep_chain": keep_chain,
                    "optimize_hydrogens": optimize_hydrogens,
                    "use_propka": use_propka,
                },
                "statistics": stats,
                "output_files": {
                    "pqr": self.format_output_path(output_pqr),
                },
                "working_path": self.format_output_path(output_dir),
            }

            if generate_pdb and output_pdb and os.path.exists(output_pdb):
                result.data["output_files"]["protonated_pdb"] = (
                    self.format_output_path(output_pdb)
                )

            if use_propka and propka_output and os.path.exists(propka_output):
                result.data["output_files"]["propka"] = (
                    self.format_output_path(propka_output)
                )

            # --- Set output file paths ---
            result.files["output"] = {
                "pqr": self.format_output_path(output_pqr),
            }

            if generate_pdb and output_pdb and os.path.exists(output_pdb):
                result.files["output"]["protonated_pdb"] = (
                    self.format_output_path(output_pdb)
                )

            if use_propka and propka_output and os.path.exists(propka_output):
                result.files["output"]["propka"] = (
                    self.format_output_path(propka_output)
                )

            result.success = True
            pdb_msg = (
                " and protonated PDB generated"
                if generate_pdb and output_pdb and os.path.exists(output_pdb)
                else ""
            )
            result.message = (
                f"PDB2PQR conversion completed for {case_name} "
                f"({stats['total_atoms']} atoms, "
                f"{stats['hydrogen_atoms']} hydrogens added){pdb_msg}"
            )

            stream_log(
                result.message, node_id=self.node_id, progress=100
            )
            return result.to_json()

        except NodeException:
            raise
        except Exception as e:
            log_message(f"Error in PDB2PQR: {str(e)}")
            raise NodeException("pdb2pqr conversion", str(e))
