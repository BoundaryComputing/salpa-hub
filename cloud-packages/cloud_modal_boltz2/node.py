"""
Cloud Modal Boltz-2 - Mode B (BoCoFlow Credits)

Protein structure prediction using Boltz-2 on Modal's H100 GPU infrastructure.

This node is a client stub that calls the BoCoFlow API Gateway,
which then routes requests to Modal endpoints with Proxy Auth tokens.

Unlike modal-user nodes (Mode A), users don't need their own Modal account.
Instead, they pay with BoCoFlow credits.

Architecture:
    1. User authenticates with Firebase (token in BOCOFLOW_CLOUD_AUTH_TOKEN)
    2. This node sends sequence/YAML to API Gateway with Firebase token
    3. API Gateway verifies token, checks credits
    4. API Gateway calls Modal with Proxy Auth tokens
    5. Boltz-2 runs on H100 GPU, returns structure prediction
    6. Result (tarball with CIF files) returns through this node
    7. Credits are deducted from user's account

Based on: https://modal.com/docs/examples/boltz_predict
"""

import base64
import os
from datetime import datetime

import requests
from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeResult
from bocoflow_core.stream_logger import post_with_progress, stream_log
from bocoflow_core.parameters import (
    BooleanParameter,
    FolderParameter,
    StringParameter,
    TextParameter,
)


class CloudModalBoltz2(Node):
    """
    Boltz-2 Protein Structure Prediction (BoCoFlow Credits - Mode B).

    This is a CLIENT STUB - actual computation happens on Modal cloud (H100 GPU).
    The workflow engine executes this node, which calls the API Gateway.

    Boltz-2 is a state-of-the-art protein structure prediction model that can
    predict 3D structures from amino acid sequences.

    Prerequisites:
    - User must be logged in with Firebase
    - User must have sufficient BoCoFlow credits

    No Modal account or 'modal setup' required!
    """

    # NOTE: Metadata (name, hashtags, num_in, num_out) comes from meta.toml.
    # NOTE: EXECUTION_STRATEGY and ENVIRONMENT are auto-detected via shared_environment in meta.toml.

    # API Gateway endpoint (unified route — not Modal directly!)
    API_ENDPOINT = (
        os.environ.get(
            "BOCOFLOW_CLOUD_API_URL",
            "https://bocoflow-api-gateway-823406908684.us-central1.run.app",
        )
        + "/api/cloud/nodes/boltz2/execute"
    )

    OPTIONS = {
        "sequence": TextParameter(
            "Protein Sequence",
            default="",
            docstring=(
                "Protein amino acid sequence (one-letter codes). "
                "Can also come from predecessor node. "
                "Example: MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"
            ),
        ),
        "yaml_config": TextParameter(
            "YAML Configuration",
            default="",
            docstring=(
                "Advanced: Provide full Boltz YAML configuration. "
                "If provided, overrides the sequence parameter. "
                "See Boltz-2 documentation for YAML format."
            ),
        ),
        "msa_mode": TextParameter(
            "MSA Mode",
            default="empty",
            docstring=(
                "Multiple Sequence Alignment mode:\n"
                "- 'empty': No MSA, fast (~45s) but reduced accuracy\n"
                "- 'server': Use MSA server, slow (10-30min) but best accuracy\n"
                "- 'provided': Use provided A3M file (via msa_a3m parameter)"
            ),
        ),
        "msa_a3m": TextParameter(
            "MSA A3M File (Base64)",
            default="",
            docstring=(
                "Optional: Pre-computed MSA in A3M format (base64 encoded). "
                "Only used when msa_mode='provided'. "
                "Generate using ColabFold or MMseqs2."
            ),
        ),
        "output_folder": FolderParameter(
            "Output Folder",
            default="",
            docstring=(
                "Folder for output files. Leave empty to use the workflow's working directory. "
                "Outputs will be saved as: {prefix}.tar.gz and {prefix}_model_0.cif"
            ),
        ),
        "output_prefix": TextParameter(
            "Output Prefix",
            default="",
            docstring=(
                "Prefix for output filenames. Leave empty to auto-generate.\n"
                "Auto-generated format: boltz_{sequence_start}_{timestamp}\n"
                "Example: 'hemoglobin_alpha' → hemoglobin_alpha.tar.gz, hemoglobin_alpha_model_0.cif"
            ),
        ),
    }

    # Cloud execution metadata for UI display
    CLOUD_CONFIG = {
        "provider": "modal",
        "credential_mode": "bocoflow",  # Mode B
        "api_endpoint": "/api/cloud/nodes/boltz2/execute",
        "requires_login": True,
        "requires_gpu": True,
        "gpu_type": "H100",
        "credits_per_call": 0.50,
        "estimated_duration": "3-10 minutes",
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute by calling the API Gateway (which calls Modal)."""
        log_message("Starting CloudModalBoltz2 execution (Mode B)")

        result = NodeResult()
        result.metadata.update(
            {
                "node_type": "CloudModalBoltz2",
                "execution_time": datetime.now().isoformat(),
                "credential_mode": "bocoflow",  # Mode B
                "gpu": "H100",
            }
        )

        # Get auth token from environment (injected by BF2 worker)
        auth_token = os.environ.get("BOCOFLOW_CLOUD_AUTH_TOKEN")

        if not auth_token:
            result.success = False
            result.message = (
                "Cloud authentication required. Please sign in to use cloud nodes.\n"
                "This node requires BoCoFlow cloud credits (Mode B).\n"
                "Unlike Mode A nodes, you don't need your own Modal account."
            )
            return result.to_json()

        # Get parameters
        sequence = flow_vars["sequence"].get_value()
        yaml_config = flow_vars["yaml_config"].get_value()
        msa_mode = flow_vars["msa_mode"].get_value() or "empty"
        msa_a3m = flow_vars["msa_a3m"].get_value() or ""
        output_folder = flow_vars["output_folder"].get_value() or ""
        output_prefix = flow_vars["output_prefix"].get_value() or ""

        # Get sequence from predecessor if not provided directly
        if not sequence and not yaml_config and predecessor_data:
            pred_data = predecessor_data[0] if predecessor_data else {}
            if isinstance(pred_data, dict):
                # Try common field names for sequence data
                sequence = (
                    pred_data.get("sequence", "")
                    or pred_data.get("protein_sequence", "")
                    or pred_data.get("fasta_sequence", "")
                    or pred_data.get("output_text", "")
                )

        # Validate input
        if not sequence and not yaml_config:
            result.success = False
            result.message = (
                "No protein sequence provided. "
                "Please provide a sequence in the 'Protein Sequence' field "
                "or connect a predecessor node that outputs a sequence."
            )
            return result.to_json()

        # Prepare request payload
        payload = {
            "node_info": {
                "node_id": getattr(self, "node_id", "unknown"),
                "node_type": "CloudModalBoltz2",
            },
            "predecessor_data": {
                "sequence": sequence,
                "yaml_config": yaml_config,
            },
            "options": {
                "sequence": sequence,
                "yaml_config": yaml_config,
                "msa_mode": msa_mode,
                "msa_a3m": msa_a3m,
            },
        }

        headers = {
            "Authorization": f"Bearer {auth_token}",  # Firebase token
            "Content-Type": "application/json",
        }

        try:
            log_message(f"Calling API Gateway: {self.API_ENDPOINT}")
            log_message(f"Sequence length: {len(sequence)} amino acids")
            log_message(f"MSA mode: {msa_mode}")

            # Boltz-2 can take several minutes - use longer timeout
            stream_log(
                "Calling Boltz-2 API... First call may take 2-3 min (cold start).",
                node_id=self.node_id,
                progress=10,
            )
            response = post_with_progress(
                url=self.API_ENDPOINT,
                json=payload,
                headers=headers,
                timeout=1800,
                node_id=self.node_id,
                service_name="Boltz-2",
                cold_start_hint="cold starts take 2-3 min",
            )

            if response.status_code == 200:
                cloud_result = response.json()
                modal_result = cloud_result.get("result", {})
                usage_info = cloud_result.get("usage", {})

                # Check if prediction succeeded
                if modal_result.get("status") == "error":
                    result.success = False
                    result.message = (
                        f"Boltz-2 prediction failed: {modal_result.get('error', 'Unknown error')}"
                    )
                    return result.to_json()

                # Save output tarball if present
                output_tarball = modal_result.get("output_tarball_base64", "")
                output_path = None
                output_size = 0
                cif_path = None  # Will be set if CIF is extracted
                file_prefix = None  # Will be set if we generate output files
                final_folder = None  # Will be set if we generate output files

                if output_tarball:
                    # Decode the tarball
                    output_bytes = base64.b64decode(output_tarball)
                    output_size = len(output_bytes)

                    # === Determine output folder ===
                    # Priority: user-specified > workflow dir > Downloads > /tmp
                    if output_folder:
                        # User specified output folder
                        if output_folder.startswith(("abs:", "rel:")):
                            # Handle path prefixes
                            workflow_dir = os.environ.get("BOCOFLOW_WORKFLOW_DIR", "")
                            if output_folder.startswith("abs:"):
                                final_folder = Path(output_folder[4:])
                            elif workflow_dir:
                                final_folder = Path(workflow_dir) / output_folder[4:]
                            else:
                                final_folder = Path(output_folder[4:])
                        else:
                            final_folder = Path(output_folder)
                    else:
                        # No folder specified - use working directory
                        workflow_dir = os.environ.get("BOCOFLOW_WORKFLOW_DIR", "")
                        if workflow_dir:
                            final_folder = Path(workflow_dir)
                        else:
                            # Fallback to user's Downloads folder
                            downloads_dir = Path.home() / "Downloads"
                            if downloads_dir.exists():
                                final_folder = downloads_dir
                            else:
                                # Last resort: /tmp
                                final_folder = Path("/tmp")

                    # === Generate filename prefix ===
                    if output_prefix:
                        # User specified prefix
                        file_prefix = output_prefix
                    else:
                        # Auto-generate: boltz_{seq_start}_{timestamp}
                        seq_short = sequence[:8] if len(sequence) >= 8 else sequence
                        # Clean up sequence (remove non-alphanumeric)
                        seq_short = "".join(c for c in seq_short if c.isalnum())
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                        file_prefix = f"boltz_{seq_short}_{timestamp}"

                    # === Build final paths ===
                    final_folder.mkdir(parents=True, exist_ok=True)
                    output_path = final_folder / f"{file_prefix}.tar.gz"

                    # Write tarball
                    output_path.write_bytes(output_bytes)
                    log_message(f"Saved output to {output_path} ({output_size} bytes)")

                    # Also extract the CIF file for convenience
                    try:
                        import io
                        import tarfile

                        cif_path = None
                        with tarfile.open(fileobj=io.BytesIO(output_bytes), mode="r:gz") as tar:
                            for member in tar.getmembers():
                                if member.name.endswith(".cif"):
                                    # Extract CIF with consistent naming
                                    # Use file_prefix + _model_0.cif (matching Boltz output)
                                    cif_path = final_folder / f"{file_prefix}_model_0.cif"
                                    with tar.extractfile(member) as f:
                                        cif_path.write_bytes(f.read())
                                    log_message(f"Extracted CIF structure to {cif_path}")
                                    break
                    except Exception as e:
                        log_message(f"Warning: Could not extract CIF file: {e}")
                else:
                    log_message("Warning: No output tarball received from Modal")

                # Build result data
                output_files_list = modal_result.get("output_files", [])

                result.success = True
                if output_path and output_path.exists():
                    msg_parts = [
                        f"Boltz-2 structure prediction completed.",
                        f"Output saved to {output_path} ({output_size} bytes).",
                    ]
                    if cif_path and cif_path.exists():
                        msg_parts.append(f"CIF structure: {cif_path}")
                    msg_parts.append(f"Duration: {usage_info.get('duration_seconds', 0):.2f}s")
                    result.message = " ".join(msg_parts)
                elif output_tarball:
                    result.message = (
                        f"Boltz-2 structure prediction completed. "
                        f"Output available ({output_size} bytes, {len(output_files_list)} files). "
                        f"Duration: {usage_info.get('duration_seconds', 0):.2f}s"
                    )
                else:
                    result.message = (
                        f"Boltz-2 structure prediction completed but no output files found. "
                        f"Duration: {usage_info.get('duration_seconds', 0):.2f}s"
                    )

                result.data = {
                    "output_file": str(output_path)
                    if output_path and output_path.exists()
                    else None,
                    "cif_file": str(cif_path) if cif_path and cif_path.exists() else None,
                    "output_folder": str(final_folder) if final_folder else None,
                    "output_prefix": file_prefix,
                    "output_file_size": output_size,
                    "output_files": output_files_list,
                    "output_tarball_available": bool(output_tarball),
                    "sequence_length": modal_result.get("sequence_length", len(sequence)),
                    "processing_time_seconds": modal_result.get("processing_time_seconds", 0),
                    "modal_metadata": modal_result.get("modal_metadata", {}),
                    "job_id": cloud_result.get("job_id"),
                    "usage": {
                        "duration_seconds": usage_info.get("duration_seconds", 0),
                        "cost_usd": usage_info.get("cost_usd", 0),
                    },
                    "status": "completed",
                    "credential_mode": "bocoflow",
                }
                result.metadata["cloud_job_id"] = cloud_result.get("job_id")

            elif response.status_code == 401:
                result.success = False
                result.message = "Authentication failed. Please sign in again."

            elif response.status_code == 402:
                result.success = False
                result.message = "Insufficient credits. Please purchase more credits."

            elif response.status_code == 503:
                result.success = False
                result.message = (
                    "Modal cloud service temporarily unavailable. "
                    "The H100 GPU may be scaling up. Please try again in a few minutes."
                )

            else:
                error_detail = ""
                try:
                    error_detail = response.json().get("detail", response.text)
                except Exception:
                    error_detail = response.text
                result.success = False
                result.message = f"API error ({response.status_code}): {error_detail}"

        except requests.Timeout:
            result.success = False
            result.message = (
                "Request timed out. Boltz-2 predictions can take several minutes. "
                "Please try again or use a shorter sequence."
            )

        except requests.RequestException as e:
            result.success = False
            result.message = f"Network error: {str(e)}"

        except Exception as e:
            log_message(f"Unexpected error in CloudModalBoltz2: {str(e)}")
            result.success = False
            result.message = f"Unexpected error: {str(e)}"

        return result.to_json()
