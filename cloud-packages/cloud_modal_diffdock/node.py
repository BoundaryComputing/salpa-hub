"""
Cloud Modal DiffDock - Mode B (BoCoFlow Credits)

Blind protein-ligand docking using DiffDock on Modal's A10G GPU infrastructure.

This node is a client stub that calls the BoCoFlow API Gateway,
which then routes requests to Modal endpoints with Proxy Auth tokens.

Unlike modal-user nodes (Mode A), users don't need their own Modal account.
Instead, they pay with BoCoFlow credits.

Architecture:
    1. User authenticates with Firebase (token in BOCOFLOW_CLOUD_AUTH_TOKEN)
    2. This node sends protein PDB + ligand to API Gateway with Firebase token
    3. API Gateway verifies token, checks credits
    4. API Gateway calls Modal with Proxy Auth tokens
    5. DiffDock runs on A10G GPU, returns ranked docking poses
    6. Result (tarball with SDF files) returns through this node
    7. Credits are deducted from user's account

IMPORTANT — Cold Start Warning:
    The first call after the GPU container scales to zero takes ~10 minutes
    (loading ESM-2 protein language model, 2.6 GB, into GPU memory).
    Subsequent calls within 2 minutes are "warm" and take only ~40 seconds.

    This node is designed for INTERACTIVE single-protein docking, NOT for
    batch processing. For batch docking (many protein-ligand pairs), use
    a local DiffDock installation or a dedicated HPC pipeline instead.

Reference: Corso et al., ICLR 2023 (MIT License)
"""

import base64
import io
import os
import tarfile
from datetime import datetime

import requests
from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeResult
from bocoflow_core.stream_logger import post_with_progress, stream_log
from bocoflow_core.parameters import (
    FileParameterEdit,
    FolderParameter,
    IntegerParameter,
    StringParameter,
    TextParameter,
)


class CloudModalDiffdock(Node):
    """
    DiffDock Blind Protein-Ligand Docking (BoCoFlow Credits - Mode B).

    This is a CLIENT STUB - actual computation happens on Modal cloud (A10G GPU).
    The workflow engine executes this node, which calls the API Gateway.

    DiffDock is a generative diffusion model that predicts how small molecules
    bind to protein targets — a critical step in drug discovery pipelines.

    Timing:
    - Cold start: ~10 minutes (ESM-2 model loading into GPU memory)
    - Warm start: ~40 seconds (container reused within 2-minute window)
    - NOT suitable for batch processing of many protein-ligand pairs

    Prerequisites:
    - User must be logged in with Firebase
    - User must have sufficient BoCoFlow credits

    No Modal account or 'modal setup' required!
    """

    # NOTE: Metadata (name, hashtags, num_in, num_out) comes from meta.toml.
    # NOTE: EXECUTION_STRATEGY and ENVIRONMENT are auto-detected from pixi.toml.

    # API Gateway endpoint (unified route — not Modal directly!)
    API_ENDPOINT = (
        os.environ.get(
            "BOCOFLOW_CLOUD_API_URL",
            "https://bocoflow-api-gateway-823406908684.us-central1.run.app",
        )
        + "/api/cloud/nodes/modal-diffdock/execute"
    )

    OPTIONS = {
        "protein_pdb_file": FileParameterEdit(
            "Protein PDB",
            docstring="Protein structure file (PDB format) for docking target",
        ),
        "ligand_smiles": StringParameter(
            "Ligand SMILES",
            default="",
            docstring=(
                "SMILES string for the ligand molecule (primary input). "
                "Example: C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1 (Erlotinib)"
            ),
        ),
        "ligand_sdf_file": FileParameterEdit(
            "Ligand SDF",
            docstring="Alternative: Ligand structure file (SDF format). Used if SMILES is empty.",
        ),
        "num_poses": IntegerParameter(
            "Number of Poses",
            default=10,
            docstring="Number of binding poses to generate (1-40)",
        ),
        "inference_steps": IntegerParameter(
            "Inference Steps",
            default=20,
            docstring="Number of denoising steps (10-40, higher = more accurate but slower)",
        ),
        "samples_per_complex": IntegerParameter(
            "Samples per Complex",
            default=10,
            docstring="Number of samples per protein-ligand complex (1-40)",
        ),
        "output_folder": FolderParameter(
            "Output Folder",
            default="",
            docstring=(
                "Folder for output files. Leave empty to use the workflow's working directory. "
                "Outputs: {prefix}.tar.gz and extracted SDF pose files"
            ),
        ),
        "output_prefix": TextParameter(
            "Output Prefix",
            default="",
            docstring=(
                "Prefix for output filenames. Leave empty to auto-generate.\n"
                "Auto-generated format: diffdock_{timestamp}\n"
                "Example: 'aspirin_dock' -> aspirin_dock.tar.gz, aspirin_dock_rank1.sdf"
            ),
        ),
    }

    # Cloud execution metadata for UI display
    CLOUD_CONFIG = {
        "provider": "modal",
        "credential_mode": "bocoflow",  # Mode B
        "api_endpoint": "/api/cloud/nodes/modal-diffdock/execute",
        "requires_login": True,
        "requires_gpu": True,
        "gpu_type": "A10G",
        "credits_per_call": 0.10,
        "estimated_duration": "40-600 seconds (cold start ~10 min, warm ~40s). Not for batch use.",
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute by calling the API Gateway (which calls Modal)."""
        stream_log(
            "Starting DiffDock cloud execution... "
            "First call may take ~10 min (cold start). Warm calls take ~40s.",
            node_id=self.node_id, progress=0,
        )

        result = NodeResult()
        result.metadata.update(
            {
                "node_type": "CloudModalDiffdock",
                "execution_time": datetime.now().isoformat(),
                "credential_mode": "bocoflow",  # Mode B
                "gpu": "A10G",
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

        # -- Read protein PDB content --
        pdb_path = flow_vars["protein_pdb_file"].get_value()
        protein_pdb = ""

        if pdb_path:
            resolved_path = self.resolve_path(pdb_path)
            if resolved_path and os.path.isfile(resolved_path):
                with open(resolved_path, "r") as f:
                    protein_pdb = f.read()

        # Try predecessor data if no PDB from file
        if not protein_pdb and predecessor_data:
            pred_data = predecessor_data[0] if predecessor_data else {}
            if isinstance(pred_data, dict):
                protein_pdb = (
                    pred_data.get("pdb_content", "")
                    or pred_data.get("protein_pdb", "")
                )
                if not protein_pdb:
                    output_file = pred_data.get("output_file", "")
                    if output_file and os.path.isfile(output_file):
                        with open(output_file, "r") as f:
                            protein_pdb = f.read()

        if not protein_pdb:
            result.success = False
            result.message = (
                "No protein PDB content provided. "
                "Please select a PDB file or connect a predecessor node."
            )
            stream_log(f"Error: {result.message}", node_id=self.node_id, level="error")
            return result.to_json()

        # -- Read ligand --
        ligand_smiles = flow_vars["ligand_smiles"].get_value() or ""
        ligand_sdf = ""

        sdf_path = flow_vars["ligand_sdf_file"].get_value()
        if sdf_path and not ligand_smiles:
            resolved_sdf = self.resolve_path(sdf_path)
            if resolved_sdf and os.path.isfile(resolved_sdf):
                with open(resolved_sdf, "r") as f:
                    ligand_sdf = f.read()

        if not ligand_smiles and not ligand_sdf:
            result.success = False
            result.message = (
                "No ligand provided. "
                "Please enter a SMILES string or select an SDF file."
            )
            return result.to_json()

        # Get parameters
        num_poses = flow_vars["num_poses"].get_value()
        inference_steps = flow_vars["inference_steps"].get_value()
        samples_per_complex = flow_vars["samples_per_complex"].get_value()
        output_folder = flow_vars["output_folder"].get_value() or ""
        output_prefix = flow_vars["output_prefix"].get_value() or ""

        # Prepare request payload
        payload = {
            "node_info": {
                "node_id": getattr(self, "node_id", "unknown"),
                "node_type": "CloudModalDiffdock",
            },
            "predecessor_data": {
                "protein_pdb": protein_pdb,
            },
            "options": {
                "protein_pdb": protein_pdb,
                "ligand_smiles": ligand_smiles,
                "ligand_sdf": ligand_sdf,
                "num_poses": num_poses,
                "inference_steps": inference_steps,
                "samples_per_complex": samples_per_complex,
            },
        }

        headers = {
            "Authorization": f"Bearer {auth_token}",  # Firebase token
            "Content-Type": "application/json",
        }

        try:
            stream_log(
                f"Calling DiffDock API (poses={num_poses}, steps={inference_steps})...",
                node_id=self.node_id, progress=10,
            )
            log_message(f"Protein PDB: {len(protein_pdb)} chars, Ligand: {ligand_smiles[:50] if ligand_smiles else 'SDF file'}")

            # DiffDock: warm ~40s, cold start ~590s (ESM-2 model loading)
            response = post_with_progress(
                url=self.API_ENDPOINT,
                json=payload,
                headers=headers,
                timeout=900,
                node_id=self.node_id,
                service_name="DiffDock",
                cold_start_hint="cold starts take up to 10 min",
            )

            if response.status_code == 200:
                cloud_result = response.json()
                modal_result = cloud_result.get("result", {})
                usage_info = cloud_result.get("usage", {})

                # Check if prediction succeeded
                if modal_result.get("status") == "error":
                    result.success = False
                    result.message = (
                        f"DiffDock prediction failed: {modal_result.get('error', 'Unknown error')}"
                    )
                    return result.to_json()

                # Extract result data
                output_tarball = modal_result.get("output_tarball_base64", "")
                actual_poses = modal_result.get("num_poses", 0)
                top_confidence = modal_result.get("top_confidence", 0.0)
                confidence_scores = modal_result.get("confidence_scores", [])
                processing_time = modal_result.get("processing_time_seconds", 0)

                output_path = None
                output_size = 0
                extracted_sdfs = []
                final_folder = None
                file_prefix = None

                if output_tarball:
                    # Decode the tarball
                    output_bytes = base64.b64decode(output_tarball)
                    output_size = len(output_bytes)

                    # Determine output folder
                    # Use working_path (set by worker) as the base for relative paths
                    working_dir = getattr(self, "working_path", "") or os.environ.get("BOCOFLOW_WORKFLOW_DIR", "")
                    if output_folder:
                        if output_folder.startswith("abs:"):
                            final_folder = output_folder[4:]
                        elif output_folder.startswith("rel:"):
                            rel_part = output_folder[4:]
                            if working_dir:
                                final_folder = os.path.join(working_dir, rel_part)
                            else:
                                final_folder = rel_part
                        else:
                            final_folder = output_folder
                    else:
                        if working_dir:
                            final_folder = working_dir
                        else:
                            downloads_dir = os.path.join(os.path.expanduser("~"), "Downloads")
                            if os.path.exists(downloads_dir):
                                final_folder = downloads_dir
                            else:
                                final_folder = "/tmp"

                    # Generate filename prefix
                    if output_prefix:
                        file_prefix = output_prefix
                    else:
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                        file_prefix = f"diffdock_{timestamp}"

                    # Save tarball
                    os.makedirs(final_folder, exist_ok=True)
                    output_path = os.path.join(final_folder, f"{file_prefix}.tar.gz")
                    with open(output_path, "wb") as f:
                        f.write(output_bytes)
                    log_message(f"Saved output to {output_path} ({output_size} bytes)")

                    # Extract SDF files for convenience
                    try:
                        with tarfile.open(fileobj=io.BytesIO(output_bytes), mode="r:gz") as tar:
                            for member in tar.getmembers():
                                if member.name.endswith(".sdf") and member.isfile():
                                    sdf_name = os.path.basename(member.name)
                                    sdf_dest = os.path.join(final_folder, f"{file_prefix}_{sdf_name}")
                                    with tar.extractfile(member) as src:
                                        with open(sdf_dest, "wb") as dst:
                                            dst.write(src.read())
                                    extracted_sdfs.append(sdf_dest)
                                    log_message(f"Extracted SDF: {sdf_dest}")
                    except Exception as e:
                        log_message(f"Warning: Could not extract SDF files: {e}")
                else:
                    log_message("Warning: No output tarball received from Modal")

                # Build result data
                stream_log(
                    f"DiffDock completed: {actual_poses} poses, top confidence={top_confidence:.4f}",
                    node_id=self.node_id, progress=90,
                )
                result.success = True
                if output_path and os.path.exists(output_path):
                    msg_parts = [
                        f"DiffDock generated {actual_poses} docking poses.",
                        f"Top confidence: {top_confidence:.4f}.",
                        f"Output: {output_path} ({output_size} bytes).",
                    ]
                    if extracted_sdfs:
                        msg_parts.append(f"Extracted {len(extracted_sdfs)} SDF files.")
                    msg_parts.append(f"Duration: {usage_info.get('duration_seconds', 0):.2f}s")
                    result.message = " ".join(msg_parts)
                else:
                    result.message = (
                        f"DiffDock generated {actual_poses} docking poses. "
                        f"Top confidence: {top_confidence:.4f}. "
                        f"Duration: {usage_info.get('duration_seconds', 0):.2f}s"
                    )

                result.data = {
                    "output_file": output_path if output_path and os.path.exists(output_path) else None,
                    "output_folder": final_folder,
                    "output_prefix": file_prefix,
                    "extracted_sdfs": extracted_sdfs,
                    "output_file_size": output_size,
                    "output_files": modal_result.get("output_files", []),
                    "output_tarball_available": bool(output_tarball),
                    "num_poses": actual_poses,
                    "top_confidence": top_confidence,
                    "confidence_scores": confidence_scores,
                    "processing_time_seconds": processing_time,
                    "modal_metadata": modal_result.get("modal_metadata", {}),
                    "job_id": cloud_result.get("job_id"),
                    "protein_pdb": protein_pdb,
                    "ligand_smiles": ligand_smiles,
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
                stream_log(f"Error: {result.message}", node_id=self.node_id, level="error")

            elif response.status_code == 402:
                result.success = False
                result.message = "Insufficient credits. Please purchase more credits."
                stream_log(f"Error: {result.message}", node_id=self.node_id, level="error")

            elif response.status_code == 503:
                result.success = False
                result.message = (
                    "Modal cloud service temporarily unavailable. "
                    "The A10G GPU may be scaling up. Please try again in a few minutes."
                )
                stream_log(f"Error: {result.message}", node_id=self.node_id, level="error")

            else:
                error_detail = ""
                try:
                    error_detail = response.json().get("detail", response.text)
                except Exception:
                    error_detail = response.text
                result.success = False
                result.message = f"API error ({response.status_code}): {error_detail}"
                stream_log(f"Error: {result.message}", node_id=self.node_id, level="error")

        except requests.Timeout:
            result.success = False
            result.message = (
                "Request timed out (15 min limit). "
                "DiffDock cold starts take ~10 min (ESM-2 model loading). "
                "Warm calls take ~40s. Please try again — the container may now be warm."
            )

        except requests.RequestException as e:
            result.success = False
            result.message = f"Network error: {str(e)}"

        except Exception as e:
            log_message(f"Unexpected error in CloudModalDiffdock: {str(e)}")
            result.success = False
            result.message = f"Unexpected error: {str(e)}"

        return result.to_json()
