"""
Cloud Modal Chai-1 - Mode B (BoCoFlow Credits)

Multi-modal structure prediction using Chai-1 on Modal's H100 GPU infrastructure.

This node is a client stub that calls the BoCoFlow API Gateway,
which then routes requests to Modal endpoints with Proxy Auth tokens.

Unlike modal-user nodes (Mode A), users don't need their own Modal account.
Instead, they pay with BoCoFlow credits.

Chai-1 can predict structures of proteins, nucleic acids, small molecules,
and their complexes from FASTA-format input with entity type annotations.

Architecture:
    1. User authenticates with Firebase (token in BOCOFLOW_CLOUD_AUTH_TOKEN)
    2. This node sends FASTA/sequence to API Gateway with Firebase token
    3. API Gateway verifies token, checks credits
    4. API Gateway calls Modal with Proxy Auth tokens
    5. Chai-1 runs on H100 GPU, returns structure predictions
    6. Result (tarball with CIF files) returns through this node
    7. Credits are deducted from user's account

Reference: https://github.com/chaidiscovery/chai-lab
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
    IntegerParameter,
    TextParameter,
)


#: How to run this node on its own -- the values `salpa smoke` feeds it. Strings
#: starting with `demo_data/` resolve relative to this directory. Running needs a
#: Salpa account with cloud access; without one the node stops at authentication,
#: which is what `salpa smoke` will report. See demo_data/README.md.
DEMO_CONFIG = {
    "fasta_input": '>protein|name=trp-cage\nNLYIQWLKDGGPSSGRPPPS',
}


class CloudModalChai1(Node):
    """
    Chai-1 Multi-Modal Structure Prediction (BoCoFlow Credits - Mode B).

    This is a CLIENT STUB - actual computation happens on Modal cloud (H100 GPU).
    The workflow engine executes this node, which calls the API Gateway.

    Chai-1 can predict 3D structures of:
    - Proteins (from amino acid sequences)
    - Protein-ligand complexes (protein + SMILES)
    - Protein-nucleic acid complexes (protein + DNA/RNA)
    - Multi-chain complexes

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
        + "/api/cloud/nodes/chai1/execute"
    )

    OPTIONS = {
        "fasta_input": TextParameter(
            "FASTA Input",
            default="",
            docstring=(
                "FASTA-format input with entity type annotations.\n"
                "Supports protein, ligand (SMILES), DNA, and RNA entities.\n"
                "If provided, overrides protein_sequence and ligand_smiles.\n\n"
                "Example:\n"
                ">protein|name=hemoglobin\n"
                "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH\n"
                ">ligand|name=aspirin\n"
                "CC(=O)OC1=CC=CC=C1C(=O)O"
            ),
        ),
        "protein_sequence": TextParameter(
            "Protein Sequence",
            default="",
            docstring=(
                "Simple protein amino acid sequence (one-letter codes). "
                "Auto-wrapped to FASTA format if fasta_input is not provided. "
                "Can also come from predecessor node. "
                "Example: MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"
            ),
        ),
        "ligand_smiles": TextParameter(
            "Ligand SMILES",
            default="",
            docstring=(
                "Optional ligand SMILES string for protein-ligand complex prediction. "
                "Only used when protein_sequence is provided (not fasta_input). "
                "Example: CC(=O)OC1=CC=CC=C1C(=O)O (aspirin)"
            ),
        ),
        "num_trunk_recycles": IntegerParameter(
            "Trunk Recycles",
            default=3,
            docstring=(
                "Number of trunk recycles for structure refinement. "
                "More recycles improve accuracy but increase runtime. Default: 3."
            ),
        ),
        "num_diffn_timesteps": IntegerParameter(
            "Diffusion Timesteps",
            default=200,
            docstring=(
                "Number of diffusion timesteps for structure generation. "
                "More timesteps improve accuracy but increase runtime. Default: 200."
            ),
        ),
        "seed": IntegerParameter(
            "Random Seed",
            default=42,
            docstring="Random seed for reproducibility. Default: 42.",
        ),
        "use_esm_embeddings": BooleanParameter(
            "Use ESM Embeddings",
            default=True,
            docstring=(
                "Use ESM protein language model embeddings for improved accuracy. "
                "Recommended for protein structure prediction. Default: True."
            ),
        ),
        "output_folder": FolderParameter(
            "Output Folder",
            default="",
            docstring=(
                "Folder for output files. Leave empty to use the workflow's working directory. "
                "Outputs: {prefix}.tar.gz (all samples) and {prefix}_best.cif (best structure)"
            ),
        ),
        "output_prefix": TextParameter(
            "Output Prefix",
            default="",
            docstring=(
                "Prefix for output filenames. Leave empty to auto-generate.\n"
                "Auto-generated format: chai1_{sequence_start}_{timestamp}\n"
                "Example: 'hemoglobin' -> hemoglobin.tar.gz, hemoglobin_best.cif"
            ),
        ),
    }

    # Cloud execution metadata for UI display
    CLOUD_CONFIG = {
        "provider": "modal",
        "credential_mode": "bocoflow",  # Mode B
        "api_endpoint": "/api/cloud/nodes/chai1/execute",
        "requires_login": True,
        "requires_gpu": True,
        "gpu_type": "H100",
        "credits_per_call": 0.50,
        "estimated_duration": "3-15 minutes",
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute by calling the API Gateway (which calls Modal)."""
        log_message("Starting CloudModalChai1 execution (Mode B)")

        result = NodeResult()
        result.metadata.update(
            {
                "node_type": "CloudModalChai1",
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
        fasta_input = flow_vars["fasta_input"].get_value()
        protein_sequence = flow_vars["protein_sequence"].get_value()
        ligand_smiles = flow_vars["ligand_smiles"].get_value()
        num_trunk_recycles = flow_vars["num_trunk_recycles"].get_value()
        num_diffn_timesteps = flow_vars["num_diffn_timesteps"].get_value()
        seed = flow_vars["seed"].get_value()
        use_esm_embeddings = flow_vars["use_esm_embeddings"].get_value()
        output_folder = flow_vars["output_folder"].get_value() or ""
        output_prefix = flow_vars["output_prefix"].get_value() or ""

        # Get sequence from predecessor if not provided directly
        if not fasta_input and not protein_sequence and predecessor_data:
            pred_data = predecessor_data[0] if predecessor_data else {}
            if isinstance(pred_data, dict):
                fasta_input = pred_data.get("fasta_input", "")
                if not fasta_input:
                    protein_sequence = (
                        pred_data.get("sequence", "")
                        or pred_data.get("protein_sequence", "")
                        or pred_data.get("fasta_sequence", "")
                        or pred_data.get("output_text", "")
                    )

        # Validate input
        if not fasta_input and not protein_sequence:
            result.success = False
            result.message = (
                "No input provided. Please provide either:\n"
                "- FASTA Input: Multi-entity FASTA with annotations\n"
                "- Protein Sequence: Simple amino acid sequence\n"
                "Or connect a predecessor node that outputs a sequence."
            )
            return result.to_json()

        # Prepare request payload
        payload = {
            "node_info": {
                "node_id": getattr(self, "node_id", "unknown"),
                "node_type": "CloudModalChai1",
            },
            "predecessor_data": {
                "fasta_input": fasta_input,
                "protein_sequence": protein_sequence,
                "ligand_smiles": ligand_smiles,
            },
            "options": {
                "fasta_input": fasta_input,
                "protein_sequence": protein_sequence,
                "ligand_smiles": ligand_smiles,
                "num_trunk_recycles": num_trunk_recycles,
                "num_diffn_timesteps": num_diffn_timesteps,
                "seed": seed,
                "use_esm_embeddings": use_esm_embeddings,
            },
        }

        headers = {
            "Authorization": f"Bearer {auth_token}",  # Firebase token
            "Content-Type": "application/json",
        }

        try:
            log_message(f"Calling API Gateway: {self.API_ENDPOINT}")
            if fasta_input:
                log_message(f"FASTA input length: {len(fasta_input)} chars")
            else:
                log_message(f"Protein sequence length: {len(protein_sequence)} amino acids")
                if ligand_smiles:
                    log_message(f"Ligand SMILES: {ligand_smiles[:50]}...")

            # Chai-1 can take several minutes - use longer timeout
            stream_log(
                "Calling Chai-1 API... First call may take 2-3 min (cold start).",
                node_id=self.node_id,
                progress=10,
            )
            response = post_with_progress(
                url=self.API_ENDPOINT,
                json=payload,
                headers=headers,
                timeout=1800,
                node_id=self.node_id,
                service_name="Chai-1",
                cold_start_hint="cold starts take 2-3 min",
            )

            stream_log("Received response from cloud", node_id=self.node_id, progress=50)

            if response.status_code == 200:
                cloud_result = response.json()
                modal_result = cloud_result.get("result", {})
                usage_info = cloud_result.get("usage", {})

                # Check if prediction succeeded
                if modal_result.get("status") == "error":
                    result.success = False
                    result.message = (
                        f"Chai-1 prediction failed: {modal_result.get('error', 'Unknown error')}"
                    )
                    return result.to_json()

                # Save output tarball if present
                output_tarball = modal_result.get("output_tarball_base64", "")
                stream_log("Processing output files...", node_id=self.node_id, progress=60)
                output_path = None
                output_size = 0
                cif_path = None
                file_prefix = None
                final_folder = None

                if output_tarball:
                    output_bytes = base64.b64decode(output_tarball)
                    output_size = len(output_bytes)

                    # === Determine output folder ===
                    if output_folder:
                        if output_folder.startswith(("abs:", "rel:")):
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
                        workflow_dir = os.environ.get("BOCOFLOW_WORKFLOW_DIR", "")
                        if workflow_dir:
                            final_folder = Path(workflow_dir)
                        else:
                            downloads_dir = Path.home() / "Downloads"
                            if downloads_dir.exists():
                                final_folder = downloads_dir
                            else:
                                final_folder = Path("/tmp")

                    # === Generate filename prefix ===
                    if output_prefix:
                        file_prefix = output_prefix
                    else:
                        # Auto-generate from sequence or fasta
                        if protein_sequence:
                            seq_short = protein_sequence[:8]
                        elif fasta_input:
                            # Extract first sequence from FASTA
                            lines = fasta_input.strip().split("\n")
                            seq_short = ""
                            for line in lines:
                                if not line.startswith(">"):
                                    seq_short = line[:8]
                                    break
                        else:
                            seq_short = "unknown"
                        seq_short = "".join(c for c in seq_short if c.isalnum())
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                        file_prefix = f"chai1_{seq_short}_{timestamp}"

                    # === Build final paths ===
                    final_folder.mkdir(parents=True, exist_ok=True)
                    output_path = final_folder / f"{file_prefix}.tar.gz"

                    # Write tarball
                    output_path.write_bytes(output_bytes)
                    log_message(f"Saved output to {output_path} ({output_size} bytes)")
                    stream_log(f"Output saved: {output_path.name} ({output_size} bytes)", node_id=self.node_id, progress=75)

                    # Extract the best CIF file for convenience
                    best_idx = modal_result.get("best_sample_idx", 0)
                    try:
                        import io
                        import tarfile

                        with tarfile.open(fileobj=io.BytesIO(output_bytes), mode="r:gz") as tar:
                            cif_members = [m for m in tar.getmembers() if m.name.endswith(".cif")]
                            if cif_members:
                                # Try to find the best sample CIF, otherwise use first
                                target_member = cif_members[0]
                                for m in cif_members:
                                    if f"pred.model_idx_{best_idx}" in m.name:
                                        target_member = m
                                        break

                                cif_path = final_folder / f"{file_prefix}_best.cif"
                                with tar.extractfile(target_member) as f:
                                    cif_path.write_bytes(f.read())
                                log_message(f"Extracted best CIF structure to {cif_path}")
                    except Exception as e:
                        log_message(f"Warning: Could not extract CIF file: {e}")
                else:
                    log_message("Warning: No output tarball received from Modal")

                # Build result data
                output_files_list = modal_result.get("output_files", [])
                num_samples = modal_result.get("num_samples", 0)
                best_score = modal_result.get("best_aggregate_score", 0.0)
                scores = modal_result.get("scores", [])

                result.success = True
                if output_path and output_path.exists():
                    msg_parts = [
                        f"Chai-1 structure prediction completed.",
                        f"{num_samples} samples generated.",
                        f"Best score: {best_score:.4f} (sample {best_idx}).",
                        f"Output saved to {output_path} ({output_size} bytes).",
                    ]
                    if cif_path and cif_path.exists():
                        msg_parts.append(f"Best CIF: {cif_path}")
                    msg_parts.append(f"Duration: {usage_info.get('duration_seconds', 0):.2f}s")
                    result.message = " ".join(msg_parts)
                elif output_tarball:
                    result.message = (
                        f"Chai-1 prediction completed. "
                        f"{num_samples} samples, best score: {best_score:.4f}. "
                        f"Output available ({output_size} bytes, {len(output_files_list)} files). "
                        f"Duration: {usage_info.get('duration_seconds', 0):.2f}s"
                    )
                else:
                    result.message = (
                        f"Chai-1 prediction completed but no output files found. "
                        f"Duration: {usage_info.get('duration_seconds', 0):.2f}s"
                    )

                result.data = {
                    "output_file": str(output_path) if output_path and output_path.exists() else None,
                    "cif_file": str(cif_path) if cif_path and cif_path.exists() else None,
                    "output_folder": str(final_folder) if final_folder else None,
                    "output_prefix": file_prefix,
                    "output_file_size": output_size,
                    "output_files": output_files_list,
                    "output_tarball_available": bool(output_tarball),
                    "num_samples": num_samples,
                    "best_sample_idx": best_idx,
                    "best_aggregate_score": best_score,
                    "scores": scores,
                    "fasta_input_length": modal_result.get("fasta_input_length", 0),
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
                "Request timed out. Chai-1 predictions can take several minutes. "
                "Please try again or use a smaller complex."
            )

        except requests.RequestException as e:
            result.success = False
            result.message = f"Network error: {str(e)}"

        except Exception as e:
            log_message(f"Unexpected error in CloudModalChai1: {str(e)}")
            result.success = False
            result.message = f"Unexpected error: {str(e)}"

        return result.to_json()
