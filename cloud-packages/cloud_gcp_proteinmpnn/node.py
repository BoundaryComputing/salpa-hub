"""
cloud-gcp-proteinmpnn — BoCoFlow cloud client stub for GCP Cloud Run.

Designs amino acid sequences from protein backbone structures using ProteinMPNN.
CPU-only cloud service — lightweight model (~7 MB), fast inference (1-30s).

Pattern: Standard Node class with pixi.toml for auto-detected PIXI_SUBPROCESS.
Same architecture as cloud-gcp-hello-world.

Reference: Dauparas et al., Science 378:49-56 (2022)
"""

import os
from datetime import datetime

import requests
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit,
    FloatParameter,
    FolderParameter,
    IntegerParameter,
    SelectParameter,
    StringParameter,
    TextParameter,
)
from bocoflow_core.stream_logger import stream_log


class CloudGcpProteinmpnn(Node):
    """
    Protein sequence design from backbone structure using ProteinMPNN.

    This is a client stub — actual computation happens on GCP Cloud Run (CPU).
    Authentication is injected via BOCOFLOW_CLOUD_AUTH_TOKEN environment variable.

    Input: PDB file with backbone atoms (N, CA, C, O)
    Output: FASTA file with designed sequences + scores
    """

    # NOTE: Metadata (name, hashtags, num_in, num_out) comes from meta.toml.
    # NOTE: force_to_run is inherited from Node.BASE_OPTIONS — do NOT add it here.
    # NOTE: EXECUTION_STRATEGY and ENVIRONMENT are auto-detected via shared_environment in meta.toml.

    # Cloud API endpoint — configurable via environment
    API_ENDPOINT = (
        os.environ.get(
            "BOCOFLOW_CLOUD_API_URL",
            "https://bocoflow-api-gateway-823406908684.us-central1.run.app",
        )
        + "/api/cloud/nodes/proteinmpnn/execute"
    )

    OPTIONS = {
        "pdb_input": FileParameterEdit(
            "Input PDB",
            docstring="Backbone structure to design sequences for (PDB format)",
        ),
        "chains_to_design": TextParameter(
            "Chains to Design",
            default="",
            docstring="Chain IDs to redesign, comma-separated (e.g., 'A,B'). Leave empty for all chains.",
        ),
        "num_sequences": IntegerParameter(
            "Number of Sequences",
            default=8,
            docstring="Number of sequence designs to generate (1-100)",
        ),
        "sampling_temperature": FloatParameter(
            "Sampling Temperature",
            default=0.1,
            docstring="Higher = more diverse sequences (0.0-1.0, recommended: 0.1-0.3)",
        ),
        "model_variant": SelectParameter(
            "Model Variant",
            default="vanilla",
            options=["vanilla", "soluble", "ca_only"],
            docstring="vanilla: general-purpose, soluble: soluble proteins, ca_only: Cα-only structures",
        ),
        "checkpoint": SelectParameter(
            "Backbone Noise Level",
            default="v_48_020",
            options=["v_48_020", "v_48_002", "v_48_010", "v_48_030"],
            docstring="Noise tolerance: 002=crystal structures, 020=default, 030=AlphaFold/RFDiffusion outputs",
        ),
        "fixed_positions": TextParameter(
            "Fixed Positions",
            default="",
            docstring="Residue positions to keep fixed, e.g., 'A:1,2,3;B:10-15'. Leave empty to redesign all.",
        ),
        "omit_amino_acids": TextParameter(
            "Omit Amino Acids",
            default="",
            docstring="Amino acids to exclude from design, e.g., 'C,M' (single-letter codes)",
        ),
        "seed": IntegerParameter(
            "Random Seed",
            default=37,
            docstring="Random seed for reproducibility (0 = random)",
        ),
        "output_folder": FolderParameter(
            "Output Folder",
            docstring="Directory for output FASTA files",
        ),
        "output_prefix": StringParameter(
            "Output Prefix",
            default="",
            docstring="Prefix for output filenames",
        ),
    }

    # Cloud execution metadata for UI display
    CLOUD_CONFIG = {
        "provider": "gcp",
        "credential_mode": "bocoflow",
        "api_endpoint": "/api/cloud/nodes/proteinmpnn/execute",
        "requires_login": True,
        "credits_per_call": 0.01,
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute the cloud API call for ProteinMPNN sequence design."""
        stream_log(
            "Starting ProteinMPNN cloud execution...",
            node_id=self.node_id,
            progress=0,
        )

        try:
            result = NodeResult()
            result.metadata.update(
                {
                    "execution_time": datetime.now().isoformat(),
                    "credential_mode": "bocoflow",
                }
            )

            # ── Auth token ───────────────────────────────────────────────
            auth_token = os.environ.get("BOCOFLOW_CLOUD_AUTH_TOKEN")
            if not auth_token:
                raise NodeException(
                    "cloud-gcp-proteinmpnn",
                    "Cloud authentication required. Please sign in to use cloud nodes.",
                )

            # ── Read PDB content ─────────────────────────────────────────
            pdb_path = flow_vars["pdb_input"].get_value()
            pdb_content = ""

            if pdb_path:
                resolved_path = self.resolve_path(pdb_path)
                if resolved_path and os.path.isfile(resolved_path):
                    with open(resolved_path, "r") as f:
                        pdb_content = f.read()

            # Try predecessor data if no PDB from file
            if not pdb_content and predecessor_data:
                pred_data = predecessor_data[0] if predecessor_data else {}
                if isinstance(pred_data, dict):
                    pdb_content = pred_data.get("pdb_content", "")
                    if not pdb_content:
                        # Try output_file from predecessor
                        output_file = pred_data.get("output_file", "")
                        if output_file and os.path.isfile(output_file):
                            with open(output_file, "r") as f:
                                pdb_content = f.read()

            if not pdb_content:
                raise NodeException(
                    "cloud-gcp-proteinmpnn",
                    "No PDB content provided. Please select an input PDB file or connect a predecessor node.",
                )

            # ── Read parameters ──────────────────────────────────────────
            num_sequences = flow_vars["num_sequences"].get_value()
            sampling_temperature = flow_vars["sampling_temperature"].get_value()
            model_variant = flow_vars["model_variant"].get_value()
            checkpoint = flow_vars["checkpoint"].get_value()
            chains_to_design = flow_vars["chains_to_design"].get_value()
            fixed_positions = flow_vars["fixed_positions"].get_value()
            omit_amino_acids = flow_vars["omit_amino_acids"].get_value()
            seed = flow_vars["seed"].get_value()
            output_folder = flow_vars["output_folder"].get_value()
            output_prefix = flow_vars["output_prefix"].get_value()

            # ── Prepare request payload ──────────────────────────────────
            payload = {
                "node_info": {
                    "node_id": getattr(self, "node_id", "unknown"),
                    "node_type": "CloudGcpProteinmpnn",
                },
                "predecessor_data": {},
                "options": {
                    "pdb_content": pdb_content,
                    "num_sequences": num_sequences,
                    "sampling_temperature": sampling_temperature,
                    "model_variant": model_variant,
                    "checkpoint": checkpoint,
                    "chains_to_design": chains_to_design,
                    "fixed_positions": fixed_positions,
                    "omit_amino_acids": omit_amino_acids,
                    "seed": seed,
                },
            }

            headers = {
                "Authorization": f"Bearer {auth_token}",
                "Content-Type": "application/json",
            }

            stream_log(
                f"Calling ProteinMPNN API ({model_variant}/{checkpoint}, "
                f"{num_sequences} sequences, T={sampling_temperature})...",
                node_id=self.node_id,
                progress=20,
            )

            # ── Call API Gateway ─────────────────────────────────────────
            response = requests.post(
                self.API_ENDPOINT,
                json=payload,
                headers=headers,
                timeout=300,
            )

            # ── Handle response ──────────────────────────────────────────
            if response.status_code == 200:
                cloud_result = response.json()
                cloud_data = cloud_result.get("result", {})
                usage_info = cloud_result.get("usage", {})

                sequences = cloud_data.get("sequences", [])
                native_sequence = cloud_data.get("native_sequence", "")
                processing_time = cloud_data.get("processing_time_seconds", 0)

                stream_log(
                    f"Received {len(sequences)} designed sequences",
                    node_id=self.node_id,
                    progress=60,
                )

                # ── Save FASTA output ────────────────────────────────────
                output_file_path = ""
                if output_folder:
                    resolved_folder = self.resolve_path(output_folder)
                    os.makedirs(resolved_folder, exist_ok=True)

                    prefix = output_prefix or "proteinmpnn"
                    fasta_filename = f"{prefix}_designed.fasta"
                    output_file_path = os.path.join(resolved_folder, fasta_filename)

                    with open(output_file_path, "w") as f:
                        # Write native sequence
                        if native_sequence:
                            f.write(f">native\n{native_sequence}\n")

                        # Write designed sequences
                        for seq_info in sequences:
                            idx = seq_info.get("sample_index", 0)
                            score = seq_info.get("score", 0)
                            recovery = seq_info.get("seq_recovery", 0)
                            f.write(
                                f">design_{idx}, "
                                f"score={score:.4f}, "
                                f"global_score={seq_info.get('global_score', 0):.4f}, "
                                f"seq_recovery={recovery:.4f}\n"
                                f"{seq_info['sequence']}\n"
                            )

                    stream_log(
                        f"Saved FASTA to {fasta_filename}",
                        node_id=self.node_id,
                        progress=80,
                    )

                # ── Build result ─────────────────────────────────────────
                result.success = True
                result.message = (
                    f"ProteinMPNN designed {len(sequences)} sequences. "
                    f"Best score: {sequences[0]['score']:.4f}. "
                    f"Duration: {processing_time:.1f}s"
                    if sequences
                    else "ProteinMPNN completed with no sequences."
                )
                result.data = {
                    "sequences": sequences,
                    "native_sequence": native_sequence,
                    "designed_chains": cloud_data.get("designed_chains", []),
                    "fixed_chains": cloud_data.get("fixed_chains", []),
                    "num_sequences": len(sequences),
                    "model_variant": cloud_data.get("model_variant", model_variant),
                    "model_name": cloud_data.get("model_name", checkpoint),
                    "processing_time_seconds": processing_time,
                    "backend": cloud_data.get("backend", "cpu"),
                    "job_id": cloud_result.get("job_id"),
                    "output_file": output_file_path,
                    "pdb_content": pdb_content,
                    "usage": {
                        "duration_seconds": usage_info.get("duration_seconds", 0),
                        "cost_usd": usage_info.get("cost_usd", 0),
                    },
                }
                result.metadata["cloud_job_id"] = cloud_result.get("job_id")

                stream_log(
                    result.message,
                    node_id=self.node_id,
                    progress=100,
                )
                return result.to_json()

            elif response.status_code == 401:
                raise NodeException(
                    "cloud-gcp-proteinmpnn",
                    "Authentication failed. Please sign in again.",
                )

            elif response.status_code == 402:
                raise NodeException(
                    "cloud-gcp-proteinmpnn",
                    "Insufficient credits. Please purchase more credits.",
                )

            elif response.status_code == 403:
                raise NodeException(
                    "cloud-gcp-proteinmpnn",
                    "Usage quota exceeded. Please upgrade your plan.",
                )

            elif response.status_code == 503:
                raise NodeException(
                    "cloud-gcp-proteinmpnn",
                    "Cloud service temporarily unavailable. Please try again.",
                )

            else:
                error_detail = ""
                try:
                    error_detail = response.json().get("detail", response.text)
                except Exception:
                    error_detail = response.text
                raise NodeException(
                    "cloud-gcp-proteinmpnn",
                    f"Cloud API error ({response.status_code}): {error_detail}",
                )

        except NodeException:
            raise
        except requests.Timeout:
            raise NodeException(
                "cloud-gcp-proteinmpnn",
                "Cloud request timed out. The service may be scaling up — please retry.",
            )
        except requests.RequestException as e:
            raise NodeException("cloud-gcp-proteinmpnn", f"Network error: {e}")
        except Exception as e:
            raise NodeException(
                "cloud-gcp-proteinmpnn", f"Unexpected error: {e}"
            )
