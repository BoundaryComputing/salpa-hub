"""
cloud-gcp-evo2 — a Salpa Compute node running on GCP Cloud Run.

DNA foundation model (Evo2 7B) for sequence scoring, embedding extraction,
and de novo DNA generation. GPU cloud service (L4 GPU).

Pattern: Standard Node class with shared_environment for auto-detected PIXI_SUBPROCESS.
The gateway address is read when the node runs (_salpa_compute.api_base).

IMPORTANT — Cold Start Warning:
    The first call after the GPU container scales to zero may take several minutes
    (loading 7B model into GPU memory). With weights pre-baked, cold starts ~2-5 min.
    Warm starts take ~10-30 seconds.

Reference: Arc Institute, 2025 (Apache 2.0)
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import requests
from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FileParameterEdit,
    FolderParameter,
    IntegerParameter,
    SelectParameter,
    StringParameter,
    TextParameter,
)
from bocoflow_core.stream_logger import post_with_progress, stream_log

try:  # loaded as a package: node_runner, the server, tests
    from . import _salpa_compute as sc
except ImportError:
    try:  # the node's folder is on sys.path
        import _salpa_compute as sc
    except ImportError:  # anything else: load the file that sits beside this one
        import importlib.util

        _spec = importlib.util.spec_from_file_location(
            "_salpa_compute_cloud_gcp_evo2",
            str(Path(__file__).with_name("_salpa_compute.py")),
        )
        sc = importlib.util.module_from_spec(_spec)
        sys.modules[_spec.name] = sc
        _spec.loader.exec_module(sc)

#: Sent with every request, so the gateway records which version of this node called it.
PACKAGE_NAME = "cloud-gcp-evo2"
PACKAGE_VERSION = "1.0.4"
SERVICE = "evo2"
#: The timeout a run needs: the gateway's 15-minute limit for Evo2 and a few minutes more.
RECOMMENDED_TIMEOUT = 1200
#: The gateway stops waiting for Evo2 at 900 s; wait a minute more for its answer.
SERVICE_MAX_WAIT = 960


#: How to run this node on its own -- the values `salpa smoke` feeds it. Strings
#: starting with `demo_data/` resolve relative to this directory. Running needs a
#: Salpa account with cloud access; without one the node stops at authentication,
#: which is what `salpa smoke` will report. See demo_data/README.md.
DEMO_CONFIG = {
    "sequence_file": "demo_data/lac_operator.fasta",
    "mode": "score",
}


class CloudGcpEvo2(Node):
    """
    DNA foundation model using Evo2 7B on Salpa Compute (GCP Cloud Run, L4 GPU).

    This is a client stub — actual computation happens in the cloud.
    Authentication is injected via BOCOFLOW_CLOUD_AUTH_TOKEN environment variable.

    Modes:
    - score: Compute per-position log-likelihoods for a DNA sequence
    - embed: Extract embeddings from a model layer (mean-pooled vector)
    - generate: Generate new DNA from a prompt sequence

    Input: DNA sequence (ACGT characters)
    Output: Log-likelihoods, embeddings, or generated sequence
    """

    OPTIONS = {
        "sequence": TextParameter(
            "DNA Sequence",
            default="",
            docstring=(
                "DNA sequence (ACGT characters only). "
                "For scoring: the full sequence to evaluate. "
                "For generation: the prompt sequence to extend."
            ),
        ),
        "sequence_file": FileParameterEdit(
            "Sequence File",
            docstring="Alternative: FASTA or plain text file containing DNA sequence",
        ),
        "mode": SelectParameter(
            "Mode",
            default="score",
            options=["score", "embed", "generate"],
            docstring=(
                "score: per-position log-likelihoods; "
                "embed: layer embeddings (mean-pooled); "
                "generate: extend prompt with new DNA"
            ),
        ),
        "n_tokens": IntegerParameter(
            "Generate Length",
            default=256,
            docstring="Number of tokens to generate (generate mode only, 1-4096)",
        ),
        "temperature": StringParameter(
            "Temperature",
            default="1.0",
            docstring="Generation temperature (generate mode, default 1.0)",
        ),
        "top_k": IntegerParameter(
            "Top-K",
            default=4,
            docstring="Top-k sampling (generate mode, default 4)",
        ),
        "output_folder": FolderParameter(
            "Output Folder",
            docstring=(
                "Directory for output files (embeddings JSON, generated sequences). "
                "Leave empty to use the workflow's folder; a relative folder is placed "
                "inside it."
            ),
        ),
        "output_prefix": StringParameter(
            "Output Prefix",
            default="",
            docstring="Prefix for output filenames (auto-generated if empty)",
        ),
    }

    # Cloud execution metadata for UI display
    CLOUD_CONFIG = {
        "provider": "gcp",
        "credential_mode": "bocoflow",
        "api_endpoint": "/api/cloud/nodes/evo2/execute",
        "requires_login": True,
        "requires_gpu": True,
        "gpu_type": "L4",
        "credits_per_call": 0.10,
        "estimated_duration": "10-60 seconds (cold start ~2-5 min, warm ~10s)",
    }

    def execute(self, predecessor_data, flow_vars):
        """Execute the cloud API call for Evo2."""
        stream_log(
            "Starting Evo2 on Salpa Compute (GCP L4)... "
            "First call may take ~2-5 min (cold start).",
            node_id=self.node_id,
            progress=0,
        )
        deadline = None

        try:
            result = NodeResult()
            result.metadata.update(
                {
                    "execution_time": datetime.now().isoformat(),
                    "credential_mode": "bocoflow",
                    "gpu": "L4",
                }
            )

            # -- Auth token --
            auth_token = os.environ.get("BOCOFLOW_CLOUD_AUTH_TOKEN")
            if not auth_token:
                raise NodeException("cloud-gcp-evo2", sc.SIGN_IN_MESSAGE)

            # -- Read DNA sequence --
            sequence = flow_vars["sequence"].get_value() or ""

            # Try file input
            if not sequence:
                seq_path = flow_vars["sequence_file"].get_value()
                if seq_path:
                    resolved_path = self.resolve_path(seq_path)
                    if resolved_path and os.path.isfile(resolved_path):
                        with open(resolved_path, "r") as f:
                            content = f.read().strip()
                        # Strip FASTA header if present
                        lines = content.split("\n")
                        seq_lines = [l for l in lines if not l.startswith(">")]
                        sequence = "".join(seq_lines).replace(" ", "").replace("\n", "")

            # Try predecessor data
            if not sequence and predecessor_data:
                pred_data = predecessor_data[0] if predecessor_data else {}
                if isinstance(pred_data, dict):
                    sequence = (
                        pred_data.get("sequence", "")
                        or pred_data.get("dna_sequence", "")
                        or pred_data.get("generated_sequence", "")
                    )

            if not sequence:
                raise NodeException(
                    "cloud-gcp-evo2",
                    "No DNA sequence provided. Enter a sequence, select a file, or connect a predecessor node.",
                )

            # -- Read parameters --
            mode = flow_vars["mode"].get_value() or "score"
            n_tokens = flow_vars["n_tokens"].get_value()
            temperature = float(flow_vars["temperature"].get_value() or "1.0")
            top_k = flow_vars["top_k"].get_value()
            output_folder = flow_vars["output_folder"].get_value()
            output_prefix = flow_vars["output_prefix"].get_value()

            deadline = sc.deadline(flow_vars, RECOMMENDED_TIMEOUT, SERVICE_MAX_WAIT)
            if deadline.warning:
                stream_log(deadline.warning, node_id=self.node_id, level="warning")

            # -- Prepare request payload --
            payload = {
                "node_info": {
                    "node_id": getattr(self, "node_id", "unknown"),
                    "node_type": "CloudGcpEvo2",
                    "package": PACKAGE_NAME,
                    "package_version": PACKAGE_VERSION,
                },
                "predecessor_data": {},
                "options": {
                    "sequence": sequence,
                    "mode": mode,
                    "n_tokens": n_tokens,
                    "temperature": temperature,
                    "top_k": top_k,
                },
                "client_request_id": sc.new_client_request_id(),
            }

            headers = {
                "Authorization": f"Bearer {auth_token}",
                "Content-Type": "application/json",
            }

            stream_log(
                f"Calling Evo2 on Salpa Compute (mode={mode}, seq_len={len(sequence)})...",
                node_id=self.node_id,
                progress=20,
            )

            # -- Call the gateway --
            response = post_with_progress(
                url=sc.execute_url(SERVICE),
                json=payload,
                headers=headers,
                timeout=deadline.post_timeout,
                node_id=self.node_id,
                service_name="Evo2",
                cold_start_hint="cold starts take 2-5 min",
            )

            # -- Handle response --
            if response.status_code != 200:
                raise NodeException("cloud-gcp-evo2", sc.http_error_message(response, "Evo2"))

            cloud_result = response.json()
            job_id = cloud_result.get("job_id")
            failure = sc.failure_message(cloud_result)
            if failure:
                raise NodeException(
                    "cloud-gcp-evo2",
                    f"Evo2 failed: {failure}" + (f" (job {job_id})" if job_id else ""),
                )

            cloud_data = cloud_result.get("result")
            if not isinstance(cloud_data, dict) or not cloud_data:
                # An empty result is never a success: older gateways relayed a failed
                # service run this way.
                raise NodeException(
                    "cloud-gcp-evo2",
                    "Evo2 returned no result" + (f" (job {job_id})." if job_id else "."),
                )
            usage_info = cloud_result.get("usage") or {}
            processing_time = cloud_data.get("processing_time_seconds", 0) or 0

            stream_log(
                f"Evo2 {mode} completed ({processing_time:.1f}s)",
                node_id=self.node_id,
                progress=60,
            )

            # -- Save output files --
            resolved_folder, folder_note = sc.resolve_output_dir(self, output_folder)
            if folder_note:
                stream_log(folder_note, node_id=self.node_id, level="warning")
            os.makedirs(resolved_folder, exist_ok=True)
            prefix = output_prefix or f"evo2_{mode}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            prefix = prefix.replace("/", "_").replace("\\", "_")

            # Save results as JSON
            output_file_path = os.path.join(str(resolved_folder), f"{prefix}.json")
            with open(output_file_path, "w") as f:
                json.dump(cloud_data, f, indent=2)
            result.files["output"]["results_json"] = self.format_output_path(output_file_path)

            # For generate mode, also save the sequence as FASTA
            if mode == "generate" and cloud_data.get("generated_sequence"):
                fasta_path = os.path.join(str(resolved_folder), f"{prefix}.fasta")
                with open(fasta_path, "w") as f:
                    f.write(f">evo2_generated len={cloud_data.get('total_length', 0)}\n")
                    seq = cloud_data["generated_sequence"]
                    # Wrap at 80 chars
                    for i in range(0, len(seq), 80):
                        f.write(seq[i : i + 80] + "\n")
                result.files["output"]["generated_fasta"] = self.format_output_path(fasta_path)

            stream_log(
                f"Saved results to {prefix}.json",
                node_id=self.node_id,
                progress=90,
            )

            # -- Build result --
            result.success = True

            if mode == "score":
                mean_ll = cloud_data.get("mean_log_likelihood", 0)
                result.message = (
                    f"Evo2 scored {len(sequence)} bp sequence. "
                    f"Mean log-likelihood: {mean_ll:.4f}. "
                    f"Duration: {processing_time:.1f}s"
                )
            elif mode == "embed":
                emb_dim = cloud_data.get("embedding_dim", 0)
                result.message = (
                    f"Evo2 extracted {emb_dim}-dim embedding for {len(sequence)} bp sequence. "
                    f"Duration: {processing_time:.1f}s"
                )
            elif mode == "generate":
                gen_len = cloud_data.get("generated_length", 0)
                result.message = (
                    f"Evo2 generated {gen_len} bp from {len(sequence)} bp prompt. "
                    f"Duration: {processing_time:.1f}s"
                )

            result.data = {
                **cloud_data,
                "output_file": output_file_path,
                "sequence_input": sequence[:100] + ("..." if len(sequence) > 100 else ""),
                "job_id": job_id,
                "usage": {
                    "duration_seconds": usage_info.get("duration_seconds", 0),
                    "cost_usd": usage_info.get("cost_usd", 0),
                },
            }
            result.metadata["cloud_job_id"] = job_id

            stream_log(
                result.message,
                node_id=self.node_id,
                progress=100,
            )
            return result.to_json()

        except NodeException:
            raise
        except requests.Timeout:
            raise NodeException(
                "cloud-gcp-evo2",
                "No answer from Salpa Compute in time. "
                + (
                    (deadline.warning if deadline else None)
                    or "Cold starts may take ~2-5 min. Please try again."
                ),
            )
        except requests.RequestException as e:
            raise NodeException("cloud-gcp-evo2", f"Network error: {e}")
        except Exception as e:
            raise NodeException("cloud-gcp-evo2", f"Unexpected error: {e}")
