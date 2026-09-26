"""
Cloud Modal Chai-1 — a Salpa Compute node.

Multi-modal structure prediction using Chai-1 on an H100 GPU.

This node is a client stub. It sends the input to the Salpa Compute gateway, which checks
the user's sign-in and quota and runs Chai-1 on Modal. No Modal account is needed.

Chai-1 can predict structures of proteins, nucleic acids, small molecules,
and their complexes from FASTA-format input with entity type annotations.

How the result comes back:
    - A result that fits in the reply (up to 16 MiB compressed, nearly every run) arrives
      whole.
    - A larger one arrives as its structure samples and scores plus a link to the full
      archive. The node downloads the archive into the workflow folder, checks its size
      and SHA-256, and then has the server copy deleted. Salpa Compute never keeps it
      longer than 24 hours.

Reference: https://github.com/chaidiscovery/chai-lab
"""

import os
import sys
import tarfile
from datetime import datetime
from pathlib import Path

import requests
from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeResult
from bocoflow_core.parameters import (
    BooleanParameter,
    FolderParameter,
    IntegerParameter,
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
            "_salpa_compute_cloud_modal_chai1",
            str(Path(__file__).with_name("_salpa_compute.py")),
        )
        sc = importlib.util.module_from_spec(_spec)
        sys.modules[_spec.name] = sc
        _spec.loader.exec_module(sc)

#: Sent with every request, so the gateway records which version of this node called it.
PACKAGE_NAME = "cloud-modal-chai1"
PACKAGE_VERSION = "1.0.7"
SERVICE = "chai1"
#: The timeout a run needs: Chai-1's 15-minute limit, a cold start, and the download.
RECOMMENDED_TIMEOUT = 1500
#: The gateway stops waiting at 1800 s; wait a minute more for its answer.
SERVICE_MAX_WAIT = 1860

#: How to run this node on its own -- the values `salpa smoke` feeds it. Strings
#: starting with `demo_data/` resolve relative to this directory. Running needs a
#: Salpa account with cloud access; without one the node stops at authentication,
#: which is what `salpa smoke` will report. See demo_data/README.md.
DEMO_CONFIG = {
    "fasta_input": ">protein|name=trp-cage\nNLYIQWLKDGGPSSGRPPPS",
}


def _main_files(names, folder, prefix, best_idx):
    """{role: (member name, destination)}: the best sample's structure and scores."""
    chosen = {}
    cifs = [n for n in names if n.endswith(".cif")]
    best = [n for n in cifs if f"pred.model_idx_{best_idx}." in n]
    structure = (best or cifs or [None])[0]
    if structure:
        chosen["structure"] = (structure, folder / f"{prefix}_best.cif")
    scores = [n for n in names if f"scores.model_idx_{best_idx}." in n and n.endswith(".npz")]
    if scores:
        chosen["scores"] = (scores[0], folder / f"{prefix}_best_scores.npz")
    return chosen


class CloudModalChai1(Node):
    """
    Chai-1 multi-modal structure prediction on Salpa Compute (H100 GPU).

    This is a client stub: the prediction runs in the cloud, and this node sends the
    input and saves the result into the workflow folder.

    Chai-1 can predict 3D structures of:
    - Proteins (from amino acid sequences)
    - Protein-ligand complexes (protein + SMILES)
    - Protein-nucleic acid complexes (protein + DNA/RNA)
    - Multi-chain complexes

    Prerequisites:
    - The user is signed in to Salpa.
    - The user's Salpa Compute quota is not used up.

    No Modal account or 'modal setup' is needed.
    """

    # NOTE: Metadata (name, hashtags, num_in, num_out) comes from meta.toml.
    # NOTE: EXECUTION_STRATEGY and ENVIRONMENT are auto-detected via shared_environment in meta.toml.
    # NOTE: The gateway address is read when the node runs (_salpa_compute.api_base), so
    #       BOCOFLOW_CLOUD_API_URL set after import still counts.

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
                "Folder for output files. Leave empty to use the workflow's folder; a "
                "relative folder is placed inside it. "
                "Outputs: {prefix}.tar.gz (all samples), {prefix}_best.cif (best structure) "
                "and {prefix}_best_scores.npz"
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
        """Send the input to Salpa Compute and save what comes back."""
        log_message("Starting CloudModalChai1 (Salpa Compute)")

        result = NodeResult()
        result.metadata.update(
            {
                "node_type": "CloudModalChai1",
                "execution_time": datetime.now().isoformat(),
                "credential_mode": "bocoflow",  # Mode B
                "gpu": "H100",
            }
        )

        # The worker passes the signed-in user's token in this variable.
        auth_token = os.environ.get("BOCOFLOW_CLOUD_AUTH_TOKEN")
        if not auth_token:
            result.success = False
            result.message = sc.SIGN_IN_MESSAGE
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

        deadline = sc.deadline(flow_vars, RECOMMENDED_TIMEOUT, SERVICE_MAX_WAIT)
        if deadline.warning:
            stream_log(deadline.warning, node_id=self.node_id, level="warning")

        payload = {
            "node_info": {
                "node_id": getattr(self, "node_id", "unknown"),
                "node_type": "CloudModalChai1",
                "package": PACKAGE_NAME,
                "package_version": PACKAGE_VERSION,
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
            # A result too large for the reply comes as its main files plus a download link.
            "result_delivery": "archive",
            "client_request_id": sc.new_client_request_id(),
        }
        inline_max = sc.requested_inline_max_bytes()
        if inline_max is not None:
            payload["inline_max_bytes"] = inline_max

        headers = {
            "Authorization": f"Bearer {auth_token}",  # Firebase token
            "Content-Type": "application/json",
        }

        try:
            url = sc.execute_url(SERVICE)
            log_message(f"Calling Salpa Compute: {url}")
            if fasta_input:
                log_message(f"FASTA input length: {len(fasta_input)} chars")
            else:
                log_message(f"Protein sequence length: {len(protein_sequence)} amino acids")

            stream_log(
                "Calling Chai-1 on Salpa Compute... First call may take 2-3 min (cold start).",
                node_id=self.node_id,
                progress=10,
            )
            # Stop in the app cancels this run on Salpa Compute.
            with sc.cancellable(payload["client_request_id"], deadline.post_timeout):
                response = post_with_progress(
                    url=url,
                    json=payload,
                    headers=headers,
                    timeout=deadline.post_timeout,
                    node_id=self.node_id,
                    service_name="Chai-1",
                    cold_start_hint="cold starts take 2-3 min",
                )

            stream_log("Received response from cloud", node_id=self.node_id, progress=50)

            if response.status_code != 200:
                result.success = False
                result.message = sc.http_error_message(response, "Chai-1")
                return result.to_json()

            cloud_result = response.json()
            job_id = cloud_result.get("job_id")
            result.metadata["cloud_job_id"] = job_id

            failure = sc.failure_message(cloud_result)
            if failure:
                tail = sc.log_tail(cloud_result)
                if tail:
                    stream_log(
                        f"Chai-1 output before it stopped:\n{tail}",
                        node_id=self.node_id,
                        level="error",
                    )
                result.success = False
                result.message = f"Chai-1 prediction failed: {failure}" + (
                    f" (job {job_id})" if job_id else ""
                )
                return result.to_json()

            modal_result = cloud_result.get("result") or {}
            usage_info = cloud_result.get("usage") or {}

            final_folder, folder_note = sc.resolve_output_dir(self, output_folder)
            if output_prefix:
                file_prefix = output_prefix.replace("/", "_").replace("\\", "_")
            else:
                # Auto-generate from sequence or fasta
                if protein_sequence:
                    seq_short = protein_sequence[:8]
                elif fasta_input:
                    seq_short = ""
                    for line in fasta_input.strip().split("\n"):
                        if not line.startswith(">"):
                            seq_short = line[:8]
                            break
                else:
                    seq_short = "unknown"
                seq_short = "".join(c for c in seq_short if c.isalnum())
                file_prefix = f"chai1_{seq_short}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

            stream_log("Saving the result...", node_id=self.node_id, progress=60)
            try:
                saved = sc.save_result(
                    cloud_result,
                    final_folder,
                    file_prefix,
                    auth_token,
                    node_id=self.node_id,
                    deadline=deadline,
                )
            except sc.ResultError as exc:
                result.success = False
                result.message = f"Chai-1 prediction completed, but nothing could be saved. {exc}"
                return result.to_json()

            warnings = list(saved.warnings)
            if folder_note:
                warnings.append(folder_note)

            best_idx = modal_result.get("best_sample_idx", 0)
            files = {}
            try:
                chosen = _main_files(sc.members(saved.source), final_folder, file_prefix, best_idx)
                written = sc.extract(saved.source, {m: dest for m, dest in chosen.values()})
                files = {
                    role: str(dest) for role, (member, dest) in chosen.items() if member in written
                }
            except (tarfile.TarError, EOFError, OSError) as exc:
                warnings.append(f"The best structure could not be taken out of the result: {exc}")

            for warning in warnings:
                stream_log(warning, node_id=self.node_id, level="warning")

            if "structure" not in files:
                # A prediction without a structure is not a result, whatever the reply said.
                result.success = False
                result.message = (
                    "Chai-1 finished, but the result holds no structure"
                    + (f" (job {job_id})" if job_id else "")
                    + f". What came back was saved to {saved.tarball}."
                )
                return result.to_json()

            num_samples = modal_result.get("num_samples", 0)
            best_score = modal_result.get("best_aggregate_score", 0.0) or 0.0
            scores = modal_result.get("scores", [])
            tarball = saved.tarball
            structure = files.get("structure")
            duration = usage_info.get("duration_seconds", 0) or 0

            parts = [
                "Chai-1 structure prediction completed.",
                f"{num_samples} samples generated.",
                f"Best score: {best_score:.4f} (sample {best_idx}).",
            ]
            if structure:
                parts.append(f"Best CIF: {structure}.")
            if saved.complete:
                parts.append(f"Full result saved to {tarball} ({saved.tarball_bytes} bytes).")
            else:
                parts.append(f"Main files saved to {tarball} ({saved.tarball_bytes} bytes).")
            if saved.downloaded and saved.server_copy_deleted:
                parts.append("The copy held by Salpa Compute was deleted.")
            parts.append(f"Duration: {duration:.2f}s")
            if warnings:
                parts.append(f"Note: {warnings[0]}")

            result.success = True
            result.message = " ".join(parts)
            result.data = {
                "output_file": str(tarball) if tarball else None,
                "cif_file": structure,
                "output_folder": str(final_folder),
                "output_prefix": file_prefix,
                "output_file_size": saved.tarball_bytes,
                "output_files": saved.files,
                "output_tarball_available": bool(tarball),
                "num_samples": num_samples,
                "best_sample_idx": best_idx,
                "best_aggregate_score": best_score,
                "scores": scores,
                "fasta_input_length": modal_result.get("fasta_input_length", 0),
                "processing_time_seconds": modal_result.get("processing_time_seconds", 0),
                "modal_metadata": modal_result.get("modal_metadata", {}),
                "job_id": job_id,
                "usage": {
                    "duration_seconds": usage_info.get("duration_seconds", 0),
                    "cost_usd": usage_info.get("cost_usd", 0),
                },
                "status": "completed",
                "credential_mode": "bocoflow",
                "structure_file": structure,
                "scores_file": files.get("scores"),
                "archive_complete": saved.complete,
                "archive_sha256": saved.sha256,
                "delivery": saved.content,
                "warnings": warnings,
            }
            if tarball:
                result.files["output"]["archive"] = self.format_output_path(str(tarball))
            for role, path in files.items():
                result.files["output"][role] = self.format_output_path(path)
            stream_log("Chai-1 result saved.", node_id=self.node_id, progress=100)

        except requests.Timeout:
            # The node gives up; the run need not go on without it.
            stopped = sc.cancel_run(payload["client_request_id"])
            result.success = False
            result.message = (
                "No answer from Salpa Compute in time"
                + ("; the run there was cancelled. " if stopped else ". ")
                + (
                    deadline.warning
                    or "Chai-1 predictions can take several minutes; try again, or use a smaller complex."
                )
            )

        except requests.RequestException as e:
            result.success = False
            result.message = f"Network error: {str(e)}"

        except Exception as e:
            log_message(f"Unexpected error in CloudModalChai1: {str(e)}")
            result.success = False
            result.message = f"Unexpected error: {str(e)}"

        return result.to_json()
