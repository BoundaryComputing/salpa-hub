"""
Cloud Modal DiffDock — a Salpa Compute node.

Blind protein-ligand docking using DiffDock on an A10G GPU.

This node is a client stub. It sends the protein and ligand to the Salpa Compute gateway,
which checks the user's sign-in and quota and runs DiffDock on Modal. No Modal account is
needed.

How the result comes back:
    - A result that fits in the reply (up to 16 MiB compressed, every run so far) arrives
      whole.
    - A larger one arrives as its SDF poses plus a link to the full archive. The node
      downloads the archive into the workflow folder, checks its size and SHA-256, and
      then has the server copy deleted. Salpa Compute never keeps it longer than 24 hours.

IMPORTANT — Cold Start Warning:
    The first call after the GPU container scales to zero takes ~10 minutes
    (loading ESM-2 protein language model, 2.6 GB, into GPU memory).
    Subsequent calls within 2 minutes are "warm" and take only ~40 seconds.

    This node is designed for INTERACTIVE single-protein docking, NOT for
    batch processing. For batch docking (many protein-ligand pairs), use
    a local DiffDock installation or a dedicated HPC pipeline instead.

Reference: Corso et al., ICLR 2023 (MIT License)
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
    FileParameterEdit,
    FolderParameter,
    IntegerParameter,
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
            "_salpa_compute_cloud_modal_diffdock",
            str(Path(__file__).with_name("_salpa_compute.py")),
        )
        sc = importlib.util.module_from_spec(_spec)
        sys.modules[_spec.name] = sc
        _spec.loader.exec_module(sc)

#: Sent with every request, so the gateway records which version of this node called it.
PACKAGE_NAME = "cloud-modal-diffdock"
PACKAGE_VERSION = "1.0.5"
SERVICE = "modal-diffdock"
#: The timeout a run needs: DiffDock's 15-minute limit, a cold start, and the download.
RECOMMENDED_TIMEOUT = 1500
#: The gateway stops waiting for DiffDock at 900 s; wait a minute more for its answer.
SERVICE_MAX_WAIT = 960


def _main_files(names, folder, prefix):
    """{member name: destination} for every SDF pose, and the top pose's destination."""
    picks = {}
    top = None
    for name in names:
        if not name.endswith(".sdf"):
            continue
        base = sc.safe_name(name)
        dest = folder / f"{prefix}_{base}"
        if dest in picks.values():
            continue  # two complexes with the same pose name: keep the first
        picks[name] = dest
        if top is None and base.startswith("rank1"):
            top = dest
        if base == "rank1.sdf":
            top = dest
    return picks, top


class CloudModalDiffdock(Node):
    """
    DiffDock blind protein-ligand docking on Salpa Compute (A10G GPU).

    This is a client stub: the docking runs in the cloud, and this node sends the input
    and saves the ranked poses into the workflow folder.

    DiffDock is a generative diffusion model that predicts how small molecules
    bind to protein targets — a critical step in drug discovery pipelines.

    Timing:
    - Cold start: ~10 minutes (ESM-2 model loading into GPU memory)
    - Warm start: ~40 seconds (container reused within 2-minute window)
    - NOT suitable for batch processing of many protein-ligand pairs

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
                "Folder for output files. Leave empty to use the workflow's folder; a "
                "relative folder is placed inside it. "
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
        """Send the input to Salpa Compute and save what comes back."""
        stream_log(
            "Starting DiffDock on Salpa Compute... "
            "First call may take ~10 min (cold start). Warm calls take ~40s.",
            node_id=self.node_id,
            progress=0,
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

        # The worker passes the signed-in user's token in this variable.
        auth_token = os.environ.get("BOCOFLOW_CLOUD_AUTH_TOKEN")
        if not auth_token:
            result.success = False
            result.message = sc.SIGN_IN_MESSAGE
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
                protein_pdb = pred_data.get("pdb_content", "") or pred_data.get("protein_pdb", "")
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
                "No ligand provided. Please enter a SMILES string or select an SDF file."
            )
            return result.to_json()

        # Get parameters
        num_poses = flow_vars["num_poses"].get_value()
        inference_steps = flow_vars["inference_steps"].get_value()
        samples_per_complex = flow_vars["samples_per_complex"].get_value()
        output_folder = flow_vars["output_folder"].get_value() or ""
        output_prefix = flow_vars["output_prefix"].get_value() or ""

        deadline = sc.deadline(flow_vars, RECOMMENDED_TIMEOUT, SERVICE_MAX_WAIT)
        if deadline.warning:
            stream_log(deadline.warning, node_id=self.node_id, level="warning")

        payload = {
            "node_info": {
                "node_id": getattr(self, "node_id", "unknown"),
                "node_type": "CloudModalDiffdock",
                "package": PACKAGE_NAME,
                "package_version": PACKAGE_VERSION,
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
            # A result too large for the reply comes as its poses plus a download link.
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
            stream_log(
                f"Calling DiffDock (poses={num_poses}, steps={inference_steps})...",
                node_id=self.node_id,
                progress=10,
            )
            log_message(
                f"Protein PDB: {len(protein_pdb)} chars, "
                f"ligand: {'SMILES' if ligand_smiles else 'SDF file'}"
            )

            # DiffDock: warm ~40s, cold start ~590s (ESM-2 model loading)
            # Stop in the app cancels this run on Salpa Compute.
            with sc.cancellable(payload["client_request_id"], deadline.post_timeout):
                response = post_with_progress(
                    url=sc.execute_url(SERVICE),
                    json=payload,
                    headers=headers,
                    timeout=deadline.post_timeout,
                    node_id=self.node_id,
                    service_name="DiffDock",
                    cold_start_hint="cold starts take up to 10 min",
                )

            if response.status_code != 200:
                result.success = False
                result.message = sc.http_error_message(response, "DiffDock")
                stream_log(f"Error: {result.message}", node_id=self.node_id, level="error")
                return result.to_json()

            cloud_result = response.json()
            job_id = cloud_result.get("job_id")
            result.metadata["cloud_job_id"] = job_id

            failure = sc.failure_message(cloud_result)
            if failure:
                tail = sc.log_tail(cloud_result)
                if tail:
                    stream_log(
                        f"DiffDock output before it stopped:\n{tail}",
                        node_id=self.node_id,
                        level="error",
                    )
                result.success = False
                result.message = f"DiffDock prediction failed: {failure}" + (
                    f" (job {job_id})" if job_id else ""
                )
                stream_log(f"Error: {result.message}", node_id=self.node_id, level="error")
                return result.to_json()

            modal_result = cloud_result.get("result") or {}
            usage_info = cloud_result.get("usage") or {}
            actual_poses = modal_result.get("num_poses", 0)
            top_confidence = modal_result.get("top_confidence", 0.0) or 0.0
            confidence_scores = modal_result.get("confidence_scores", [])
            processing_time = modal_result.get("processing_time_seconds", 0)

            final_folder, folder_note = sc.resolve_output_dir(self, output_folder)
            if output_prefix:
                file_prefix = output_prefix.replace("/", "_").replace("\\", "_")
            else:
                file_prefix = f"diffdock_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

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
                result.message = f"DiffDock completed, but nothing could be saved. {exc}"
                return result.to_json()

            warnings = list(saved.warnings)
            if folder_note:
                warnings.append(folder_note)

            extracted_sdfs = []
            top_pose = None
            try:
                picks, top = _main_files(sc.members(saved.source), final_folder, file_prefix)
                written = sc.extract(saved.source, picks)
                extracted_sdfs = [str(picks[m]) for m in picks if m in written]
                top_pose = str(top) if top is not None and str(top) in extracted_sdfs else None
            except (tarfile.TarError, EOFError, OSError) as exc:
                warnings.append(f"The poses could not be taken out of the result: {exc}")

            for warning in warnings:
                stream_log(warning, node_id=self.node_id, level="warning")

            tarball = saved.tarball
            if not extracted_sdfs:
                # Docking without a pose is not a result, whatever the reply said.
                result.success = False
                result.message = (
                    "DiffDock finished, but the result holds no pose"
                    + (f" (job {job_id})" if job_id else "")
                    + f". What came back was saved to {tarball}."
                )
                stream_log(f"Error: {result.message}", node_id=self.node_id, level="error")
                return result.to_json()
            duration = usage_info.get("duration_seconds", 0) or 0
            stream_log(
                f"DiffDock completed: {actual_poses} poses, top confidence={top_confidence:.4f}",
                node_id=self.node_id,
                progress=90,
            )
            parts = [
                f"DiffDock generated {actual_poses} docking poses.",
                f"Top confidence: {top_confidence:.4f}.",
            ]
            if saved.complete:
                parts.append(f"Output: {tarball} ({saved.tarball_bytes} bytes).")
            else:
                parts.append(f"Main files: {tarball} ({saved.tarball_bytes} bytes).")
            if extracted_sdfs:
                parts.append(f"Extracted {len(extracted_sdfs)} SDF files.")
            if saved.downloaded and saved.server_copy_deleted:
                parts.append("The copy held by Salpa Compute was deleted.")
            parts.append(f"Duration: {duration:.2f}s")
            if warnings:
                parts.append(f"Note: {warnings[0]}")

            result.success = True
            result.message = " ".join(parts)
            result.data = {
                "output_file": str(tarball) if tarball else None,
                "output_folder": str(final_folder),
                "output_prefix": file_prefix,
                "extracted_sdfs": extracted_sdfs,
                "output_file_size": saved.tarball_bytes,
                "output_files": saved.files,
                "output_tarball_available": bool(tarball),
                "num_poses": actual_poses,
                "top_confidence": top_confidence,
                "confidence_scores": confidence_scores,
                "processing_time_seconds": processing_time,
                "modal_metadata": modal_result.get("modal_metadata", {}),
                "job_id": job_id,
                "protein_pdb": protein_pdb,
                "ligand_smiles": ligand_smiles,
                "usage": {
                    "duration_seconds": usage_info.get("duration_seconds", 0),
                    "cost_usd": usage_info.get("cost_usd", 0),
                },
                "status": "completed",
                "credential_mode": "bocoflow",
                "top_pose_file": top_pose,
                "archive_complete": saved.complete,
                "archive_sha256": saved.sha256,
                "delivery": saved.content,
                "warnings": warnings,
            }
            if tarball:
                result.files["output"]["archive"] = self.format_output_path(str(tarball))
            if top_pose:
                result.files["output"]["top_pose"] = self.format_output_path(top_pose)
            for index, path in enumerate(extracted_sdfs):
                result.files["output"][f"pose_{index + 1}"] = self.format_output_path(path)

        except requests.Timeout:
            # The node gives up; the run need not go on without it.
            stopped = sc.cancel_run(payload["client_request_id"])
            result.success = False
            result.message = (
                "No answer from Salpa Compute in time"
                + ("; the run there was cancelled. " if stopped else ". ")
                + (
                    deadline.warning
                    or "DiffDock cold starts take ~10 min (ESM-2 model loading); warm calls take "
                    "~40s. Try again: the container may now be warm."
                )
            )

        except requests.RequestException as e:
            result.success = False
            result.message = f"Network error: {str(e)}"

        except Exception as e:
            log_message(f"Unexpected error in CloudModalDiffdock: {str(e)}")
            result.success = False
            result.message = f"Unexpected error: {str(e)}"

        return result.to_json()
