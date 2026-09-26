"""
Cloud Modal Boltz-2 — a Salpa Compute node.

Protein structure prediction using Boltz-2 on an H100 GPU.

This node is a client stub. It sends the input to the Salpa Compute gateway, which checks
the user's sign-in and quota and runs Boltz-2 on Modal. No Modal account is needed.

How the result comes back:
    - A result that fits in the reply (up to 16 MiB compressed) arrives whole.
    - A larger one (MSA-server runs, large complexes) arrives as its main files -- the
      structure, confidence and pLDDT -- plus a link to the full archive. The node
      downloads the archive into the workflow folder, checks its size and SHA-256, and
      then has the server copy deleted. Salpa Compute never keeps it longer than 24 hours.

Based on: https://modal.com/docs/examples/boltz_predict
"""

import os
import re
import sys
import tarfile
from datetime import datetime
from pathlib import Path

import requests
from bocoflow_core.logger import log_message
from bocoflow_core.node import Node, NodeResult
from bocoflow_core.parameters import FolderParameter, TextParameter
from bocoflow_core.stream_logger import post_with_progress, stream_log

try:  # loaded as a package: node_runner, the server, tests
    from . import _salpa_compute as sc
except ImportError:
    try:  # the node's folder is on sys.path
        import _salpa_compute as sc
    except ImportError:  # anything else: load the file that sits beside this one
        import importlib.util

        _spec = importlib.util.spec_from_file_location(
            "_salpa_compute_cloud_modal_boltz2",
            str(Path(__file__).with_name("_salpa_compute.py")),
        )
        sc = importlib.util.module_from_spec(_spec)
        sys.modules[_spec.name] = sc
        _spec.loader.exec_module(sc)

#: Sent with every request, so the gateway records which version of this node called it.
PACKAGE_NAME = "cloud-modal-boltz2"
PACKAGE_VERSION = "1.0.7"
SERVICE = "boltz2"
#: The timeout a run needs: Boltz-2's 30-minute limit, a cold start, and the download.
RECOMMENDED_TIMEOUT = 2400
#: The gateway stops waiting for Boltz-2 at 1800 s; wait a minute more for its answer.
SERVICE_MAX_WAIT = 1860

#: How to run this node on its own -- the values `salpa smoke` feeds it. Strings
#: starting with `demo_data/` resolve relative to this directory. Running needs a
#: Salpa account with cloud access; without one the node stops at authentication,
#: which is what `salpa smoke` will report. See demo_data/README.md.
DEMO_CONFIG = {
    "sequence": "NLYIQWLKDGGPSSGRPPPS",
    "msa_mode": "empty",
}


def _pick(names, *patterns):
    """The first member name matching the first pattern that matches anything."""
    for pattern in patterns:
        rx = re.compile(pattern)
        for name in names:
            if rx.search(name):
                return name
    return None


def _main_files(names, folder, prefix):
    """{role: (member name, destination)} for the files a user opens first."""
    chosen = {}
    structure = _pick(
        names,
        r"(^|/)predictions/.*_model_0\.(cif|pdb)$",
        r"(^|/)predictions/.*\.(cif|pdb)$",
        r"\.(cif|pdb)$",
    )
    if structure:
        chosen["structure"] = (structure, folder / f"{prefix}_model_0{Path(structure).suffix}")
    confidence = _pick(
        names, r"(^|/)confidence_[^/]*_model_0\.json$", r"(^|/)confidence_[^/]*\.json$"
    )
    if confidence:
        chosen["confidence"] = (confidence, folder / f"{prefix}_confidence_model_0.json")
    plddt = _pick(names, r"(^|/)plddt_[^/]*_model_0\.npz$", r"(^|/)plddt_[^/]*\.npz$")
    if plddt:
        chosen["plddt"] = (plddt, folder / f"{prefix}_plddt_model_0.npz")
    affinity = _pick(names, r"(^|/)affinity_[^/]*\.json$")
    if affinity:
        chosen["affinity"] = (affinity, folder / f"{prefix}_affinity.json")
    return chosen


class CloudModalBoltz2(Node):
    """
    Boltz-2 protein structure prediction on Salpa Compute (H100 GPU).

    This is a client stub: the prediction runs in the cloud, and this node sends the
    input and saves the result into the workflow folder.

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
                "Folder for output files. Leave empty to use the workflow's folder; a "
                "relative folder is placed inside it. "
                "Outputs: {prefix}.tar.gz (the full result), {prefix}_model_0.cif, "
                "{prefix}_confidence_model_0.json and {prefix}_plddt_model_0.npz"
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
        """Send the input to Salpa Compute and save what comes back."""
        log_message("Starting CloudModalBoltz2 (Salpa Compute)")

        result = NodeResult()
        result.metadata.update(
            {
                "node_type": "CloudModalBoltz2",
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

        deadline = sc.deadline(flow_vars, RECOMMENDED_TIMEOUT, SERVICE_MAX_WAIT)
        if deadline.warning:
            stream_log(deadline.warning, node_id=self.node_id, level="warning")

        payload = {
            "node_info": {
                "node_id": getattr(self, "node_id", "unknown"),
                "node_type": "CloudModalBoltz2",
                "package": PACKAGE_NAME,
                "package_version": PACKAGE_VERSION,
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
            log_message(f"Sequence length: {len(sequence)} amino acids")
            log_message(f"MSA mode: {msa_mode}")

            stream_log(
                "Calling Boltz-2 on Salpa Compute... First call may take 2-3 min (cold start).",
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
                    service_name="Boltz-2",
                    cold_start_hint="cold starts take 2-3 min",
                )

            if response.status_code != 200:
                result.success = False
                result.message = sc.http_error_message(response, "Boltz-2")
                return result.to_json()

            cloud_result = response.json()
            job_id = cloud_result.get("job_id")
            result.metadata["cloud_job_id"] = job_id

            failure = sc.failure_message(cloud_result)
            if failure:
                tail = sc.log_tail(cloud_result)
                if tail:
                    stream_log(
                        f"Boltz-2 output before it stopped:\n{tail}",
                        node_id=self.node_id,
                        level="error",
                    )
                result.success = False
                result.message = f"Boltz-2 prediction failed: {failure}" + (
                    f" (job {job_id})" if job_id else ""
                )
                return result.to_json()

            modal_result = cloud_result.get("result") or {}
            usage_info = cloud_result.get("usage") or {}

            final_folder, folder_note = sc.resolve_output_dir(self, output_folder)
            if output_prefix:
                file_prefix = output_prefix.replace("/", "_").replace("\\", "_")
            else:
                seq_short = "".join(c for c in (sequence or "")[:8] if c.isalnum())
                file_prefix = f"boltz_{seq_short}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

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
                result.message = f"Boltz-2 prediction completed, but nothing could be saved. {exc}"
                return result.to_json()

            warnings = list(saved.warnings)
            if folder_note:
                warnings.append(folder_note)

            files = {}
            try:
                chosen = _main_files(sc.members(saved.source), final_folder, file_prefix)
                written = sc.extract(saved.source, {m: dest for m, dest in chosen.values()})
                files = {
                    role: str(dest) for role, (member, dest) in chosen.items() if member in written
                }
            except (tarfile.TarError, EOFError, OSError) as exc:
                warnings.append(f"The main files could not be taken out of the result: {exc}")

            for warning in warnings:
                stream_log(warning, node_id=self.node_id, level="warning")

            tarball = saved.tarball
            if "structure" not in files:
                # A prediction without a structure is not a result, whatever the reply said.
                result.success = False
                result.message = (
                    "Boltz-2 finished, but the result holds no structure"
                    + (f" (job {job_id})" if job_id else "")
                    + f". What came back was saved to {tarball}."
                )
                return result.to_json()
            structure = files.get("structure")
            duration = usage_info.get("duration_seconds", 0) or 0
            parts = ["Boltz-2 structure prediction completed."]
            if structure:
                parts.append(f"Structure: {structure}.")
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
                "cif_file": structure if structure and structure.endswith(".cif") else None,
                "output_folder": str(final_folder),
                "output_prefix": file_prefix,
                "output_file_size": saved.tarball_bytes,
                "output_files": saved.files,
                "output_tarball_available": bool(tarball),
                "sequence_length": modal_result.get("sequence_length", len(sequence)),
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
                "confidence_file": files.get("confidence"),
                "plddt_file": files.get("plddt"),
                "affinity_file": files.get("affinity"),
                "archive_complete": saved.complete,
                "archive_sha256": saved.sha256,
                "delivery": saved.content,
                "warnings": warnings,
            }
            if tarball:
                result.files["output"]["archive"] = self.format_output_path(str(tarball))
            for role, path in files.items():
                result.files["output"][role] = self.format_output_path(path)
            stream_log("Boltz-2 result saved.", node_id=self.node_id, progress=100)

        except requests.Timeout:
            # The node gives up; the run need not go on without it.
            stopped = sc.cancel_run(payload["client_request_id"])
            result.success = False
            result.message = (
                "No answer from Salpa Compute in time"
                + ("; the run there was cancelled. " if stopped else ". ")
                + (
                    deadline.warning
                    or "Boltz-2 predictions can take up to 30 minutes; try again, or use a shorter sequence."
                )
            )

        except requests.RequestException as e:
            result.success = False
            result.message = f"Network error: {str(e)}"

        except Exception as e:
            log_message(f"Unexpected error in CloudModalBoltz2: {str(e)}")
            result.success = False
            result.message = f"Unexpected error: {str(e)}"

        return result.to_json()
