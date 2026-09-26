"""Salpa Compute client helper, shared by every Salpa Compute stub.

Each stub package carries an identical copy of this file; ``HELPER_VERSION`` says which
revision it is. It needs only the standard library and ``requests``, and runs on
Python 3.9, so it works against any core the app has shipped.

It answers the questions every stub has to get right:

- where to send the job: ``api_base()`` reads ``BOCOFLOW_CLOUD_API_URL`` on every call and
  otherwise uses https://compute.salpa.app;
- where the results go: ``resolve_output_dir()`` puts them in the workflow folder, never
  in the folder the node process happens to run in;
- how long to wait: ``deadline()`` keeps the wait inside the node's own timeout, so the stub
  stops cleanly before the app stops it;
- what the gateway answered: ``failure_message()``, ``inline_content()``,
  ``inline_tarball()``;
- how to receive a result too large to come back in the reply: the reply then carries the
  main files and a signed link to the full archive, and ``save_result()`` downloads it
  (``download_archive()``: resumable, size and SHA-256 checked, never sending the sign-in
  token to storage), then asks for the server copy to be deleted (``delete_archive()``).
  Salpa Compute deletes an archive nobody fetched within 24 hours.
"""

import base64
import hashlib
import io
import logging
import os
import re
import shutil
import tarfile
import tempfile
import time
import uuid
from pathlib import Path

import requests

HELPER_VERSION = 3

DEFAULT_API_URL = "https://compute.salpa.app"
API_URL_ENV = "BOCOFLOW_CLOUD_API_URL"
#: Tests and end-to-end checks lower the inline limit with this, to force the archive path.
INLINE_MAX_ENV = "BOCOFLOW_CLOUD_INLINE_MAX_BYTES"
INLINE_MAX_BYTES = 16 * 1024 * 1024

#: What the app stops a node at when its settings carry no timeout.
EXECUTOR_DEFAULT_TIMEOUT = 3600.0
#: Kept back from the node's timeout, so the stub gives up before the app kills it.
SAFETY_MARGIN_SECONDS = 30.0
MIN_POST_TIMEOUT = 5.0

#: A deadline assumes downloads run at least this fast, plus a fixed overhead.
ASSUMED_DOWNLOAD_RATE = 512 * 1024
DOWNLOAD_OVERHEAD_SECONDS = 30.0

DOWNLOAD_ATTEMPTS = 5
RETRY_BASE_SECONDS = 2.0
PROGRESS_INTERVAL_SECONDS = 5.0
CONNECT_TIMEOUT = 15
READ_TIMEOUT = 120
CHUNK_BYTES = 1024 * 1024
DISK_MARGIN_BYTES = 64 * 1024 * 1024
DELETE_TIMEOUT = 30

_MIB = 1024 * 1024
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


# ---------------------------------------------------------------------------
# Where to send the job
# ---------------------------------------------------------------------------


def api_base():
    """The gateway's base URL, read from the environment on every call."""
    value = (os.environ.get(API_URL_ENV) or "").strip()
    return (value or DEFAULT_API_URL).rstrip("/")


def execute_url(service):
    """The route that runs ``service`` (for example ``boltz2``)."""
    return f"{api_base()}/api/cloud/nodes/{service}/execute"


def archive_url(job_id):
    """The route that re-signs (GET) or deletes (DELETE) a job's archive."""
    return f"{api_base()}/api/cloud/jobs/{job_id}/archive"


def new_client_request_id():
    """An identifier for this request, stored with the job by the gateway."""
    return uuid.uuid4().hex


def requested_inline_max_bytes():
    """The inline limit to ask for, when ``BOCOFLOW_CLOUD_INLINE_MAX_BYTES`` sets one."""
    raw = (os.environ.get(INLINE_MAX_ENV) or "").strip()
    if not raw:
        return None
    try:
        value = int(raw)
    except ValueError:
        return None
    return max(0, min(value, INLINE_MAX_BYTES))


# ---------------------------------------------------------------------------
# Where the results go
# ---------------------------------------------------------------------------


def _fallback_dir():
    downloads = Path.home() / "Downloads"
    if downloads.is_dir():
        return downloads
    return Path(tempfile.gettempdir())


def workflow_dir(node):
    """The workflow folder of this run, or None when the run has none."""
    try:
        resolved = node.resolve_path("rel:.")
    except Exception:
        return None
    if not resolved or str(resolved).startswith("rel:"):
        return None
    # Without a working path the core falls back to the custom-nodes folder, which is
    # not a place for results.
    custom = getattr(node, "custom_nodes_dir", None)
    if custom and os.path.normpath(str(resolved)) == os.path.normpath(str(custom)):
        return None
    return Path(os.path.normpath(str(resolved)))


def resolve_output_dir(node, folder):
    """Where to save results: ``(path, note)``.

    An empty folder means the workflow folder. A relative one (``rel:x`` or ``x``) is
    placed inside the workflow folder. Only a run with no workflow folder at all falls
    back to ``~/Downloads``, and ``note`` then says so, for the node to show the user.
    """
    text = str(folder).strip() if folder else ""
    wdir = workflow_dir(node)
    if text:
        try:
            resolved = str(node.resolve_path(text) or text)
        except Exception:
            resolved = text
        for prefix in ("rel:", "abs:", "node:"):
            if resolved.startswith(prefix):
                resolved = resolved[len(prefix) :]
                break
        path = Path(os.path.expanduser(resolved))
        if path.is_absolute():
            return path, None
        if wdir is not None:
            return wdir / path, None
        base = _fallback_dir()
        return base / path, (
            f"This run has no workflow folder, so the output folder was placed under {base}."
        )
    if wdir is not None:
        return wdir, None
    base = _fallback_dir()
    return base, f"This run has no workflow folder, so the results were saved to {base}."


# ---------------------------------------------------------------------------
# How long to wait
# ---------------------------------------------------------------------------


def _minutes(seconds):
    minutes = float(seconds) / 60.0
    if minutes >= 1 and abs(minutes - round(minutes)) < 0.05:
        return f"{int(round(minutes))} min"
    if minutes >= 1:
        return f"{minutes:.1f} min"
    return f"{int(round(float(seconds)))} s"


def _timeout_from(flow_vars):
    if not flow_vars:
        return None
    try:
        param = flow_vars.get("_timeout_seconds")
    except AttributeError:
        return None
    if param is None:
        return None
    value = param.get_value() if hasattr(param, "get_value") else param
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if value > 0 else None


class Deadline:
    """The node's own timeout, counted from when the stub started.

    ``post_timeout`` is how long to wait for the gateway: the service's own limit, but
    never past the node's timeout (less a margin). ``warning`` is set when the node's
    timeout is shorter than a run of this service can take.
    """

    def __init__(
        self, timeout_seconds, recommended_seconds, service_max_seconds, clock=time.monotonic
    ):
        self._clock = clock
        self.started = clock()
        self.configured = timeout_seconds
        self.timeout_seconds = (
            float(timeout_seconds) if timeout_seconds else EXECUTOR_DEFAULT_TIMEOUT
        )
        self.recommended_seconds = float(recommended_seconds)
        self.service_max_seconds = float(service_max_seconds)
        self.warning = None
        if timeout_seconds and float(timeout_seconds) < self.recommended_seconds:
            self.warning = (
                f"This node's timeout is {_minutes(timeout_seconds)}, and a run can take up to "
                f"{_minutes(recommended_seconds)}. Set Advanced Options > Timeout (seconds) to "
                f"{int(self.recommended_seconds)}, or a long run will be stopped before its "
                f"result arrives."
            )

    def elapsed(self):
        return self._clock() - self.started

    def remaining(self):
        return self.timeout_seconds - self.elapsed()

    @property
    def post_timeout(self):
        usable = self.remaining() - SAFETY_MARGIN_SECONDS
        return max(MIN_POST_TIMEOUT, min(self.service_max_seconds, usable))

    def can_download(self, nbytes):
        needed = float(nbytes or 0) / ASSUMED_DOWNLOAD_RATE + DOWNLOAD_OVERHEAD_SECONDS
        return self.remaining() - SAFETY_MARGIN_SECONDS >= needed

    def expired(self):
        """True when the stub must stop now to finish before the node's timeout."""
        return self.remaining() <= SAFETY_MARGIN_SECONDS / 2


def deadline(flow_vars, recommended_seconds, service_max_seconds):
    """A :class:`Deadline` from the node's ``_timeout_seconds`` setting."""
    return Deadline(_timeout_from(flow_vars), recommended_seconds, service_max_seconds)


# ---------------------------------------------------------------------------
# Reading the gateway's answer
# ---------------------------------------------------------------------------


def _text(value):
    return str(value).strip() if value else ""


def failure_message(cloud_result):
    """The reason a job failed, or None when it succeeded."""
    if not isinstance(cloud_result, dict):
        return "Salpa Compute sent a reply this node cannot read."
    result = cloud_result.get("result")
    status = _text(cloud_result.get("status")).lower()
    if status in ("failed", "error", "cancelled"):
        reason = _text(cloud_result.get("error"))
        if not reason and isinstance(result, dict):
            reason = _text(result.get("error"))
        return reason or "The job failed, and Salpa Compute gave no reason."
    if isinstance(result, dict) and _text(result.get("status")).lower() == "error":
        return _text(result.get("error")) or "The service reported an error without a reason."
    return None


def log_tail(cloud_result):
    """The last lines of the tool's output that came with a failure, if any."""
    if not isinstance(cloud_result, dict):
        return ""
    result = cloud_result.get("result")
    if isinstance(result, dict):
        return _text(result.get("log_tail"))
    return ""


def inline_content(result):
    """What the reply itself carries: ``"full"``, ``"essentials"`` or ``"none"``."""
    if not isinstance(result, dict):
        return "none"
    content = result.get("inline_content")
    if content in ("full", "essentials", "none"):
        return content
    # A reply from before this field existed carries the whole result, or nothing.
    return "full" if result.get("output_tarball_base64") else "none"


def inline_tarball(result):
    """The tar.gz the reply carries, as bytes (empty when there is none)."""
    if not isinstance(result, dict):
        return b""
    encoded = result.get("output_tarball_base64") or ""
    return base64.b64decode(encoded) if encoded else b""


SIGN_IN_MESSAGE = (
    "This node runs on Salpa Compute. Sign in to Salpa to use it; "
    "no Modal or Google Cloud account is needed."
)


def http_error_message(response, label):
    """A sentence for a gateway answer other than HTTP 200."""
    status = response.status_code
    detail = ""
    try:
        body = response.json()
        if isinstance(body, dict):
            detail = _text(body.get("detail") or body.get("error"))
    except ValueError:
        detail = _text(response.text)[:500]
    if status == 401:
        return "Your Salpa sign-in has expired or was refused. Sign in again, then rerun the node."
    if status in (402, 403):
        return detail or "Your Salpa Compute quota is used up for this period."
    if status == 404:
        return detail or f"Salpa Compute at {api_base()} does not offer {label}."
    if status == 429:
        # The gateway says why: most often two GPU runs are already in progress.
        return detail or "Salpa Compute is busy. Try again in a minute."
    if status == 503:
        return (
            f"{label} is temporarily unavailable on Salpa Compute; its GPU may be starting. "
            f"Try again in a few minutes." + (f" ({detail})" if detail else "")
        )
    if status == 504:
        return f"Salpa Compute stopped waiting for {label}." + (f" ({detail})" if detail else "")
    return f"Salpa Compute answered HTTP {status}" + (f": {detail}" if detail else ".")


# ---------------------------------------------------------------------------
# Taking files out of a tar.gz
# ---------------------------------------------------------------------------


def _open_tar(source):
    if isinstance(source, (bytes, bytearray)):
        return tarfile.open(fileobj=io.BytesIO(bytes(source)), mode="r:gz")
    return tarfile.open(name=str(source), mode="r:gz")


def members(source):
    """The names of the regular files in a tar.gz (bytes or a path), in archive order."""
    names = []
    with _open_tar(source) as tar:
        for member in tar:
            if member.isfile():
                names.append(member.name)
    return names


def safe_name(name, default="file"):
    """A file name taken from an archive member, safe to join onto a folder."""
    base = re.split(r"[\\/]", str(name))[-1]
    base = re.sub(r"[^A-Za-z0-9._-]", "_", base).lstrip(".")
    return base or default


def _write_stream(dest, stream):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f".{dest.name}.{os.getpid()}.tmp")
    try:
        with open(tmp, "wb") as out:
            shutil.copyfileobj(stream, out, CHUNK_BYTES)
        os.replace(tmp, dest)
    except BaseException:
        _unlink(tmp)
        raise


def write_atomic(dest, data):
    """Write bytes to ``dest`` through a temporary file, so no half-written file remains."""
    _write_stream(dest, io.BytesIO(data))


def extract(source, picks):
    """Copy chosen regular files out of a tar.gz.

    ``picks`` maps an archive member name to the path it is written to. Only regular files
    are taken, each through a temporary file and ``os.replace``; nothing else in the
    archive is written anywhere. Returns ``{member name: path written}``.
    """
    wanted = {str(name): Path(dest) for name, dest in dict(picks).items()}
    written = {}
    if not wanted:
        return written
    with _open_tar(source) as tar:
        for member in tar:
            if not member.isfile() or member.name in written:
                continue
            dest = wanted.get(member.name)
            if dest is None:
                continue
            handle = tar.extractfile(member)
            if handle is None:
                continue
            with handle:
                _write_stream(dest, handle)
            written[member.name] = dest
    return written


# ---------------------------------------------------------------------------
# Receiving an archive
# ---------------------------------------------------------------------------


class DownloadError(Exception):
    """The archive could not be downloaded; ``gone`` means it is no longer held."""

    def __init__(self, message, gone=False):
        super().__init__(message)
        self.gone = gone


def _unlink(path):
    try:
        os.unlink(str(path))
    except OSError:
        pass


def _mib(nbytes):
    return f"{float(nbytes or 0) / _MIB:.1f} MiB"


def _describe(exc):
    # A requests exception names the URL, and a storage URL carries its signature, so only
    # the kind of failure is ever reported.
    if isinstance(exc, requests.Timeout):
        return "the connection timed out"
    if isinstance(exc, requests.ConnectionError):
        return "the connection failed"
    return f"a network error ({type(exc).__name__})"


class _QuietHttpLogs:
    """urllib3 logs request lines, signed URLs included, at DEBUG. Silence it meanwhile."""

    def __init__(self):
        self._logger = logging.getLogger("urllib3")
        self._level = self._logger.level

    def __enter__(self):
        if self._logger.getEffectiveLevel() < logging.WARNING:
            self._logger.setLevel(logging.WARNING)
        return self

    def __exit__(self, *exc):
        self._logger.setLevel(self._level)
        return False


def _stream_progress(message, node_id, progress=None, level="info"):
    try:
        from bocoflow_core.stream_logger import stream_log
    except Exception:
        return
    try:
        stream_log(message, node_id=node_id, progress=progress, level=level)
    except Exception:
        pass


class _Progress:
    def __init__(self, node_id, total, start=70, end=90):
        self.node_id = node_id
        self.total = max(int(total or 0), 1)
        self.start = start
        self.end = end
        self._last = 0.0

    def update(self, received, force=False):
        now = time.monotonic()
        if not force and now - self._last < PROGRESS_INTERVAL_SECONDS:
            return
        self._last = now
        fraction = min(1.0, float(received) / self.total)
        _stream_progress(
            f"Downloading the full result: {_mib(received)} of {_mib(self.total)}",
            self.node_id,
            progress=round(self.start + (self.end - self.start) * fraction),
        )


def _archive_info(archive):
    if not isinstance(archive, dict):
        raise DownloadError("the reply carried no usable download link")
    url = _text(archive.get("url"))
    sha256 = _text(archive.get("sha256")).lower()
    try:
        size = int(archive.get("bytes"))
    except (TypeError, ValueError):
        size = 0
    if not url or size <= 0 or not _SHA256.match(sha256):
        raise DownloadError("the reply carried an incomplete download link")
    job_id = _text(archive.get("job_id"))
    refresh_path = _text(archive.get("refresh_path"))
    if not refresh_path and job_id:
        refresh_path = f"/api/cloud/jobs/{job_id}/archive"
    return {
        "url": url,
        "bytes": size,
        "sha256": sha256,
        "job_id": job_id,
        "refresh_path": refresh_path,
    }


def _check_disk(folder, size):
    try:
        free = shutil.disk_usage(str(folder)).free
    except OSError:
        return
    if free < size + DISK_MARGIN_BYTES:
        raise DownloadError(
            f"there is not enough free disk space in {folder} "
            f"({_mib(free)} free, {_mib(size + DISK_MARGIN_BYTES)} needed)"
        )


def _refresh(info, auth_token):
    """A freshly signed link from the gateway, or None to keep trying the old one."""
    path = info.get("refresh_path")
    if not path:
        return None
    if not path.startswith("/"):
        path = "/" + path
    try:
        response = requests.get(
            api_base() + path,
            headers={"Authorization": f"Bearer {auth_token}"},
            timeout=(CONNECT_TIMEOUT, 60),
        )
    except requests.RequestException:
        return None
    with response:
        if response.status_code in (404, 410):
            raise DownloadError(
                "Salpa Compute no longer holds this result (it expired or was deleted)",
                gone=True,
            )
        if response.status_code != 200:
            return None
        try:
            data = response.json()
        except ValueError:
            return None
    if not isinstance(data, dict):
        return None
    new_sha = _text(data.get("sha256")).lower()
    if new_sha and new_sha != info["sha256"]:
        raise DownloadError("the refreshed link points at a different archive")
    return _text(data.get("url")) or None


def _range_start(response):
    match = re.match(r"bytes\s+(\d+)-", response.headers.get("Content-Range", ""))
    return int(match.group(1)) if match else None


def download_archive(archive, dest, auth_token, node_id=None, deadline=None):
    """Download a job's archive to ``dest`` and check it.

    The file is streamed into ``<dest>.part``; an interrupted transfer resumes with a
    ``Range`` request, and a link the storage refuses (HTTP 400/401/403/404) is replaced by
    a fresh one from the gateway. The storage request never carries the sign-in token. The
    result must match the announced size and SHA-256 before it is renamed to ``dest``; a
    mismatch is downloaded once more from scratch, then given up. At most
    ``DOWNLOAD_ATTEMPTS`` requests are made. Returns ``{"path", "bytes", "sha256"}``, or
    raises :class:`DownloadError` (after removing the partial file).
    """
    info = _archive_info(archive)
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    total = info["bytes"]
    _check_disk(dest.parent, total)

    part = dest.with_name(dest.name + ".part")
    _unlink(part)  # a leftover from an earlier run is not this archive
    url = info["url"]
    hasher = hashlib.sha256()
    received = 0
    mismatches = 0
    last_error = "no attempt was made"
    progress = _Progress(node_id, total)

    def restart():
        nonlocal hasher, received
        hasher = hashlib.sha256()
        received = 0
        _unlink(part)

    def out_of_time():
        return deadline is not None and deadline.expired()

    def read_timeout():
        if deadline is None:
            return READ_TIMEOUT
        return max(5.0, min(READ_TIMEOUT, deadline.remaining() - SAFETY_MARGIN_SECONDS / 2))

    with _QuietHttpLogs():
        for attempt in range(1, DOWNLOAD_ATTEMPTS + 1):
            if attempt > 1:
                time.sleep(RETRY_BASE_SECONDS * (2 ** (attempt - 2)))
            if out_of_time():
                last_error = "this node's timeout was reached"
                break
            headers = {}
            if 0 < received < total:
                headers["Range"] = f"bytes={received}-"
            try:
                response = requests.get(
                    url, headers=headers, stream=True, timeout=(CONNECT_TIMEOUT, read_timeout())
                )
            except requests.RequestException as exc:
                last_error = _describe(exc)
                continue

            with response:
                status = response.status_code
                if status in (400, 401, 403, 404):
                    last_error = f"the storage refused the link (HTTP {status})"
                    try:
                        fresh = _refresh(info, auth_token)
                    except DownloadError:
                        _unlink(part)
                        raise
                    if fresh:
                        url = fresh
                    continue
                if status == 416:
                    last_error = "the storage refused the resume request (HTTP 416)"
                    restart()
                    continue
                if status == 200:
                    restart()
                    mode = "wb"
                elif status == 206 and received > 0 and _range_start(response) == received:
                    mode = "ab"
                elif status == 206:
                    last_error = "the storage sent a part that was not asked for"
                    restart()
                    continue
                else:
                    last_error = f"the storage answered HTTP {status}"
                    continue

                try:
                    with open(part, mode) as out:
                        for chunk in response.iter_content(CHUNK_BYTES):
                            if not chunk:
                                continue
                            out.write(chunk)
                            hasher.update(chunk)
                            received += len(chunk)
                            if received > total or out_of_time():
                                break
                            progress.update(received)
                except requests.RequestException as exc:
                    last_error = f"the transfer stopped at {_mib(received)}: {_describe(exc)}"
                    continue

            if received < total and out_of_time():
                last_error = f"this node's timeout was reached at {_mib(received)} of {_mib(total)}"
                break
            if received < total:
                last_error = f"the transfer stopped at {_mib(received)} of {_mib(total)}"
                continue
            if received > total or hasher.hexdigest() != info["sha256"]:
                mismatches += 1
                restart()
                if mismatches >= 2:
                    raise DownloadError("the download does not match its announced checksum")
                last_error = "the download did not match its announced checksum"
                continue

            os.replace(part, dest)
            progress.update(total, force=True)
            return {"path": str(dest), "bytes": total, "sha256": info["sha256"]}

    _unlink(part)
    if out_of_time():
        raise DownloadError(last_error)
    raise DownloadError(f"it failed after {DOWNLOAD_ATTEMPTS} attempts; {last_error}")


def delete_archive(job_id, auth_token):
    """Ask the gateway to delete a job's archive. True when it confirmed."""
    if not job_id:
        return False
    try:
        response = requests.delete(
            archive_url(job_id),
            headers={"Authorization": f"Bearer {auth_token}"},
            timeout=DELETE_TIMEOUT,
        )
    except requests.RequestException:
        return False
    with response:
        return response.status_code in (200, 202, 204)


# ---------------------------------------------------------------------------
# Saving a result
# ---------------------------------------------------------------------------


class ResultError(Exception):
    """Nothing usable could be saved; the message is for the user."""


class SavedResult:
    """What :func:`save_result` put on disk.

    ``tarball`` is ``<prefix>.tar.gz``. ``complete`` says whether it holds the whole result
    or only the main files. ``source`` is where to take the main files from (bytes, or the
    tarball's path). ``warnings`` are sentences for the user.
    """

    def __init__(self):
        self.tarball = None
        self.tarball_bytes = 0
        self.complete = False
        self.sha256 = None
        self.source = None
        self.content = "none"
        self.files = []
        self.warnings = []
        self.downloaded = False
        self.server_copy_deleted = False


def save_result(cloud_result, folder, prefix, auth_token, node_id=None, deadline=None):
    """Save what the gateway delivered as ``<folder>/<prefix>.tar.gz``.

    - The whole result came in the reply: it is written as is.
    - The main files came in the reply and the full archive by link: the archive is
      downloaded, checked, and the server copy deleted. If the download cannot be done,
      the main files are saved instead, with a warning.
    - Only a link: the archive is downloaded; if that fails, :class:`ResultError`.
    - Only the main files and no link (an older gateway, or a failed upload): they are
      saved, with the service's notice as a warning.
    """
    if not isinstance(cloud_result, dict):
        raise ResultError("Salpa Compute sent a reply this node cannot read.")
    result = cloud_result.get("result")
    result = result if isinstance(result, dict) else {}
    job_id = _text(cloud_result.get("job_id"))
    archive = cloud_result.get("archive")
    archive = archive if isinstance(archive, dict) and archive.get("url") else None

    saved = SavedResult()
    saved.warnings.extend(w for w in (cloud_result.get("warnings") or []) if _text(w))
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    tar_path = folder / f"{prefix}.tar.gz"

    content = inline_content(result)
    saved.content = content
    inline = inline_tarball(result) if content != "none" else b""
    all_files = list(result.get("all_output_files") or result.get("output_files") or [])

    if content == "full":
        write_atomic(tar_path, inline)
        saved.tarball = tar_path
        saved.tarball_bytes = len(inline)
        saved.complete = True
        saved.sha256 = hashlib.sha256(inline).hexdigest()
        saved.source = inline
        saved.files = all_files
        if archive:
            # A copy the reply already made redundant.
            saved.server_copy_deleted = delete_archive(
                job_id or _text(archive.get("job_id")), auth_token
            )
        return saved

    if archive:
        size = archive.get("bytes") or 0
        if deadline is not None and not deadline.can_download(size):
            reason = f"too little time is left before this node's timeout to download {_mib(size)}"
        else:
            try:
                got = download_archive(archive, tar_path, auth_token, node_id, deadline)
            except DownloadError as exc:
                reason = str(exc)
            else:
                saved.tarball = tar_path
                saved.tarball_bytes = got["bytes"]
                saved.complete = True
                saved.sha256 = got["sha256"]
                saved.downloaded = True
                saved.source = inline if (content == "essentials" and inline) else tar_path
                saved.files = all_files
                saved.server_copy_deleted = delete_archive(
                    job_id or _text(archive.get("job_id")), auth_token
                )
                if not saved.server_copy_deleted:
                    saved.warnings.append(
                        "The copy held by Salpa Compute could not be deleted just now; "
                        "it is deleted automatically within 24 hours."
                    )
                return saved
        if content != "essentials":
            raise ResultError(
                f"The result could not be downloaded: {reason}. "
                f"Salpa Compute deletes its copy within 24 hours (job {job_id or 'unknown'})."
            )
        saved.warnings.append(
            f"The main files were saved, but the full result ({_mib(size)}) could not be "
            f"downloaded: {reason}. Salpa Compute deletes its copy within 24 hours "
            f"(job {job_id or 'unknown'})."
        )

    if content == "essentials" and not inline:
        raise ResultError(
            f"Salpa Compute announced the main files but sent none (job {job_id or 'unknown'})."
        )
    if content == "essentials":
        write_atomic(tar_path, inline)
        saved.tarball = tar_path
        saved.tarball_bytes = len(inline)
        saved.complete = False
        saved.source = inline
        saved.files = list(result.get("output_files") or [])
        if not archive:
            notice = _text(result.get("notice"))
            if notice:
                saved.warnings.append(notice)
            elif not saved.warnings:
                saved.warnings.append(
                    "Only the main files came back; the rest of the result was not delivered."
                )
        return saved

    raise ResultError(f"Salpa Compute returned no result files (job {job_id or 'unknown'}).")


# ---------------------------------------------------------------------------
# Stopping a run
# ---------------------------------------------------------------------------

CANCEL_TIMEOUT = (5, 10)
#: How long the record for Stop outlives the node's own wait.
CANCEL_RECORD_MARGIN_SECONDS = 120


def cancel_url():
    return f"{api_base()}/api/cloud/jobs/cancel"


class _Cancellable:
    """Records, for the app's Stop, which Salpa Compute run this node is waiting on.

    Stop kills a node outright, so the node cannot cancel its run itself. The app's
    server reads this record and sends the cancel, then stops the node. It needs a Salpa
    that knows how (bocoflow_core.stop_flag.register_cloud_request); with an older one
    this does nothing, and a stopped run goes on to its end as before.
    """

    def __init__(self, client_request_id, seconds):
        self.client_request_id = client_request_id
        self.seconds = seconds
        self.registered = False

    def __enter__(self):
        try:
            from bocoflow_core.stop_flag import register_cloud_request
        except ImportError:
            return self
        try:
            self.registered = bool(
                register_cloud_request(
                    self.client_request_id,
                    api_base(),
                    int(self.seconds) + CANCEL_RECORD_MARGIN_SECONDS,
                )
            )
        except Exception:
            self.registered = False
        return self

    def __exit__(self, *exc):
        if self.registered:
            try:
                from bocoflow_core.stop_flag import clear_cloud_request

                clear_cloud_request()
            except Exception:
                pass
        return False


def cancellable(client_request_id, seconds):
    """``with cancellable(request_id, wait_seconds): <send the request>``."""
    return _Cancellable(client_request_id, seconds)


def cancel_run(client_request_id):
    """Ask Salpa Compute to stop a run this node stopped waiting for. Best effort.

    Returns whether Salpa Compute accepted the cancel. The request id is the credential,
    so no sign-in token is needed.
    """
    try:
        response = requests.post(
            cancel_url(),
            json={"client_request_id": client_request_id},
            timeout=CANCEL_TIMEOUT,
        )
        response.close()
        return response.status_code == 202
    except requests.RequestException:
        return False
