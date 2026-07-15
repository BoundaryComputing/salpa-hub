"""Charge injection for the EasyParm pipeline — the step the node refactor
originally dropped.

Antechamber (in ep_mol2_generation) emits a MOL2 with **zero** partial
charges; EasyParm's `01_easyPARM.sh` then injects QM-derived charges as a
separate final step. This module reimplements that step cleanly so it can be
wrapped as the `ep_charges` node and unit-tested.

Three methods, mirroring EasyParm's ORCA charge menu:

* ``orca_resp``   — ORCA 6 native ``!RESP`` block from the ``.out`` (recommended,
                    AMBER-standard; restrained + symmetry-equivalenced by ORCA)
* ``orca_chelpg`` — ORCA CHELPG block from the ``.out`` (plain ESP fit)
* ``resp_vpot``   — EasyParm's classic path: ORCA ``.vpot`` ESP grid →
                    ``RESP_ORCA.py`` writes ``resp.in`` → AmberTools ``resp``
                    → injected (needs ``similar.dat`` for equivalencing)

All three end at :func:`inject_charges_into_mol2`, the clean re-implementation
of the bundled ``Retrieve_RESP_Charges.py``.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

# ORCA prints a per-geometry block; for a single-point there is exactly one.
# We always take the LAST block so an OPT/FREQ .out yields the converged charges.
_BLOCK_HEADERS = {
    "orca_resp": "RESP Charges",
    "orca_chelpg": "CHELPG Charges",
}
# matches e.g.  "  0   N   :      -0.616603"
_CHARGE_LINE = re.compile(r"^\s*\d+\s+[A-Za-z]{1,3}\s*:\s*(-?\d+\.\d+)\s*$")


def parse_orca_charges(out_path: str | Path, kind: str) -> list[float]:
    """Return the charges from the LAST ``<kind>`` block of an ORCA ``.out``.

    ``kind`` is ``orca_resp`` or ``orca_chelpg``. Raises ValueError if no
    such block is found (e.g. the keyword was not in the ORCA input).
    """
    header = _BLOCK_HEADERS.get(kind)
    if header is None:
        raise ValueError(f"parse_orca_charges: unknown kind {kind!r}")
    lines = Path(out_path).read_text(errors="replace").splitlines()

    blocks: list[list[float]] = []
    i = 0
    n = len(lines)
    while i < n:
        if lines[i].strip().startswith(header):
            # skip the header and the dashed separator that follows
            j = i + 1
            while j < n and set(lines[j].strip()) <= {"-"} and lines[j].strip():
                j += 1
            charges: list[float] = []
            while j < n:
                m = _CHARGE_LINE.match(lines[j])
                if m:
                    charges.append(float(m.group(1)))
                    j += 1
                    continue
                # stop at the first non-charge line (blank, dashes, "Total charge:")
                if lines[j].strip() == "" or lines[j].lstrip().startswith("-") \
                        or "Total charge" in lines[j]:
                    break
                break
            if charges:
                blocks.append(charges)
            i = j
        else:
            i += 1

    if not blocks:
        raise ValueError(
            f"No '{header}' block found in {out_path}. "
            f"Was the ORCA input run with the matching keyword?"
        )
    return blocks[-1]


def inject_charges_into_mol2(mol2_path: str | Path, charges: list[float],
                             out_path: str | Path | None = None) -> str:
    """Rewrite the charge column of a TRIPOS MOL2 ATOM block, in order.

    Clean re-implementation of the bundled ``Retrieve_RESP_Charges.py``.
    The number of charges must match the number of ATOM records.
    """
    mol2_path = Path(mol2_path)
    out_path = Path(out_path) if out_path else mol2_path
    atom_re = re.compile(
        r"(\s*\d+\s+)(\S+\s+)"
        r"(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+)(\S+\s+)(\d+\s+\S+\s+)(-?\d+\.\d+)"
    )

    out_lines: list[str] = []
    in_atoms = False
    idx = 0           # charges consumed
    n_records = 0     # ATOM records seen
    for line in mol2_path.read_text().splitlines(keepends=True):
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            idx = 0
            n_records = 0
            out_lines.append(line)
            continue
        if line.startswith("@<TRIPOS>"):
            in_atoms = False
            out_lines.append(line)
            continue
        if in_atoms and line.strip():
            m = atom_re.match(line)
            if m:
                n_records += 1
                if idx < len(charges):
                    aid, name, coords, atype, post, _old = m.groups()
                    out_lines.append(f"{aid}{name}{coords}{atype}{post}{charges[idx]:8.6f}\n")
                    idx += 1
                    continue
        out_lines.append(line)

    if n_records == 0:
        raise ValueError(f"No ATOM records found in {mol2_path}")
    if n_records != len(charges):
        raise ValueError(
            f"Charge count mismatch: {len(charges)} charges vs {n_records} ATOM "
            f"records in {mol2_path}"
        )
    out_path.write_text("".join(out_lines))
    return str(out_path)


def count_mol2_atoms(mol2_path: str | Path) -> int:
    """Number of ATOM records in a MOL2 file (for sanity-checking charge count)."""
    n = 0
    in_atoms = False
    for line in Path(mol2_path).read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_atoms = False
            continue
        if in_atoms and line.strip():
            n += 1
    return n


def _run_resp_vpot(mol2_path: Path, vpot_path: Path, similar_path: Path,
                   hess_or_geom: Path, charge: int, work_dir: Path,
                   scripts_dir: Path) -> list[float]:
    """EasyParm's classic ORCA-RESP path: .vpot grid → AmberTools `resp`.

    Mirrors lines 887-901 of 01_easyPARM.sh. Returns the fitted charges.
    """
    esp_in = work_dir / "esp.in"
    text = vpot_path.read_text().splitlines()
    # EasyParm strips the spaces from the first (count) line
    if text:
        first = re.sub(r"^(\s+)(\d+)(\s+)(\d+)", r"\1\2\4", text[0])
        text[0] = first
    esp_in.write_text("\n".join(text) + "\n")

    subprocess.run(
        ["python3", str(scripts_dir / "RESP_ORCA.py"),
         str(hess_or_geom), str(charge), str(similar_path)],
        cwd=work_dir, check=True, capture_output=True,
    )
    subprocess.run(
        ["resp", "-O", "-i", "resp.in", "-o", "resp.out", "-e", "esp.in", "-t", "esp.chg"],
        cwd=work_dir, check=True, capture_output=True,
    )
    esp_chg = work_dir / "esp.chg"
    toks = esp_chg.read_text().split()
    return [float(t) for t in toks]


def run_charges(mol2_file: str | Path, method: str,
                orca_out: str | Path | None = None,
                vpot_file: str | Path | None = None,
                similar_file: str | Path | None = None,
                geom_file: str | Path | None = None,
                charge: int = 0,
                work_dir: str | Path | None = None,
                scripts_dir: str | Path | None = None) -> dict:
    """Inject QM-derived charges into ``mol2_file`` and return a summary.

    Returns ``{method, n_atoms, net_charge, max_abs_charge, output_mol2}``.
    """
    mol2_file = Path(mol2_file)
    work_dir = Path(work_dir) if work_dir else mol2_file.parent
    n_atoms = count_mol2_atoms(mol2_file)

    if method in ("orca_resp", "orca_chelpg"):
        if not orca_out:
            raise ValueError(f"method {method!r} requires orca_out")
        charges = parse_orca_charges(orca_out, method)
    elif method == "resp_vpot":
        for label, val in (("vpot_file", vpot_file), ("similar_file", similar_file),
                           ("geom_file", geom_file)):
            if not val:
                raise ValueError(f"method 'resp_vpot' requires {label}")
        scripts_dir = Path(scripts_dir or Path(__file__).parent / "scripts")
        charges = _run_resp_vpot(mol2_file, Path(vpot_file), Path(similar_file),
                                 Path(geom_file), charge, work_dir, scripts_dir)
    else:
        raise ValueError(f"run_charges: unknown method {method!r}")

    if len(charges) != n_atoms:
        raise ValueError(
            f"{method}: parsed {len(charges)} charges but MOL2 has {n_atoms} atoms"
        )

    inject_charges_into_mol2(mol2_file, charges)
    net = sum(charges)
    return {
        "method": method,
        "n_atoms": n_atoms,
        "net_charge": net,
        "max_abs_charge": max(abs(c) for c in charges),
        "output_mol2": str(mol2_file),
    }
