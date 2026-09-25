"""
The GenAli node's sequence check, run on the package's own demo structure (3LZ0).

Until 1.2.4 the node said "Sequence agreement: all match" having compared nothing. It
looked for `*_chain_*.fasta`, which pdb_fasta_biopython stopped writing on 2026-03-31,
and those files hold only the residues that have coordinates, so on a structure with
gaps they would have reported a mismatch that is not one. The reference is the
deposited sequence, `{PDB_ID}_rcsb.fasta`.

Run from the package root: pixi run test
"""

import json
import shutil
import sys
from pathlib import Path

import pytest

try:
    from bocoflow_core.parameters import (
        BooleanParameter,
        FolderParameter,
        StringParameter,
    )
except ImportError:
    pytest.skip(
        "bocoflow_core not installed: these tests run the node through BoCoFlow",
        allow_module_level=True,
    )

_pkg_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_pkg_root))

from gen_ali.node import GenAli  # noqa: E402

DEMO = _pkg_root / "pdb_fasta_biopython" / "demo_data"


def _run(folder):
    node = GenAli(
        {"node_id": "test-gen-ali", "package_name": "pdbmdauto", "name": "Gen Ali"}
    )
    flow_vars = {
        "case_name": StringParameter("Case Name"),
        "output_dir": FolderParameter("Output Directory"),
        "append_end_in_seq": BooleanParameter("Append End Marker", default=True),
        "force_to_run": BooleanParameter("Force to Run", default=False),
    }
    flow_vars["case_name"].set_value("3lz0")
    flow_vars["output_dir"].set_value(f"abs:{folder}")
    flow_vars["append_end_in_seq"].set_value(True)
    upstream = {"case_name": "3lz0", "working_path": f"abs:{folder}", "chain_info": {}}
    return json.loads(node.execute([upstream], flow_vars))


def _folder(tmp_path, reference=None):
    shutil.copy(DEMO / "3LZ0.pdb", tmp_path / "3LZ0.pdb")
    if reference is not None:
        (tmp_path / "3LZ0_rcsb.fasta").write_text(reference)
    return tmp_path


def test_the_demo_matches_its_deposited_sequence(tmp_path):
    result = _run(_folder(tmp_path, (DEMO / "3LZ0.fasta").read_text()))
    assert result["success"] is True
    assert "Sequence check against RCSB: all 10 chain(s) match." in result["message"]
    assert result["data"]["seq_agree_pdb_fasta"] is True
    assert result["data"]["sequence_check"]["reference"] == "3LZ0_rcsb.fasta"


def test_a_difference_is_a_warning_naming_the_chains(tmp_path):
    lines = (DEMO / "3LZ0.fasta").read_text().splitlines()
    assert lines[0].startswith(">3LZ0_1|Chains A, E|")
    lines[1] = "G" + lines[1][1:]  # histone H3's first residue, A in the deposit
    result = _run(_folder(tmp_path, "\n".join(lines) + "\n"))
    assert result["success"] is True
    assert "MISMATCH in chain(s) A, E" in result["message"]
    assert result["data"]["seq_agree_pdb_fasta"] is False
    assert result["data"]["sequence_check"]["mismatched"] == ["A", "E"]


def test_without_a_deposited_sequence_nothing_is_claimed(tmp_path):
    # A structure opened from a local file has no RCSB sequence beside it.
    result = _run(_folder(tmp_path))
    assert "all match" not in result["message"]
    assert "Sequences not checked" in result["message"]
    assert result["data"]["seq_agree_pdb_fasta"] is None
    assert result["data"]["sequence_check"]["not_compared"] == list("ABCDEFGHIJ")


def test_per_chain_files_are_not_a_reference_under_either_name(tmp_path):
    # pdb_fasta_biopython's per-chain files: residues with coordinates only, 97 of
    # chain A's 135 in 3LZ0. Compared with the rebuild, they would call a correct
    # rebuild a mismatch.
    folder = _folder(tmp_path)
    observed = (
        ">3lz0_chain_A Chain A from 3LZ0.pdb - protein - 97 residues\n"
        "PHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVALFEDTNLC\n"
        "AIHAKRVTIMPKDIQLARRIRGER\n"
    )
    (folder / "chain_A.fasta").write_text(observed)
    (folder / "3lz0_chain_A.fasta").write_text(observed)
    result = _run(folder)
    assert "MISMATCH" not in result["message"]
    assert result["data"]["seq_agree_pdb_fasta"] is None
