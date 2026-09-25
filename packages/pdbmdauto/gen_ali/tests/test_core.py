"""
gen_ali's sequence check, in its pure-Python core.

gen_ali rebuilds each chain's full sequence from the residues that have coordinates
plus the ones REMARK 465 lists as missing, and ProMod3 models the chain from it. The
check compares that rebuild with the sequence the depositors gave, which
pdb_fasta_biopython saves as {PDB_ID}_rcsb.fasta when it fetches a structure. Until
1.2.4 it compared nothing and reported "all match": it looked for per-chain FASTA
files under a name no longer written, and those hold only the residues with
coordinates, so they could not have checked the rebuild anyway.

Run from the package root: pixi run test
"""

import sys
from pathlib import Path

_pkg_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_pkg_root))

from gen_ali.core import (  # noqa: E402
    find_reference_fasta,
    process_all_chains,
    read_modres,
    read_reference_sequences,
)

RCSB = """>1ABC_1|Chains A, B|Some protein|Homo sapiens (9606)
MAGKLW
>1ABC_2|Chain C[auth D]|A peptide|synthetic construct (32630)
QEEW
ETVM
"""

REMARK_465 = """REMARK 465
REMARK 465 MISSING RESIDUES
REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE
REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN
REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)
REMARK 465
REMARK 465   M RES C SSSEQI
REMARK 465     LYS A     4
REMARK 465     LEU A     5"""


def _atom(record, serial, resname, resid):
    return (
        "{:<6}{:>5} {:<4}{:1}{:>3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}"
        "{:>6.2f}{:>6.2f}          {:>2}"
    ).format(
        record,
        serial,
        " CA ",
        "",
        resname,
        "A",
        resid,
        "",
        float(resid),
        0.0,
        0.0,
        1.0,
        0.0,
        "C",
    )


MODRES_MSE = "MODRES 1ABC MSE A    2  MET  SELENOMETHIONINE"


def _structure(tmp_path, second=("ATOM", "ALA"), header=()):
    """Chain A: M A G present, K L missing (REMARK 465), W present -> MAGKLW."""
    residues = [
        ("ATOM", "MET", 1),
        (*second, 2),
        ("ATOM", "GLY", 3),
        ("ATOM", "TRP", 6),
    ]
    lines = (
        list(header)
        + [REMARK_465]
        + [
            _atom(record, i + 1, name, resid)
            for i, (record, name, resid) in enumerate(residues)
        ]
    )
    path = tmp_path / "1ABC.pdb"
    path.write_text("\n".join(lines + ["END"]) + "\n")
    return str(path)


class TestReadModres:
    def test_the_columns_of_a_real_record(self, tmp_path):
        # 3I2V's own line.
        path = tmp_path / "x.pdb"
        path.write_text(
            "MODRES 3I2V MSE A  113  MET  SELENOMETHIONINE                       \n"
        )
        assert read_modres(str(path)) == {("A", 113): ("MSE", "MET")}


class TestReadReferenceSequences:
    def test_chains_are_keyed_by_the_names_the_pdb_file_uses(self, tmp_path):
        path = tmp_path / "1ABC_rcsb.fasta"
        path.write_text(RCSB)
        assert read_reference_sequences(str(path)) == {
            "A": "MAGKLW",
            "B": "MAGKLW",
            "D": "QEEWETVM",
        }

    def test_a_per_chain_file_names_no_chains(self, tmp_path):
        path = tmp_path / "chain_A.fasta"
        path.write_text(
            ">1abc_chain_A Chain A from 1ABC.pdb - protein - 4 residues\nMAGW\n"
        )
        assert read_reference_sequences(str(path)) == {}


class TestFindReferenceFasta:
    def test_the_structures_own_file_comes_first(self, tmp_path):
        for name in ("1ABC_rcsb.fasta", "2XYZ_rcsb.fasta", "chain_A.fasta"):
            (tmp_path / name).write_text(RCSB)
        found = find_reference_fasta(str(tmp_path), str(tmp_path / "1ABC.pdb"))
        assert found == str(tmp_path / "1ABC_rcsb.fasta")

    def test_otherwise_the_only_one_even_in_a_folder_named_like_a_glob(self, tmp_path):
        folder = tmp_path / "case [1]"
        folder.mkdir()
        (folder / "1ABC_rcsb.fasta").write_text(RCSB)
        (folder / "chain_A.fasta").write_text(">x\nMAGW\n")
        found = find_reference_fasta(str(folder), str(folder / "fixed.pdb"))
        assert found == str(folder / "1ABC_rcsb.fasta")

    def test_none_when_there_is_none_or_several_to_choose_from(self, tmp_path):
        (tmp_path / "chain_A.fasta").write_text(">x\nMAGW\n")
        assert find_reference_fasta(str(tmp_path), str(tmp_path / "1ABC.pdb")) is None
        (tmp_path / "2XYZ_rcsb.fasta").write_text(RCSB)
        (tmp_path / "3DEF_rcsb.fasta").write_text(RCSB)
        assert find_reference_fasta(str(tmp_path), str(tmp_path / "1ABC.pdb")) is None


class TestTheCheck:
    def _check(self, tmp_path, reference, second=("ATOM", "ALA"), header=()):
        pdb = _structure(tmp_path, second, header)
        return process_all_chains(
            pdb, str(tmp_path / "out"), "t", append_end=True, reference_seqs=reference
        )

    def test_a_rebuild_that_matches_the_deposited_sequence(self, tmp_path):
        out = self._check(tmp_path, {"A": "MAGKLW"})
        assert out["seq_agree_all"] is True
        assert out["mismatched_chains"] == []
        assert out["not_compared_chains"] == []

    def test_a_rebuild_that_differs_is_a_mismatch_that_says_where(self, tmp_path):
        out = self._check(tmp_path, {"A": "MAGRLW"})
        assert out["seq_agree_all"] is False
        assert out["mismatched_chains"] == ["A"]
        assert out["chain_results"]["A"].seq_mismatch == (
            "reconstructed 6 residues, deposited 6; "
            "first difference at position 4: K here, R deposited"
        )

    def test_no_reference_is_not_compared_rather_than_a_match(self, tmp_path):
        out = self._check(tmp_path, {})
        assert out["seq_agree_all"] is None
        assert out["not_compared_chains"] == ["A"]
        assert out["chain_results"]["A"].seq_agree is None

    def test_a_selenomethionine_modres_declares_is_rebuilt_as_methionine(
        self, tmp_path
    ):
        # 3I2V chain A: MSE is a HETATM record, and MODRES names MET as its parent.
        # Left out, the chain came out one residue short and ProMod3 stopped
        # ("Alignment-structure mismatch"), because OpenStructure reads MSE as M.
        out = self._check(
            tmp_path, {"A": "MMGKLW"}, second=("HETATM", "MSE"), header=[MODRES_MSE]
        )
        assert out["seq_agree_all"] is True
        assert out["chain_results"]["A"].template_seq.startswith("MMG--W")

    def test_an_undeclared_hetatm_is_still_left_out_and_caught(self, tmp_path):
        out = self._check(tmp_path, {"A": "MMGKLW"}, second=("HETATM", "MSE"))
        assert out["mismatched_chains"] == ["A"]
        assert out["chain_results"]["A"].seq_mismatch.startswith(
            "reconstructed 5 residues, deposited 6"
        )
