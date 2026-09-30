"""A chain holding a residue this preparation cannot build: set aside when short, refused when long.

6LU7's N3 inhibitor is chain C, six residues, three of them non-standard with no MODRES record.
1.2.6 dropped those three and kept the other three (Ala-Val-Leu) as a free tripeptide, which was
solvated and simulated beside the protease (found inside Salpa, 2026-09-30). The rule since:
a chain of 30 residues or fewer is set aside as a ligand, with a warning; a longer one stops the
step. The residues judged are the chain's SEQRES, because a bound ligand carries a chain letter
without being part of the polymer.

Run from the package root: pixi run test
"""

import sys
from pathlib import Path

import pytest

_pkg_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_pkg_root))

from gen_ali.core import (  # noqa: E402
    PEPTIDE_LIGAND_MAX_RESIDUES,
    process_all_chains,
    read_seqres,
    unbuildable_residues,
)

EXCERPT_6LU7 = Path(__file__).parent / "data" / "6lu7_excerpt.pdb"


def _seqres(tmp_path, chains, extra=()):
    """A PDB file holding only SEQRES records (and any extra lines)."""
    lines = []
    for chain, names in chains.items():
        for i in range(0, len(names), 13):
            lines.append(
                f"SEQRES {i // 13 + 1:>3} {chain} {len(names):>4}  " + " ".join(names[i:i + 13])
            )
    path = tmp_path / "x.pdb"
    path.write_text("\n".join(lines + list(extra) + ["END"]) + "\n")
    return str(path)


class TestTheRule:
    def test_6lu7s_inhibitor_chain_is_found(self):
        assert unbuildable_residues(str(EXCERPT_6LU7)) == {
            "C": {"length": 6, "residues": [(1, "02J"), (5, "PJE"), (6, "010")]}
        }

    def test_seqres_is_read_per_chain(self):
        assert read_seqres(str(EXCERPT_6LU7)) == {
            "A": ["SER", "GLY", "PHE"],
            "C": ["02J", "ALA", "VAL", "LEU", "PJE", "010"],
        }

    def test_terminal_caps_are_left_out_and_the_chain_kept(self, tmp_path):
        path = _seqres(tmp_path, {"P": ["ACE", "ALA", "GLY", "NH2"]})
        assert unbuildable_residues(path) == {}

    def test_a_cap_name_inside_the_chain_is_not_a_cap(self, tmp_path):
        path = _seqres(tmp_path, {"P": ["ALA", "ACE", "GLY"]})
        assert unbuildable_residues(path)["P"]["residues"] == [(2, "ACE")]

    def test_a_modified_residue_modres_declares_is_built(self, tmp_path):
        path = _seqres(
            tmp_path, {"A": ["GLY", "MSE", "ALA"]},
            extra=["MODRES 1ABC MSE A    2  MET  SELENOMETHIONINE"],
        )
        assert unbuildable_residues(path) == {}

    def test_a_nucleic_acid_chain_is_not_judged(self, tmp_path):
        # DNA and RNA are set aside before modelling in any case.
        path = _seqres(tmp_path, {"I": ["DA", "5CM", "DG", "DT"]})
        assert unbuildable_residues(path) == {}

    def test_a_ligand_with_the_proteins_chain_letter_is_not_judged(self, tmp_path):
        # A bound ligand or ion is not in SEQRES, however its chain is lettered.
        path = _seqres(
            tmp_path, {"A": ["GLY", "ALA"]},
            extra=["HETATM    1 NI    NI A 201      10.000  10.000  10.000  1.00 10.00          NI"],
        )
        assert unbuildable_residues(path) == {}

    def test_a_file_without_seqres_judges_nothing(self, tmp_path):
        path = tmp_path / "model.pdb"
        path.write_text(
            "ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00 10.00           C\nEND\n"
        )
        assert unbuildable_residues(str(path)) == {}


class TestWhatTheAlignmentStepDoes:
    def test_a_short_chain_is_set_aside_and_gets_no_alignment(self, tmp_path):
        out = process_all_chains(str(EXCERPT_6LU7), str(tmp_path), "case")
        assert sorted(out["chain_results"]) == ["A"]
        assert sorted(out["set_aside_chains"]) == ["C"]
        assert out["cannot_build_chains"] == {}
        assert not (tmp_path / "C").exists()

    def test_a_longer_chain_is_refused(self, tmp_path):
        n = PEPTIDE_LIGAND_MAX_RESIDUES + 1
        names = ["ALA"] * (n - 1)
        names.insert(10, "DAL")  # a D-alanine inside a 31-residue chain
        path = _seqres(tmp_path, {"B": names})
        out = process_all_chains(path, str(tmp_path / "out"), "case")
        assert out["cannot_build_chains"] == {"B": {"length": n, "residues": [(11, "DAL")]}}
        assert out["set_aside_chains"] == {}

    @pytest.mark.parametrize("length", [PEPTIDE_LIGAND_MAX_RESIDUES])
    def test_thirty_residues_is_still_a_ligand(self, tmp_path, length):
        names = ["ALA"] * (length - 1) + ["PJE"]
        out = process_all_chains(_seqres(tmp_path, {"L": names}), str(tmp_path / "o"), "c")
        assert sorted(out["set_aside_chains"]) == ["L"]
