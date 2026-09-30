"""
merge_pdb_chains' residue filter, in its pure-Python core.

The merged structure is ProMod3's template, and ProMod3 requires it to agree residue for
residue with the alignment Generate Alignment writes. Since 1.2.5 Generate Alignment
rebuilds a modified residue the file's MODRES records declare as its standard parent
(selenomethionine, MSE, as M). This step dropped every HETATM residue, so on 3I2V chain A
the template had no residue 113 where the alignment had M, and Fix Missing Residues
stopped with "Alignment-structure mismatch at pos 113 in chain A, alignment is 'M'
structure residue is 'A'". Found running pdbmdauto 1.2.6 inside Salpa on 3I2V.

Run from the package root: pixi run test
"""

import sys
from pathlib import Path

from Bio.PDB import PDBParser

_pkg_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_pkg_root))

from merge_pdb_chains.core import merge_chains, process_merge, read_modres  # noqa: E402

# 3I2V chain A, residues 112-114, as deposited (RCSB, 2026-09-30).
MODRES_3I2V = "MODRES 3I2V MSE A  113  MET  SELENOMETHIONINE"
RESIDUES_3I2V = """\
ATOM    864  N   LEU A 112      31.545  -1.974   1.698  1.00  7.48           N
ATOM    865  CA  LEU A 112      31.567  -0.895   0.733  1.00  7.68           C
ATOM    866  C   LEU A 112      32.543  -1.173  -0.410  1.00  8.43           C
ATOM    867  O   LEU A 112      32.364  -0.604  -1.499  1.00  8.27           O
ATOM    868  CB  LEU A 112      31.943   0.430   1.454  1.00  9.41           C
ATOM    869  CG  LEU A 112      30.817   0.957   2.336  1.00  8.96           C
ATOM    870  CD1 LEU A 112      31.351   1.939   3.351  1.00 13.50           C
ATOM    871  CD2 LEU A 112      29.740   1.618   1.495  1.00 10.76           C
HETATM  872  N   MSE A 113      33.547  -2.013  -0.217  1.00  8.05           N
HETATM  873  CA  MSE A 113      34.423  -2.342  -1.345  1.00  8.94           C
HETATM  874  C   MSE A 113      33.695  -3.212  -2.360  1.00  7.71           C
HETATM  875  O   MSE A 113      33.907  -3.047  -3.554  1.00  9.06           O
HETATM  876  CB  MSE A 113      35.750  -2.942  -0.879  1.00 10.04           C
HETATM  877  CG  MSE A 113      36.635  -1.929  -0.067  1.00 11.57           C
HETATM  878 SE   MSE A 113      37.057  -0.283  -1.032  0.76 11.11          SE
HETATM  879  CE  MSE A 113      35.683   0.822  -0.343  1.00 11.32           C
ATOM    880  N   ALA A 114      32.793  -4.100  -1.897  1.00  7.94           N
ATOM    881  CA  ALA A 114      31.929  -4.833  -2.830  1.00  8.56           C
ATOM    882  C   ALA A 114      30.952  -3.894  -3.524  1.00  7.53           C
ATOM    883  O   ALA A 114      30.668  -4.033  -4.716  1.00  9.49           O
ATOM    884  CB  ALA A 114      31.224  -5.989  -2.159  1.00 10.51           C
"""
LIGAND_AND_WATER = """\
HETATM  900 NI    NI A 201      10.000  10.000  10.000  1.00 10.00          NI
HETATM  901  O   HOH A 301      12.000  12.000  12.000  1.00 10.00           O
"""


def _write(tmp_path, *parts):
    path = tmp_path / "in.pdb"
    path.write_text("".join(p if p.endswith("\n") else p + "\n" for p in parts) + "END\n")
    return str(path)


def _merge(tmp_path, *parts):
    out = tmp_path / "Merge" / "merge.pdb"
    merge_chains(_write(tmp_path, *parts), str(out), ["A"])
    return out


def _residues(pdb):
    structure = PDBParser(QUIET=True).get_structure("m", str(pdb))
    return [(r.id[0], r.id[1], r.get_resname()) for r in structure[0]["A"]]


def _atom_lines(pdb, resid):
    return [
        line
        for line in Path(pdb).read_text().splitlines()
        if line.startswith(("ATOM", "HETATM")) and int(line[22:26]) == resid
    ]


class TestReadModres:
    def test_the_columns_of_a_real_record(self, tmp_path):
        assert read_modres(_write(tmp_path, MODRES_3I2V)) == {("A", 113): ("MSE", "MET")}


class TestDeclaredModifiedResidues:
    def test_a_declared_selenomethionine_is_kept_as_methionine(self, tmp_path):
        merged = _merge(tmp_path, MODRES_3I2V, RESIDUES_3I2V)
        assert _residues(merged) == [(" ", 112, "LEU"), (" ", 113, "MET"), (" ", 114, "ALA")]

    def test_its_selenium_becomes_methionines_sulfur(self, tmp_path):
        lines = _atom_lines(_merge(tmp_path, MODRES_3I2V, RESIDUES_3I2V), 113)
        assert [line[:6] for line in lines] == ["ATOM  "] * 8
        assert [line[12:16] for line in lines] == [
            " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " SD ", " CE "
        ]
        sd = next(line for line in lines if line[12:16] == " SD ")
        assert sd[76:78] == " S"
        # The coordinates, occupancy and B-factor are the deposited selenium's.
        assert sd[30:66] == "  37.057  -0.283  -1.032  0.76 11.11"

    def test_atoms_the_parent_does_not_have_are_dropped(self, tmp_path):
        # Phosphoserine written as serine: the phosphate goes, OG stays.
        sep = """\
MODRES 1ABC SEP A    5  SER  PHOSPHOSERINE
HETATM    1  N   SEP A   5       1.000   0.000   0.000  1.00 10.00           N
HETATM    2  CA  SEP A   5       2.000   0.000   0.000  1.00 10.00           C
HETATM    3  C   SEP A   5       3.000   0.000   0.000  1.00 10.00           C
HETATM    4  O   SEP A   5       4.000   0.000   0.000  1.00 10.00           O
HETATM    5  CB  SEP A   5       2.000   1.000   0.000  1.00 10.00           C
HETATM    6  OG  SEP A   5       2.000   2.000   0.000  1.00 10.00           O
HETATM    7  P   SEP A   5       2.000   3.000   0.000  1.00 10.00           P
HETATM    8  O1P SEP A   5       2.000   4.000   0.000  1.00 10.00           O
HETATM    9  O2P SEP A   5       3.000   3.000   0.000  1.00 10.00           O
HETATM   10  O3P SEP A   5       1.000   3.000   0.000  1.00 10.00           O
"""
        lines = _atom_lines(_merge(tmp_path, sep), 5)
        assert [line[12:16].strip() for line in lines] == ["N", "CA", "C", "O", "CB", "OG"]
        assert {line[17:20] for line in lines} == {"SER"}


class TestEverythingElseIsStillLeftOut:
    def test_an_undeclared_hetatm_residue_is_left_out(self, tmp_path):
        # No MODRES: the residue stays out, as Generate Alignment also leaves it out.
        merged = _merge(tmp_path, RESIDUES_3I2V)
        assert _residues(merged) == [(" ", 112, "LEU"), (" ", 114, "ALA")]

    def test_ligands_ions_and_waters_are_left_out(self, tmp_path):
        merged = _merge(tmp_path, MODRES_3I2V, RESIDUES_3I2V, LIGAND_AND_WATER)
        assert [r[2] for r in _residues(merged)] == ["LEU", "MET", "ALA"]

    def test_a_modres_naming_another_residue_there_is_not_trusted(self, tmp_path):
        # MODRES says residue 113 is a modified CYS; the file has MSE there.
        wrong = "MODRES 3I2V CSO A  113  CYS  S-HYDROXYCYSTEINE"
        merged = _merge(tmp_path, wrong, RESIDUES_3I2V)
        assert _residues(merged) == [(" ", 112, "LEU"), (" ", 114, "ALA")]


EXCERPT_6LU7 = _pkg_root / "gen_ali" / "tests" / "data" / "6lu7_excerpt.pdb"


class TestWhichChainsAllMeans:
    """ "all" is every chain the alignment covers, when the alignment step said which.

    Generate Alignment sets aside a short chain it cannot build (6LU7's inhibitor, chain C),
    so the chain has no alignment. This step took "all" from the PDB file and merged the
    chain's three standard residues anyway, and 1.2.6 simulated them as a free tripeptide.
    """

    def test_a_chain_without_an_alignment_is_left_out(self, tmp_path):
        out = process_merge(str(EXCERPT_6LU7), str(tmp_path), "case", "all", {"A": "P1"})
        assert out.selected_chains == ["A"]
        merged = PDBParser(QUIET=True).get_structure("m", out.output_pdb)[0]
        assert [c.id for c in merged] == ["A"]

    def test_without_chain_types_all_is_still_every_chain(self, tmp_path):
        out = process_merge(str(EXCERPT_6LU7), str(tmp_path), "case", "all", None)
        assert out.selected_chains == ["A", "C"]

    def test_an_explicit_selection_is_kept(self, tmp_path):
        out = process_merge(str(EXCERPT_6LU7), str(tmp_path), "case", "A,C", {"A": "P1"})
        assert out.selected_chains == ["A", "C"]

