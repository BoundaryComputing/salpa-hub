"""The alignment Generate Alignment writes agrees with the template Merge PDB Chains writes.

ProMod3 reads the two together and stops at the first residue where they differ
("Alignment-structure mismatch at pos N"). Each node was tested alone, and 1.2.5 changed
only one of them: Generate Alignment began rebuilding a MODRES-declared residue as its
parent (MSE as M) while Merge PDB Chains still dropped every HETATM residue. The pair
disagreed on 3I2V chain A at residue 113 and Fix Missing Residues failed, which only a
run of the whole template inside Salpa showed (2026-09-30). This checks the pair.

Run from the package root: pixi run test
"""

import sys
from pathlib import Path

import pytest
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG))

from gen_ali.core import process_all_chains  # noqa: E402
from merge_pdb_chains.core import merge_chains, process_merge  # noqa: E402

MODRES = "MODRES 1ABC MSE A    3  MET  SELENOMETHIONINE"
REMARK_465 = """\
REMARK 465 MISSING RESIDUES
REMARK 465   M RES C SSSEQI
REMARK 465     LYS A     5"""


def _line(record, serial, name, resname, resid, x, element):
    field = f" {name:<3}" if len(name) < 4 and len(element) == 1 else f"{name:<4}"
    return (
        f"{record:<6}{serial:>5} {field} {resname:>3} A{resid:>4}    "
        f"{x:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{1.0:>6.2f}{10.0:>6.2f}          {element:>2}"
    )


def _structure(tmp_path, header):
    """Chain A: G L (MSE) A present, K missing (REMARK 465), W present."""
    residues = [
        ("ATOM", "GLY", 1, [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]),
        ("ATOM", "LEU", 2, [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")]),
        ("HETATM", "MSE", 3, [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C"),
                              ("CG", "C"), ("SE", "SE"), ("CE", "C")]),
        ("ATOM", "ALA", 4, [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")]),
        ("ATOM", "TRP", 6, [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")]),
        ("HETATM", "HOH", 101, [("O", "O")]),
    ]
    lines, serial = list(header) + [REMARK_465], 0
    for record, resname, resid, atoms in residues:
        for name, element in atoms:
            serial += 1
            lines.append(_line(record, serial, name, resname, resid, float(serial), element))
    path = tmp_path / "1ABC.pdb"
    path.write_text("\n".join(lines + ["END"]) + "\n")
    return str(path)


def _pair(tmp_path, header):
    pdb = _structure(tmp_path, header)
    alignment = process_all_chains(pdb, str(tmp_path / "ali"), "case")["chain_results"]["A"]
    merged = tmp_path / "Merge" / "merge.pdb"
    merge_chains(pdb, str(merged), ["A"])
    chain = PDBParser(QUIET=True).get_structure("m", str(merged))[0]["A"]
    template = "".join(seq1(r.get_resname()) for r in chain)
    track = alignment.template_seq.replace("-", "").rstrip("*")
    return track, template


@pytest.mark.parametrize(
    "header", [[MODRES], []], ids=["MSE declared by MODRES", "MSE not declared"]
)
def test_the_alignments_structure_track_is_the_templates_sequence(tmp_path, header):
    track, template = _pair(tmp_path, header)
    assert track == template


def test_a_declared_selenomethionine_is_in_both_as_methionine(tmp_path):
    track, template = _pair(tmp_path, [MODRES])
    assert template == "GLMAW"


def test_a_chain_the_alignment_sets_aside_is_not_in_the_template(tmp_path):
    # 6LU7: Generate Alignment sets chain C (the N3 inhibitor) aside, and Merge PDB Chains
    # must leave it out too, or its three standard residues reach ProMod3's template.
    excerpt = PKG / "gen_ali" / "tests" / "data" / "6lu7_excerpt.pdb"
    results = process_all_chains(str(excerpt), str(tmp_path / "ali"), "case")["chain_results"]
    chain_types = {c: r.chain_type for c, r in results.items()}
    merged = process_merge(str(excerpt), str(tmp_path), "case", "all", chain_types)
    assert merged.selected_chains == sorted(results) == ["A"]
    chain = PDBParser(QUIET=True).get_structure("m", merged.output_pdb)[0]["A"]
    assert "".join(seq1(r.get_resname()) for r in chain) == \
        results["A"].template_seq.replace("-", "").rstrip("*")

