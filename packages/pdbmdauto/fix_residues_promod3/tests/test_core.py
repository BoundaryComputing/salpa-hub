"""
Tests for Fix Missing Residues' record of what it rebuilt, and for failing loudly.

ProMod3 numbers the model by position in the target sequence, from 1, so the residues it
rebuilt are the target positions where the template has no residue. The node writes them to
rebuilt_residues.json beside fixed.pdb, which Original Atom Groups reads.

A protein chain ProMod3 cannot model used to be dropped from fixed.pdb while the step still
reported success, so every later step ran on a structure missing a chain. It now fails.

No ProMod3 here: the record is built from files, and the modelling call is replaced.
"""

import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

PACKAGE = Path(__file__).resolve().parent.parent.parent
_spec = importlib.util.spec_from_file_location(
    "fix_residues_core", PACKAGE / "fix_residues_promod3" / "core.py"
)
core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(core)


def write_fasta(path, target, template):
    path.write_text(f">target\n{target}\n>merge.A\n{template}\n")
    return path


def write_model_pdb(path, chains):
    """A fixed.pdb as ProMod3 writes it: chains in order, each numbered as given."""
    lines, serial = [], 0
    for chain, numbers in chains:
        for resid in numbers:
            for name in ("N", "CA", "C", "O"):
                serial += 1
                lines.append(
                    f"ATOM  {serial:5d}  {name:<3s} ALA {chain}{resid:4d}    "
                    f"   0.000   0.000   0.000  1.00  0.00           {name[0]}\n"
                )
        lines.append("TER\n")
    path.write_text("".join(lines) + "END\n")
    return path


def write_ali(path, template, target):
    """A MODELLER .ali: template block first, target second, as Generate Alignment writes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f">P1;caseA\nstructure:Merge:1:A:+{len(template.replace('-', ''))}:A::::\n{template}*\n\n"
        f">P1;caseA_full\nsequence:::::::::\n{target}*\n"
    )
    return path


def test_rebuilt_positions_are_where_the_template_has_no_residue(tmp_path):
    # 4Z8J chain A in miniature: a five-residue tag missing before the first crystal residue.
    tag = write_fasta(tmp_path / "a.fasta", "GSHGGSPRV", "-----SPRV")
    assert core.rebuilt_positions(str(tag)) == [1, 2, 3, 4, 5]
    loop = write_fasta(tmp_path / "b.fasta", "ACDEFGH", "AC---GH")
    assert core.rebuilt_positions(str(loop)) == [3, 4, 5]


def test_a_template_residue_absent_from_the_target_shifts_nothing(tmp_path):
    # A column where the target has a gap is not a target position at all.
    extra = write_fasta(tmp_path / "c.fasta", "AC-DEF", "ACXD-F")
    assert core.rebuilt_positions(str(extra)) == [4]


def test_the_record_lists_each_chain_in_model_order_with_its_rebuilt_residues(tmp_path):
    fixed = write_model_pdb(
        tmp_path / "fixed.pdb", [("A", range(1, 10)), ("B", range(1, 4))]
    )
    fastas = {
        "A": str(write_fasta(tmp_path / "alignment_A.fasta", "GSHGGSPRV", "-----SPRV")),
        "B": str(write_fasta(tmp_path / "alignment_B.fasta", "QLV", "-LV")),
    }
    out = core.write_rebuilt_map(
        str(fixed), fastas, str(tmp_path / "rebuilt_residues.json")
    )
    record = json.loads(Path(out).read_text())
    assert record["numbering"] == "model"
    assert record["chains"] == [
        {"chain": "A", "residues": list(range(1, 10)), "rebuilt": [1, 2, 3, 4, 5]},
        {"chain": "B", "residues": [1, 2, 3], "rebuilt": [1]},
    ]


def test_positions_the_model_does_not_contain_are_not_listed_as_rebuilt(tmp_path):
    # With terminal modelling off, the tag is simply absent from the model: nothing rebuilt.
    fixed = write_model_pdb(tmp_path / "fixed.pdb", [("A", range(6, 10))])
    fastas = {
        "A": str(write_fasta(tmp_path / "alignment_A.fasta", "GSHGGSPRV", "-----SPRV"))
    }
    out = core.write_rebuilt_map(
        str(fixed), fastas, str(tmp_path / "rebuilt_residues.json")
    )
    record = json.loads(Path(out).read_text())
    assert record["chains"] == [{"chain": "A", "residues": [6, 7, 8, 9], "rebuilt": []}]


def _case_with_alignments(tmp_path, chains):
    case = tmp_path / "case"
    for chain in chains:
        write_ali(case / chain / "homology.ali", "--ACDE", "GGACDE")
    return case


def test_a_chain_promod3_fails_on_fails_the_step(tmp_path):
    case = _case_with_alignments(tmp_path, ["A", "B"])

    def fake_fix(pdb_path, fasta_path, chain_id):
        if chain_id == "B":
            raise RuntimeError("no fragments for the gap")
        return {
            "model": object(),
            "residues_before": 4,
            "residues_after": 6,
            "gaps_initial": 1,
            "gaps_remaining": 0,
        }

    with patch.object(core, "fix_chain_residues", side_effect=fake_fix):
        result = core.process_fix_residues(
            pdb_path=str(tmp_path / "merge.pdb"),
            ali_dir=str(case),
            output_dir=str(case),
            case_name="case",
            chain_ids=["A", "B"],
            protein_chains=["A", "B"],
        )

    assert not result.success
    assert "Chain B: FAILED" in result.promod3_log
    assert "no fragments for the gap" in result.promod3_log
    assert not (case / "Merge" / "fixed.pdb").exists()


def test_a_protein_chain_without_an_alignment_fails_the_step(tmp_path):
    case = _case_with_alignments(tmp_path, ["A"])

    with patch.object(core, "fix_chain_residues") as fake_fix:
        fake_fix.return_value = {
            "model": object(),
            "residues_before": 4,
            "residues_after": 6,
            "gaps_initial": 1,
            "gaps_remaining": 0,
        }
        result = core.process_fix_residues(
            pdb_path=str(tmp_path / "merge.pdb"),
            ali_dir=str(case),
            output_dir=str(case),
            case_name="case",
            chain_ids=["A", "B"],
            protein_chains=["A", "B"],
        )

    assert not result.success
    assert "Chain B" in result.promod3_log
    assert not (case / "Merge" / "fixed.pdb").exists()
