"""
Tests for the Original Atom Groups.

OriHeavy and OriBackBone must hold the crystal residues and leave out the ones Fix Missing
Residues rebuilt, because step 8 freezes them while the rebuilt parts relax. Up to 1.2.5 the
builder matched the deposited numbers of missing_residues_chain_*.csv against a model that
ProMod3 numbers from 1, and matched by residue number alone in a .gro. On 4Z8J that left out
crystal residues A33-A37 and froze all six rebuilt ones.

The fixture is a small copy of that situation. Chain A was deposited from residue 7 with a
two-residue tag (5-6) missing; chain B was deposited from 587 with 586 missing. ProMod3's model
numbers them A 1-12 and B 1-3, so the rebuilt residues are A1, A2 and B1, while the CSVs still
say A5, A6 and B586.

gen_gmx_ndx carries the same code, so every test runs against both nodes.
"""

import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

import pytest

PACKAGE = Path(__file__).resolve().parent.parent.parent


def _load(node):
    spec = importlib.util.spec_from_file_location(
        f"{node}_core", PACKAGE / node / "core.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(params=["ori_ndx_builder", "gen_gmx_ndx"])
def core(request):
    return _load(request.param)


# The model as Fix Missing Residues writes it: chain order, residue numbers, rebuilt ones.
MODEL = [
    ("A", ["GLY", "SER"] + ["ALA"] * 10, {1, 2}),
    ("B", ["GLN", "LEU", "VAL"], {1}),
]
REBUILT = {("A", 1), ("A", 2), ("B", 1)}
CRYSTAL = {(c, r) for c, names, _ in MODEL for r in range(1, len(names) + 1)} - REBUILT


def _atoms_of(resname, last):
    """Atom names of one residue: backbone, one hydrogen, a CB unless glycine."""
    names = ["N", "H", "CA"] + ([] if resname == "GLY" else ["CB"]) + ["C"]
    return names + (["OC1", "OC2"] if last else ["O"])


def _residues(model, renumber_continuously=False):
    """(chain, resid, resname, atom names) in file order."""
    out, offset = [], 0
    for chain, names, _ in model:
        for i, resname in enumerate(names, start=1):
            resid = i + offset if renumber_continuously else i
            out.append((chain, resid, resname, _atoms_of(resname, i == len(names))))
        if renumber_continuously:
            offset += len(names)
    return out


def write_gro(path, model, renumber_continuously=False, solvent=2):
    lines, index = [], 0
    for _, resid, resname, names in _residues(model, renumber_continuously):
        for name in names:
            index += 1
            lines.append(
                f"{resid:5d}{resname:<5s}{name:>5s}{index:5d}   0.000   0.000   0.000\n"
            )
    last_resid = resid
    for w in range(solvent):
        for name in ("OW", "HW1", "HW2"):
            index += 1
            lines.append(
                f"{last_resid + 1 + w:5d}{'SOL':<5s}{name:>5s}{index:5d}   0.000   0.000   0.000\n"
            )
    path.write_text(
        f"fixture\n{index:5d}\n" + "".join(lines) + "   5.00000   5.00000   5.00000\n"
    )
    return path


def write_pdb(path, residues):
    lines, index = [], 0
    for chain, resid, resname, names in residues:
        for name in names:
            index += 1
            element = name[0]
            lines.append(
                f"ATOM  {index:5d} {name:<4s} {resname:>3s} {chain}{resid:4d}    "
                f"   0.000   0.000   0.000  1.00  0.00          {element:>2s}\n"
            )
    path.write_text("".join(lines) + "END\n")
    return path


def write_csvs(case_dir):
    """missing_residues_chain_*.csv in the DEPOSITED numbering, as Generate Alignment writes."""
    (case_dir / "missing_residues_chain_A.csv").write_text(
        "chain,res_name,one_letter,ssseq,model,insertion\n"
        "A,GLY,G,5,None,\nA,SER,S,6,None,\n"
    )
    (case_dir / "missing_residues_chain_B.csv").write_text(
        "chain,res_name,one_letter,ssseq,model,insertion\nB,GLN,Q,586,None,\n"
    )


def write_map(path, model):
    chains = [
        {
            "chain": c,
            "residues": list(range(1, len(names) + 1)),
            "rebuilt": sorted(rebuilt),
        }
        for c, names, rebuilt in model
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"numbering": "model", "chains": chains}))
    return path


def run(core, structure, ndx, case_dir, rebuilt_map=""):
    """generate_ori_ndx with `gmx make_ndx` replaced by an empty base index."""

    def fake_make_ndx(cmd, **kwargs):
        Path(cmd[cmd.index("-o") + 1]).write_text("[ System ]\n1\n")

        class Done:
            returncode = 0

        return Done()

    with patch.object(core.subprocess, "run", side_effect=fake_make_ndx):
        return core.generate_ori_ndx(
            structure_path=str(structure),
            ndx_path=str(ndx),
            missing_csv_dir=str(case_dir),
            rebuilt_map_path=str(rebuilt_map) if rebuilt_map else "",
        )


def groups(ndx):
    out, current = {}, None
    for line in Path(ndx).read_text().splitlines():
        if line.startswith("["):
            current = line.strip("[] ")
            out[current] = []
        elif current and line.strip():
            out[current] += [int(x) for x in line.split()]
    return out


def residues_of(indices, structure_atoms):
    by_index = {a["index"]: a for a in structure_atoms}
    return {(by_index[i]["chain"] or "?", by_index[i]["resid"]) for i in indices}


def gro_chain_of(core, gro):
    """Label the .gro's atoms with the chain each belongs to, for readable assertions."""
    atoms = core.read_structure_atoms(str(gro))
    chain_by_position = [c for c, names, _ in MODEL for _ in names]
    residue, previous = -1, None
    for a in atoms:
        if a["resname"] == "SOL":
            continue
        if a["resid"] != previous:
            residue, previous = residue + 1, a["resid"]
        a["chain"] = chain_by_position[residue]
    return atoms


def test_gro_groups_leave_out_exactly_the_rebuilt_residues(core, tmp_path):
    case = tmp_path / "case"
    (case / "gmx").mkdir(parents=True)
    write_csvs(case)
    gro = write_gro(case / "gmx" / "pdb2gmx.gro", MODEL)
    record = write_map(case / "Merge" / "rebuilt_residues.json", MODEL)

    result = run(core, gro, case / "gmx" / "index.ndx", case, record)

    assert result.success, result.log
    atoms = gro_chain_of(core, gro)
    g = groups(case / "gmx" / "index.ndx")
    assert residues_of(g["OriHeavy"], atoms) == CRYSTAL
    assert residues_of(g["OriBackBone"], atoms) == CRYSTAL


def test_backbone_is_n_ca_c_o_and_c_terminal_oxygens_are_not_in_it(core, tmp_path):
    # 4Z8J's 414 = 102 residues x 4 + 2 C-termini x 3: OC1 and OC2 are not backbone atoms.
    case = tmp_path / "case"
    (case / "gmx").mkdir(parents=True)
    gro = write_gro(case / "gmx" / "pdb2gmx.gro", MODEL)
    record = write_map(case / "Merge" / "rebuilt_residues.json", MODEL)

    result = run(core, gro, case / "gmx" / "index.ndx", case, record)

    assert result.success, result.log
    # 10 crystal residues in A (A12 a C-terminus), 2 in B (B3 a C-terminus).
    assert result.n_ori_backbone == (10 * 4 - 1) + (2 * 4 - 1)


def test_the_record_is_found_in_the_case_folder_without_being_named(core, tmp_path):
    case = tmp_path / "case"
    (case / "gmx").mkdir(parents=True)
    write_map(case / "Merge" / "rebuilt_residues.json", MODEL)
    assert core.locate_rebuilt_map("", str(case)) == str(
        case / "Merge" / "rebuilt_residues.json"
    )
    assert core.locate_rebuilt_map("", str(tmp_path / "elsewhere")) == ""
    named = write_map(tmp_path / "mine.json", MODEL)
    assert core.locate_rebuilt_map(str(named), str(case)) == str(named)


def test_a_structure_that_is_not_the_recorded_model_fails_instead_of_guessing(
    core, tmp_path
):
    # pdb2gmx renumbering the residues continuously would put B1 at 13: matching the record
    # by number would then free a rebuilt residue. The builder must refuse, not guess.
    case = tmp_path / "case"
    (case / "gmx").mkdir(parents=True)
    gro = write_gro(case / "gmx" / "pdb2gmx.gro", MODEL, renumber_continuously=True)
    record = write_map(case / "Merge" / "rebuilt_residues.json", MODEL)

    result = run(core, gro, case / "gmx" / "index.ndx", case, record)

    assert not result.success
    assert "Fix Missing Residues" in result.log


def test_a_structure_with_a_residue_missing_fails(core, tmp_path):
    case = tmp_path / "case"
    (case / "gmx").mkdir(parents=True)
    shorter = [MODEL[0], ("B", ["GLN", "LEU"], {1})]
    gro = write_gro(case / "gmx" / "pdb2gmx.gro", shorter)
    record = write_map(case / "Merge" / "rebuilt_residues.json", MODEL)

    result = run(core, gro, case / "gmx" / "index.ndx", case, record)

    assert not result.success
    assert "residues" in result.log


def test_a_gro_without_the_record_fails_and_says_what_is_missing(core, tmp_path):
    # A .gro has no chain IDs and carries ProMod3's numbers, so the deposited numbers in the
    # CSVs cannot place the rebuilt residues. 1.2.5 guessed, and guessed wrong.
    case = tmp_path / "case"
    (case / "gmx").mkdir(parents=True)
    write_csvs(case)
    gro = write_gro(case / "gmx" / "pdb2gmx.gro", MODEL)

    result = run(core, gro, case / "gmx" / "index.ndx", case)

    assert not result.success
    assert "rebuilt_residues.json" in result.log


def test_a_pdb_in_the_deposited_numbering_still_works_from_the_csvs(core, tmp_path):
    # Without a record, a PDB that keeps the deposited numbers and its chain IDs is matched
    # against the CSVs per chain, as before. Chain B's missing 586 must not touch chain A.
    case = tmp_path / "case"
    case.mkdir()
    write_csvs(case)
    deposited = [
        ("A", resid, "ALA", ["N", "CA", "CB", "C", "O"]) for resid in range(5, 17)
    ] + [("B", resid, "LEU", ["N", "CA", "CB", "C", "O"]) for resid in (586, 587, 588)]
    pdb = write_pdb(case / "structure.pdb", deposited)

    result = run(core, pdb, case / "index.ndx", case)

    assert result.success, result.log
    atoms = core.read_structure_atoms(str(pdb))
    expected = {(c, r) for c, r, _, _ in deposited} - {("A", 5), ("A", 6), ("B", 586)}
    assert residues_of(groups(case / "index.ndx")["OriBackBone"], atoms) == expected
