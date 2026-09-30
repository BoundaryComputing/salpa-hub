"""What the Generate Alignment node says about a chain it cannot build, run through BoCoFlow.

A short one (6LU7's N3 inhibitor, chain C) is set aside with a warning and named in the step's
message; a longer one stops the step with the chain, the residues and what to do. See
test_unbuildable.py for the rule itself.

Run from the package root: pixi run test
"""

import json
import shutil
import sys
from pathlib import Path

import pytest

try:
    from bocoflow_core.node import NodeException
    from bocoflow_core.parameters import BooleanParameter, FolderParameter, StringParameter
except ImportError:
    pytest.skip(
        "bocoflow_core not installed: these tests run the node through BoCoFlow",
        allow_module_level=True,
    )

_pkg_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_pkg_root))

import bocoflow_core.stream_logger as stream_logger  # noqa: E402

from gen_ali.node import GenAli  # noqa: E402

EXCERPT_6LU7 = Path(__file__).parent / "data" / "6lu7_excerpt.pdb"


def _run(folder, monkeypatch):
    streamed = []
    monkeypatch.setattr(stream_logger, "stream_log", lambda msg, **kw: streamed.append((msg, kw)))
    node = GenAli({"node_id": "test-gen-ali", "package_name": "pdbmdauto", "name": "Gen Ali"})
    flow_vars = {
        "case_name": StringParameter("Case Name"),
        "output_dir": FolderParameter("Output Directory"),
        "append_end_in_seq": BooleanParameter("Append End Marker", default=True),
        "force_to_run": BooleanParameter("Force to Run", default=False),
    }
    flow_vars["case_name"].set_value("case")
    flow_vars["output_dir"].set_value(f"abs:{folder}")
    flow_vars["append_end_in_seq"].set_value(True)
    upstream = {"case_name": "case", "working_path": f"abs:{folder}", "chain_info": {}}
    return json.loads(node.execute([upstream], flow_vars)), streamed


def test_a_short_chain_is_set_aside_with_a_warning(tmp_path, monkeypatch):
    shutil.copy(EXCERPT_6LU7, tmp_path / "6LU7.pdb")
    result, streamed = _run(tmp_path, monkeypatch)
    assert result["success"] is True
    assert result["data"]["pdb_chain_list"] == ["A"]
    assert result["data"]["set_aside_chains"] == {
        "C": {"length": 6, "residues": ["02J", "PJE", "010"]}
    }
    assert "Set aside as ligands: chain(s) C." in result["message"]
    warnings = [m for m, kw in streamed if kw.get("level") == "warning"]
    assert any("Chain C (6 residues) is set aside as a ligand" in m and "02J, PJE, 010" in m
               for m in warnings)


def test_a_longer_chain_stops_the_step_and_says_what_to_do(tmp_path, monkeypatch):
    names = ["ALA"] * 31
    names[10] = "DAL"
    seqres = [f"SEQRES {i // 13 + 1:>3} B {len(names):>4}  " + " ".join(names[i:i + 13])
              for i in range(0, len(names), 13)]
    atom = "ATOM      1  CA  ALA B   1       1.000   0.000   0.000  1.00 10.00           C"
    (tmp_path / "1ABC.pdb").write_text("\n".join(seqres + [atom, "END"]) + "\n")
    with pytest.raises(NodeException) as raised:
        _run(tmp_path, monkeypatch)
    text = str(raised.value)
    assert "chain B (31 residues): DAL at 11" in text
    assert "local_file" in text
