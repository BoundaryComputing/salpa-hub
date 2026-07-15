"""
Tests for the PdbFastaBiopython node (BoCoFlow wrapper).

Requires bocoflow_core to be installed.

Run: pytest tests/test_node.py -v
"""

import json
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

# ---------------------------------------------------------------------------
# Check bocoflow_core availability
# ---------------------------------------------------------------------------
try:
    from bocoflow_core.node import Node, NodeException, NodeResult
    from bocoflow_core.parameters import (
        BooleanParameter,
        FileParameterEdit,
        FolderParameter,
        SelectParameter,
        StringParameter,
    )

    BOCOFLOW_AVAILABLE = True
except ImportError:
    BOCOFLOW_AVAILABLE = False

if not BOCOFLOW_AVAILABLE:
    pytest.skip(
        "bocoflow_core not installed — these tests require BoCoFlow. "
        "Install via: pip install -e ../../../bocoflow-core",
        allow_module_level=True,
    )

# ---------------------------------------------------------------------------
# Import the node class (use package import for proper relative imports)
# ---------------------------------------------------------------------------
from bocoflow_nodes.pdb_fasta_biopython.node import PdbFastaBiopython  # noqa: E402

# Demo data
_node_dir = Path(__file__).parent.parent
DEMO_DATA_DIR = _node_dir / "demo_data"
DEMO_PDB = DEMO_DATA_DIR / "3LZ0.pdb"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def node():
    """Create a PdbFastaBiopython node instance."""
    node_info = {
        "node_id": "test-pdb-fasta-001",
        "package_name": "pdb-fasta-biopython",
        "name": "PDB FASTA Parser",
    }
    return PdbFastaBiopython(node_info)


@pytest.fixture
def output_dir():
    """Provide a temporary output directory."""
    with tempfile.TemporaryDirectory(prefix="bocoflow-pdb-test-") as tmpdir:
        yield tmpdir


def make_flow_vars(
    input_mode="local_file",
    pdb_id="",
    pdb_file="",
    case_name="test-case",
    output_dir="/tmp",
    split_chains=True,
    include_hetatm=False,
    check_pdb_header=False,
):
    """Helper: create flow_vars dict mimicking server behavior."""
    fv = {
        "input_mode": SelectParameter(
            "Input Mode", options=["pdb_id", "local_file"], default="pdb_id"
        ),
        "pdb_id": StringParameter("PDB ID"),
        "pdb_file": FileParameterEdit("PDB File"),
        "case_name": StringParameter("Case Name"),
        "output_dir": FolderParameter("Output Directory"),
        "split_chains": BooleanParameter("Split by Chain", default=True),
        "include_hetatm": BooleanParameter("Include HETATM", default=False),
        "check_pdb_header": BooleanParameter("Check PDB Header", default=False),
        "force_to_run": BooleanParameter("Force to Run", default=False),
    }
    fv["input_mode"].set_value(input_mode)
    fv["pdb_id"].set_value(pdb_id)
    fv["pdb_file"].set_value(f"abs:{pdb_file}" if pdb_file else "")
    fv["case_name"].set_value(case_name)
    fv["output_dir"].set_value(f"abs:{output_dir}")
    fv["split_chains"].set_value(split_chains)
    fv["include_hetatm"].set_value(include_hetatm)
    fv["check_pdb_header"].set_value(check_pdb_header)
    return fv


# ---------------------------------------------------------------------------
# Tests: Node Metadata
# ---------------------------------------------------------------------------

class TestNodeMetadata:
    """Verify node class definition."""

    def test_inherits_from_node(self):
        assert issubclass(PdbFastaBiopython, Node)

    def test_options_defined(self):
        opts = PdbFastaBiopython.OPTIONS
        assert "input_mode" in opts
        assert "pdb_id" in opts
        assert "pdb_file" in opts
        assert "case_name" in opts
        assert "output_dir" in opts
        assert "split_chains" in opts
        assert "include_hetatm" in opts
        assert "check_pdb_header" in opts

    def test_force_to_run_not_in_options(self):
        assert "force_to_run" not in PdbFastaBiopython.OPTIONS

    def test_option_types(self):
        opts = PdbFastaBiopython.OPTIONS
        assert isinstance(opts["input_mode"], SelectParameter)
        assert isinstance(opts["pdb_id"], StringParameter)
        assert isinstance(opts["pdb_file"], FileParameterEdit)
        assert isinstance(opts["case_name"], StringParameter)
        assert isinstance(opts["output_dir"], FolderParameter)
        assert isinstance(opts["split_chains"], BooleanParameter)
        assert isinstance(opts["include_hetatm"], BooleanParameter)
        assert isinstance(opts["check_pdb_header"], BooleanParameter)

    def test_select_parameter_options(self):
        assert PdbFastaBiopython.OPTIONS["input_mode"].options == [
            "pdb_id", "local_file"
        ]

    def test_node_instantiation(self, node):
        assert node.name == "PDB FASTA Parser"
        assert node.node_id == "test-pdb-fasta-001"


# ---------------------------------------------------------------------------
# Tests: Local file execution
# ---------------------------------------------------------------------------

class TestLocalFileExecution:
    """Test execute() with local PDB file."""

    @pytest.fixture(autouse=True)
    def _check_demo(self):
        if not DEMO_PDB.exists():
            pytest.skip("Demo PDB file 3LZ0.pdb not found")

    def test_basic_execution(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file=str(DEMO_PDB),
            case_name="3lz0",
            output_dir=output_dir,
        )
        result_json = node.execute([], flow_vars)
        result = json.loads(result_json)

        assert result["success"] is True
        assert result["data"]["case_name"] == "3lz0"
        assert result["data"]["num_chains"] >= 8

    def test_split_chains_output(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file=str(DEMO_PDB),
            case_name="3lz0",
            output_dir=output_dir,
            split_chains=True,
        )
        result = json.loads(node.execute([], flow_vars))

        # Should have per-chain FASTA files in output
        output_files = result["files"]["output"]
        assert "chain_A" in output_files

        # Verify files exist
        for label, path in output_files.items():
            # Strip abs: prefix
            actual_path = path.replace("abs:", "")
            assert os.path.exists(actual_path), f"Missing: {actual_path}"

    def test_combined_fasta(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file=str(DEMO_PDB),
            case_name="3lz0",
            output_dir=output_dir,
            split_chains=False,
        )
        result = json.loads(node.execute([], flow_vars))

        assert "combined" in result["files"]["output"]

    def test_missing_residues_check(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file=str(DEMO_PDB),
            case_name="3lz0",
            output_dir=output_dir,
            check_pdb_header=True,
        )
        result = json.loads(node.execute([], flow_vars))

        assert result["success"] is True

    def test_chain_info_in_result(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file=str(DEMO_PDB),
            case_name="3lz0",
            output_dir=output_dir,
        )
        result = json.loads(node.execute([], flow_vars))

        chain_info = result["data"]["chain_info"]
        assert "A" in chain_info
        assert chain_info["A"]["type"] == "protein"
        assert chain_info["A"]["length"] > 0

    def test_include_hetatm(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file=str(DEMO_PDB),
            case_name="3lz0",
            output_dir=output_dir,
            include_hetatm=True,
        )
        result = json.loads(node.execute([], flow_vars))
        assert result["success"] is True

    def test_missing_file_raises(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file="/nonexistent/file.pdb",
            case_name="test",
            output_dir=output_dir,
        )
        with pytest.raises(Exception):
            node.execute([], flow_vars)


# ---------------------------------------------------------------------------
# Tests: Case name resolution
# ---------------------------------------------------------------------------

class TestCaseNameResolution:
    """Test case_name fallback logic."""

    @pytest.fixture(autouse=True)
    def _check_demo(self):
        if not DEMO_PDB.exists():
            pytest.skip("Demo PDB file 3LZ0.pdb not found")

    def test_explicit_case_name(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file=str(DEMO_PDB),
            case_name="my-protein",
            output_dir=output_dir,
        )
        result = json.loads(node.execute([], flow_vars))
        assert result["data"]["case_name"] == "my-protein"

    def test_empty_case_name_local_defaults(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file=str(DEMO_PDB),
            case_name="",
            output_dir=output_dir,
        )
        result = json.loads(node.execute([], flow_vars))
        # Should default to "protein"
        assert result["data"]["case_name"] == "protein"

    def test_predecessor_case_name(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="local_file",
            pdb_file=str(DEMO_PDB),
            case_name="",
            output_dir=output_dir,
        )
        predecessor_data = [{"case_name": "upstream-case"}]
        result = json.loads(node.execute(predecessor_data, flow_vars))
        assert result["data"]["case_name"] == "upstream-case"


# ---------------------------------------------------------------------------
# Tests: PDB ID mode (mocked API)
# ---------------------------------------------------------------------------

class TestPdbIdMode:
    """Test pdb_id input mode with mocked RCSB API."""

    @patch("bocoflow_nodes.pdb_fasta_biopython.core.requests.get")
    def test_pdb_id_execution(self, mock_get, node, output_dir):
        if not DEMO_PDB.exists():
            pytest.skip("Demo PDB not found")

        # Mock API to return demo PDB content
        pdb_content = DEMO_PDB.read_text()
        fasta_content = ">3LZ0_1|Chain A|Protein\nMKWVTF\n"

        def side_effect(url, **kwargs):
            resp = type("Response", (), {"status_code": 200})()
            if "download" in url:
                resp.text = pdb_content
            else:
                resp.text = fasta_content
            return resp

        mock_get.side_effect = side_effect

        flow_vars = make_flow_vars(
            input_mode="pdb_id",
            pdb_id="3LZ0",
            case_name="",
            output_dir=output_dir,
        )
        result = json.loads(node.execute([], flow_vars))

        assert result["success"] is True
        assert result["data"]["case_name"] == "3lz0"
        assert result["data"]["input_mode"] == "pdb_id"
        assert result["data"]["num_chains"] >= 8

        # PDB file should be saved
        pdb_path = os.path.join(output_dir, "3LZ0.pdb")
        assert os.path.exists(pdb_path)

    def test_pdb_id_missing_raises(self, node, output_dir):
        flow_vars = make_flow_vars(
            input_mode="pdb_id",
            pdb_id="",
            output_dir=output_dir,
        )
        with pytest.raises(Exception):
            node.execute([], flow_vars)
