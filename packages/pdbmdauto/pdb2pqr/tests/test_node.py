"""
Tests for the PDB2PQR node (BoCoFlow wrapper).

Requires bocoflow_core to be installed.

Run: pytest tests/test_node.py -v
"""

import json
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Check bocoflow_core availability
# ---------------------------------------------------------------------------
try:
    from bocoflow_core.node import Node, NodeException, NodeResult
    from bocoflow_core.parameters import (
        BooleanParameter,
        FileParameterEdit,
        FloatParameter,
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
from bocoflow_nodes.pdb2pqr.node import PDB2PQR  # noqa: E402

# Demo data
_node_dir = Path(__file__).parent.parent
DEMO_DATA_DIR = _node_dir / "demo_data"
DEMO_PDB = DEMO_DATA_DIR / "mini.pdb"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def node():
    """Create a PDB2PQR node instance."""
    node_info = {
        "node_id": "test-pdb2pqr-001",
        "package_name": "pdb2pqr",
        "name": "PDB to PQR Converter",
    }
    return PDB2PQR(node_info)


@pytest.fixture
def output_dir():
    """Provide a temporary output directory."""
    with tempfile.TemporaryDirectory(prefix="bocoflow-pdb2pqr-test-") as tmpdir:
        yield tmpdir


def make_flow_vars(
    case_name="test-case",
    input_pdb="",
    output_dir="/tmp",
    force_field="AMBER",
    ph=7.0,
    keep_chain=True,
    optimize_hydrogens=True,
    include_header=True,
    use_propka=True,
    log_level="INFO",
    custom_pdb2pqr_path="",
    generate_pdb=True,
):
    """Helper: create flow_vars dict mimicking server behavior."""
    fv = {
        "case_name": StringParameter("Case Name"),
        "input_pdb": FileParameterEdit("Input PDB File"),
        "output_dir": FolderParameter("Output Directory"),
        "force_field": SelectParameter(
            "Force Field",
            default="AMBER",
            options=["AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON"],
        ),
        "ph": FloatParameter("pH Value", default=7.0),
        "keep_chain": BooleanParameter("Keep Chain IDs", default=True),
        "optimize_hydrogens": BooleanParameter("Optimize Hydrogens", default=True),
        "include_header": BooleanParameter("Include Header", default=True),
        "use_propka": BooleanParameter("Use PROPKA", default=True),
        "log_level": SelectParameter(
            "Log Level",
            default="INFO",
            options=["DEBUG", "INFO", "WARNING", "ERROR"],
        ),
        "custom_pdb2pqr_path": StringParameter("Custom PDB2PQR Path", default=""),
        "generate_pdb": BooleanParameter("Generate PDB from PQR", default=True),
        "force_to_run": BooleanParameter("Force to Run", default=False),
    }
    fv["case_name"].set_value(case_name)
    fv["input_pdb"].set_value(f"abs:{input_pdb}" if input_pdb else "")
    fv["output_dir"].set_value(f"abs:{output_dir}")
    fv["force_field"].set_value(force_field)
    fv["ph"].set_value(ph)
    fv["keep_chain"].set_value(keep_chain)
    fv["optimize_hydrogens"].set_value(optimize_hydrogens)
    fv["include_header"].set_value(include_header)
    fv["use_propka"].set_value(use_propka)
    fv["log_level"].set_value(log_level)
    fv["custom_pdb2pqr_path"].set_value(custom_pdb2pqr_path)
    fv["generate_pdb"].set_value(generate_pdb)
    return fv


# ---------------------------------------------------------------------------
# Tests: Node Metadata
# ---------------------------------------------------------------------------


class TestNodeMetadata:
    """Verify node class definition."""

    def test_inherits_from_node(self):
        assert issubclass(PDB2PQR, Node)

    def test_class_attributes(self):
        assert PDB2PQR.name == "PDB to PQR Converter"
        assert PDB2PQR.node_key == "PDB2PQR"
        assert PDB2PQR.category == "io"

    def test_options_present(self):
        opts = PDB2PQR.OPTIONS
        expected_keys = [
            "case_name",
            "input_pdb",
            "output_dir",
            "force_field",
            "ph",
            "keep_chain",
            "optimize_hydrogens",
            "include_header",
            "use_propka",
            "log_level",
            "custom_pdb2pqr_path",
            "generate_pdb",
        ]
        for key in expected_keys:
            assert key in opts, f"Missing option: {key}"

    def test_options_count(self):
        assert len(PDB2PQR.OPTIONS) == 12

    def test_force_to_run_not_in_options(self):
        assert "force_to_run" not in PDB2PQR.OPTIONS

    def test_option_types(self):
        opts = PDB2PQR.OPTIONS
        assert isinstance(opts["case_name"], StringParameter)
        assert isinstance(opts["input_pdb"], FileParameterEdit)
        assert isinstance(opts["output_dir"], FolderParameter)
        assert isinstance(opts["force_field"], SelectParameter)
        assert isinstance(opts["ph"], FloatParameter)
        assert isinstance(opts["keep_chain"], BooleanParameter)
        assert isinstance(opts["optimize_hydrogens"], BooleanParameter)
        assert isinstance(opts["include_header"], BooleanParameter)
        assert isinstance(opts["use_propka"], BooleanParameter)
        assert isinstance(opts["log_level"], SelectParameter)
        assert isinstance(opts["custom_pdb2pqr_path"], StringParameter)
        assert isinstance(opts["generate_pdb"], BooleanParameter)

    def test_option_defaults(self):
        opts = PDB2PQR.OPTIONS
        assert opts["force_field"].default == "AMBER"
        assert opts["ph"].default == 7.0
        assert opts["keep_chain"].default is True
        assert opts["optimize_hydrogens"].default is True
        assert opts["use_propka"].default is True
        assert opts["generate_pdb"].default is True
        assert opts["log_level"].default == "INFO"

    def test_force_field_options(self):
        ff = PDB2PQR.OPTIONS["force_field"]
        assert "AMBER" in ff.options
        assert "CHARMM" in ff.options
        assert len(ff.options) == 6

    def test_node_instantiation(self, node):
        assert node.name == "PDB to PQR Converter"
        assert node.node_id == "test-pdb2pqr-001"


# ---------------------------------------------------------------------------
# Tests: Execution (mocked core functions)
# ---------------------------------------------------------------------------


class TestExecution:
    """Test execute() with mocked core functions."""

    @patch("bocoflow_nodes.pdb2pqr.node.convert_pqr_to_pdb")
    @patch("bocoflow_nodes.pdb2pqr.node.extract_pqr_statistics")
    @patch("bocoflow_nodes.pdb2pqr.node.run_pdb2pqr")
    @patch("bocoflow_nodes.pdb2pqr.node.find_pdb2pqr_executable")
    @patch("bocoflow_nodes.pdb2pqr.node.stream_log")
    def test_execute_success(
        self,
        mock_stream_log,
        mock_find_exe,
        mock_run,
        mock_stats,
        mock_convert,
        node,
        output_dir,
    ):
        mock_find_exe.return_value = "/usr/bin/pdb2pqr"
        mock_run.return_value = ("Success\n", 0)
        mock_stats.return_value = {
            "total_atoms": 100,
            "hydrogen_atoms": 40,
            "non_hydrogen_atoms": 60,
        }

        # Create fake PQR output file so the check passes
        pqr_path = os.path.join(output_dir, "mycase_structure.pqr")
        with open(pqr_path, "w") as f:
            f.write("ATOM      1  N   ALA A   1\nEND\n")

        # Create fake PDB output file
        pdb_path = os.path.join(output_dir, "mycase_protonated.pdb")
        with open(pdb_path, "w") as f:
            f.write("ATOM      1  N   ALA A   1\nEND\n")

        # Need a real input file too
        input_pdb = os.path.join(output_dir, "input.pdb")
        with open(input_pdb, "w") as f:
            f.write("ATOM      1  N   ALA A   1\nEND\n")

        flow_vars = make_flow_vars(
            case_name="mycase",
            input_pdb=input_pdb,
            output_dir=output_dir,
        )

        result_json = node.execute([], flow_vars)
        result = json.loads(result_json)

        assert result["success"] is True
        assert result["data"]["case_name"] == "mycase"
        assert result["data"]["statistics"]["total_atoms"] == 100
        assert result["data"]["statistics"]["hydrogen_atoms"] == 40
        mock_find_exe.assert_called_once()
        mock_run.assert_called_once()
        mock_stats.assert_called_once()
        mock_convert.assert_called_once()

    @patch("bocoflow_nodes.pdb2pqr.node.stream_log")
    def test_execute_missing_input(self, mock_stream_log, node, output_dir):
        flow_vars = make_flow_vars(
            case_name="test",
            input_pdb="/nonexistent/file.pdb",
            output_dir=output_dir,
        )
        with pytest.raises(NodeException):
            node.execute([], flow_vars)

    @patch("bocoflow_nodes.pdb2pqr.node.find_pdb2pqr_executable")
    @patch("bocoflow_nodes.pdb2pqr.node.run_pdb2pqr")
    @patch("bocoflow_nodes.pdb2pqr.node.stream_log")
    def test_execute_pdb2pqr_failure(
        self, mock_stream_log, mock_run, mock_find_exe, node, output_dir
    ):
        mock_find_exe.return_value = "pdb2pqr"
        mock_run.return_value = ("Error: bad input\n", 1)

        input_pdb = os.path.join(output_dir, "input.pdb")
        with open(input_pdb, "w") as f:
            f.write("ATOM      1  N   ALA A   1\nEND\n")

        flow_vars = make_flow_vars(
            case_name="test",
            input_pdb=input_pdb,
            output_dir=output_dir,
        )
        with pytest.raises(NodeException, match="pdb2pqr failed"):
            node.execute([], flow_vars)

    @patch("bocoflow_nodes.pdb2pqr.node.extract_pqr_statistics")
    @patch("bocoflow_nodes.pdb2pqr.node.run_pdb2pqr")
    @patch("bocoflow_nodes.pdb2pqr.node.find_pdb2pqr_executable")
    @patch("bocoflow_nodes.pdb2pqr.node.stream_log")
    def test_execute_without_generate_pdb(
        self,
        mock_stream_log,
        mock_find_exe,
        mock_run,
        mock_stats,
        node,
        output_dir,
    ):
        mock_find_exe.return_value = "pdb2pqr"
        mock_run.return_value = ("Success\n", 0)
        mock_stats.return_value = {
            "total_atoms": 50,
            "hydrogen_atoms": 20,
            "non_hydrogen_atoms": 30,
        }

        # Create fake PQR output
        pqr_path = os.path.join(output_dir, "noconv_structure.pqr")
        with open(pqr_path, "w") as f:
            f.write("ATOM      1  N   ALA A   1\nEND\n")

        input_pdb = os.path.join(output_dir, "input.pdb")
        with open(input_pdb, "w") as f:
            f.write("ATOM      1  N   ALA A   1\nEND\n")

        flow_vars = make_flow_vars(
            case_name="noconv",
            input_pdb=input_pdb,
            output_dir=output_dir,
            generate_pdb=False,
        )

        result = json.loads(node.execute([], flow_vars))

        assert result["success"] is True
        assert "protonated_pdb" not in result["data"]["output_files"]
        assert "protonated_pdb" not in result["files"]["output"]

    @patch("bocoflow_nodes.pdb2pqr.node.convert_pqr_to_pdb")
    @patch("bocoflow_nodes.pdb2pqr.node.extract_pqr_statistics")
    @patch("bocoflow_nodes.pdb2pqr.node.run_pdb2pqr")
    @patch("bocoflow_nodes.pdb2pqr.node.find_pdb2pqr_executable")
    @patch("bocoflow_nodes.pdb2pqr.node.stream_log")
    def test_execute_with_propka(
        self,
        mock_stream_log,
        mock_find_exe,
        mock_run,
        mock_stats,
        mock_convert,
        node,
        output_dir,
    ):
        mock_find_exe.return_value = "pdb2pqr"
        mock_run.return_value = ("Success\n", 0)
        mock_stats.return_value = {
            "total_atoms": 50,
            "hydrogen_atoms": 20,
            "non_hydrogen_atoms": 30,
        }

        # Create fake outputs
        pqr_path = os.path.join(output_dir, "pk_structure.pqr")
        pdb_path = os.path.join(output_dir, "pk_protonated.pdb")
        propka_path = os.path.join(output_dir, "pk_propka.out")
        for path in [pqr_path, pdb_path, propka_path]:
            with open(path, "w") as f:
                f.write("fake content\n")

        input_pdb = os.path.join(output_dir, "input.pdb")
        with open(input_pdb, "w") as f:
            f.write("ATOM      1  N   ALA A   1\nEND\n")

        flow_vars = make_flow_vars(
            case_name="pk",
            input_pdb=input_pdb,
            output_dir=output_dir,
            use_propka=True,
        )

        result = json.loads(node.execute([], flow_vars))

        assert result["success"] is True
        assert "propka" in result["data"]["output_files"]
        assert "propka" in result["files"]["output"]


# ---------------------------------------------------------------------------
# Tests: Case name resolution
# ---------------------------------------------------------------------------


class TestCaseNameResolution:
    """Test case_name fallback logic."""

    @patch("bocoflow_nodes.pdb2pqr.node.convert_pqr_to_pdb")
    @patch("bocoflow_nodes.pdb2pqr.node.extract_pqr_statistics")
    @patch("bocoflow_nodes.pdb2pqr.node.run_pdb2pqr")
    @patch("bocoflow_nodes.pdb2pqr.node.find_pdb2pqr_executable")
    @patch("bocoflow_nodes.pdb2pqr.node.stream_log")
    def test_case_name_from_predecessor(
        self,
        mock_stream_log,
        mock_find_exe,
        mock_run,
        mock_stats,
        mock_convert,
        node,
        output_dir,
    ):
        mock_find_exe.return_value = "pdb2pqr"
        mock_run.return_value = ("Success\n", 0)
        mock_stats.return_value = {
            "total_atoms": 10,
            "hydrogen_atoms": 5,
            "non_hydrogen_atoms": 5,
        }

        # Create fake output files using the predecessor case_name
        pqr_path = os.path.join(output_dir, "upstream_structure.pqr")
        pdb_path = os.path.join(output_dir, "upstream_protonated.pdb")
        for path in [pqr_path, pdb_path]:
            with open(path, "w") as f:
                f.write("ATOM      1  N   ALA A   1\nEND\n")

        input_pdb = os.path.join(output_dir, "input.pdb")
        with open(input_pdb, "w") as f:
            f.write("ATOM      1  N   ALA A   1\nEND\n")

        flow_vars = make_flow_vars(
            case_name="",
            input_pdb=input_pdb,
            output_dir=output_dir,
        )

        predecessor_data = [{"case_name": "upstream"}]
        result = json.loads(node.execute(predecessor_data, flow_vars))

        assert result["data"]["case_name"] == "upstream"
