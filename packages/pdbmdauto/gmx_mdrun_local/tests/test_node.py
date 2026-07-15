"""
Level 2 Tests: Node Integration (BoCoFlow Dependent)

These tests verify the BoCoFlow node wrapper in node.py, testing:
- Node metadata and OPTIONS definition
- Parameter extraction from flow_vars
- Integration with core functions
- Result formatting for BoCoFlow

Run with: pytest test_node.py -v

Note: These tests require bocoflow_core to be installed.
For pure algorithm testing without BoCoFlow, see test_core.py.
"""

import sys
from pathlib import Path

import pytest

# ========================================================================
# IMPORTANT: Check bocoflow_core availability BEFORE any other imports
# that might trigger node.py loading (which imports bocoflow_core)
# ========================================================================
try:
    from bocoflow_core.node import NodeException

    BOCOFLOW_AVAILABLE = True
except ImportError:
    BOCOFLOW_AVAILABLE = False

# Skip entire module at collection time if bocoflow_core not available
if not BOCOFLOW_AVAILABLE:
    pytest.skip(
        "bocoflow_core not installed - Level 2 tests require BoCoFlow. "
        "Install via: pip install bocoflow-core or add to pixi.toml",
        allow_module_level=True,
    )

# ========================================================================
# Now safe to import everything else (bocoflow_core is available)
# ========================================================================
import importlib.util
import json
import os
import tempfile
import types
from unittest.mock import MagicMock, patch

# Set up module hierarchy for relative imports
_module_dir = Path(__file__).parent.parent

# Create a fake parent package to enable relative imports
_fake_package_name = "gmx_mdrun_local_pkg"
_fake_package = types.ModuleType(_fake_package_name)
_fake_package.__path__ = [str(_module_dir)]
sys.modules[_fake_package_name] = _fake_package

# Import core module
_core_spec = importlib.util.spec_from_file_location(
    f"{_fake_package_name}.core", _module_dir / "core.py"
)
core = importlib.util.module_from_spec(_core_spec)
sys.modules[f"{_fake_package_name}.core"] = core
_core_spec.loader.exec_module(core)

# Patch the relative import in node.py by pre-loading the core module
# This is necessary because node.py uses "from .core import ..."
_fake_package.core = core

# Now import node module with patched relative import
_node_spec = importlib.util.spec_from_file_location(
    f"{_fake_package_name}.node", _module_dir / "node.py"
)
node_module = importlib.util.module_from_spec(_node_spec)
node_module.__package__ = _fake_package_name
sys.modules[f"{_fake_package_name}.node"] = node_module
_node_spec.loader.exec_module(node_module)

GmxMdRunLocal = node_module.GmxMdRunLocal


def create_mock_node_info():
    """Create a minimal node_info dictionary for testing."""
    return {
        "node_id": "test_node_001",
        "node_type": "simulation",
        "node_key": "GmxMdRunLocal",
        "name": "Test GROMACS Node",
    }


# =============================================================================
# Test: Node Metadata (Class-level attributes)
# =============================================================================


class TestNodeMetadata:
    """Test node class metadata and configuration."""

    def test_node_name(self):
        """Test node has correct name."""
        assert GmxMdRunLocal.name == "GROMACS MD Run (Local)"

    def test_node_key(self):
        """Test node has correct key."""
        assert GmxMdRunLocal.node_key == "GmxMdRunLocal"

    def test_node_category(self):
        """Test node category."""
        assert GmxMdRunLocal.category == "simulation"

    def test_node_tags(self):
        """Test node has appropriate tags."""
        assert "gromacs" in GmxMdRunLocal.tags
        assert "molecular-dynamics" in GmxMdRunLocal.tags
        assert "local" in GmxMdRunLocal.tags

    def test_node_ports(self):
        """Test node has correct number of ports."""
        assert GmxMdRunLocal.num_in == 1
        assert GmxMdRunLocal.num_out == 1


# =============================================================================
# Test: Node OPTIONS
# =============================================================================


class TestNodeOptions:
    """Test node OPTIONS definition."""

    def test_required_options_exist(self):
        """Test that all required options are defined."""
        required = [
            "case_name",
            "run_label",
            "input_top_file",
            "input_gro_file",
            "input_mdp_file",
            "input_ndx_file",
            "num_threads",
            "max_warnings",
            "verbose",
            "force_to_run",
        ]
        for opt in required:
            assert opt in GmxMdRunLocal.OPTIONS, f"Missing option: {opt}"

    def test_run_label_default(self):
        """Test run_label has correct default value."""
        run_label = GmxMdRunLocal.OPTIONS["run_label"]
        assert run_label.default == "md"

    def test_num_threads_default(self):
        """Test num_threads default is 0 (auto)."""
        num_threads = GmxMdRunLocal.OPTIONS["num_threads"]
        assert num_threads.default == 0

    def test_max_warnings_default(self):
        """Test max_warnings default is 10."""
        max_warnings = GmxMdRunLocal.OPTIONS["max_warnings"]
        assert max_warnings.default == 10

    def test_verbose_default(self):
        """Test verbose default is True."""
        verbose = GmxMdRunLocal.OPTIONS["verbose"]
        assert verbose.default is True

    def test_file_parameters_have_docstrings(self):
        """Test that file parameters have documentation."""
        file_params = ["input_top_file", "input_gro_file", "input_mdp_file"]
        for param_name in file_params:
            param = GmxMdRunLocal.OPTIONS[param_name]
            assert hasattr(param, "docstring")
            assert param.docstring is not None


# =============================================================================
# Test: Node Instance Creation
# =============================================================================


class TestNodeInstance:
    """Test node instance creation."""

    def test_create_instance(self):
        """Test creating a node instance."""
        node_info = create_mock_node_info()
        node = GmxMdRunLocal(node_info)
        assert node is not None

    def test_node_has_node_id(self):
        """Test node has node_id from node_info."""
        node_info = create_mock_node_info()
        node = GmxMdRunLocal(node_info)
        assert node.node_id == "test_node_001"


# =============================================================================
# Test: Execute Method - Parameter Extraction
# =============================================================================


class TestExecuteParameterExtraction:
    """Test parameter extraction in execute method."""

    @pytest.fixture
    def mock_flow_vars(self):
        """Create mock flow_vars dictionary."""

        def make_var(value):
            mock = MagicMock()
            mock.get_value.return_value = value
            return mock

        return {
            "case_name": make_var("test_case"),
            "run_label": make_var("nvt"),
            "input_top_file": make_var("/path/to/topol.top"),
            "input_gro_file": make_var("/path/to/conf.gro"),
            "input_mdp_file": make_var("/path/to/nvt.mdp"),
            "input_ndx_file": make_var(None),
            "num_threads": make_var(4),
            "max_warnings": make_var(5),
            "verbose": make_var(True),
            "force_to_run": make_var(False),
        }

    @patch.object(GmxMdRunLocal, "resolve_path", side_effect=lambda x: x)
    def test_execute_extracts_run_label(self, mock_resolve, mock_flow_vars):
        """Test that execute extracts run_label correctly."""
        # Patch the core functions at module level
        with patch.object(node_module, "check_gromacs_available", return_value=True), patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=True,
                message="OK",
                tpr_file=None,
                gro_file=None,
                xtc_file=None,
                edr_file=None,
                log_file=None,
            ),
        ) as mock_sim:
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)
            node.format_output_path = lambda x: x

            node.execute([], mock_flow_vars)

            # Verify run_md_simulation was called with correct run_label
            call_kwargs = mock_sim.call_args[1]
            assert call_kwargs["run_label"] == "nvt"

    @patch.object(GmxMdRunLocal, "resolve_path", side_effect=lambda x: x)
    def test_execute_extracts_num_threads(self, mock_resolve, mock_flow_vars):
        """Test that execute extracts num_threads correctly."""
        with patch.object(node_module, "check_gromacs_available", return_value=True), patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=True,
                message="OK",
                tpr_file=None,
                gro_file=None,
                xtc_file=None,
                edr_file=None,
                log_file=None,
            ),
        ) as mock_sim:
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)
            node.format_output_path = lambda x: x

            node.execute([], mock_flow_vars)

            call_kwargs = mock_sim.call_args[1]
            assert call_kwargs["num_threads"] == 4


# =============================================================================
# Test: Execute Method - Error Handling
# =============================================================================


class TestExecuteErrorHandling:
    """Test error handling in execute method."""

    @pytest.fixture
    def mock_flow_vars(self):
        """Create mock flow_vars dictionary."""

        def make_var(value):
            mock = MagicMock()
            mock.get_value.return_value = value
            return mock

        return {
            "case_name": make_var(""),
            "run_label": make_var("md"),
            "input_top_file": make_var("/path/to/topol.top"),
            "input_gro_file": make_var("/path/to/conf.gro"),
            "input_mdp_file": make_var("/path/to/md.mdp"),
            "input_ndx_file": make_var(None),
            "num_threads": make_var(0),
            "max_warnings": make_var(10),
            "verbose": make_var(True),
            "force_to_run": make_var(False),
        }

    def test_gromacs_not_available_raises_exception(self, mock_flow_vars):
        """Test that missing GROMACS raises NodeException."""
        with patch.object(node_module, "check_gromacs_available", return_value=False):
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)

            with pytest.raises(NodeException) as exc_info:
                node.execute([], mock_flow_vars)

            assert "GROMACS not found" in str(exc_info.value)

    @patch.object(GmxMdRunLocal, "resolve_path", side_effect=lambda x: x)
    def test_simulation_failure_raises_exception(self, mock_resolve, mock_flow_vars):
        """Test that simulation failure raises NodeException."""
        with patch.object(node_module, "check_gromacs_available", return_value=True), patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=False,
                message="grompp failed: Error in topology",
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)
            node.format_output_path = lambda x: x

            with pytest.raises(NodeException) as exc_info:
                node.execute([], mock_flow_vars)

            assert "grompp failed" in str(exc_info.value)


# =============================================================================
# Test: Execute Method - Result Formatting
# =============================================================================


class TestExecuteResultFormatting:
    """Test result formatting in execute method."""

    @pytest.fixture
    def mock_flow_vars(self):
        """Create mock flow_vars dictionary."""

        def make_var(value):
            mock = MagicMock()
            mock.get_value.return_value = value
            return mock

        return {
            "case_name": make_var("protein_sim"),
            "run_label": make_var("production"),
            "input_top_file": make_var("/work/topol.top"),
            "input_gro_file": make_var("/work/conf.gro"),
            "input_mdp_file": make_var("/work/md.mdp"),
            "input_ndx_file": make_var(None),
            "num_threads": make_var(8),
            "max_warnings": make_var(10),
            "verbose": make_var(True),
            "force_to_run": make_var(False),
        }

    @patch.object(GmxMdRunLocal, "resolve_path", side_effect=lambda x: x)
    def test_result_is_valid_json(self, mock_resolve, mock_flow_vars):
        """Test that execute returns valid JSON."""
        with patch.object(node_module, "check_gromacs_available", return_value=True), patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=True,
                message="Completed",
                tpr_file="/work/production.tpr",
                gro_file="/work/production.gro",
                xtc_file=None,
                edr_file=None,
                log_file="/work/production.log",
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)
            node.format_output_path = lambda x: f"abs:{x}"

            result_json = node.execute([], mock_flow_vars)

            # Should be valid JSON
            result = json.loads(result_json)
            assert isinstance(result, dict)

    @patch.object(GmxMdRunLocal, "resolve_path", side_effect=lambda x: x)
    def test_result_contains_required_fields(self, mock_resolve, mock_flow_vars):
        """Test that result contains required fields."""
        with patch.object(node_module, "check_gromacs_available", return_value=True), patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=True,
                message="Completed",
                tpr_file="/work/production.tpr",
                gro_file="/work/production.gro",
                xtc_file=None,
                edr_file=None,
                log_file=None,
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)
            node.format_output_path = lambda x: f"abs:{x}"

            result_json = node.execute([], mock_flow_vars)
            result = json.loads(result_json)

            # Check required fields
            assert "success" in result
            assert "message" in result
            assert "data" in result
            assert result["success"] is True

    @patch.object(GmxMdRunLocal, "resolve_path", side_effect=lambda x: x)
    def test_result_data_contains_case_name(self, mock_resolve, mock_flow_vars):
        """Test that result data contains case_name."""
        with patch.object(node_module, "check_gromacs_available", return_value=True), patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=True,
                message="OK",
                tpr_file=None,
                gro_file=None,
                xtc_file=None,
                edr_file=None,
                log_file=None,
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)
            node.format_output_path = lambda x: x

            result_json = node.execute([], mock_flow_vars)
            result = json.loads(result_json)

            assert result["data"]["case_name"] == "protein_sim"
            assert result["data"]["run_label"] == "production"


# =============================================================================
# Test: Predecessor Data Handling
# =============================================================================


class TestPredecessorDataHandling:
    """Test handling of predecessor_data."""

    @pytest.fixture
    def mock_flow_vars_no_case(self):
        """Create mock flow_vars with empty case_name."""

        def make_var(value):
            mock = MagicMock()
            mock.get_value.return_value = value
            return mock

        return {
            "case_name": make_var(""),  # Empty - should use predecessor
            "run_label": make_var("md"),
            "input_top_file": make_var("/path/to/topol.top"),
            "input_gro_file": make_var("/path/to/conf.gro"),
            "input_mdp_file": make_var("/path/to/md.mdp"),
            "input_ndx_file": make_var(None),
            "num_threads": make_var(0),
            "max_warnings": make_var(10),
            "verbose": make_var(True),
            "force_to_run": make_var(False),
        }

    @patch.object(GmxMdRunLocal, "resolve_path", side_effect=lambda x: x)
    def test_case_name_from_predecessor(self, mock_resolve, mock_flow_vars_no_case):
        """Test that case_name falls back to predecessor data."""
        with patch.object(node_module, "check_gromacs_available", return_value=True), patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=True,
                message="OK",
                tpr_file=None,
                gro_file=None,
                xtc_file=None,
                edr_file=None,
                log_file=None,
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)
            node.format_output_path = lambda x: x

            predecessor_data = [{"case_name": "from_predecessor"}]

            result_json = node.execute(predecessor_data, mock_flow_vars_no_case)
            result = json.loads(result_json)

            assert result["data"]["case_name"] == "from_predecessor"

    @patch.object(GmxMdRunLocal, "resolve_path", side_effect=lambda x: x)
    def test_default_case_name_when_no_predecessor(self, mock_resolve, mock_flow_vars_no_case):
        """Test default case_name when no predecessor data."""
        with patch.object(node_module, "check_gromacs_available", return_value=True), patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=True,
                message="OK",
                tpr_file=None,
                gro_file=None,
                xtc_file=None,
                edr_file=None,
                log_file=None,
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)
            node.format_output_path = lambda x: x

            result_json = node.execute([], mock_flow_vars_no_case)
            result = json.loads(result_json)

            assert result["data"]["case_name"] == "simulation"  # default


# =============================================================================
# Test: Integration with Core Functions
# =============================================================================


class TestCoreIntegration:
    """Test integration between node wrapper and core functions."""

    @pytest.fixture
    def mock_flow_vars(self):
        """Create mock flow_vars dictionary."""

        def make_var(value):
            mock = MagicMock()
            mock.get_value.return_value = value
            return mock

        return {
            "case_name": make_var("integration_test"),
            "run_label": make_var("em"),
            "input_top_file": make_var("/sim/system.top"),
            "input_gro_file": make_var("/sim/start.gro"),
            "input_mdp_file": make_var("/sim/em.mdp"),
            "input_ndx_file": make_var("/sim/index.ndx"),
            "num_threads": make_var(2),
            "max_warnings": make_var(3),
            "verbose": make_var(False),
            "force_to_run": make_var(False),
        }

    @patch.object(GmxMdRunLocal, "resolve_path", side_effect=lambda x: x)
    def test_core_function_called_with_correct_args(self, mock_resolve, mock_flow_vars):
        """Test that core function is called with all correct arguments."""
        with patch.object(node_module, "check_gromacs_available", return_value=True), patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=True,
                message="OK",
                tpr_file=None,
                gro_file=None,
                xtc_file=None,
                edr_file=None,
                log_file=None,
            ),
        ) as mock_sim:
            node_info = create_mock_node_info()
            node = GmxMdRunLocal(node_info)
            node.format_output_path = lambda x: x

            node.execute([], mock_flow_vars)

            # Verify run_md_simulation was called with correct arguments
            mock_sim.assert_called_once()
            call_kwargs = mock_sim.call_args[1]

            assert call_kwargs["top_file"] == "/sim/system.top"
            assert call_kwargs["gro_file"] == "/sim/start.gro"
            assert call_kwargs["mdp_file"] == "/sim/em.mdp"
            assert call_kwargs["ndx_file"] == "/sim/index.ndx"
            assert call_kwargs["run_label"] == "em"
            assert call_kwargs["num_threads"] == 2
            assert call_kwargs["max_warnings"] == 3
            assert call_kwargs["verbose"] is False
