"""
Level 2 Tests: Node Integration (BoCoFlow Dependent)

These tests verify the BoCoFlow node wrapper in node.py, testing:
- Node metadata and OPTIONS definition
- HPCNodeBase integration and abstract method implementations
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
    from bocoflow_core.hpc_node import HPCNodeBase
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
_fake_package_name = "gmx_mdrun_pkg"
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

GmxMdRun = node_module.GmxMdRun


def create_mock_node_info():
    """Create a minimal node_info dictionary for testing."""
    return {
        "node_id": "test_node_001",
        "node_type": "simulation",
        "node_key": "GmxMdRun",
        "name": "Test GROMACS Node",
    }


# =============================================================================
# Test: Node Metadata (Class-level attributes)
# =============================================================================


class TestNodeMetadata:
    """Test node class metadata and configuration."""

    def test_node_name(self):
        """Test node has correct name."""
        assert GmxMdRun.name == "GROMACS MD Run"

    def test_node_key(self):
        """Test node has correct key."""
        assert GmxMdRun.node_key == "GmxMdRun"

    def test_node_category(self):
        """Test node category."""
        assert GmxMdRun.category == "simulation"

    def test_node_tags(self):
        """Test node has appropriate tags."""
        assert "gromacs" in GmxMdRun.tags
        assert "molecular-dynamics" in GmxMdRun.tags
        assert "hpc" in GmxMdRun.tags
        assert "slurm" in GmxMdRun.tags

    def test_node_ports(self):
        """Test node has correct number of ports."""
        assert GmxMdRun.num_in == 1
        assert GmxMdRun.num_out == 1

    def test_inherits_from_hpcnodebase(self):
        """Test node inherits from HPCNodeBase."""
        assert issubclass(GmxMdRun, HPCNodeBase)


# =============================================================================
# Test: Node OPTIONS
# =============================================================================


class TestNodeOptions:
    """Test node OPTIONS definition."""

    def test_hpc_options_included(self):
        """Test that HPC options from HPCNodeBase are included."""
        hpc_options = ["execution_mode", "hpc_profile", "slurm_script", "force_resubmit"]
        for opt in hpc_options:
            assert opt in GmxMdRun.OPTIONS, f"Missing HPC option: {opt}"

    def test_node_specific_options_exist(self):
        """Test that node-specific options are defined."""
        node_options = [
            "case_name",
            "run_label",
            "input_top_file",
            "input_gro_file",
            "input_mdp_file",
            "input_ndx_file",
            "num_threads",
            "max_warnings",
            "force_to_run",
        ]
        for opt in node_options:
            assert opt in GmxMdRun.OPTIONS, f"Missing option: {opt}"

    def test_run_label_default(self):
        """Test run_label has correct default value."""
        run_label = GmxMdRun.OPTIONS["run_label"]
        assert run_label.default == "md"

    def test_execution_mode_default(self):
        """Test execution_mode defaults to local."""
        execution_mode = GmxMdRun.OPTIONS["execution_mode"]
        assert execution_mode.default == "local"

    def test_file_parameters_have_docstrings(self):
        """Test that file parameters have documentation."""
        file_params = ["input_top_file", "input_gro_file", "input_mdp_file"]
        for param_name in file_params:
            param = GmxMdRun.OPTIONS[param_name]
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
        node = GmxMdRun(node_info)
        assert node is not None

    def test_node_has_node_id(self):
        """Test node has node_id from node_info."""
        node_info = create_mock_node_info()
        node = GmxMdRun(node_info)
        assert node.node_id == "test_node_001"


# =============================================================================
# Test: HPCNodeBase Abstract Methods Implementation
# =============================================================================


class TestHPCAbstractMethods:
    """Test HPCNodeBase abstract method implementations."""

    @pytest.fixture
    def mock_flow_vars(self):
        """Create mock flow_vars dictionary."""

        def make_var(value):
            mock = MagicMock()
            mock.get_value.return_value = value
            return mock

        return {
            "execution_mode": make_var("local"),
            "hpc_profile": make_var(""),
            "slurm_script": make_var(""),
            "force_resubmit": make_var(False),
            "case_name": make_var("test_case"),
            "run_label": make_var("nvt"),
            "input_top_file": make_var("/path/to/topol.top"),
            "input_gro_file": make_var("/path/to/conf.gro"),
            "input_mdp_file": make_var("/path/to/nvt.mdp"),
            "input_ndx_file": make_var("/path/to/index.ndx"),
            "num_threads": make_var("4"),
            "max_warnings": make_var("5"),
            "force_to_run": make_var(False),
        }

    def test_get_input_files(self, mock_flow_vars):
        """Test get_input_files returns correct file list."""
        node_info = create_mock_node_info()
        node = GmxMdRun(node_info)
        node.resolve_path = lambda x: x  # Identity function

        files = node.get_input_files(mock_flow_vars)

        assert "/path/to/topol.top" in files
        assert "/path/to/conf.gro" in files
        assert "/path/to/nvt.mdp" in files
        assert "/path/to/index.ndx" in files

    def test_get_output_files(self, mock_flow_vars):
        """Test get_output_files returns correct patterns."""
        node_info = create_mock_node_info()
        node = GmxMdRun(node_info)

        files = node.get_output_files(mock_flow_vars)

        assert "nvt.tpr" in files
        assert "nvt.gro" in files
        assert "nvt.xtc" in files
        assert "nvt.edr" in files
        assert "nvt.log" in files

    def test_get_template_variables(self, mock_flow_vars):
        """Test get_template_variables returns correct dict."""
        node_info = create_mock_node_info()
        node = GmxMdRun(node_info)

        vars = node.get_template_variables(mock_flow_vars)

        assert vars["RUN_LABEL"] == "nvt"
        assert vars["INPUT_TOP_FILE"] == "topol.top"
        assert vars["INPUT_GRO_FILE"] == "conf.gro"
        assert vars["INPUT_MDP_FILE"] == "nvt.mdp"
        assert vars["INPUT_NDX_FILE"] == "index.ndx"


# =============================================================================
# Test: run_local Method
# =============================================================================


class TestRunLocalMethod:
    """Test run_local method (local execution)."""

    @pytest.fixture
    def mock_flow_vars(self):
        """Create mock flow_vars dictionary."""

        def make_var(value):
            mock = MagicMock()
            mock.get_value.return_value = value
            return mock

        return {
            "execution_mode": make_var("local"),
            "hpc_profile": make_var(""),
            "slurm_script": make_var(""),
            "force_resubmit": make_var(False),
            "case_name": make_var("protein_sim"),
            "run_label": make_var("production"),
            "input_top_file": make_var("/work/topol.top"),
            "input_gro_file": make_var("/work/conf.gro"),
            "input_mdp_file": make_var("/work/md.mdp"),
            "input_ndx_file": make_var(None),
            "num_threads": make_var("8"),
            "max_warnings": make_var("10"),
            "force_to_run": make_var(False),
        }

    @patch.object(GmxMdRun, "resolve_path", side_effect=lambda x: x)
    def test_run_local_calls_core_function(self, mock_resolve, mock_flow_vars):
        """Test that run_local calls core.run_md_simulation."""
        with patch.object(
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
                grompp_returncode=0,
                mdrun_returncode=0,
            ),
        ) as mock_sim:
            node_info = create_mock_node_info()
            node = GmxMdRun(node_info)
            node.format_output_path = lambda x: f"abs:{x}"

            result = node.run_local([], mock_flow_vars)

            mock_sim.assert_called_once()
            call_kwargs = mock_sim.call_args[1]
            assert call_kwargs["run_label"] == "production"
            assert call_kwargs["num_threads"] == 8
            assert call_kwargs["max_warnings"] == 10

    @patch.object(GmxMdRun, "resolve_path", side_effect=lambda x: x)
    def test_run_local_returns_dict(self, mock_resolve, mock_flow_vars):
        """Test that run_local returns a dict (not JSON string)."""
        with patch.object(
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
                grompp_returncode=0,
                mdrun_returncode=0,
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRun(node_info)
            node.format_output_path = lambda x: x

            result = node.run_local([], mock_flow_vars)

            assert isinstance(result, dict)
            assert "success" in result

    @patch.object(GmxMdRun, "resolve_path", side_effect=lambda x: x)
    def test_run_local_failure_raises_exception(self, mock_resolve, mock_flow_vars):
        """Test that simulation failure raises NodeException."""
        with patch.object(
            node_module,
            "run_md_simulation",
            return_value=MagicMock(
                success=False,
                message="grompp failed: Error in topology",
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRun(node_info)
            node.format_output_path = lambda x: x

            with pytest.raises(NodeException) as exc_info:
                node.run_local([], mock_flow_vars)

            assert "grompp failed" in str(exc_info.value)


# =============================================================================
# Test: Execute Method - Mode Switching
# =============================================================================


class TestExecuteModeSwitching:
    """Test execute method switches between local and remote modes."""

    @pytest.fixture
    def mock_flow_vars_local(self):
        """Create mock flow_vars for local execution."""

        def make_var(value):
            mock = MagicMock()
            mock.get_value.return_value = value
            return mock

        return {
            "execution_mode": make_var("local"),
            "hpc_profile": make_var(""),
            "slurm_script": make_var(""),
            "force_resubmit": make_var(False),
            "case_name": make_var(""),
            "run_label": make_var("md"),
            "input_top_file": make_var("/path/to/topol.top"),
            "input_gro_file": make_var("/path/to/conf.gro"),
            "input_mdp_file": make_var("/path/to/md.mdp"),
            "input_ndx_file": make_var(None),
            "num_threads": make_var("0"),
            "max_warnings": make_var("10"),
            "force_to_run": make_var(False),
        }

    @pytest.fixture
    def mock_flow_vars_remote(self):
        """Create mock flow_vars for remote execution."""

        def make_var(value):
            mock = MagicMock()
            mock.get_value.return_value = value
            return mock

        return {
            "execution_mode": make_var("remote"),
            "hpc_profile": make_var("my_cluster"),
            "slurm_script": make_var("#!/bin/bash\n#SBATCH --job-name=test"),
            "force_resubmit": make_var(False),
            "case_name": make_var(""),
            "run_label": make_var("md"),
            "input_top_file": make_var("/path/to/topol.top"),
            "input_gro_file": make_var("/path/to/conf.gro"),
            "input_mdp_file": make_var("/path/to/md.mdp"),
            "input_ndx_file": make_var(None),
            "num_threads": make_var("0"),
            "max_warnings": make_var("10"),
            "force_to_run": make_var(False),
        }

    @patch.object(GmxMdRun, "resolve_path", side_effect=lambda x: x)
    @patch.object(GmxMdRun, "run_local")
    def test_execute_calls_run_local_for_local_mode(
        self, mock_run_local, mock_resolve, mock_flow_vars_local
    ):
        """Test that execute calls run_local when execution_mode is 'local'."""
        mock_run_local.return_value = {"success": True, "message": "OK"}

        node_info = create_mock_node_info()
        node = GmxMdRun(node_info)

        node.execute([], mock_flow_vars_local)

        mock_run_local.assert_called_once()

    def test_execute_remote_requires_hpc_profile(self, mock_flow_vars_remote):
        """Test that remote execution requires HPC profile."""
        # Remove hpc_profile
        mock_flow_vars_remote["hpc_profile"].get_value.return_value = ""

        node_info = create_mock_node_info()
        node = GmxMdRun(node_info)

        with pytest.raises(NodeException) as exc_info:
            node.execute([], mock_flow_vars_remote)

        assert "HPC Profile is required" in str(exc_info.value)

    def test_execute_remote_requires_slurm_script(self, mock_flow_vars_remote):
        """Test that remote execution requires SLURM script."""
        # Remove slurm_script
        mock_flow_vars_remote["slurm_script"].get_value.return_value = ""

        node_info = create_mock_node_info()
        node = GmxMdRun(node_info)

        with pytest.raises(NodeException) as exc_info:
            node.execute([], mock_flow_vars_remote)

        assert "SLURM job script is required" in str(exc_info.value)


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
            "execution_mode": make_var("local"),
            "hpc_profile": make_var(""),
            "slurm_script": make_var(""),
            "force_resubmit": make_var(False),
            "case_name": make_var(""),  # Empty - should use predecessor
            "run_label": make_var("md"),
            "input_top_file": make_var("/path/to/topol.top"),
            "input_gro_file": make_var("/path/to/conf.gro"),
            "input_mdp_file": make_var("/path/to/md.mdp"),
            "input_ndx_file": make_var(None),
            "num_threads": make_var("0"),
            "max_warnings": make_var("10"),
            "force_to_run": make_var(False),
        }

    @patch.object(GmxMdRun, "resolve_path", side_effect=lambda x: x)
    def test_case_name_from_predecessor(self, mock_resolve, mock_flow_vars_no_case):
        """Test that case_name falls back to predecessor data."""
        with patch.object(
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
                grompp_returncode=0,
                mdrun_returncode=0,
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRun(node_info)
            node.format_output_path = lambda x: x

            predecessor_data = [{"case_name": "from_predecessor"}]

            result = node.run_local(predecessor_data, mock_flow_vars_no_case)

            assert result["data"]["case_name"] == "from_predecessor"

    @patch.object(GmxMdRun, "resolve_path", side_effect=lambda x: x)
    def test_default_case_name_when_no_predecessor(
        self, mock_resolve, mock_flow_vars_no_case
    ):
        """Test default case_name when no predecessor data."""
        with patch.object(
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
                grompp_returncode=0,
                mdrun_returncode=0,
            ),
        ):
            node_info = create_mock_node_info()
            node = GmxMdRun(node_info)
            node.format_output_path = lambda x: x

            result = node.run_local([], mock_flow_vars_no_case)

            assert result["data"]["case_name"] == "protein"  # default


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
            "execution_mode": make_var("local"),
            "hpc_profile": make_var(""),
            "slurm_script": make_var(""),
            "force_resubmit": make_var(False),
            "case_name": make_var("integration_test"),
            "run_label": make_var("em"),
            "input_top_file": make_var("/sim/system.top"),
            "input_gro_file": make_var("/sim/start.gro"),
            "input_mdp_file": make_var("/sim/em.mdp"),
            "input_ndx_file": make_var("/sim/index.ndx"),
            "num_threads": make_var("2"),
            "max_warnings": make_var("3"),
            "force_to_run": make_var(False),
        }

    @patch.object(GmxMdRun, "resolve_path", side_effect=lambda x: x)
    def test_core_function_called_with_correct_args(self, mock_resolve, mock_flow_vars):
        """Test that core function is called with all correct arguments."""
        with patch.object(
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
                grompp_returncode=0,
                mdrun_returncode=0,
            ),
        ) as mock_sim:
            node_info = create_mock_node_info()
            node = GmxMdRun(node_info)
            node.format_output_path = lambda x: x

            node.run_local([], mock_flow_vars)

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
            assert call_kwargs["verbose"] is True
