"""
Level 1 Tests: Core Functions (BoCoFlow Independent)

These tests verify the pure Python functions in core.py without any
BoCoFlow dependencies. They can be run quickly and independently.

Run with: pytest test_core.py -v

Tests are organized by function:
1. test_build_grompp_command_* - Command building for grompp
2. test_build_mdrun_command_* - Command building for mdrun
3. test_dataclasses_* - Configuration dataclasses
4. test_run_* - Execution functions (mocked subprocess)
5. test_check_gromacs_* - GROMACS availability checks
"""

import importlib.util
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Import core module using importlib (handles hyphenated directory name)
_module_dir = Path(__file__).parent.parent
_core_spec = importlib.util.spec_from_file_location("core", _module_dir / "core.py")
core = importlib.util.module_from_spec(_core_spec)
_core_spec.loader.exec_module(core)

# Extract functions and classes from core module
GromppConfig = core.GromppConfig
MdrunConfig = core.MdrunConfig
SimulationResult = core.SimulationResult
build_grompp_command = core.build_grompp_command
build_mdrun_command = core.build_mdrun_command
check_gromacs_available = core.check_gromacs_available
get_gromacs_version = core.get_gromacs_version
run_grompp = core.run_grompp
run_md_simulation = core.run_md_simulation
run_mdrun = core.run_mdrun


# =============================================================================
# Test: GromppConfig Dataclass
# =============================================================================


class TestGromppConfig:
    """Test GromppConfig dataclass."""

    def test_basic_config(self):
        """Test creating a basic GromppConfig."""
        config = GromppConfig(
            mdp_file="md.mdp",
            gro_file="conf.gro",
            top_file="topol.top",
            output_tpr="md.tpr",
        )
        assert config.mdp_file == "md.mdp"
        assert config.gro_file == "conf.gro"
        assert config.top_file == "topol.top"
        assert config.output_tpr == "md.tpr"
        assert config.ndx_file is None
        assert config.maxwarn == 10  # default
        assert config.restraint_file is None

    def test_full_config(self):
        """Test GromppConfig with all options."""
        config = GromppConfig(
            mdp_file="nvt.mdp",
            gro_file="em.gro",
            top_file="system.top",
            output_tpr="nvt.tpr",
            ndx_file="index.ndx",
            maxwarn=5,
            restraint_file="posre.gro",
        )
        assert config.ndx_file == "index.ndx"
        assert config.maxwarn == 5
        assert config.restraint_file == "posre.gro"


class TestMdrunConfig:
    """Test MdrunConfig dataclass."""

    def test_basic_config(self):
        """Test creating a basic MdrunConfig."""
        config = MdrunConfig(deffnm="md")
        assert config.deffnm == "md"
        assert config.num_threads == 0  # default (auto)
        assert config.verbose is True  # default
        assert config.gpu_ids is None
        assert config.extra_args == []

    def test_full_config(self):
        """Test MdrunConfig with all options."""
        config = MdrunConfig(
            deffnm="production",
            num_threads=8,
            verbose=False,
            gpu_ids="0,1",
            extra_args=["-nsteps", "1000"],
        )
        assert config.num_threads == 8
        assert config.verbose is False
        assert config.gpu_ids == "0,1"
        assert config.extra_args == ["-nsteps", "1000"]


class TestSimulationResult:
    """Test SimulationResult dataclass."""

    def test_success_result(self):
        """Test creating a successful SimulationResult."""
        result = SimulationResult(
            success=True,
            message="Completed",
            tpr_file="/path/to/md.tpr",
            gro_file="/path/to/md.gro",
        )
        assert result.success is True
        assert result.tpr_file == "/path/to/md.tpr"
        assert result.gro_file == "/path/to/md.gro"
        assert result.xtc_file is None  # not all files generated

    def test_failure_result(self):
        """Test creating a failed SimulationResult."""
        result = SimulationResult(
            success=False,
            message="grompp failed",
            grompp_returncode=1,
            grompp_stderr="Error in topology",
        )
        assert result.success is False
        assert result.grompp_returncode == 1
        assert "topology" in result.grompp_stderr


# =============================================================================
# Test: build_grompp_command
# =============================================================================


class TestBuildGromppCommand:
    """Test build_grompp_command function."""

    def test_basic_command(self):
        """Test building a basic grompp command."""
        config = GromppConfig(
            mdp_file="md.mdp",
            gro_file="conf.gro",
            top_file="topol.top",
            output_tpr="md.tpr",
        )
        cmd = build_grompp_command(config)

        assert cmd[0:2] == ["gmx", "grompp"]
        assert "-f" in cmd and "md.mdp" in cmd
        assert "-c" in cmd and "conf.gro" in cmd
        assert "-p" in cmd and "topol.top" in cmd
        assert "-o" in cmd and "md.tpr" in cmd
        assert "-maxwarn" in cmd and "10" in cmd
        # -r should default to gro_file
        r_idx = cmd.index("-r")
        assert cmd[r_idx + 1] == "conf.gro"

    def test_command_with_index_file(self):
        """Test grompp command includes index file when specified."""
        config = GromppConfig(
            mdp_file="md.mdp",
            gro_file="conf.gro",
            top_file="topol.top",
            output_tpr="md.tpr",
            ndx_file="index.ndx",
        )
        cmd = build_grompp_command(config)

        assert "-n" in cmd
        n_idx = cmd.index("-n")
        assert cmd[n_idx + 1] == "index.ndx"

    def test_command_with_restraint_file(self):
        """Test grompp command uses restraint file for -r flag."""
        config = GromppConfig(
            mdp_file="md.mdp",
            gro_file="conf.gro",
            top_file="topol.top",
            output_tpr="md.tpr",
            restraint_file="posre.gro",
        )
        cmd = build_grompp_command(config)

        r_idx = cmd.index("-r")
        assert cmd[r_idx + 1] == "posre.gro"

    def test_command_with_custom_maxwarn(self):
        """Test grompp command with custom maxwarn."""
        config = GromppConfig(
            mdp_file="md.mdp",
            gro_file="conf.gro",
            top_file="topol.top",
            output_tpr="md.tpr",
            maxwarn=0,
        )
        cmd = build_grompp_command(config)

        maxwarn_idx = cmd.index("-maxwarn")
        assert cmd[maxwarn_idx + 1] == "0"


# =============================================================================
# Test: build_mdrun_command
# =============================================================================


class TestBuildMdrunCommand:
    """Test build_mdrun_command function."""

    def test_basic_command(self):
        """Test building a basic mdrun command."""
        config = MdrunConfig(deffnm="md")
        cmd = build_mdrun_command(config)

        assert cmd[0:2] == ["gmx", "mdrun"]
        assert "-deffnm" in cmd and "md" in cmd
        assert "-v" in cmd  # verbose by default

    def test_command_with_threads(self):
        """Test mdrun command with thread specification."""
        config = MdrunConfig(deffnm="md", num_threads=4)
        cmd = build_mdrun_command(config)

        assert "-nt" in cmd
        nt_idx = cmd.index("-nt")
        assert cmd[nt_idx + 1] == "4"

    def test_command_without_verbose(self):
        """Test mdrun command without verbose flag."""
        config = MdrunConfig(deffnm="md", verbose=False)
        cmd = build_mdrun_command(config)

        assert "-v" not in cmd

    def test_command_with_gpu(self):
        """Test mdrun command with GPU specification."""
        config = MdrunConfig(deffnm="md", gpu_ids="0")
        cmd = build_mdrun_command(config)

        assert "-gpu_id" in cmd
        gpu_idx = cmd.index("-gpu_id")
        assert cmd[gpu_idx + 1] == "0"

    def test_command_with_extra_args(self):
        """Test mdrun command with extra arguments."""
        config = MdrunConfig(
            deffnm="md",
            extra_args=["-nsteps", "1000", "-cpt", "15"],
        )
        cmd = build_mdrun_command(config)

        assert "-nsteps" in cmd and "1000" in cmd
        assert "-cpt" in cmd and "15" in cmd

    def test_zero_threads_not_added(self):
        """Test that 0 threads (auto) doesn't add -nt flag."""
        config = MdrunConfig(deffnm="md", num_threads=0)
        cmd = build_mdrun_command(config)

        assert "-nt" not in cmd


# =============================================================================
# Test: run_grompp (mocked subprocess)
# =============================================================================


class TestRunGrompp:
    """Test run_grompp function with mocked subprocess."""

    @pytest.fixture
    def temp_dir_with_files(self):
        """Create temp directory with mock input files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create mock input files
            (Path(tmpdir) / "md.mdp").touch()
            (Path(tmpdir) / "conf.gro").touch()
            (Path(tmpdir) / "topol.top").touch()
            yield tmpdir

    def test_missing_mdp_file(self, temp_dir_with_files):
        """Test that missing MDP file raises FileNotFoundError."""
        config = GromppConfig(
            mdp_file="nonexistent.mdp",
            gro_file="conf.gro",
            top_file="topol.top",
            output_tpr="md.tpr",
        )

        with pytest.raises(FileNotFoundError) as exc_info:
            run_grompp(config, temp_dir_with_files)

        assert "MDP" in str(exc_info.value)

    def test_missing_gro_file(self, temp_dir_with_files):
        """Test that missing GRO file raises FileNotFoundError."""
        config = GromppConfig(
            mdp_file="md.mdp",
            gro_file="nonexistent.gro",
            top_file="topol.top",
            output_tpr="md.tpr",
        )

        with pytest.raises(FileNotFoundError) as exc_info:
            run_grompp(config, temp_dir_with_files)

        assert "GRO" in str(exc_info.value)

    def test_missing_ndx_file(self, temp_dir_with_files):
        """Test that missing NDX file raises FileNotFoundError."""
        config = GromppConfig(
            mdp_file="md.mdp",
            gro_file="conf.gro",
            top_file="topol.top",
            output_tpr="md.tpr",
            ndx_file="nonexistent.ndx",
        )

        with pytest.raises(FileNotFoundError) as exc_info:
            run_grompp(config, temp_dir_with_files)

        assert "NDX" in str(exc_info.value)

    @patch("subprocess.run")
    def test_successful_grompp(self, mock_run, temp_dir_with_files):
        """Test successful grompp execution."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        config = GromppConfig(
            mdp_file="md.mdp",
            gro_file="conf.gro",
            top_file="topol.top",
            output_tpr="md.tpr",
        )

        result = run_grompp(config, temp_dir_with_files)

        assert result.returncode == 0
        mock_run.assert_called_once()
        # Verify correct command was built
        call_args = mock_run.call_args
        cmd = call_args[0][0]
        assert cmd[0:2] == ["gmx", "grompp"]


# =============================================================================
# Test: run_mdrun (mocked subprocess)
# =============================================================================


class TestRunMdrun:
    """Test run_mdrun function with mocked subprocess."""

    @pytest.fixture
    def temp_dir_with_tpr(self):
        """Create temp directory with mock TPR file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "md.tpr").touch()
            yield tmpdir

    def test_missing_tpr_file(self):
        """Test that missing TPR file raises FileNotFoundError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = MdrunConfig(deffnm="md")

            with pytest.raises(FileNotFoundError) as exc_info:
                run_mdrun(config, tmpdir)

            assert "TPR" in str(exc_info.value)

    @patch("subprocess.run")
    def test_successful_mdrun(self, mock_run, temp_dir_with_tpr):
        """Test successful mdrun execution."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        config = MdrunConfig(deffnm="md")

        result = run_mdrun(config, temp_dir_with_tpr)

        assert result.returncode == 0
        mock_run.assert_called_once()


# =============================================================================
# Test: run_md_simulation (mocked subprocess)
# =============================================================================


class TestRunMdSimulation:
    """Test run_md_simulation orchestrator function."""

    @pytest.fixture
    def temp_simulation_dir(self):
        """Create temp directory with all required input files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "topol.top").touch()
            (Path(tmpdir) / "conf.gro").touch()
            (Path(tmpdir) / "md.mdp").touch()
            yield tmpdir

    def test_missing_input_files(self):
        """Test handling of missing input files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_md_simulation(
                top_file="nonexistent.top",
                gro_file="nonexistent.gro",
                mdp_file="nonexistent.mdp",
                working_dir=tmpdir,
            )

            assert result.success is False
            assert "not found" in result.message.lower()

    @patch("subprocess.run")
    def test_grompp_failure(self, mock_run, temp_simulation_dir):
        """Test handling of grompp failure."""
        mock_run.return_value = MagicMock(
            returncode=1,
            stdout="",
            stderr="Error in topology file",
        )

        result = run_md_simulation(
            top_file="topol.top",
            gro_file="conf.gro",
            mdp_file="md.mdp",
            working_dir=temp_simulation_dir,
        )

        assert result.success is False
        assert result.grompp_returncode == 1
        assert "grompp failed" in result.message

    @patch("subprocess.run")
    def test_mdrun_failure(self, mock_run, temp_simulation_dir):
        """Test handling of mdrun failure."""

        # First call (grompp) succeeds, second call (mdrun) fails
        def side_effect(*args, **kwargs):
            cmd = args[0]
            if "grompp" in cmd:
                # Create TPR file for grompp success
                (Path(temp_simulation_dir) / "md.tpr").touch()
                return MagicMock(returncode=0, stdout="", stderr="")
            else:
                return MagicMock(returncode=1, stdout="", stderr="mdrun error")

        mock_run.side_effect = side_effect

        result = run_md_simulation(
            top_file="topol.top",
            gro_file="conf.gro",
            mdp_file="md.mdp",
            working_dir=temp_simulation_dir,
        )

        assert result.success is False
        assert result.mdrun_returncode == 1
        assert "mdrun failed" in result.message
        # TPR should still be recorded
        assert result.tpr_file is not None

    @patch("subprocess.run")
    def test_successful_simulation(self, mock_run, temp_simulation_dir):
        """Test successful simulation run."""

        def side_effect(*args, **kwargs):
            cmd = args[0]
            if "grompp" in cmd:
                # Create TPR file
                (Path(temp_simulation_dir) / "md.tpr").touch()
            else:
                # Create output files
                (Path(temp_simulation_dir) / "md.gro").touch()
                (Path(temp_simulation_dir) / "md.log").touch()
            return MagicMock(returncode=0, stdout="", stderr="")

        mock_run.side_effect = side_effect

        result = run_md_simulation(
            top_file="topol.top",
            gro_file="conf.gro",
            mdp_file="md.mdp",
            working_dir=temp_simulation_dir,
            run_label="md",
        )

        assert result.success is True
        assert "completed successfully" in result.message
        assert result.tpr_file is not None
        assert result.gro_file is not None
        assert result.log_file is not None

    @patch("subprocess.run")
    def test_custom_parameters(self, mock_run, temp_simulation_dir):
        """Test simulation with custom parameters."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")
        # Create TPR file
        (Path(temp_simulation_dir) / "nvt.tpr").touch()

        result = run_md_simulation(
            top_file="topol.top",
            gro_file="conf.gro",
            mdp_file="md.mdp",
            working_dir=temp_simulation_dir,
            run_label="nvt",
            num_threads=4,
            max_warnings=5,
            verbose=False,
        )

        # Verify grompp was called with maxwarn=5
        grompp_call = mock_run.call_args_list[0]
        grompp_cmd = grompp_call[0][0]
        maxwarn_idx = grompp_cmd.index("-maxwarn")
        assert grompp_cmd[maxwarn_idx + 1] == "5"

        # Verify mdrun was called with -nt 4 and without -v
        mdrun_call = mock_run.call_args_list[1]
        mdrun_cmd = mdrun_call[0][0]
        assert "-nt" in mdrun_cmd
        assert "-v" not in mdrun_cmd


# =============================================================================
# Test: GROMACS availability checks
# =============================================================================


class TestGromacsAvailability:
    """Test GROMACS availability check functions."""

    @patch("subprocess.run")
    def test_gromacs_available(self, mock_run):
        """Test when GROMACS is available."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="GROMACS version 2023.3",
            stderr="",
        )

        assert check_gromacs_available() is True

    @patch("subprocess.run")
    def test_gromacs_not_available_returncode(self, mock_run):
        """Test when gmx returns non-zero."""
        mock_run.return_value = MagicMock(returncode=1)

        assert check_gromacs_available() is False

    @patch("subprocess.run")
    def test_gromacs_not_available_not_found(self, mock_run):
        """Test when gmx command not found."""
        mock_run.side_effect = FileNotFoundError()

        assert check_gromacs_available() is False

    @patch("subprocess.run")
    def test_get_version(self, mock_run):
        """Test getting GROMACS version."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="GROMACS version 2023.3\nSome other info",
            stderr="",
        )

        version = get_gromacs_version()
        assert version == "2023.3"

    @patch("subprocess.run")
    def test_get_version_not_available(self, mock_run):
        """Test version when GROMACS not available."""
        mock_run.side_effect = FileNotFoundError()

        version = get_gromacs_version()
        assert version is None
