"""
Tests for pdb-tools-clean node

Run with: pixi run test
"""

import os
import tempfile
from pathlib import Path

import pytest

# Mock bocoflow_core for testing without full installation
try:
    from bocoflow_core.node import NodeInputParams, NodeOutputParams
except ImportError:
    # Create mock classes for standalone testing
    class NodeInputParams:
        def __init__(self, input_data=None, parameters=None, context=None):
            self.input_data = input_data
            self.parameters = parameters or {}
            self.context = context or {}

    class NodeOutputParams:
        def __init__(self, output_data=None, status="success", message=""):
            self.output_data = output_data
            self.status = status
            self.message = message

# Import the node after mocking
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from node import PDBToolsClean


# Sample PDB content for testing
SAMPLE_PDB = """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.986  -0.760  -1.217  1.00  0.00           C
ATOM      6  N   GLY A   2       3.320   1.540   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       3.970   2.850   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2       5.480   2.720   0.000  1.00  0.00           C
ATOM      9  O   GLY A   2       6.030   1.610   0.000  1.00  0.00           O
ATOM     10  N   ALA B   1      10.000   0.000   0.000  1.00  0.00           N
ATOM     11  CA  ALA B   1      11.458   0.000   0.000  1.00  0.00           C
ATOM     12  C   ALA B   1      12.009   1.420   0.000  1.00  0.00           C
ATOM     13  O   ALA B   1      11.246   2.390   0.000  1.00  0.00           O
HETATM   14  O   HOH A 100       5.000   5.000   5.000  1.00  0.00           O
HETATM   15  O   HOH A 101       6.000   6.000   6.000  1.00  0.00           O
END
"""


class TestPDBToolsClean:
    """Test suite for PDBToolsClean node."""

    def setup_method(self):
        """Set up test fixtures."""
        self.node = PDBToolsClean()
        self.temp_dir = tempfile.mkdtemp()
        self.test_pdb = Path(self.temp_dir) / "test.pdb"
        with open(self.test_pdb, 'w') as f:
            f.write(SAMPLE_PDB)

    def teardown_method(self):
        """Clean up test files."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_node_instantiation(self):
        """Test that node can be instantiated."""
        assert self.node is not None

    def test_get_parameters_schema(self):
        """Test that parameters schema is valid."""
        schema = self.node.get_parameters_schema()
        assert isinstance(schema, dict)
        assert "type" in schema
        assert schema["type"] == "object"
        assert "input_pdb" in schema["properties"]

    def test_count_atoms(self):
        """Test atom counting function."""
        count = self.node._count_atoms(self.test_pdb)
        # 13 ATOM + 2 HETATM = 15 total
        assert count == 15

    def test_get_chains(self):
        """Test chain extraction function."""
        chains = self.node._get_chains(self.test_pdb)
        assert sorted(chains) == ["A", "B"]

    def test_resolve_path_absolute(self):
        """Test path resolution with absolute path."""
        path = self.node._resolve_path("/absolute/path.pdb", "/working")
        assert str(path) == "/absolute/path.pdb"

    def test_resolve_path_relative(self):
        """Test path resolution with relative path."""
        path = self.node._resolve_path("relative/path.pdb", "/working")
        assert str(path) == "/working/relative/path.pdb"

    def test_generate_output_path(self):
        """Test output path generation."""
        output = self.node._generate_output_path(Path("/test/input.pdb"), "_clean")
        assert str(output) == "/test/input_clean.pdb"

    def test_describe_operations(self):
        """Test operation description."""
        ops = self.node._describe_operations(
            chains="A,B",
            remove_hetatm=True,
            remove_hydrogens=True,
            renumber_residues=1
        )
        assert "Selected chains: A,B" in ops
        assert "Removed HETATM records" in ops
        assert "Removed hydrogen atoms" in ops
        assert "Renumbered residues starting from 1" in ops
        assert "Tidied PDB format" in ops

    def test_build_pipeline_all_options(self):
        """Test pipeline building with all options."""
        pipeline = self.node._build_pipeline(
            input_path=self.test_pdb,
            chains="A,B",
            remove_hetatm=True,
            remove_hydrogens=True,
            renumber_residues=1
        )
        # Should have: selchain, delhetatm, delelem, reres, tidy
        assert len(pipeline) == 5

    def test_build_pipeline_minimal(self):
        """Test pipeline building with minimal options."""
        pipeline = self.node._build_pipeline(
            input_path=self.test_pdb,
            chains="",
            remove_hetatm=False,
            remove_hydrogens=False,
            renumber_residues=0
        )
        # Should have only: tidy
        assert len(pipeline) == 1
        assert pipeline[0] == ["pdb_tidy"]

    def test_execute_missing_input(self):
        """Test execution with missing input file."""
        input_params = NodeInputParams(
            parameters={"input_pdb": ""},
            context={"working_path": self.temp_dir}
        )
        result = self.node.execute(input_params)
        assert result.status == "error"
        assert "No input PDB file specified" in result.message

    def test_execute_file_not_found(self):
        """Test execution with non-existent file."""
        input_params = NodeInputParams(
            parameters={"input_pdb": "nonexistent.pdb"},
            context={"working_path": self.temp_dir}
        )
        result = self.node.execute(input_params)
        assert result.status == "error"
        assert "not found" in result.message.lower()


# Integration tests (require pdb-tools installed)
@pytest.mark.skipif(
    os.system("pdb_tidy --help > /dev/null 2>&1") != 0,
    reason="pdb-tools not installed"
)
class TestPDBToolsCleanIntegration:
    """Integration tests requiring pdb-tools."""

    def setup_method(self):
        """Set up test fixtures."""
        self.node = PDBToolsClean()
        self.temp_dir = tempfile.mkdtemp()
        self.test_pdb = Path(self.temp_dir) / "test.pdb"
        with open(self.test_pdb, 'w') as f:
            f.write(SAMPLE_PDB)

    def teardown_method(self):
        """Clean up test files."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_execute_basic_clean(self):
        """Test basic cleaning execution."""
        input_params = NodeInputParams(
            parameters={
                "input_pdb": str(self.test_pdb),
                "remove_hetatm": True,
                "output_suffix": "_clean"
            },
            context={"working_path": self.temp_dir}
        )
        result = self.node.execute(input_params)

        assert result.status == "success"
        assert "output_pdb" in result.output_data

        output_path = Path(result.output_data["output_pdb"])
        assert output_path.exists()

        # Check HETATM removed
        assert result.output_data["atoms_removed"] >= 2

    def test_execute_chain_selection(self):
        """Test chain selection."""
        input_params = NodeInputParams(
            parameters={
                "input_pdb": str(self.test_pdb),
                "chains": "A",
                "remove_hetatm": True,
                "output_suffix": "_chainA"
            },
            context={"working_path": self.temp_dir}
        )
        result = self.node.execute(input_params)

        assert result.status == "success"
        assert result.output_data["chains_selected"] == ["A"]
