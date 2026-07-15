"""
Tests for pdb-fasta-biopython core.py — pure Python logic.

Tests cover:
  - RCSB API fetching (mocked)
  - PDB structure parsing (real, using demo_data/3LZ0.pdb)
  - Sequence extraction
  - Missing residues extraction
  - FASTA file writing
  - Missing residues CSV writing

Run: pytest tests/test_core.py -v
"""

import csv
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Import via package path (avoids bare "core" collision across packages)
from bocoflow_nodes.pdb_fasta_biopython.core import (  # noqa: E402
    ChainInfo,
    MissingResidue,
    extract_missing_residues,
    extract_sequences,
    fetch_fasta_from_rcsb,
    fetch_pdb_from_rcsb,
    parse_pdb_structure,
    write_fasta_files,
    write_missing_residues_csv,
)

# Path to demo data
_pkg_dir = Path(__file__).parent.parent
DEMO_DATA_DIR = _pkg_dir / "demo_data"
DEMO_PDB = DEMO_DATA_DIR / "3LZ0.pdb"
DEMO_FASTA = DEMO_DATA_DIR / "3LZ0.fasta"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def output_dir():
    """Provide a temporary output directory."""
    with tempfile.TemporaryDirectory(prefix="bocoflow-pdb-test-") as tmpdir:
        yield tmpdir


@pytest.fixture
def structure_3lz0():
    """Parse the 3LZ0 demo PDB into a BioPython Structure."""
    if not DEMO_PDB.exists():
        pytest.skip("Demo PDB file 3LZ0.pdb not found")
    return parse_pdb_structure(str(DEMO_PDB), "3LZ0", is_file=True)


@pytest.fixture
def chains_3lz0(structure_3lz0):
    """Extract chains from 3LZ0."""
    return extract_sequences(structure_3lz0, include_hetatm=False)


# ---------------------------------------------------------------------------
# Tests: RCSB API (mocked)
# ---------------------------------------------------------------------------

class TestRcsbApi:
    """Test RCSB API fetching with mocked HTTP requests."""

    @patch("bocoflow_nodes.pdb_fasta_biopython.core.requests.get")
    def test_fetch_pdb_success(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.text = "HEADER  MOCK PDB\nATOM      1  N   ALA A   1\nEND\n"
        mock_get.return_value = mock_resp

        result = fetch_pdb_from_rcsb("3LZ0")

        assert "HEADER" in result
        mock_get.assert_called_once()
        assert "3LZ0" in mock_get.call_args[0][0]

    @patch("bocoflow_nodes.pdb_fasta_biopython.core.requests.get")
    def test_fetch_pdb_not_found(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.status_code = 404
        mock_get.return_value = mock_resp

        with pytest.raises(ValueError, match="Failed to fetch PDB"):
            fetch_pdb_from_rcsb("ZZZZ")

    @patch("bocoflow_nodes.pdb_fasta_biopython.core.requests.get")
    def test_fetch_fasta_success(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.text = ">3LZ0_1|Chain A|Protein\nMKWVTFISLL\n"
        mock_get.return_value = mock_resp

        result = fetch_fasta_from_rcsb("3LZ0")

        assert ">3LZ0" in result
        mock_get.assert_called_once()

    @patch("bocoflow_nodes.pdb_fasta_biopython.core.requests.get")
    def test_fetch_fasta_not_found(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.status_code = 404
        mock_get.return_value = mock_resp

        with pytest.raises(ValueError, match="Failed to fetch FASTA"):
            fetch_fasta_from_rcsb("ZZZZ")

    @patch("bocoflow_nodes.pdb_fasta_biopython.core.requests.get")
    def test_fetch_pdb_case_insensitive(self, mock_get):
        """PDB ID should be uppercased in the URL."""
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.text = "HEADER\nEND\n"
        mock_get.return_value = mock_resp

        fetch_pdb_from_rcsb("3lz0")

        url = mock_get.call_args[0][0]
        assert "3LZ0" in url


# ---------------------------------------------------------------------------
# Tests: PDB parsing (real data)
# ---------------------------------------------------------------------------

class TestPdbParsing:
    """Test PDB parsing with real 3LZ0 demo data."""

    def test_parse_pdb_from_file(self, structure_3lz0):
        assert structure_3lz0 is not None
        assert structure_3lz0.id == "3LZ0"

    def test_parse_pdb_from_text(self, output_dir):
        if not DEMO_PDB.exists():
            pytest.skip("Demo PDB not found")
        pdb_text = DEMO_PDB.read_text()
        structure = parse_pdb_structure(pdb_text, "test", is_file=False)
        assert structure is not None

    def test_extract_sequences_chain_count(self, chains_3lz0):
        # 3LZ0 has chains A-H (protein) + I,J (DNA) = 10 chains
        assert len(chains_3lz0) >= 8  # At least the 8 protein chains

    def test_extract_sequences_chain_types(self, chains_3lz0):
        # Chains A-H should be protein
        for cid in ["A", "B", "C", "D", "E", "F", "G", "H"]:
            if cid in chains_3lz0:
                assert chains_3lz0[cid].chain_type == "protein"

        # Chains I,J should be nucleic_acid
        for cid in ["I", "J"]:
            if cid in chains_3lz0:
                assert chains_3lz0[cid].chain_type == "nucleic_acid"

    def test_chain_info_fields(self, chains_3lz0):
        chain_a = chains_3lz0.get("A")
        assert chain_a is not None
        assert chain_a.chain_id == "A"
        assert chain_a.length > 0
        assert chain_a.protein_residues > 0
        assert len(chain_a.sequence) == chain_a.length

    def test_include_hetatm(self, structure_3lz0):
        chains_no_het = extract_sequences(structure_3lz0, include_hetatm=False)
        chains_with_het = extract_sequences(structure_3lz0, include_hetatm=True)

        # With HETATM may have more residues or same
        for cid in chains_no_het:
            if cid in chains_with_het:
                assert chains_with_het[cid].length >= chains_no_het[cid].length

    def test_sequences_are_valid_amino_acids(self, chains_3lz0):
        valid_aa = set("ACDEFGHIKLMNPQRSTVWYX")
        valid_na = set("ACGTUX")
        for cid, info in chains_3lz0.items():
            if info.chain_type == "protein":
                assert all(c in valid_aa for c in info.sequence), (
                    f"Chain {cid} has invalid AA: {info.sequence}"
                )
            else:
                assert all(c in valid_na for c in info.sequence), (
                    f"Chain {cid} has invalid NA: {info.sequence}"
                )


# ---------------------------------------------------------------------------
# Tests: Missing residues
# ---------------------------------------------------------------------------

class TestMissingResidues:
    """Test missing residues extraction."""

    def test_extract_missing_residues(self, structure_3lz0):
        missing = extract_missing_residues(structure_3lz0)
        # 3LZ0 does have missing residues, result should be a dict
        assert isinstance(missing, dict)

    def test_missing_residue_fields(self, structure_3lz0):
        missing = extract_missing_residues(structure_3lz0)
        for chain_id, residues in missing.items():
            for mr in residues:
                assert isinstance(mr, MissingResidue)
                assert mr.chain == chain_id
                assert len(mr.res_name) == 3
                assert len(mr.one_letter) == 1


# ---------------------------------------------------------------------------
# Tests: FASTA writing
# ---------------------------------------------------------------------------

class TestFastaWriting:
    """Test FASTA file output."""

    def test_write_split_chains(self, chains_3lz0, output_dir):
        files = write_fasta_files(
            chains_3lz0, output_dir, "3lz0", split_chains=True
        )
        # Should have one file per chain
        assert len(files) == len(chains_3lz0)
        for label, path in files.items():
            assert os.path.exists(path)
            assert path.endswith(".fasta")
            content = open(path).read()
            assert content.startswith(">")

    def test_write_combined(self, chains_3lz0, output_dir):
        files = write_fasta_files(
            chains_3lz0, output_dir, "3lz0", split_chains=False
        )
        assert "combined" in files
        combined_path = files["combined"]
        assert os.path.exists(combined_path)

        content = open(combined_path).read()
        # Should have one header per chain
        headers = [line for line in content.splitlines() if line.startswith(">")]
        assert len(headers) == len(chains_3lz0)

    def test_fasta_file_naming(self, output_dir):
        chains = {
            "A": ChainInfo("A", "MKWVTF", 6, 6, 0, "protein"),
        }
        files = write_fasta_files(chains, output_dir, "test_case", split_chains=True)
        assert "chain_A" in files
        assert "test_case_chain_A.fasta" in files["chain_A"]


# ---------------------------------------------------------------------------
# Tests: Missing residues CSV writing
# ---------------------------------------------------------------------------

class TestMissingResiduesCsv:
    """Test CSV output for missing residues."""

    def test_write_csv(self, output_dir):
        missing = {
            "A": [
                MissingResidue("A", "ALA", "A", 1),
                MissingResidue("A", "GLY", "G", 2),
            ],
            "B": [
                MissingResidue("B", "LEU", "L", 10),
            ],
        }
        csv_files = write_missing_residues_csv(missing, output_dir)

        assert len(csv_files) == 2
        assert "A" in csv_files
        assert "B" in csv_files

        # Verify CSV content
        with open(csv_files["A"]) as f:
            reader = csv.reader(f)
            rows = list(reader)
        # header + 2 data rows
        assert len(rows) == 3
        assert rows[0][0] == "chain"
        assert rows[1][2] == "A"  # one_letter for ALA

    def test_write_csv_empty(self, output_dir):
        csv_files = write_missing_residues_csv({}, output_dir)
        assert csv_files == {}

    def test_csv_with_real_data(self, structure_3lz0, output_dir):
        missing = extract_missing_residues(structure_3lz0)
        if not missing:
            pytest.skip("No missing residues in 3LZ0")
        csv_files = write_missing_residues_csv(missing, output_dir)
        for chain_id, path in csv_files.items():
            assert os.path.exists(path)
            assert os.path.getsize(path) > 0
