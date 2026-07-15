"""Module 2: Ligand (metal complex) preparation.

Canonicalize the XYZ structure with OpenBabel, then build a NetworkX
molecular graph from the output.
"""

import logging
import subprocess
from pathlib import Path

import networkx as nx

from .xyz2graph import build_molecular_graph

logger = logging.getLogger(__name__)


def canonicalize_xyz(
    input_xyz: Path,
    output_xyz: Path,
    obabel_path: str = "obabel",
    use_python_api: bool = False,
) -> Path:
    """Canonicalize a ligand XYZ file using OpenBabel.

    Tries the ``obabel`` CLI first. If that fails (e.g. dylib conflict from
    mgltools), falls back to the OpenBabel Python API automatically.
    Set *use_python_api=True* to skip the CLI attempt.

    Args:
        input_xyz: Source .xyz file.
        output_xyz: Destination for the canonicalized .xyz.
        obabel_path: Path or command name for the obabel executable.
        use_python_api: If True, use Python API directly instead of CLI.

    Returns:
        The output_xyz path.
    """
    if output_xyz.exists():
        logger.info("Canonicalized XYZ already exists: %s", output_xyz)
        return output_xyz

    output_xyz.parent.mkdir(parents=True, exist_ok=True)

    if not use_python_api:
        try:
            cmd = [
                obabel_path,
                "-ixyz", str(input_xyz),
                "-oxyz",
                "-O", str(output_xyz),
                "--canonical",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                logger.info("Canonicalized XYZ (CLI) → %s", output_xyz)
                return output_xyz
            logger.warning("obabel CLI failed (rc=%d), trying Python API", result.returncode)
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            logger.warning("obabel CLI unavailable (%s), trying Python API", e)

    # Python API fallback
    return _canonicalize_xyz_python_api(input_xyz, output_xyz)


def _canonicalize_xyz_python_api(input_xyz: Path, output_xyz: Path) -> Path:
    """Canonicalize using OpenBabel Python bindings."""
    from openbabel import openbabel as ob

    conv = ob.OBConversion()
    conv.SetInAndOutFormats("xyz", "xyz")
    conv.AddOption("canonical", ob.OBConversion.OUTOPTIONS)

    mol = ob.OBMol()
    conv.ReadFile(mol, str(input_xyz))
    conv.WriteFile(mol, str(output_xyz))
    conv.CloseOutFile()

    logger.info("Canonicalized XYZ (Python API) → %s", output_xyz)
    return output_xyz


def build_graph_from_xyz(xyz_path: Path) -> nx.Graph:
    """Build a NetworkX molecular graph from an XYZ file.

    Node attributes: ``element``, ``xyz``.
    Edge attributes: ``length``.

    Args:
        xyz_path: Path to the .xyz file.

    Returns:
        A NetworkX Graph.
    """
    G = build_molecular_graph(xyz_path)
    logger.info(
        "Built molecular graph: %d atoms, %d bonds",
        G.number_of_nodes(), G.number_of_edges(),
    )
    return G


def prepare_ligand(
    xyz_path: Path,
    output_dir: Path,
    obabel_path: str = "obabel",
) -> dict:
    """Full ligand preparation: canonicalize → build graph.

    Args:
        xyz_path: Input .xyz file for the metal complex.
        output_dir: Working directory.
        obabel_path: obabel executable.

    Returns:
        Dict with keys: ``canonical_xyz``, ``graph``.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = xyz_path.stem

    canonical_xyz = output_dir / f"{stem}_c.xyz"
    canonicalize_xyz(xyz_path, canonical_xyz, obabel_path=obabel_path)

    graph = build_graph_from_xyz(canonical_xyz)

    return {
        "canonical_xyz": canonical_xyz,
        "graph": graph,
    }
