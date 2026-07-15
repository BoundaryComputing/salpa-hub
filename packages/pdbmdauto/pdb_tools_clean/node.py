"""
pdb-tools-clean - PDB File Cleaning Node

Uses pdb-tools (https://github.com/haddocking/pdb-tools) to:
- Select specific chains from a PDB file
- Remove HETATM records (heteroatoms like water, ligands)
- Produce a valid, tidy PDB file

This is useful for preparing PDB structures for molecular simulations
or other computational analyses that require clean protein structures.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from bocoflow_core.node import NodeBase, NodeInputParams, NodeOutputParams


class PDBToolsClean(NodeBase):
    """
    Clean PDB files using pdb-tools.

    This node performs a common PDB cleaning pipeline:
    1. Select specific chain(s) (optional)
    2. Remove HETATM records (water, ligands, etc.)
    3. Tidy the PDB file to ensure valid format

    Uses the pdb-tools suite: https://www.bonvinlab.org/pdb-tools/
    """

    def execute(self, input: NodeInputParams) -> NodeOutputParams:
        """
        Execute the PDB cleaning pipeline.

        Parameters (from UI):
            input_pdb: Path to input PDB file
            chains: Comma-separated chain IDs to select (e.g., "A,B"). Empty = all chains.
            remove_hetatm: Whether to remove HETATM records (default: True)
            remove_hydrogens: Whether to remove hydrogen atoms (default: False)
            renumber_residues: Starting residue number for renumbering (0 = no renumbering)
            output_suffix: Suffix for output file (default: "_clean")

        Returns:
            output_pdb: Path to the cleaned PDB file
            chains_selected: List of chains in output
            atoms_removed: Number of atoms removed
        """
        params = input.parameters
        context = input.context or {}

        # Get parameters
        input_pdb = params.get("input_pdb", "")
        chains = params.get("chains", "").strip()
        remove_hetatm = params.get("remove_hetatm", True)
        remove_hydrogens = params.get("remove_hydrogens", False)
        renumber_residues = params.get("renumber_residues", 0)
        output_suffix = params.get("output_suffix", "_clean")

        # Validate input
        if not input_pdb:
            return NodeOutputParams(
                output_data={"error": "No input PDB file specified"},
                status="error",
                message="Please specify an input PDB file"
            )

        # Resolve input path
        working_path = context.get("working_path", ".")
        input_path = self._resolve_path(input_pdb, working_path)

        if not input_path.exists():
            return NodeOutputParams(
                output_data={"error": f"Input file not found: {input_path}"},
                status="error",
                message=f"File not found: {input_pdb}"
            )

        # Generate output path
        output_pdb = self._generate_output_path(input_path, output_suffix)

        try:
            # Build the pdb-tools pipeline
            pipeline = self._build_pipeline(
                input_path=input_path,
                chains=chains,
                remove_hetatm=remove_hetatm,
                remove_hydrogens=remove_hydrogens,
                renumber_residues=renumber_residues
            )

            # Count atoms before
            atoms_before = self._count_atoms(input_path)

            # Execute the pipeline
            self._execute_pipeline(pipeline, input_path, output_pdb)

            # Count atoms after
            atoms_after = self._count_atoms(output_pdb)
            atoms_removed = atoms_before - atoms_after

            # Get chains in output
            chains_in_output = self._get_chains(output_pdb)

            result = {
                "output_pdb": str(output_pdb),
                "input_pdb": str(input_path),
                "chains_selected": chains_in_output,
                "atoms_before": atoms_before,
                "atoms_after": atoms_after,
                "atoms_removed": atoms_removed,
                "operations": self._describe_operations(
                    chains, remove_hetatm, remove_hydrogens, renumber_residues
                )
            }

            return NodeOutputParams(
                output_data=result,
                status="success",
                message=f"Cleaned PDB saved to {output_pdb.name} ({atoms_removed} atoms removed)"
            )

        except subprocess.CalledProcessError as e:
            return NodeOutputParams(
                output_data={"error": str(e), "stderr": e.stderr},
                status="error",
                message=f"pdb-tools error: {e.stderr}"
            )
        except Exception as e:
            return NodeOutputParams(
                output_data={"error": str(e)},
                status="error",
                message=f"Error: {e}"
            )

    def _resolve_path(self, path_str: str, working_path: str) -> Path:
        """Resolve path relative to working directory."""
        path = Path(path_str)
        if path.is_absolute():
            return path
        return Path(working_path) / path

    def _generate_output_path(self, input_path: Path, suffix: str) -> Path:
        """Generate output path with suffix."""
        stem = input_path.stem
        return input_path.parent / f"{stem}{suffix}.pdb"

    def _build_pipeline(
        self,
        input_path: Path,
        chains: str,
        remove_hetatm: bool,
        remove_hydrogens: bool,
        renumber_residues: int
    ) -> list:
        """Build the pdb-tools command pipeline."""
        commands = []

        # Select chains (if specified)
        if chains:
            chain_arg = chains.replace(" ", "")  # Remove spaces
            commands.append(["pdb_selchain", f"-{chain_arg}"])

        # Remove HETATM records
        if remove_hetatm:
            commands.append(["pdb_delhetatm"])

        # Remove hydrogens
        if remove_hydrogens:
            commands.append(["pdb_delelem", "-H"])

        # Renumber residues
        if renumber_residues != 0:
            commands.append(["pdb_reres", str(renumber_residues)])

        # Always tidy at the end for valid PDB
        commands.append(["pdb_tidy"])

        return commands

    def _execute_pipeline(self, pipeline: list, input_path: Path, output_path: Path):
        """Execute the pdb-tools pipeline using shell pipes."""
        if not pipeline:
            # Just copy the file if no operations
            import shutil
            shutil.copy(input_path, output_path)
            return

        # Build shell command with pipes
        # e.g., "pdb_selchain -A,D input.pdb | pdb_delhetatm | pdb_tidy > output.pdb"
        first_cmd = pipeline[0]
        first_cmd_str = " ".join(first_cmd) + f" {input_path}"

        if len(pipeline) > 1:
            rest_cmds = [" ".join(cmd) for cmd in pipeline[1:]]
            full_cmd = first_cmd_str + " | " + " | ".join(rest_cmds)
        else:
            full_cmd = first_cmd_str

        full_cmd += f" > {output_path}"

        # Execute
        result = subprocess.run(
            full_cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )

    def _count_atoms(self, pdb_path: Path) -> int:
        """Count ATOM and HETATM records in a PDB file."""
        count = 0
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        count += 1
        except Exception:
            pass
        return count

    def _get_chains(self, pdb_path: Path) -> list:
        """Get list of chain IDs in a PDB file."""
        chains = set()
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") and len(line) > 21:
                        chain_id = line[21]
                        if chain_id.strip():
                            chains.add(chain_id)
        except Exception:
            pass
        return sorted(list(chains))

    def _describe_operations(
        self,
        chains: str,
        remove_hetatm: bool,
        remove_hydrogens: bool,
        renumber_residues: int
    ) -> list:
        """Describe the operations performed."""
        ops = []
        if chains:
            ops.append(f"Selected chains: {chains}")
        if remove_hetatm:
            ops.append("Removed HETATM records")
        if remove_hydrogens:
            ops.append("Removed hydrogen atoms")
        if renumber_residues != 0:
            ops.append(f"Renumbered residues starting from {renumber_residues}")
        ops.append("Tidied PDB format")
        return ops

    def get_parameters_schema(self) -> dict:
        """Define the node parameters schema."""
        return {
            "type": "object",
            "properties": {
                "input_pdb": {
                    "type": "string",
                    "title": "Input PDB File",
                    "description": "Path to the input PDB file",
                    "format": "file-path"
                },
                "chains": {
                    "type": "string",
                    "title": "Select Chains",
                    "description": "Comma-separated chain IDs to keep (e.g., 'A,B'). Leave empty to keep all chains.",
                    "default": ""
                },
                "remove_hetatm": {
                    "type": "boolean",
                    "title": "Remove HETATM",
                    "description": "Remove heteroatoms (water, ligands, ions, etc.)",
                    "default": True
                },
                "remove_hydrogens": {
                    "type": "boolean",
                    "title": "Remove Hydrogens",
                    "description": "Remove hydrogen atoms from the structure",
                    "default": False
                },
                "renumber_residues": {
                    "type": "integer",
                    "title": "Renumber Residues",
                    "description": "Starting number for residue renumbering (0 = no renumbering)",
                    "default": 0
                },
                "output_suffix": {
                    "type": "string",
                    "title": "Output Suffix",
                    "description": "Suffix to append to output filename",
                    "default": "_clean"
                }
            },
            "required": ["input_pdb"]
        }
