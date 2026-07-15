# PDB FASTA Parser (pdb-fasta-biopython)

Parse PDB protein structures and extract FASTA sequences using BioPython.

## Features

- **RCSB API fetch**: Download PDB + FASTA from RCSB by PDB ID (e.g. "3LZ0")
- **Local file processing**: Parse local PDB files
- **Chain splitting**: Create separate FASTA files for each chain
- **HETATM handling**: Optionally include ligands and modified residues
- **Missing residues**: Extract structural gaps from PDB header
- **Protein + DNA/RNA**: Handles both macromolecule types

## Usage

1. Add the **PDB FASTA Parser** node to your workflow
2. Choose input mode:
   - **PDB ID**: Enter a 4-character PDB identifier (fetches from RCSB)
   - **Local file**: Select a local PDB file
3. Configure options (chain splitting, HETATM, missing residues)
4. Set the output directory
5. Execute

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| Input Mode | Select | `pdb_id` | `pdb_id` (RCSB fetch) or `local_file` |
| PDB ID | String | - | 4-char RCSB PDB ID (e.g. "3LZ0") |
| PDB File | File | - | Path to local PDB file |
| Case Name | String | auto | Identifier (falls back to PDB ID or "protein") |
| Output Directory | Folder | - | Where to write output files |
| Split by Chain | Boolean | true | Create per-chain FASTA files |
| Include HETATM | Boolean | false | Include ligands/modified residues |
| Check PDB Header | Boolean | false | Extract missing residues info |

## Output

### PDB ID mode
- `{PDB_ID}.pdb` - Downloaded PDB file
- `{PDB_ID}_rcsb.fasta` - RCSB reference FASTA
- `{case}_chain_{X}.fasta` - Per-chain FASTA (if split)
- `missing_residues_chain_{X}.csv` - Missing residues (if enabled)

### Local file mode
- `{case}_chain_{X}.fasta` - Per-chain FASTA (if split)
- `{case}.fasta` - Combined FASTA (if not split)
- `missing_residues_chain_{X}.csv` - Missing residues (if enabled)

## Demo Data

The `demo_data/` directory contains structure 3LZ0 (a multi-chain protein) for testing.

## Technical Details

| Property | Value |
|----------|-------|
| Execution Strategy | PIXI_SUBPROCESS |
| Dependencies | biopython, requests |
| Base Class | Node |
| Ports | 1 in, 1 out |
