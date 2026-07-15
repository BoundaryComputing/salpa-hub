# PDB Tools Clean

A BoCoFlow node for cleaning and preparing PDB files using [pdb-tools](https://github.com/haddocking/pdb-tools).

## Overview

This node performs common PDB file cleaning operations:
- **Select chains** - Keep only specific chains from multi-chain structures
- **Remove HETATM** - Remove heteroatoms (water, ligands, ions, cofactors)
- **Remove hydrogens** - Strip hydrogen atoms from the structure
- **Renumber residues** - Renumber residues starting from a specific number
- **Tidy PDB** - Ensure valid PDB format with proper records

## Use Cases

- Prepare protein structures for molecular dynamics simulations
- Clean experimental structures for computational analysis
- Extract specific chains from protein complexes
- Remove crystallographic waters and ligands

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_pdb` | string | (required) | Path to input PDB file |
| `chains` | string | "" | Comma-separated chain IDs to keep (e.g., "A,B"). Empty = all chains |
| `remove_hetatm` | boolean | true | Remove HETATM records (water, ligands, etc.) |
| `remove_hydrogens` | boolean | false | Remove hydrogen atoms |
| `renumber_residues` | integer | 0 | Starting residue number (0 = no renumbering) |
| `output_suffix` | string | "_clean" | Suffix for output filename |

## Output

| Field | Type | Description |
|-------|------|-------------|
| `output_pdb` | string | Path to cleaned PDB file |
| `chains_selected` | list | Chain IDs in the output file |
| `atoms_before` | integer | Number of atoms in input |
| `atoms_after` | integer | Number of atoms in output |
| `atoms_removed` | integer | Number of atoms removed |
| `operations` | list | List of operations performed |

## Example Usage

### Basic Cleaning
Remove all heteroatoms and tidy the PDB:
```
Input PDB: 1brs.pdb
Remove HETATM: true
Output: 1brs_clean.pdb
```

### Select Specific Chains
Keep only chains A and D:
```
Input PDB: 1brs.pdb
Chains: A,D
Remove HETATM: true
Output: 1brs_clean.pdb (contains only chains A and D)
```

### Full Cleanup Pipeline
```
Input PDB: structure.pdb
Chains: A
Remove HETATM: true
Remove Hydrogens: true
Renumber Residues: 1
Output: structure_clean.pdb
```

## About pdb-tools

[pdb-tools](https://github.com/haddocking/pdb-tools) is a collection of Python scripts for manipulating PDB files, developed by the [Bonvin Lab](https://www.bonvinlab.org/) at Utrecht University.

### Citation

If you use this node in your research, please cite:

> Rodrigues, J. P. G. L. M., Teixeira, J. M. C., Trellet, M. & Bonvin, A. M. J. J.
> pdb-tools: a swiss army knife for molecular structures.
> F1000Research 7, 1961 (2018).
> DOI: [10.12688/f1000research.17456.1](https://doi.org/10.12688/f1000research.17456.1)

## Requirements

- Python >= 3.9
- pdb-tools >= 2.5.0 (installed via conda-forge)

## Installation

Install from the BoCoFlow Marketplace:

1. Open BoCoFlow
2. Navigate to **Marketplace**
3. Search for "pdb-tools-clean"
4. Click **Install**

Or manually copy this folder into your BoCoFlow installed-nodes directory.

## License

MIT License
