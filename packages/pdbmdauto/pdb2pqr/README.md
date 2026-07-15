# PDB2PQR Node

**Version**: 2.0.0
**Author**: BoCoFlow Development Team
**License**: MIT

## Overview

Converts PDB files to PQR format with protonation states, hydrogen addition, and partial charge assignment using PDB2PQR. Supports PROPKA for pH-dependent pKa predictions and multiple force fields (AMBER, CHARMM, PARSE, etc.).

## Features

- **Hydrogen Addition**: Automatically adds missing hydrogen atoms
- **Protonation States**: pH-dependent protonation via PROPKA
- **Multiple Force Fields**: AMBER, CHARMM, PARSE, TYL06, PEOEPB, SWANSON
- **Charge Assignment**: Partial charges for all atoms
- **PDB Generation**: Optionally converts PQR to protonated PDB via MDAnalysis

## Parameters

### Required

- **Case Name**: Identifier for the conversion job
- **Input PDB**: PDB file to process
- **Output Directory**: Where to save output files

### Optional

- **Force Field**: AMBER (default), CHARMM, PARSE, TYL06, PEOEPB, SWANSON
- **pH Value**: pH for protonation (default: 7.0)
- **Keep Chain IDs**: Preserve chain identifiers (default: true)
- **Optimize Hydrogens**: Optimize hydrogen positions (default: true)
- **Use PROPKA**: Enable pKa predictions (default: true)
- **Generate PDB**: Create protonated PDB from PQR (default: true)

## Output Files

- `{case_name}_structure.pqr` — PQR file with charges and radii
- `{case_name}_protonated.pdb` — Protonated PDB file (optional)
- `{case_name}_propka.out` — PROPKA pKa predictions (if enabled)

## Example Workflow

```
ESMFold Prediction --> PDB2PQR --> OpenMM MD Simulation
```

## Demo Data

The `demo_data/` directory contains `mini.pdb`, a small peptide structure for testing.

## Dependencies

- Python >= 3.9
- MDAnalysis >= 2.4.0
- pdb2pqr >= 3.0.0

## Citation

Dolinsky TJ et al. PDB2PQR: An Automated Pipeline for the Setup of
Poisson-Boltzmann Electrostatics Calculations. *Nucleic Acids Research*, 2004.
DOI: [10.1093/nar/gkm276](https://doi.org/10.1093/nar/gkm276)
