# PDB2PQR Demo Data

## mini.pdb

A minimal peptide structure (Ala-Ala-Gly, 15 atoms) for testing PDB2PQR conversion.

- **Residues**: ALA-ALA-GLY (3 residues, chain A)
- **Atoms**: 15 heavy atoms
- **Purpose**: Quick functional testing of pdb2pqr protonation and PQR output

This file is used by the unit tests and E2E tests. It runs fast through pdb2pqr
(seconds, not minutes) and produces a valid PQR output with added hydrogens.
