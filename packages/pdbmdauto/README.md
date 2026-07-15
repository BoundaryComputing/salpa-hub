# pdbmdauto — automated protein structure preparation for molecular dynamics

Automated preparation of protein structures for molecular dynamics (MD) simulation: homology-based reconstruction of missing residues, pH-dependent protonation and force-field topology generation, solvation and ionization, and staged energy minimization and equilibration. The complete preparation is provided as a single, reproducible workflow.

## Overview

`pdbmdauto` converts an experimentally determined protein structure — supplied either as a Protein Data Bank (PDB) identifier or as a local coordinate file — into a solvated, equilibrated system suitable for production molecular dynamics. It integrates established structural-bioinformatics and simulation tools (Biopython, ProMod3, PDB2PQR/PROPKA, GROMACS) into a package of fourteen nodes. The default workflow chains eleven of these steps in a fixed sequence; the remaining nodes provide alternatives for specific scenarios. Each node executes in an isolated, version-pinned environment (via pixi), and intermediate products are written to disk for inspection at every stage.

## Scope and intended use

The package addresses the structure-preparation stage that precedes molecular dynamics, which is typically the most labor-intensive and error-prone part of an MD study. Experimentally determined structures are frequently incomplete (unresolved residues and loops), lack explicit hydrogen atoms, and carry no description of protonation state or solvent environment — all of which must be resolved before simulation. `pdbmdauto` performs these steps in a fixed, reproducible order and retains the intermediate files at each stage, supporting inspection and method reporting.

## Methodology

The pipeline is organized into four stages.

**1. Reconstruction of missing residues (homology modeling).** The deposited structure is retrieved and its per-chain sequences are extracted. Residues absent from the coordinates are identified against the full reference sequence and reconstructed by homology modeling (ProMod3), yielding a chemically complete and continuous model prior to parameterization.

**2. Protonation and topology generation.** Protonation states of titratable residues are predicted at a specified pH (PROPKA); hydrogen atoms and partial charges are assigned (PDB2PQR); and a force-field topology is generated (GROMACS `pdb2gmx`), followed by staged energy minimization. Because crystallographic data seldom resolve hydrogen positions, and because protonation is pH-dependent, this stage substantially determines the electrostatics of the simulated system.

**3. Solvation, ionization, and equilibration.** Atoms originating from the experimental structure are distinguished from reconstructed atoms so that the former may be positionally restrained while the modeled regions relax. The system is solvated, neutralized with counter-ions at physiological ionic strength, and equilibrated through minimization and restrained molecular dynamics.

**4. Production simulation.** A production MD run is performed locally (GROMACS `mdrun`); an alternative node provides HPC/SLURM execution for extended simulations.

## Pipeline nodes

The default workflow comprises eleven steps in the order below. The "GMX MD Relaxation" node is applied twice — once before solvation (on the reconstructed regions, under restraints) and once after solvation (on the full system).

| # | Node | Function | Primary tool |
|---|------|----------|--------------|
| 1 | PDB FASTA Parser | Retrieve structure; extract per-chain sequences; identify missing residues | Biopython, RCSB PDB |
| 2 | Generate Alignment | Align the resolved sequence to the full reference sequence, per chain | — |
| 3 | Multi-Chain Alignment | Consolidate per-chain alignments for multi-chain structures | — |
| 4 | Merge PDB Chains | Combine selected chains into a single model | Biopython |
| 5 | Fix Missing Residues | Reconstruct missing residues and loops by homology modeling | ProMod3 |
| 6 | pKa + GROMACS EM | Predict protonation at target pH; assign hydrogens and charges; generate topology; energy-minimize | PROPKA, PDB2PQR, GROMACS |
| 7 | Original Atom Groups | Define index groups distinguishing experimental from reconstructed atoms | GROMACS |
| 8 | GMX MD Relaxation (restrained) | Minimize and equilibrate reconstructed regions under position restraints on experimental atoms | GROMACS |
| 9 | Solvate & Ionize | Add solvent box and neutralizing counter-ions | GROMACS |
| 10 | GMX MD Relaxation (solvated) | Minimize and equilibrate the solvated system | GROMACS |
| 11 | GROMACS MD Run (Local) | Execute production molecular dynamics | GROMACS |

Nodes included in the package but not part of the default pipeline: a standalone PDB2PQR node (protonation only), PDB Clean (chain selection and removal of heteroatoms/solvent), an alternative index-group builder, and a GROMACS MD Run node with HPC/SLURM support.

## Example workflow (4Z8J)

The bundled template `workflows/pdbmdauto-pipeline.json` prepares PDB entry 4Z8J (the SNX27 PDZ domain in complex with a parathyroid hormone receptor C-terminal peptide; two chains; five unresolved residues) at pH 7. The workflow requires no manual configuration — the structure is retrieved automatically — and completes within a few minutes for this system on a typical workstation. To prepare a different target, set the input identifier (or a local file) and the pH on the first node and execute.

## Outputs and interpretation

The pipeline produces a solvated, equilibrated system and a production trajectory in standard GROMACS formats (`.gro` coordinates, `.xtc` trajectory, `.edr` energies), together with the corresponding topology. Recommended validation:

- Confirm continuity of the reconstructed regions and the absence of chain breaks.
- Verify that minimization and equilibration energies decrease and stabilize; monotonic divergence indicates steric clashes or an ill-formed system.
- Confirm net-zero system charge following ionization.
- The reported missing-residue records identify the reconstructed regions; these are model-derived and should be interpreted with greater caution than experimentally resolved regions.

Reconstructed loops are predictions rather than measurements, force-field models are approximate, and the short demonstration run is intended to exercise the pipeline rather than to provide converged sampling; production studies require substantially longer simulation.

## Requirements and platform support

Supported on Linux and macOS (Apple Silicon under Rosetta). Native Windows is not currently supported, as GROMACS and ProMod3 have no Windows conda packages; WSL2 or Docker is recommended in the interim (issue #48). Each node provisions its own version-pinned environment via [pixi](https://pixi.sh); no manual dependency installation is required.

## References

Cite each tool according to its own guidance.

- **GROMACS** — Abraham, M. J., *et al.* (2015). *GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers.* SoftwareX 1–2, 19–25. https://doi.org/10.1016/j.softx.2015.06.001
- **ProMod3 / SWISS-MODEL** — https://swissmodel.expasy.org
- **PDB2PQR** — Dolinsky, T. J., *et al.* (2004). *PDB2PQR: an automated pipeline for the setup of Poisson–Boltzmann electrostatics calculations.* Nucleic Acids Research. https://www.poissonboltzmann.org
- **PROPKA** — Olsson, M. H. M., *et al.* (2011). *PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions.* Journal of Chemical Theory and Computation.
- **Biopython** — Cock, P. J. A., *et al.* (2009). *Biopython: freely available Python tools for computational molecular biology and bioinformatics.* Bioinformatics.

---

_Part of Salpa Hub. Licensing for this package is being finalized; see the Hub's `LICENSING.md`._
