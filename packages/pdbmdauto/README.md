# pdbmdauto — automated protein structure preparation for molecular dynamics

> **New to molecular dynamics?** Start with [`TUTORIAL.md`](TUTORIAL.md): the same pipeline in
> plain language, with a worked example. This README is the formal reference — scope,
> methodology, node table, outputs, references — and the two are kept apart on purpose.

Automated preparation of protein structures for molecular dynamics (MD) simulation: homology-based
reconstruction of unresolved residues, pH-dependent protonation and force-field topology generation,
solvation and ionization, and staged energy minimization and equilibration. The complete preparation
is provided as a single, reproducible workflow.

## Overview

`pdbmdauto` converts an experimentally determined protein structure — supplied either as a Protein
Data Bank (PDB) identifier or as a local coordinate file — into a solvated, energy-relaxed system
suitable for production molecular dynamics. It integrates established structural-bioinformatics and
simulation tools (Biopython, ProMod3/OpenStructure, PDB2PQR/PROPKA, GROMACS) into a package of
thirteen nodes. The bundled workflow chains eleven of these steps in a fixed sequence; the remaining
three provide alternatives for specific scenarios (see *Pipeline nodes*). All nodes run in one
shared, version-pinned environment provisioned by [pixi](https://pixi.sh), and every intermediate
product is written to disk for inspection.

## Scope and intended use

The package addresses the structure-preparation stage that precedes molecular dynamics, which is
typically the most labor-intensive and error-prone part of an MD study. Experimentally determined
structures are frequently incomplete (unresolved residues and loops), lack explicit hydrogen atoms,
and carry no description of protonation state or solvent environment — all of which must be resolved
before simulation. `pdbmdauto` performs these steps in a fixed, reproducible order and retains the
intermediate files at each stage, supporting inspection and method reporting.

Two scope limits apply. Water, ligands, cofactors and other hetero residues present in the
experimental structure are removed when the chains are merged (step 4); a system requiring a bound
ligand needs its own parameterization outside this package. And the production run is deliberately
short (2 ps), an exercise of the prepared system rather than sampling.

## Methodology

The pipeline is organized into four stages.

**1. Reconstruction of unresolved residues (homology modeling).** The deposited structure is
retrieved and its per-chain sequences extracted. Residues present in the sequence but absent from
the coordinates — the entry's `REMARK 465` records — are identified, an alignment of the resolved
residues to the full sequence is written per chain, and the missing residues are built by ProMod3:
internal gaps from its fragment database, terminal extensions by its terminus modeling, followed by
sidechain reconstruction and a brief minimization. The result is a chemically complete, continuous
model. Reconstructed residues are predictions, not measurements, and should be treated as such in
any analysis.

**2. Protonation and topology generation.** Protonation states of titratable residues are predicted
at the requested pH (PROPKA); hydrogen atoms and partial charges are assigned (PDB2PQR, AMBER
naming); and a force-field topology is generated (GROMACS `pdb2gmx`, amber99sb, TIP3P water),
followed by two steepest-descent minimizations in vacuo (unconstrained, then with H-bond
constraints). Because crystallographic data seldom resolve hydrogen positions, and because
protonation is pH-dependent, this stage substantially determines the electrostatics of the
simulated system.

**3. Restrained relaxation, solvation and ionization.** Atoms originating from the experimental
structure are distinguished from reconstructed atoms by two index groups (heavy atoms and backbone
atoms of the resolved residues). The reconstructed regions are relaxed under position restraints on
the experimental atoms (two short NVT runs at 300 K, then two conjugate-gradient minimizations),
the system is placed in a 5 nm cubic box, solvated, neutralized and brought to 0.15 M NaCl, and
minimized once more in water.

**4. Production simulation.** A short production MD run is performed locally (GROMACS `mdrun`,
1000 steps of 2 fs at 300 K, velocity rescaling, no pressure coupling); an alternative node submits
the same step to an HPC scheduler (SLURM) for extended simulation.

## Pipeline nodes

The bundled workflow comprises eleven steps in the order below. The "GMX MD Relaxation" node is
applied twice — once before solvation (on the reconstructed regions, under restraints, in vacuo)
and once after solvation (on the full system, minimization only).

| # | Node | Function | Primary tool |
|---|------|----------|--------------|
| 1 | PDB FASTA Parser | Retrieve structure; extract per-chain sequences; list unresolved residues | Biopython, RCSB PDB |
| 2 | Generate Alignment | Align the resolved sequence to the full reference sequence, per chain | — |
| 3 | Multi-Chain Alignment | Consolidate per-chain alignments for multi-chain structures | — |
| 4 | Merge PDB Chains | Combine selected chains into a single model; strip water and hetero residues | Biopython |
| 5 | Fix Missing Residues | Reconstruct unresolved residues and termini by homology modeling | ProMod3 / OpenStructure |
| 6 | pKa + GROMACS EM | Predict protonation at target pH; assign hydrogens and charges; generate topology; minimize in vacuo | PROPKA, PDB2PQR, GROMACS |
| 7 | Original Atom Groups | Define index groups distinguishing experimental from reconstructed atoms | GROMACS |
| 8 | GMX MD Relaxation (restrained) | Relax reconstructed regions under position restraints on experimental atoms | GROMACS |
| 9 | GMX Solvate & Ionize | Cubic box, solvent, neutralizing counter-ions at 0.15 M | GROMACS |
| 10 | GMX MD Relaxation (solvated) | Minimize the solvated system | GROMACS |
| 11 | GROMACS MD Run (Local) | Execute production molecular dynamics | GROMACS |

Nodes included in the package but not part of the bundled workflow: a standalone **PDB2PQR** node
(protonation only, without topology generation), **Generate GROMACS Index**, an alternative
index-group builder, and **GROMACS MD Run** with HPC/SLURM submission.

## Example workflow (4Z8J)

The bundled template `workflows/pdbmdauto-pipeline.json` prepares PDB entry 4Z8J — the SNX27 PDZ
domain in complex with the C-terminal PDZ-binding motif of the parathyroid hormone receptor, at
0.95 Å — at pH 7. The entry has two chains (the PDZ domain, 101 residues in sequence, and an
8-residue peptide) and six residues absent from the coordinates: Gly33–Gly37 of the domain and
Gln586 of the peptide, all at N-termini. The workflow rebuilds them, protonates the complex, and
delivers a solvated system of 12,193 atoms (3,509 waters, 12 Na⁺ and 11 Cl⁻, net charge zero) in a
5 nm cube, with a 2 ps production trajectory.

The workflow requires no manual configuration: the structure is retrieved from the RCSB PDB at run
time (network access is required for this example), and the only input at load is a working
directory. To prepare a different target, set the PDB identifier (or a local file) and the pH on
the first node and execute. Its walkthrough, `workflows/pdbmdauto-pipeline.md`, describes every
step's inputs, outputs and parameters, with figures from a reference run.

### Measured run times

| Machine | Run | Time |
|---|---|---|
| Intel MacBook Pro (Core i7-9750H, 6 cores, macOS 15.7), Salpa 0.3.1 | full pipeline, three runs | **3 to 5½ min** (179 s, 313 s, 329 s) |
| same | first install, empty package cache (3.1 GB environment) | 3½ min |
| Apple Silicon Mac Mini (M-series, native arm64) | full pipeline | ≈ 2 min |
| Windows 11, WSL2 (virtual machine) | first install | ≈ 20 min |

The spread on the Intel laptop is thermal: its GROMACS stages ran 2.5× slower on a back-to-back run
than on a cold machine. The restrained relaxation (step 8) is the dominant cost.

## Outputs and interpretation

All files are written under `<working directory>/pdbmdauto-e2e-full/e2e_4z8j/` for the bundled
example (the directory is named after the workflow's case name). The files to examine, by stage:

- `4Z8J.pdb`, `chain_A.fasta`, `chain_B.fasta`, `missing_residues_chain_*.csv` — the retrieved
  structure, its sequences, and the unresolved residues (step 1).
- `Merge/merge.pdb` (the experimental chains, hetero atoms removed) and `Merge/fixed.pdb` (the
  completed model; ProMod3 renumbers residues from 1) — steps 4–5. Confirm continuity of the
  rebuilt regions and the absence of chain breaks.
- `gmx/propka.pqr`, `gmx/protonated.pdb`, `gmx/pdb2gmx.top`, `gmx/em_hbonds.log` — protonation,
  topology and the in-vacuo minimization (step 6). The log reports the final potential energy and
  maximum force; in the reference run the minimization converged to Fmax < 1000 kJ mol⁻¹ nm⁻¹ in
  15 steps.
- `gmx/index.ndx` — the `OriHeavy` and `OriBackBone` groups (step 7); their sizes should equal the
  heavy-atom and backbone-atom counts of the experimentally resolved residues.
- `gmx/ion.gro`, `gmx/topol.top` — the solvated, neutralized system (step 9). Confirm net-zero
  charge in the topology's molecule counts.
- `gmx/em.log`, `gmx/em.edr` — the solvated minimization (step 10). With the shipped `em.mdp`
  (500 steps, emtol 500) the run stops at its step limit with Fmax ≈ 3 × 10³ kJ mol⁻¹ nm⁻¹; this
  is sufficient to relieve the worst contacts before the short demonstration run, and a production
  study should raise `nsteps` until the tolerance is met.
- `gmx/md.gro`, `gmx/md.trr`, `gmx/md.edr`, `gmx/md.log` — the production coordinates, trajectory
  (11 frames at 100-step intervals), energies and log (step 11).

Recommended validation: minimization energies decrease and stabilize (monotonic divergence indicates
steric clashes or an ill-formed system); the system is neutral after ionization; the rebuilt
residues are continuous with their chains. The reconstructed residues are model-derived and should
be interpreted with more caution than experimentally resolved regions; force-field models are
approximate; and the 2 ps run exercises the pipeline rather than providing converged sampling —
production studies require substantially longer simulation.

## Requirements and platform support

Supported on Linux (x86-64 and ARM64) and on macOS on both Intel and Apple Silicon processors,
with native builds of GROMACS 2026 and ProMod3 3.6 (no Rosetta). Native Windows is not supported —
GROMACS and ProMod3 have no Windows conda packages — so on Windows the application provisions the
package inside WSL2 automatically; the first installation there is correspondingly longer
(approximately 20 minutes on a virtual machine), after which the package runs like any other. The
package installs one shared, version-pinned environment of about 3.1 GB via pixi; no manual
dependency installation is required. The bundled example needs network access to retrieve 4Z8J.

## References

Cite each tool according to its own guidance.

- **GROMACS** — Abraham, M. J., *et al.* (2015). *GROMACS: High performance molecular simulations
  through multi-level parallelism from laptops to supercomputers.* SoftwareX 1–2, 19–25.
  https://doi.org/10.1016/j.softx.2015.06.001
- **ProMod3** — Studer, G., Tauriello, G., Bienert, S., Biasini, M., Johner, N., Schwede, T. (2021).
  *ProMod3—A versatile homology modelling toolbox.* PLoS Computational Biology 17(1), e1008667.
  https://doi.org/10.1371/journal.pcbi.1008667
- **OpenStructure** — Biasini, M., *et al.* (2013). *OpenStructure: an integrated software framework
  for computational structural biology.* Acta Crystallographica D 69, 701–709.
  https://doi.org/10.1107/S0907444913007051
- **PDB2PQR** — Dolinsky, T. J., *et al.* (2004). *PDB2PQR: an automated pipeline for the setup of
  Poisson–Boltzmann electrostatics calculations.* Nucleic Acids Research 32, W665–W667.
  https://doi.org/10.1093/nar/gkh381 — and Jurrus, E., *et al.* (2018). *Improvements to the APBS
  biomolecular solvation software suite.* Protein Science 27, 112–128.
- **PROPKA** — Olsson, M. H. M., Søndergaard, C. R., Rostkowski, M., Jensen, J. H. (2011). *PROPKA3:
  consistent treatment of internal and surface residues in empirical pKa predictions.* Journal of
  Chemical Theory and Computation 7, 525–537. https://doi.org/10.1021/ct100578z
- **Biopython** — Cock, P. J. A., *et al.* (2009). *Biopython: freely available Python tools for
  computational molecular biology and bioinformatics.* Bioinformatics 25, 1422–1423.
  https://doi.org/10.1093/bioinformatics/btp163
- **4Z8J** — Clairfeuille, T., Teasdale, R. D., Collins, B. M., Pavlos, N. *Crystal structure of the
  SNX27 PDZ domain bound to the C-terminal PTHR PDZ binding motif.* RCSB Protein Data Bank entry
  4Z8J (2015). https://www.rcsb.org/structure/4Z8J

---

_Part of the Salpa Hub. Licensed under MIT; the third-party tools it invokes are installed, not
redistributed — see `NOTICE` and the Hub's `LICENSING.md`._
