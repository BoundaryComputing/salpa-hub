# HSA + Ferrocene (Sudlow site I)

Docks ferrocene into human serum albumin and reports where it binds.

| | |
|---|---|
| Receptor | **1AO6**, human serum albumin, chain A |
| Ligand | **ferrocene**, Fe(C₅H₅)₂ |
| Site | **Sudlow site I**, subdomain IIA |
| Box | centre `34.64, 32.92, 36.18`, 20 Å cube |
| Charges | GFN1-xTB, neutral, closed shell |
| Runtime | 2–3 minutes |

Everything it needs ships with the package. Install `metaldock-vwf`, open the
template, set a working directory, run. Where each number on this page comes
from is listed at the end.

## The structures

![Ferrocene: an iron atom sandwiched between two cyclopentadienyl rings](figures/metaldock-hsa-ferrocene-ligand.jpg)

Ferrocene as the template ships it. The ten Fe–C contacts fall between 2.064 and
2.066 Å; gas-phase electron diffraction measures 2.064 Å. Nothing declares that
bonding — the graph builder infers it from interatomic distance against covalent
radii, so η⁵ coordination falls out of the geometry.

![The docked pose inside Sudlow site I, contacting residues as sticks, albumin in cartoon](figures/metaldock-hsa-ferrocene-pocket.jpg)

The best-scoring pose. Closest approaches, hydrogens included: ILE290 1.5 Å,
LEU260 1.5 Å, LEU238 1.9 Å, SER287 1.9 Å, ARG257 2.0 Å.

## The pipeline

Six nodes in a line. Each writes into the working directory and hands what it
made to the next, so a parameter left empty is *inherited*, not unset.

### 1 · Protein Prep — `pdb2pqr30`, `prepare_receptor4`

Strips every HETATM record, adds hydrogens at a chosen pH, writes a PDBQT with
partial charges and atom types.

| | |
|---|---|
| takes | `1ao6_A.pdb` — albumin, chain A, 578 residues |
| writes | `clean_1ao6_A.pdb`, `1ao6_A_protonated.pdb`, `clean_1ao6_A.pdbqt` |
| set here | `ph 7.4` · `clean true` |

Stripping all HETATM removes ligands, cofactors and waters alike — so a heme
would go too, which matters if you dock *into* a metalloprotein.

### 2 · Ligand Prep — OpenBabel

Builds the molecular graph: which atoms exist, which are bonded.

| | |
|---|---|
| takes | `ferrocene.xyz` |
| writes | `ferrocene_c.xyz`, `mol_graph.json` — 21 atoms, 30 bonds |
| set here | `metal_symbol Fe` |

The 30 bonds are 10 C–H, 10 C–C and 10 Fe–C.

### 3 · QM Charges — GFN1-xTB

Computes per-atom partial charges and bond orders, and writes them onto the
graph. Docking scores electrostatics, and a metal centre needs charges no force
field supplies.

| | |
|---|---|
| takes | `mol_graph.json`, `ferrocene_c.xyz` — both inherited |
| writes | `enriched_graph.json` — CM5 charges, Wiberg bond orders |
| set here | `engine xtb` · `charge 0` · `spin 0` |

CM5 charges require **GFN1**. GFN2 returns Mulliken only.

### 4 · Ligand PDBQT

Converts the enriched graph into AutoDock's ligand format and decides which
bonds may rotate. Bonds inside the coordination sphere are frozen so the complex
cannot pull itself apart during the search.

| | |
|---|---|
| takes | `enriched_graph.json` |
| writes | `hsa_fe_ligand.pdbqt` |
| set here | `freeze_coordination_sphere true` · `vacant_site false` |

A vacant coordination site is marked with a `DD` dummy atom, which AutoDock
scores as a hydrogen-bond donor. Ferrocene is saturated, so the template leaves
`vacant_site` off. Before 0.4.2 it was on, and every run carried a dummy 1.0 Å
from the iron; turning it off moved the best pose by 0.02 Å and its score from
−3.01 to −3.00 kcal/mol.

### 5 · AutoDock Run — `autogrid4`, `autodock4`

Builds the energy grids, then runs the genetic-algorithm search.

| | |
|---|---|
| takes | `hsa_fe_ligand.pdbqt`, `clean_1ao6_A.pdbqt`, `enriched_graph.json` |
| writes | `*.map`, `hsa_fe_ligand_clean_1ao6_A.dlg`, `hsa_fe_ligand_1…10.pdbqt` |
| set here | `box_center 34.64,32.92,36.18` · `box_size 20,20,20` · `num_poses 10` |

**The box centre is the choice that decides the result.** Left empty it centres
on the ligand's own metal atom, which is only meaningful if the ligand's
coordinates came out of the receptor's crystal. Here it is Sudlow site I, taken
from R-warfarin co-crystallised in 2BXD and transferred into the 1AO6 frame by
superposition over all 578 Cα atoms (RMSD 0.88 Å). Centring instead on the
centroid of Trp214 — the residue that lines the site — lands 9.7 Å away, puts the
pocket on the box boundary, and returns a weaker pose (−2.63 kcal/mol, residues
328–354) in the wrong subdomain without any error.

**On the iron.** Fe, Zn and Mn are the three metals MetalDock supplies no fitted
parameters for; they fall back to AutoDock 4's stock `atom_par` line, while every
other metal gets four re-fitted pairwise well depths. Fe's stock line is
`Rii 1.30 Å, epsii 0.010 kcal/mol` — 87x shallower than Mn's 0.875, 15x shallower
than aliphatic carbon's 0.150 — so the iron adds almost no dispersion and reaches
the score through its xTB charge instead. Here
that is fine: the iron sits between two Cp rings and barely touches protein, so
the carbons do the binding. Read the number with more caution for an Fe complex
whose metal is exposed and coordinating protein donors.

### 6 · Results Analysis

Reads the docking log: binding energy per pose, ligand efficiency, and which
residues line each pose.

| | |
|---|---|
| takes | the `.dlg`, `clean_1ao6_A.pdb` |
| writes | `hsa_fe_analysis.json` |
| set here | `cutoff 4.0` · `reference_xyz` empty |

`reference_xyz` is empty because there is no crystallographic pose to compare
against, so no RMSD is computed.

## Looking at the result

Two files, already in the same coordinate frame:

```bash
vmd protein/clean_1ao6_A.pdb docking/hsa_fe_ligand_1.xyz
```

`_1` … `_10` are numbered in the order AutoDock ran them, not by score;
`binding_energies` in `hsa_fe_analysis.json` lists them in the same order, so the
best is the one with the most negative value. Here all ten land in the same place
(within 0.08 Å of each other), so `_1` shows the result. The `.pdbqt` poses hold
identical coordinates and also open, but VMD needs telling they are PDB
(`mol new … type pdb`).

Of the rest of `docking/`: the `.dlg` is the log, `.gpf`/`.dpf` are inputs, and
the `.map` grids will **not** open in VMD — they are AutoGrid format, not CCP4.
`clean_1ao6_A.maps.xyz` holds the box extents, which read back as the centre and
size set above.

## Expected result

```
poses 10   best −3.00   median −3.00   spread under 0.01 kcal/mol

contacts 11–12:
  TYR150  ARG222  LEU238  ARG257  LEU260  ALA261
  ILE264  LYS286  SER287  HIS288  ILE290  ALA291
```

The iron lands 3.6 Å from the box centre, which is where warfarin sits in 2BXD:
ferrocene docks into the warfarin pocket. The energy is modest because ferrocene
is small and neutral; the pose is the result worth reading.

The exact contact count varies between runs — AutoDock seeds its search from
`pid time`, and HIS288 is the one that comes and goes. The binding energy, the
convergence and the pocket do not vary; judge a re-run on those.

## Regenerating the figures

`workflows/vmd/` holds the TCL. Both renders come from a real run of this
template, using VMD's Tachyon renderer with no display attached.

## Where these numbers come from

- **Runtime, energies, contacts, the iron's place and the spread of the poses:**
  three runs of this template (metaldock-vwf 0.4.2) on an Intel Mac on
  2026-09-24, each node run by `bocoflow_core.node_runner` in the package's own
  environment, as the app runs it. They took 148–155 s, 135–142 s of it in
  AutoDock. On an earlier version of the template, an Apple Silicon Mac mini
  running it under Rosetta took 155 s and 191 s (2026-09-04), and a Windows PC
  through WSL2 about two minutes (2026-09-05).
- **Closest approaches:** the run of 2026-08-25 that rendered the figure; the
  2026-09-24 runs give the same five to within 0.07 Å.
- **Fe–C distances:** the shipped `ferrocene.xyz`, a GFN1-xTB geometry. The
  measured value is gas-phase electron diffraction: Haaland, A.; Nilsson, J. E.
  *Acta Chem. Scand.* 1968, 22, 2653.
- **Box centre, the 0.88 Å superposition and the Trp214 attempt:** our own site
  selection on 2026-08-25, recorded in `CHANGELOG.md` 0.4.0. 2BXD is from
  Ghuman, J. et al. *J. Mol. Biol.* 2005, 353, 38–52; 1AO6 is from Sugio, S.
  et al. *Protein Eng.* 1999, 12, 439–446.
- **Well depths:** AutoDock 4's parameter library, the `atom_par` lines in
  `_vendor/metaldock_modules/data/metal_dock.dat`.
