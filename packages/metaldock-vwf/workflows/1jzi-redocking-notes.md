<!-- Deliberately NOT named metaldock-1jzi-re-pipeline.md.
     scripts/sync-workflows.mjs on the website generates a page for
     workflows/<template-id>.md, so a doc under any other name stays in the
     package and off salpa.app. The 1JZI case is a redocking kept for
     reference; the HSA template is the one users are pointed at. -->

# 1JZI redocking — Re(phen)(CO)₃ back into azurin

A redocking: the rhenium complex is taken out of the crystal structure 1JZI and
docked back into it, so the right answer is known and each pose can be scored as
RMSD against it. It checks the method on a solved structure; it predicts nothing.

- **Protein:** 1JZI (*Pseudomonas aeruginosa* azurin)
- **Ligand:** Re(phen)(CO)₃(His83) — 29 atoms, metal = **Re**, one vacant
  coordination site
- **Where the case comes from:** the inputs are the ORCA example in MetalDock's
  repository (`examples/example_runs/vacancy_coordination_sphere/ORCA/`). The
  MetalDock paper (Hakkennes et al., *J. Chem. Inf. Model.* 2023) does not report
  1JZI, so there is no published value to reproduce. Every number on this page
  comes from our own runs, and says which.

## Running it

`workflows/metaldock-1jzi-re-pipeline.json` stays in this package but is not in
the template library: `package.toml` lists only the HSA + Ferrocene template.
To run it, open **Load Workflow → Upload File** in Salpa, choose this file, and
pick a working folder. Nothing to download: the inputs ship inside the nodes'
`demo_data/` directories and the charges come from GFN1-xTB, which installs with
the package.

Everything below describes what that workflow contains, what it gives, and how
to run the same case with DFT charges instead.

## Inputs

| File | Where it lives | Feeds |
|---|---|---|
| `1jzi.pdb` | `mdock_protein_prep/demo_data/` | Protein Prep — *Protein PDB* |
| `1jzi_D_REP.xyz` | `mdock_ligand_prep/demo_data/` | Ligand Prep — *Ligand XYZ* |
| `1jzi_D_REP.xyz` | `mdock_results_analysis/demo_data/` | Results Analysis — *Reference XYZ* (for RMSD) |

Both originate from the MetalDock repository's
`examples/example_runs/vacancy_coordination_sphere/ORCA/`.

## Graph topology

A single linear chain. Each node forwards its predecessor's `data` keys and
appends its own, so the receptor PDBQT made in step 1 is still available to the
docking node in step 5 without a separate edge:

```
[Protein Prep] → [Ligand Prep] → [QM Charges]
        → [Ligand PDBQT] → [AutoDock Run] → [Results Analysis]
```

## Node configuration

| # | Node | Key parameters |
|---|------|----------------|
| 1 | **Protein Prep** | `case_name=1jzi_re`, `pdb_file=node:demo_data/1jzi.pdb`, `output_dir=rel:protein`, `ph=7.0`, `clean=true` (MGLTools auto-detected) |
| 2 | **Ligand Prep** | `xyz_file=node:demo_data/1jzi_D_REP.xyz`, `metal_symbol=Re`, `output_dir=rel:ligand` |
| 3 | **QM Charges** | `engine=xtb`, `geom_opt=false`, `charge=1`, `spin=0`, `ncpu=4`, `output_dir=rel:qm` |
| 4 | **Ligand PDBQT** | `vacant_site=true`, `max_torsions=32`, `freeze_coordination_sphere=true` (metal symbol + graph inherited) |
| 5 | **AutoDock Run** | `num_poses=10`, `box_center=1.65,-7.803,27.176`, `box_size=20,20,20` (ligand/receptor PDBQT + graph inherited) |
| 6 | **Results Analysis** | `cutoff=4.0`, `reference_xyz=node:demo_data/1jzi_D_REP.xyz` (protein PDB + poses + heavy-atom count inherited) |

Leave the file-input fields on nodes 3–5 **empty** — they auto-discover
`graph_json`, `canonical_xyz`, `ligand_pdbqt`, `receptor_pdbqt`, `dlg_path`,
`pose_xyz_paths`, `cleaned_pdb` and `n_heavy_atoms` from upstream `data`.

## What it gives

Measured on an Intel Mac (osx-64) on 2026-09-24 with metaldock-vwf 0.4.2: three
runs of this workflow, each node run by `bocoflow_core.node_runner` in the
package's own environment, as the app runs it. Ten poses per run, thirty in all.

| Poses | ΔG (kcal/mol) | RMSD from the crystal pose | Residues within 4 Å |
|---|---|---|---|
| 21 of 30 | −4.64 to −4.65 | 5.44–5.48 Å | 12 |
| 9 of 30 | −4.57 to −4.58 | 4.56–4.59 Å | 9 |

Each run took 169–172 s, 159–162 s of it in AutoDock.

No pose is near the crystal pose, and the best-scoring ones are not the closest.
The twelve residues of the better-scoring poses include His83, the residue the
rhenium binds in the crystal, so the complex lands beside its site, not in it.
The box is centred on the crystal position of the metal, so landing nearby is
partly the box's doing.

Poses are numbered in the order AutoDock ran them, not by score: `_1` is the
first run's result, and the energies in `analysis/1jzi_re_analysis.json` follow
the same order. `best_pose_index` in that file names the best-scoring pose,
counting from 0; in the third run it was 1, pose `_2`, while pose `_1` scored
−4.58 kcal/mol.

**RMSD values before 0.4.2 are not comparable.** Results Analysis used to pair
the atoms of a pose with those of the reference by their line in the two files,
and the two files list them in different orders. It reported 5.50–5.91 Å for the
same thirty poses, and scored the crystal pose itself at 3.30 Å from its own
coordinates. It now pairs atoms by chemistry. See `CHANGELOG.md`, 0.4.2.

## Running it with ORCA instead

The default engine is semi-empirical. For DFT charges, download ORCA (free for
academic use) from [orcaforum.kofo.mpg.de](https://orcaforum.kofo.mpg.de), then
on **QM Charges** set:

- `engine` → `orca`
- `orca_path` → the extracted directory (e.g. `external/orca_6_1_1_macosx_intel_openmpi411`)
- `orcasimpleinput` → the method line. The workflow carries `B3LYP def2-SVP`, the
  line our reference charges were computed with (below). MetalDock's own ORCA
  example for this case uses `B3LYP D3BJ def2-TZVP`.

On macOS, clear the download quarantine first: `xattr -dr com.apple.quarantine external/orca_*/`.

Our runs with ORCA 6.1.1 charges, three poses each, re-scored with the corrected
RMSD on 2026-09-24:

| ORCA method line | Runs | ΔG (kcal/mol) | RMSD from the crystal pose |
|---|---|---|---|
| `B3LYP def2-SVP` | 4 (2026-03-18, 2026-06-01, 2026-06-01, 2026-06-02) | −5.55 | 5.48–5.49 Å |
| `B3LYP D3BJ def2-TZVP` | 1 (2026-03-18) | −5.53 to −5.54 | 5.49 Å |

These are the "−5.54 kcal/mol, RMSD 5.91 Å" that earlier versions of this page
gave as the paper's reference values. They were our def2-TZVP run of 2026-03-18,
not the paper's, and 5.91 Å was the mis-paired RMSD. The MetalDock paper computes
its charges with ADF (AMS 2021, TZP/B3LYP/COSMO with D3-BJ and ZORA), not ORCA.

## Expected outputs

- **QM Charges** → `qm/enriched_graph.json` with CM5 charges. With xtb the Re
  charge is **+0.749**; our ORCA calculation (ORCA 6.1.1, B3LYP/def2-SVP) gives
  **+0.704**, and its charges for all 29 atoms ship as
  `mdock_qm_charges/demo_data/1jzi_re_orca_reference_graph.json`. The mean
  absolute deviation between the two is 0.065 e.
- **Ligand PDBQT** → `pdbqt/1jzi_re_ligand.pdbqt` with `ROOT`/`ENDROOT`, an `Re`
  atom, and a `DD` dummy atom at the vacant site.
- **AutoDock Run** → `docking/*.dlg` plus the affinity maps (including
  `clean_1jzi.Re.map`, the metal-specific one) and per-pose `.xyz`/`.pdbqt`.
- **Results Analysis** → `analysis/1jzi_re_analysis.json` with binding energies,
  ligand efficiencies, interacting residues and RMSD.

## Notes

- **Single point, not an optimisation:** these settings compute the charges at
  the crystal geometry (`geom_opt=false`), as MetalDock's own example for this
  case does. `geom_opt=true` optimises the complex first, which is much slower,
  and with ORCA considerably so.
- **Targeted box:** `box_center` is given explicitly. Left empty, the docking box
  centres on the metal atom's coordinates instead — nearly the same thing here,
  but not reproducible across inputs.
- **Genetic algorithm:** AutoDock Run uses its own fixed settings (population 150,
  2,500,000 evaluations). MetalDock's example for this case overrides them
  (population 200, 30,000 generations, elitism 2, mutation rate 0.2).
- **Scores are rankings.** AutoDock4 binding energies are approximate. Treat them
  as hypotheses to test, not as measured affinities.
