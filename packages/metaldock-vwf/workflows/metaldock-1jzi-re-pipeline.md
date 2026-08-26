# Demo workflow — 1JZI Re complex docking

Reproduces the **1JZI Re-complex** case from the MetalDock paper as a Salpa
workflow built from the six `metaldock-vwf` nodes.

- **Protein:** 1JZI (*Pseudomonas aeruginosa* azurin)
- **Ligand:** Re(phen)(CO)₃(His83) — 29 atoms, metal = **Re**, one vacant
  coordination site
- **Published reference (ORCA/DFT charges):** ΔG ≈ −5.54 kcal/mol, ~12
  interacting residues, RMSD ≈ 5.9 Å

## Just run it

`workflows/metaldock-1jzi-re-pipeline.json` is the installable template — load
**MetalDock 1JZI Re Pipeline** from the template library, set a working
directory, and execute. Nothing to download: the inputs ship inside the nodes'
`demo_data/` directories and the charges come from GFN1-xTB, which installs with
the package.

Everything below describes what that template contains, and how to run the same
case with DFT charges instead.

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

## Running it with ORCA instead

The default engine is semi-empirical. To reproduce the paper's numbers, download
ORCA (free for academic use) from [orcaforum.kofo.mpg.de](https://orcaforum.kofo.mpg.de),
then on **QM Charges** set:

- `engine` → `orca`
- `orca_path` → the extracted directory (e.g. `external/orca_6_1_1_macosx_intel_openmpi411`)
- `orcasimpleinput` → `B3LYP def2-SVP` (what the validated runs used)

On macOS, clear the download quarantine first: `xattr -dr com.apple.quarantine external/orca_*/`.

## Expected outputs

- **QM Charges** → `qm/enriched_graph.json` with CM5 charges. With xtb the Re
  charge is ≈ **+0.749**; with ORCA it is ≈ **+0.704**. The full ORCA reference
  for all 29 atoms is in `mdock_qm_charges/demo_data/1jzi_re_orca_reference_graph.json`
  — mean absolute deviation between the two is 0.065 e.
- **Ligand PDBQT** → `pdbqt/1jzi_re_ligand.pdbqt` with `ROOT`/`ENDROOT`, an `Re`
  atom, and a `DD` dummy atom at the vacant site.
- **AutoDock Run** → `docking/*.dlg` plus the affinity maps (including
  `clean_1jzi.Re.map`, the metal-specific one) and per-pose `.xyz`/`.pdbqt`.
- **Results Analysis** → `analysis/1jzi_re_analysis.json` with binding energies,
  ligand efficiencies, interacting residues and RMSD.

## Notes

- **Single-point vs geometry optimization:** these settings use a single point
  (`geom_opt=false`) for speed. Set `geom_opt=true` to match the paper exactly
  (much slower, and with ORCA considerably so).
- **Targeted box:** `box_center` is given explicitly. Left empty, the docking box
  centres on the metal atom's coordinates instead — nearly the same thing here,
  but not reproducible across inputs.
- **Scores are rankings.** AutoDock4 binding energies are approximate. Treat them
  as hypotheses to test, not as measured affinities.
