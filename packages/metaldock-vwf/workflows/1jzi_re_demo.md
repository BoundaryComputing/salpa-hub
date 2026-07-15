# Demo workflow — 1JZI Re complex docking

Reproduces the **1JZI Re-complex** case from the MetalDock paper as a BoCoFlow
workflow built from the six `metaldock-vwf` nodes.

- **Protein:** 1JZI (*Pseudomonas aeruginosa* azurin)
- **Ligand:** Re(phen)(CO)₃(His83) — 29 atoms, metal = **Re**
- **Reference result:** ΔG ≈ −5.54 kcal/mol, ~12 interacting residues, RMSD ≈ 5.9 Å

## Inputs (from the cloned MetalDock examples)

```
collect/MetalDock/examples/example_runs/vacancy_coordination_sphere/ORCA/
  ├── 1jzi.pdb          → mdock_protein_prep  (Protein PDB)
  └── 1jzi_D_REP.xyz    → mdock_ligand_prep   (Ligand XYZ)
```

ORCA is downloaded separately (free for academic use). Point `mdock_qm_charges`
at it via the **ORCA Directory** parameter (e.g.
`external/orca_6_1_1_macosx_intel_openmpi411`).

## Graph topology

The pipeline is a single linear chain — each node forwards its predecessor's
`data` keys and appends its own, so the receptor PDBQT made in step 1 is still
available to the docking node in step 5 without a separate edge:

```
[mdock_protein_prep] → [mdock_ligand_prep] → [mdock_qm_charges]
        → [mdock_ligand_pdbqt] → [mdock_autodock_run] → [mdock_results_analysis]
```

## Node configuration

| # | Node | Key parameters |
|---|------|----------------|
| 1 | **Protein Prep** | `case_name=1jzi_re`, `pdb_file=1jzi.pdb`, `output_dir=<run>`, `ph=7.0`, `clean=true` (MGLTools auto-detected) |
| 2 | **Ligand Prep** | `xyz_file=1jzi_D_REP.xyz`, `metal_symbol=Re`, `output_dir=<run>` (set **Force OpenBabel Python API** if mgltools shadows the obabel binary) |
| 3 | **QM Charges** | `engine=orca`, `geom_opt=false`, `charge=1`, `spin=0`, `ncpu=4`, `orca_path=<orca dir>`, `orcasimpleinput=B3LYP def2-SVP` |
| 4 | **Ligand PDBQT** | `vacant_site=true`, `max_torsions=32`, `freeze_coordination_sphere=true` (metal symbol + graph inherited from predecessors) |
| 5 | **AutoDock Run** | `num_poses=10`, `box_center=1.65,-7.803,27.176`, `box_size=20,20,20` (ligand/receptor PDBQT + graph inherited) |
| 6 | **Results Analysis** | `cutoff=4.0`, optionally `reference_xyz=1jzi_D_REP.xyz` for RMSD (protein PDB + poses + heavy-atom count inherited) |

Leave the file-input fields on nodes 3–6 **empty** — they auto-discover
`graph_json`, `canonical_xyz`, `ligand_pdbqt`, `receptor_pdbqt`, `dlg_path`,
`pose_xyz_paths`, `cleaned_pdb`, and `n_heavy_atoms` from upstream `data`.

## Expected outputs

- **mdock_qm_charges** → `enriched_graph.json` with CM5 charges (Re ≈ +0.7 e).
- **mdock_ligand_pdbqt** → `1jzi_re_ligand.pdbqt` with `ROOT`/`ENDROOT`, a `Re`
  atom, and a `DD` dummy atom at the vacant site.
- **mdock_autodock_run** → `*.dlg` + per-pose `.xyz`/`.pdbqt`.
- **mdock_results_analysis** → `1jzi_re_analysis.json` with binding energies,
  ligand efficiencies, interacting residues (and RMSD if a reference was given).

## Notes

- **Single-point vs geometry optimization:** the validated runs use single-point
  (`geom_opt=false`) for speed; set `geom_opt=true` to match the paper exactly
  (much slower).
- **Targeted box:** `box_center` is given explicitly here. With it empty, the
  docking box centers on the metal atom's coordinates instead.
- Once wired in the GUI, export the workflow to a `.bcflow` file and drop it next
  to this doc so it can be re-imported directly.
