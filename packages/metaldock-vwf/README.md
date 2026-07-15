# metaldock-vwf

Visual workflow nodes for **metal-protein docking** — the refactored
[MetalDock](https://github.com/MatthijsHak/MetalDock) pipeline wired as chainable
BoCoFlow nodes.

Each node is a thin wrapper around one function in
`src/metaldock_modules/` (this repo). The science is **not** rewritten here — the
modules are already deep-refactored (pure functions, explicit tool paths, no global
state). These nodes only handle BoCoFlow's I/O contract: read parameters, resolve
file paths, call the module, forward outputs to downstream nodes.

## Pipeline

```
mdock_protein_prep ─→ mdock_ligand_prep ─→ mdock_qm_charges ─→ mdock_ligand_pdbqt ─→ mdock_autodock_run ─→ mdock_results_analysis
   (PDB → PDBQT)        (XYZ → graph)       (graph + CM5)        (graph → PDBQT)        (autogrid+autodock)     (energies, RMSD)
```

The chain is **linear**: every node forwards all of its predecessor's `data` keys
and appends its own, so the receptor PDBQT produced at the start is still available
to `mdock_autodock_run` near the end without re-wiring. Molecular graphs are passed
between nodes as JSON files (`utils.save_graph` / `utils.load_graph`), since BoCoFlow
serializes node outputs to disk.

## Nodes

| Node | Wraps | Key outputs (forwarded as `data` keys) |
|------|-------|----------------------------------------|
| **mdock_protein_prep** | `protein_prep.prepare_protein` | `receptor_pdbqt`, `cleaned_pdb`, `protonated_pdb` |
| **mdock_ligand_prep** | `ligand_prep.prepare_ligand` | `canonical_xyz`, `graph_json`, `n_heavy_atoms`, `metal_symbol` |
| **mdock_qm_charges** | `qm_charges.run_qm_and_enrich_graph` | `graph_json` (enriched), `qm_energy`, `qm_run_type` |
| **mdock_ligand_pdbqt** | `ligand_pdbqt.create_ligand_pdbqt` | `ligand_pdbqt` |
| **mdock_autodock_run** | `autodock_run.run_autodock` | `dlg_path`, `pose_xyz_paths`, `pose_pdbqt_paths` |
| **mdock_results_analysis** | `results_analysis.analyze_docking_results` | `binding_energies`, `interacting_residues`, `rmsd_values` |

## Locating the `metaldock_modules` package

Each node resolves the `metaldock_modules` import at runtime, in this order:

1. `METALDOCK_SRC` environment variable (a directory containing `metaldock_modules/`)
2. the node-bundled `scripts/` directory (populated at package-build/release time)
3. the source tree fallback: `<package>/../../src/` (works for in-repo dev)

In development (running from this repo) no setup is needed — fallback (3) finds
the repo's `src/metaldock_modules`.

## External tools

The shared `metaldock_vwf` pixi environment provides openbabel, pdb2pqr, autogrid4,
autodock4, mgltools, and ase. MGLTools' `prepare_receptor4.py` / `prepare_gpf4.py`
run under the bundled `pythonsh`; nodes auto-detect the MGLTools directory from the
environment but expose an override.

**ORCA** (the recommended QM engine) is not a conda package — download it separately
(free for academic use) and pass its directory to `mdock_qm_charges` via the
`orca_path` parameter.

## Example

See [`workflows/1jzi_re_demo.md`](workflows/1jzi_re_demo.md) for the full 1JZI
Re-complex docking workflow.
