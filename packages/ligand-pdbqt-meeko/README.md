# Ligand PDBQT (Meeko) (Salpa node)

Docking-ready PDBQT ligands from an SDF. Wraps [Meeko](https://github.com/forlilab/Meeko)
(LGPL-2.1-only), the ligand preparation tool of the AutoDock family, with RDKit underneath.

## What it does

AutoDock Vina and AutoDock4 read ligands as PDBQT: atoms typed for the force field, Gasteiger
partial charges, non-polar hydrogens merged into their heavy atoms, and the rotatable bonds
arranged as a tree of branches. This node takes an SDF of 3D ligands — what `smiles-to-3d`
writes — and prepares every record with Meeko's `MoleculePreparation`, writing one
`<name>.pdbqt` per record into the output directory and forwarding the list for a docking
node. It is the step between conformer generation and docking, and the way off MGLTools'
`prepare_ligand4.py`.

## Inputs

- **Ligands (SDF)** — a 3D SDF, hydrogens preferably explicit, **or**
- the SDF carried from a predecessor node under `data["sdf_file"]` or `data["output_file"]`
  (what `smiles-to-3d` forwards).

A record without 3D coordinates fails the node, naming the record: a 2D sketch is not a
ligand pose. A record without explicit hydrogens gets them added with coordinates, and the
result says so per ligand (`hydrogens_added`), because invented hydrogen positions are a
modelling decision the user should know about.

## Parameters

| Option | Default | Meaning |
|--------|---------|---------|
| Rigid macrocycles | off | Keep macrocycle rings rigid (AutoDock4 needs this; Vina handles Meeko's flexible glue-atom rings) |
| Hydrated docking | off | Add Meeko's explicit water sites |
| Flexible amides | off | Let amide C–N bonds rotate, as AutoDock does not by default |
| Output Directory | | Where the `.pdbqt` files go |

## Output

- `<name>.pdbqt` per record, named by the SDF title (made file-safe; duplicates numbered)
- `data["ligands"]` — per ligand: name, file, PDBQT atom count, input atom count, torsions
  (`TORSDOF`), total Gasteiger charge, whether hydrogens were added
- `data["pdbqt_files"]` — the paths, in order, for a docking node
- `data["output_dir"]`

A record Meeko cannot prepare fails the node with an `input` error naming it. Nothing is
skipped silently.

## Deployment

**Local (PIXI_SUBPROCESS)** — CPU only. `meeko` (pure Python) and `rdkit` from conda-forge;
both build for every platform this package declares.

## Status

Built with the `salpa-node` skill and verified on 2026-09-05:

- `salpa validate` — 0 errors, 0 warnings, import checks **ran**
- `salpa smoke` — runs on `demo_data/ligands.sdf`, refuses a missing input, deterministic,
  idempotent, path-independent
- 9 tests, with PDBQT atom and torsion counts hand-worked from AutoDock's rules (ethanol 4
  atoms / 1 torsion, aspirin 14 / 4, caffeine 14 / 0) and a same-result check for an input
  whose hydrogens had to be added
- Chained: with the Ligands parameter blank it takes the SDF `smiles-to-3d` forwards

Not verified: running on a Salpa canvas downstream of `smiles-to-3d`; docking the output.
