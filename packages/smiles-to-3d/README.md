# SMILES to 3D (Salpa node)

3D conformers from SMILES, written as SDF. Wraps [RDKit](https://www.rdkit.org/) (BSD-3-Clause):
ETKDGv3 embedding, then MMFF94 optimisation, one SDF record per conformer.

## What it does

A SMILES string is a 2D description; docking, MD and most property tools want 3D coordinates.
This node adds hydrogens, embeds one or more 3D conformers with a fixed random seed (so a run
repeats exactly), optimises each with the MMFF94 force field, and writes them all to
`conformers.sdf` with the MMFF94 energy of each as a property. It sits at the start of a
ligand-preparation pipeline — before protonation, charge assignment, or docking — and is the
usual first step when a ligand exists only as a name or a drawing.

## Inputs

- **SMILES (one per line)** — paste molecules directly, `SMILES name` per line, **or**
- **SMILES file** — a `.smi`/`.txt` file in the same convention (`#` comments allowed), **or**
- SMILES carried from a predecessor node under `data["smiles"]`.

## Parameters

| Option | Default | Meaning |
|--------|---------|---------|
| Conformers per molecule | 1 | How many conformers to embed and optimise; all are written |
| Random seed | 42 | ETKDG seed; the same seed gives the same coordinates |
| Output Directory | | Where `conformers.sdf` goes |

## Output

- `conformers.sdf` — one record per conformer, titled with the molecule's name, hydrogens
  kept, with the properties `SMILES`, `conformer_id` and `MMFF94_energy_kcal_mol`
- `data["molecules"]` — per molecule: name, SMILES, atom count (with H), heavy-atom count,
  conformer count, lowest MMFF94 energy in kcal/mol
- `data["smiles"]` — the SMILES list, forwarded so a downstream node can use it
- `data["output_file"]` — the SDF path

An unparseable SMILES fails the node with an `input` error naming it. Nothing is skipped
silently: a docking run on three of four requested ligands is a wrong answer, not a partial one.

## Deployment

**Local (PIXI_SUBPROCESS)** — CPU only. One dependency, `rdkit` from conda-forge, which builds
for every platform this package declares: linux-64, linux-aarch64, osx-64, osx-arm64, win-64.

## Status

Built with the `salpa-node` skill and verified on 2026-09-05:

- `salpa validate` — 0 errors, 0 warnings, import checks **ran**
- `salpa smoke` — runs on `demo_data/ligands.smi`, refuses a missing input, deterministic,
  idempotent, path-independent
- `salpa platforms --solve` — the declared list matches the evidence (linux-aarch64 was added
  by it)
- 9 tests, including atom counts hand-worked from the formulas (ethanol 9, aspirin 21,
  caffeine 24 atoms with hydrogens) and a same-seed determinism check
- On the demo data: ethanol MMFF94 −1.34 kcal/mol, aspirin +18.91, caffeine −122.53; the
  ethanol C–O bond comes out at 1.420 Å against a reference of about 1.43 Å

Not verified: running on a Salpa canvas with an upstream node feeding `data["smiles"]`.
