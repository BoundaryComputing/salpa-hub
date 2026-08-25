# metaldock-vwf

Visual workflow nodes for **metal–protein docking** — the
[MetalDock](https://github.com/MatthijsHak/MetalDock) pipeline refactored into six
chainable Salpa nodes.

> **New to molecular docking?** Start with [`TUTORIAL.md`](TUTORIAL.md) — the same
> pipeline explained without jargon, for a bench audience, in about five minutes of
> compute. This README is the formal reference.

## Overview

Molecular docking predicts how a ligand binds a protein and ranks candidate poses
by an approximate energy function. Standard docking parameter sets do not describe
transition-metal centres: the charge distribution around the metal, and the
coordination geometry that follows from it, are outside what the empirical
force field was fitted to. MetalDock addresses this by deriving partial charges
for the metal complex from a quantum-mechanical calculation and feeding those into
an otherwise conventional AutoDock4 run.

This package packages that pipeline as workflow nodes. Each node wraps one
function in `metaldock_modules` (vendored under `_vendor/`); the scientific code
is a pure-function refactor of MetalDock, not a rewrite. The nodes handle only the
Salpa I/O contract: read parameters, resolve paths, call the module, forward
outputs downstream.

## Scope and intended use

Appropriate for: predicting binding modes of metal-organic complexes against a
prepared receptor, and ranking those poses relative to one another.

Not appropriate for: absolute binding free energies. AutoDock4 scores are
approximate rankings. Treat them as hypotheses to test experimentally, not as
measured affinities.

The supported metals are those with fitted Lennard-Jones parameters
(V, Cr, Co, Ni, Cu, Mo, Ru, Rh, Pd, Re, Os, Pt) plus Fe, Zn and Mn, which
AutoDock4 parameterises internally.

## Methodology

| Stage | What happens |
|---|---|
| Receptor preparation | HETATM records are optionally stripped, the structure is protonated at a chosen pH with **pdb2pqr** (PROPKA), and converted to PDBQT by AutoDockTools' `prepare_receptor4.py`. |
| Ligand preparation | The metal-complex geometry is canonicalized with **OpenBabel** and converted to a molecular graph — atoms with elements and coordinates, edges from covalent-radius adjacency. |
| Charge derivation | A quantum calculation on the complex yields **CM5 partial charges** and bond orders, which are written onto the graph. Bonds the calculation does not support are removed, so the torsion tree that follows is built on quantum connectivity rather than distance alone. |
| Ligand PDBQT | The enriched graph becomes an AutoDock PDBQT with a ROOT/BRANCH torsion tree. Bonds in the metal's coordination sphere can be frozen, and a `DD` dummy atom marks a vacant coordination site. |
| Docking | `autogrid4` computes affinity maps over a box; `autodock4` runs a Lamarckian genetic-algorithm search and writes a DLG. Poses are extracted to `.xyz`/`.pdbqt`. |
| Analysis | Binding energies, ligand efficiencies, residues within a cutoff of each pose, and RMSD against a reference geometry if one is supplied. |

### Charge engines

CM5 charges are the input this pipeline is built around, and four engines can
produce them.

| Engine | Method | Obtained how | Use for |
|---|---|---|---|
| **xtb** (default) | GFN1-xTB, semi-empirical | Installed with the package | Getting a pipeline working; screening |
| ORCA | DFT | Downloaded separately; free for academic use | Numbers you intend to publish |
| Gaussian | DFT | Commercial licence | Existing Gaussian workflows |
| ADF | DFT | Commercial licence | Existing ADF workflows |

xtb is the default because it is the only one of the four that requires no
user-supplied binary, which is what makes the bundled workflow runnable
unattended. It is semi-empirical and therefore approximate. On the 1JZI Re
reference case it reproduces the ORCA metal charge to 0.046 e (+0.749 against
+0.704) with a mean absolute deviation of 0.065 e across all 29 atoms, in under a
second rather than minutes. The largest deviations are on the carbonyl and
pyridyl heteroatoms, where xtb is systematically more polar.

Only **GFN1**-xTB reports CM5 charges; GFN2 reports Mulliken charges, which are
not interchangeable. The engine rejects any other parametrisation rather than
returning the wrong quantity under the right name.

## Pipeline nodes

The chain is linear. Every node forwards all of its predecessor's `data` keys and
appends its own, so the receptor PDBQT produced in step 1 is still available to
the docking node in step 5 without a separate edge. Molecular graphs travel
between nodes as JSON files, since node outputs are serialized to disk.

```
Protein Prep → Ligand Prep → QM Charges → Ligand PDBQT → AutoDock Run → Results Analysis
```

| Node | Wraps | Tool | Key outputs |
|------|-------|------|-------------|
| **mdock_protein_prep** | `protein_prep.prepare_protein` | pdb2pqr, MGLTools | `receptor_pdbqt`, `cleaned_pdb`, `protonated_pdb` |
| **mdock_ligand_prep** | `ligand_prep.prepare_ligand` | OpenBabel | `canonical_xyz`, `graph_json`, `n_heavy_atoms`, `metal_symbol` |
| **mdock_qm_charges** | `qm_charges.run_qm_and_enrich_graph` | xtb / ORCA / Gaussian / ADF | `graph_json` (enriched), `qm_energy`, `qm_run_type` |
| **mdock_ligand_pdbqt** | `ligand_pdbqt.create_ligand_pdbqt` | — | `ligand_pdbqt` |
| **mdock_autodock_run** | `autodock_run.run_autodock` | autogrid4, autodock4, MGLTools | `dlg_path`, `pose_xyz_paths`, `pose_pdbqt_paths` |
| **mdock_results_analysis** | `results_analysis.analyze_docking_results` | — | `binding_energies`, `interacting_residues`, `rmsd_values` |

File-input fields on nodes 3–6 may be left empty; each auto-discovers what it
needs (`graph_json`, `canonical_xyz`, `ligand_pdbqt`, `receptor_pdbqt`,
`dlg_path`, `pose_xyz_paths`, `cleaned_pdb`, `n_heavy_atoms`) from upstream
`data`.

## Example workflow

The package ships **two** installable templates. They answer different questions,
and the difference is the point.

### 1. MetalDock 1JZI Re Pipeline — a redocking

`workflows/metaldock-1jzi-re-pipeline.json`, reproducing the 1JZI case from the
MetalDock paper.

- **Receptor:** 1JZI, *Pseudomonas aeruginosa* azurin
- **Ligand:** Re(phen)(CO)₃(His83) — 29 atoms, vacant coordination site
- **Inputs:** bundled in each node's `demo_data/`; nothing to download
- **Charges:** GFN1-xTB single point, charge +1, closed shell

Load it from the template library and run it. Every input path is a `node:`
reference into a node's own `demo_data/`, and every output path is `rel:` under
the workflow's working directory, so it runs wherever it is opened.

`workflows/1jzi_re_demo.md` documents the same case in more detail, including how
to run it against ORCA instead.

### 2. HSA + Ferrocene (Sudlow site I) — a prediction

`workflows/metaldock-hsa-ferrocene.json`. The receptor has no metal complex in
it, so there is no crystallographic pose to recover and no RMSD to report.

- **Receptor:** 1AO6, human serum albumin, chain A (HSA is monomeric in solution)
- **Ligand:** ferrocene, Fe(C₅H₅)₂ — an independent GFN1-xTB geometry
- **Site:** Sudlow site I, subdomain IIA, centred on the position of R-warfarin
  co-crystallised in 2BXD and superposed into the 1AO6 frame
- **Charges:** GFN1-xTB, neutral, closed shell

This is the case that exercises **Fe**, which takes AutoDock's internal
parameters rather than the tabulated LJ terms Re uses.

It also demonstrates the step the redocking hides: **choosing the binding site.**
Centring the box on the centroid of Trp214 — the residue that lines Sudlow I, and
the obvious guess — lands 9.7 Å from the pocket and returns a weaker pose in the
wrong subdomain, without erroring. `workflows/hsa_ferrocene_demo.md` sets both
attempts side by side.

## Outputs and interpretation

- **QM Charges** → `enriched_graph.json`. On this case the Re charge should be
  near +0.7 e; `mdock_qm_charges/demo_data/1jzi_re_orca_reference_graph.json`
  holds the ORCA/DFT values for every atom to compare against.
- **Ligand PDBQT** → a PDBQT with `ROOT`/`ENDROOT`, the metal atom, and a `DD`
  dummy atom at the vacant coordination site.
- **AutoDock Run** → a `.dlg` plus per-pose `.xyz` and `.pdbqt`.
- **Results Analysis** → `<case>_analysis.json` with binding energies, ligand
  efficiencies, interacting residues, and RMSD when a reference was supplied.

Published reference values for this case with ORCA/DFT charges: ΔG ≈ −5.54
kcal/mol, roughly 12 interacting residues, RMSD ≈ 5.9 Å against the
crystallographic pose. Semi-empirical charges will not reproduce these exactly.

For the HSA case there is **no** reference pose, so `reference_xyz` is left empty
and `rmsd_values` comes back empty with it. That is correct: a prediction has no
answer to be scored against, and reporting an RMSD anyway would be theatre.
Judge it on binding energy and on which residues line the pose — here −3.00
kcal/mol with 11 of 12 contacts in subdomain IIA.

## Requirements and platforms

**linux-64 and osx-64 only.** MGLTools and AutoDock4 have no osx-arm64 or win-64
conda builds; Apple Silicon runs the environment under Rosetta, and Windows needs
WSL2 or Docker.

The package builds **two** pixi environments rather than one. The nodes run in
`default` (Python 3, with OpenBabel, pdb2pqr, AutoDock4/AutoGrid4, ASE, xtb,
pandas, NetworkX, SciPy). MGLTools gets its own environment because its conda
package replaces `bin/python` with Python 2.7 — in a shared prefix that stops
every node from launching, since nodes start as `python -m
bocoflow_core.node_runner`. The split lets MGLTools keep the Python 2.7 its
AutoDockTools scripts need without imposing it on anything else. Nodes locate it
as a sibling of their own prefix; see `_find_mgltools_dir()`.

ORCA is not a conda package. To use it, download it from
[the ORCA forum](https://orcaforum.kofo.mpg.de) (free for academic use) and pass
its directory to **mdock_qm_charges** via the *ORCA Directory* parameter, or put
`orca` on `PATH` / set `ASE_ORCA_COMMAND`.

### Licensing

Academic / non-commercial use only, because MGLTools is distributed by Scripps
under a non-commercial licence. A planned migration to Meeko (Apache-2.0) will
remove that constraint and restore osx-arm64 and win-64. See `NOTICE`.

### Locating `metaldock_modules`

Each node resolves the import at runtime, in order:

1. `import metaldock_modules` (already installed / on path)
2. the `METALDOCK_SRC` environment variable
3. the node-bundled `scripts/` directory
4. `<package>/_vendor/` — how an installed package finds it
5. `<repo>/src/` — the in-repo development fallback

## References

- Hakkennes, M. et al. *MetalDock: An Open-Source Docking Tool for Metal-Organic
  Compounds.* J. Chem. Inf. Model. 2023. doi:10.1021/acs.jcim.3c01582
- Morris, G. M. et al. *AutoDock4 and AutoDockTools4.* J. Comput. Chem. 2009.
  doi:10.1002/jcc.21256
- Bannwarth, C. et al. *Extended tight-binding quantum chemistry methods (xtb).*
  WIREs Comput. Mol. Sci. 2021. doi:10.1002/wcms.1493
- Marenich, A. V. et al. *Charge Model 5.* J. Chem. Theory Comput. 2012.
  doi:10.1021/ct200866d
- Dolinsky, T. J. et al. *PDB2PQR.* Nucleic Acids Res. 2004.
  doi:10.1093/nar/gkh381
