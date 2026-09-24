# Changelog

All notable changes to the `metaldock-vwf` package are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/); versioning is
[SemVer](https://semver.org/).

## [0.4.2] — 2026-09-24

### Fixed

- **RMSD pairs each atom of a pose with the same atom of the reference.** Results
  Analysis paired them by their line in the two files. AutoDock Run writes poses in
  the ligand PDBQT's torsion-tree order (Re, C, O, … for 1JZI), and a reference XYZ
  keeps its own order (O, C, Re, …), so every pair was two different atoms and every
  RMSD reported against a reference was wrong. Written out as a pose, the 1JZI
  crystal pose scored 3.30 Å against its own coordinates. A reference with another
  atom count was cut to fit, with a warning in the log.

  Both structures now become bond graphs, and the RMSD is the lowest over every
  pairing that keeps elements and bonds: symmetry-corrected, taken where the pose
  lies without superposition, over heavy atoms. A reference that is not the same
  molecule is refused rather than scored. The thirty poses of three 1JZI runs go
  from 5.50–5.91 Å to 4.56–5.48 Å; the best-scoring ones sit at 5.44–5.48 Å.

  Upstream MetalDock pairs each atom with the nearest atom of its element instead,
  which allows pairings no bond pattern does and reads a misplaced pose as closer
  than it is. On a pose far from the crystal this node therefore reports more than
  MetalDock would.

  On 2BZH, the one case from the MetalDock paper we have docked with these modules
  (its inputs are in MetalDock's repository, not in this package), our run of
  2026-03-18 re-scores from 3.34 Å to 0.45 Å. The paper reports 0.475 Å. That case
  was reproduced all along; the old RMSD hid it.

- **Results Analysis reports the best-scoring pose, and says which it is.** Its
  message paired the lowest binding energy with the contact count of pose 1, but
  AutoDock numbers poses by GA run, not by score. In one of three 1JZI runs the
  message read `Best ΔG -4.65 kcal/mol, 9 contacting residues (pose 1)`: pose 1
  had scored −4.58 with 9 contacts, the best pose −4.65 with 12. It now reads
  `Best ΔG -4.65 kcal/mol (pose 2), 12 contacting residues`, and `best_pose_index`,
  counting from 0 (the first of equals), is in the node's output and in
  `<case>_analysis.json`. The choice is made among the poses analysed, since the
  energy pattern also matches the log's cluster summary.

- **AutoDock Run stops when autogrid4 or autodock4 fails.** Both report failure only
  through their exit code, and neither was checked. A failed autogrid4 let
  autodock4 run on missing maps, and the node reported "Docked 0 pose(s)" as a
  success; the failure surfaced one node later as "No pose_xyz_paths found in
  predecessor data", with nothing to say why. It now stops at the failing step and
  quotes the tool's own last lines.

- **The protein-prep pH now changes the protonation.** `mdock_protein_prep` called
  pdb2pqr with `--with-ph` but no `--titration-state-method`, and pdb2pqr applies
  the pH only to a titration method's results, so PROPKA never ran: the 1JZI
  receptor was identical at pH 4 and pH 10. It now passes
  `--titration-state-method propka`. Upstream MetalDock passes the same flags as
  before, so this is a deliberate departure from it. At the 1JZI workflow's pH 7.0
  the receptor is unchanged. At the HSA template's pH 7.4 five residues change:
  Lys106 and Lys199 lose a proton, and Glu244, His247 and His288 gain one. The HSA
  result does not move: −3.00 kcal/mol and the same twelve contacts.

- **Running protein prep again rebuilds its outputs.** A protonated PDB or receptor
  PDBQT already in the output folder was returned untouched, so a changed pH had no
  effect on a second run. A failed pdb2pqr is now an error that quotes pdb2pqr's
  reason, instead of a warning followed by a missing file.

- **The HSA template no longer puts a vacant-site dummy on ferrocene.** It set
  `vacant_site=true`, so every run carried a `DD` dummy atom 1.0 Å from the iron,
  scored as a hydrogen-bond donor, although ferrocene is saturated and the
  walkthrough said none was added. It is now off. The best pose moved by 0.02 Å and
  its score from −3.01 to −3.00 kcal/mol.

### Documentation

- **The 1JZI "paper reference" was our own run.** The README, tutorial and 1JZI notes
  gave ΔG −5.54 kcal/mol and RMSD 5.91 Å as the MetalDock paper's ORCA/DFT values
  for 1JZI. The paper does not report 1JZI, and it computes its charges with ADF,
  not ORCA. The numbers were our validation run of 2026-03-18 (ORCA 6.1.1,
  B3LYP-D3BJ/def2-TZVP), and 5.91 Å was the mis-paired RMSD; corrected, it is
  5.49 Å. The notes also gave `B3LYP def2-SVP` as the line "the validated runs
  used". That line made the shipped reference charges and our other ORCA runs; the
  run quoted used upstream's `B3LYP D3BJ def2-TZVP`. Every number in the README, the
  tutorial and both workflow docs now says where it comes from, and the 1JZI
  workflow's own description no longer calls it the paper's case.
- **The tutorial walks through the HSA template.** It sent readers to the 1JZI
  template, which the Hub has not offered since 2026-08-26.
- **Runtimes are measured.** The HSA template takes 2–3 minutes (148–155 s in three
  runs on an Intel Mac), not "about 11 minutes".
- **Poses are numbered in the order AutoDock ran them, not by score.** The
  walkthrough said `_1` … `_10` were ranked best first.
- **NOTICE** lists the HSA template's inputs, gives the ORCA level of the reference
  charges, and names Meeko's licence correctly: LGPL-2.1, not Apache-2.0, as its
  LICENSE, conda-forge and PyPI all state.

## [0.4.1] — 2026-09-14

### Fixed

- **Every parameter in `workflows/metaldock-1jzi-re-pipeline.json` now holds one value.** A
  template records each parameter in four places. Loading reads one of them, and that one was
  right, but `option_types`, the copy stored with each parameter's type, still held other values:
  the nodes' defaults for `case_name` (`complex`, against `1jzi_re` on all six nodes) and for the
  ORCA line (`PBE def2-TZVP CPCM(Water)`, against `B3LYP def2-SVP`), and an empty `metal_symbol`
  and `box_center` (against `Re` and `1.65,-7.803,27.176`). Nothing that runs changed, and the HSA
  template already agreed with itself. The file is now what the app writes when the template is
  loaded and exported again.

## [0.4.0] — 2026-08-25

### Added

- **A second worked example: ferrocene into human serum albumin.** The 1JZI case
  is a *redocking* — the complex comes out of the same crystal it goes back into,
  so the binding site is known before you start and success is measured as RMSD
  against the answer. That is a validation exercise, and it quietly teaches a
  workflow nobody can reproduce on a new target.

  This one is a prediction. The receptor is apo (PDB 1AO6, chain A — HSA is
  monomeric in solution, so the crystallographic second copy is dropped), the
  ligand is an independent GFN1-xTB geometry rather than something extracted from
  the protein, and there is no reference pose, so `reference_xyz` is deliberately
  empty: RMSD-to-truth does not exist for a prediction and reporting one would be
  theatre.

  Suggested by Sylvestre Bonnet, who noted the demo should dock into a protein
  that does not already contain a metal complex, and named albumin as the system
  to do it with.

- **Fe is exercised for the first time.** Fe sits in `INTERNAL_PARAM_METALS`, so
  it takes a different path from Re: `get_lj_params` returns `None`, no
  `nbp_r_eps` lines are written, and docking relies on `atom_par Fe` from the
  parameter library. That branch had never been run. It works — xtb converges CM5
  charges on Fe(II), the graph builder resolves eta-5 sandwich bonding (30 bonds
  for 21 atoms = 10 C-H + 10 C-C + 10 Fe-C), and AutoDock docks it.

### The part worth reading

Choosing the site is the step the 1JZI demo hides, and getting it wrong does not
raise an error — it returns a worse answer somewhere else.

The first attempt centred the box on the **centroid of Trp214**, the residue that
lines Sudlow site I. That is 9.7 A from the pocket: with a 20 A box it put the
real site on the boundary. Result: -2.63 kcal/mol, poses spread over 0.45
kcal/mol, and contacts in residues 328-354 — subdomain IIB, the wrong place.

The shipped template takes the centre from **R-warfarin co-crystallised in 2BXD**,
transferred into the 1AO6 frame by Kabsch superposition on all 578 CA atoms
(RMSD 0.88 A). Result: -3.00 kcal/mol, all poses within 0.01 kcal/mol, and 11 of
12 contacting residues in subdomain IIA — TYR150, ARG222, LEU238, ARG257, LEU260,
ALA261, ILE264, LYS286, SER287, HIS288, ILE290, ALA291.

A residue centroid is not a pocket centroid. A co-crystallised ligand is.

## [0.3.2] — 2026-08-23

### Fixed

- **Docking still failed when the *working* directory contained a space.** 0.3.1
  fixed the AutoDock parameter file, which was one site; this is a second and
  independent one. Running the 1JZI example from `~/Salpa Runs/metaldock-1jzi`
  died in the first node with:

      AssertionError: /Users/…/Salpa does't exist

  The cause is MGLTools' own `pythonsh` wrapper, whose final line is

      exec $python $pyflags $@

  with `$@` **unquoted**. The shell re-splits every argument on whitespace after
  Python's `subprocess` has already passed a correct argv, so no amount of care
  on the calling side prevents it — and passing a list rather than a string does
  not help, which is why this survived review.

  Both affected calls (`prepare_receptor4` and `prepare_gpf4`) now invoke the
  interpreter that `pythonsh` wraps, supplying the `PYTHONHOME` and `PYTHONPATH`
  it would have set. No shell is involved at any point, so this holds for a space
  anywhere in the path — the working directory, the script location, or the
  user's home — rather than only where we thought to look.

  Verified against the failing input: the receptor PDBQT is written with 2384
  atoms, while the same call through `pythonsh` still fails identically.

## [0.3.1] — 2026-08-22

### Fixed

- **Docking failed in the packaged app because the install path contains a
  space.** Salpa installs nodes under `~/Library/Application Support/`, and the
  node wrote that absolute path into the GPF and DPF as `parameter_file`.
  AutoGrid and AutoDock read those files as whitespace-delimited text and the
  format has no quoting syntax, so the path was truncated at the space:
  `autodock4: FATAL ERROR: Sorry, I can't find or open /Users/…/Application`.
  The parameter library is now staged into the run directory and referenced by
  bare filename, which is how every other path in these files already works and
  is robust even when the user's own working directory contains a space.
  Development runs were unaffected because the development home contains no
  space, which is why every earlier run of this package passed.
- **A failed `prepare_gpf4` no longer passes silently.** A non-zero exit was
  logged as a warning and execution continued on a truncated grid parameter
  file, so the run failed ten minutes later inside AutoDock instead of
  immediately at the step that actually broke.
- **`ligand_types` is located by keyword rather than by line number.** The
  parser took line index 5 unconditionally; when the grid parameter file was
  malformed it harvested numbers from the `nbp_r_eps` block as atom types and
  requested maps such as `clean_1jzi.0.2966.map`. A well-formed file one line
  out of position produced the same class of error. A missing `ligand_types`
  line is now an error rather than silent corruption.

## [0.3.0] — 2026-08-20

The package now runs. Everything below the first heading is a prerequisite for
that sentence being true.

### Fixed

- **MGLTools no longer stops every node from starting.** The `mgltools 1.5.7`
  conda package ships its own `bin/python` and clobbers the one conda's `python`
  package installs, so `$PREFIX/bin/python` in the shared environment was
  **Python 2.7**. BoCoFlow launches every node as
  `pixi run python -m bocoflow_core.node_runner`, and installs `bocoflow_core`
  into the environment with the same `python` — so the install silently failed
  and no node could launch. MGLTools now gets its own pixi environment
  (`[feature.mgltools]` + `no-default-feature`), leaving `default` a clean
  Python 3 while AutoDockTools keeps the Python 2.7 its scripts require. The
  nodes locate it as a sibling of their own prefix.

  This was previously patched after the fact by `operations/fix_metaldock_env.sh`,
  which edited the *built* environment and therefore reverted on every rebuild
  and could never ship to a user. The fix is now in the manifest.

- **`pandas` was missing from the environment.** `metaldock_modules/qm_charges.py`
  imports it at module level for the Hirshfeld→CM5 conversion. It is declared in
  the project's own `pyproject.toml` and was dropped when the shared environment
  was written by hand, so every QM run in the app would have died on ImportError.
  Neither structural validation nor an import check catches this; only running
  the node does.

- **OpenBabel stopped printing a dlopen error on every call.** MGLTools ships a
  libcairo built in 2012 against a libpng12 at an absolute path on the machine
  that built it. Removing MGLTools from the nodes' environment removed the noise.

### Added

- **`xtb` charge engine, and it is now the default.** GFN1-xTB is the only
  supported engine that needs no user-supplied binary: it installs from
  conda-forge with the package and returns CM5 charges and Wiberg bond orders —
  precisely the two quantities the pipeline consumes — in under a second. That
  is what makes the bundled workflow runnable unattended.

  It is semi-empirical, not DFT. On the 1JZI Re case it reproduces the ORCA
  metal charge to 0.046 e (+0.749 against +0.704), mean absolute deviation
  0.065 e over 29 atoms. ORCA remains the accuracy reference; `orca`,
  `gaussian` and `adf` are unchanged.

  **Only GFN1 prints CM5 charges.** GFN2 prints Mulliken charges, which are not
  interchangeable, and the engine refuses any other parametrisation rather than
  returning the wrong quantity under the right name.

- **An installable workflow template**, `workflows/metaldock-1jzi-pipeline.json`,
  registered in `package.toml`. Produced by the app's own
  `/api/workflow/{id}/export-template` from a workflow that was imported,
  configured and **executed end to end** — not hand-written JSON, whose edges do
  not paint (see 0.2.1). Every input is a `node:` path into a node's own
  `demo_data/` and every output a `rel:` path under the working directory, so it
  carries no trace of the machine that produced it.

- **Bundled demo data.** `1jzi.pdb` and `1jzi_D_REP.xyz` now ship inside the node
  directories that consume them, so the template needs no downloads. They
  previously lived in a gitignored clone of the MetalDock repository, which is
  why no example could be shipped at all. `mdock_qm_charges/demo_data/` also
  carries the ORCA reference graph to compare xtb against. Provenance is
  recorded in `NOTICE`.

- **`DEMO_CONFIG` in all six nodes** — the values that make each node run against
  its own `demo_data/`, declared once and read by both `salpa smoke` and the
  shipped template.

- Tests for the xtb readers, including that the CM5 column is read rather than
  the Mulliken column beside it, and that `wbo`'s 1-based indices are shifted to
  the graph's 0-based ones.

- **`TUTORIAL.md`** — the plain-language lesson the Hub quality gate asks for,
  written for a wet-lab audience with no simulation background. Kept strictly
  separate from the README, which stays formal reference documentation: the two
  registers serve different readers and mixing them serves neither. Covers what
  docking is, why metals break the default scoring, the fast-vs-careful engine
  choice with the measured numbers, how to read the outputs, what "wrong" looks
  like, and the honest limits — rigid protein, no solvent, and scores that are
  rankings rather than affinities.

### Changed

- `mdock_qm_charges` gains `xtb_path`, `xtb_solvent` and `xtb_accuracy`; the
  engine list is now `xtb | orca | gaussian | adf` and defaults to `xtb`
  (previously `orca`, which no user could run without downloading it first).
- README rewritten as formal reference documentation.
- All six node versions → 0.2.0.

### Removed

- `workflows/1jzi-re-docking/workflow.bcflow`, superseded by the registered
  template. It was a plain export rather than a template — never listed in
  `[package.workflows]`, so the app never offered it — and it carried `abs:`
  input paths pointing at one machine, which made it a file that imports and
  then cannot run. Its slide moved to `workflows/1jzi-re-pipeline-slide.html`.
  `operations/export-showcase.spec.ts`, which regenerated it, is obsolete.

### Known limitations

- Still linux-64 and osx-64 only, and still academic / non-commercial, both
  because of MGLTools. The Meeko migration that lifts each of those is planned
  and not done — see `dev-notes/mgltools-to-meeko-migration.md`.

## [0.2.1] — 2026-06-23

### Fixed
- **Showcase workflow edges now render.** `workflows/1jzi-re-docking/workflow.bcflow`
  is now a **real BoCoFlow UI export** (`exportedBy: BoCoFlow BF2`) instead of
  hand-generated JSON. The old `gen_bcflow.py` output imported fine but its edges
  never painted on the canvas: it wired links into each node's flow-control port
  (`flow-in`) and never created the data `in-0` port, so at render time the link's
  last point resolved to `(0,0)` and `VPLinkModel.getSVGPath` suppressed the line.
  A genuine app export wires `out-0 → in-0` with render-compatible geometry.
- Replaced the Python generator with `operations/export-showcase.spec.ts`, which
  assembles + configures the workflow through the live UI, saves, and exports it —
  the always-correct, render-verified generator. Removed `workflows/gen_bcflow.py`.

## [0.2.0] — 2026-06-02

### Changed
- **Slide-style node names** — dropped the `MetalDock:` prefix from all six
  `display_name`s (`"MetalDock: Protein Prep"` → `"Protein Prep"`, etc.). Folder
  names (`mdock_*`) and `class_name`s (`Mdock*`) are unchanged, so node identity
  and existing `.bcflow` imports are unaffected — the rename is purely cosmetic.
  Mirrors `metalparm-vwf` v1.25.0. Regenerated `registry.json` and
  `workflows/1jzi-re-docking/workflow.bcflow` (only the `name` fields changed);
  updated `gen_bcflow.py`, `workflows/1jzi_re_demo.md`, and
  `operations/metaldock-pipeline.spec.ts`.

## [0.1.0] — 2026-06-01

Initial release. Wraps the 6 refactored MetalDock pipeline modules
(`src/metaldock_modules/`) as chainable BoCoFlow nodes — the first
BoCoFlow integration of the MetalDock pipeline.

### Added
- `mdock_protein_prep` — clean PDB → pdb2pqr protonate → prepare_receptor4 PDBQT.
- `mdock_ligand_prep` — canonicalize XYZ (OpenBabel) → build molecular graph (saved as JSON).
- `mdock_qm_charges` — ORCA/Gaussian/ADF DFT → CM5 charges + Mayer bond orders → enriched graph.
- `mdock_ligand_pdbqt` — enriched graph → ROOT/BRANCH PDBQT with metal-aware torsion freezing.
- `mdock_autodock_run` — GPF/DPF generation → autogrid4 → autodock4 → pose extraction.
- `mdock_results_analysis` — binding energies, ligand efficiency, interacting residues, RMSD.
- Shared `metaldock_vwf` pixi environment (openbabel, pdb2pqr, autogrid4, autodock4, mgltools, ase).
- Linear data-flow convention: each node forwards predecessor `data` keys and appends its own.
- `metaldock_modules` runtime resolution (METALDOCK_SRC env → bundled `scripts/` → repo `src/` fallback).
- `workflows/1jzi_re_demo.md` — full 1JZI Re-complex docking walkthrough.
- `workflows/1jzi-re-docking/workflow.bcflow` — importable showcase workflow, plus
  `workflows/gen_bcflow.py` (deterministic generator) — roadmap Phase 2 (G4).
- `operations/metaldock-pipeline.spec.ts` — Playwright E2E driving the 6-node pipeline
  through the BoCoFlow GUI (build/wire/configure + ORCA-free prep execution), with
  `playwright.config.ts`, `run.sh`, and `README.md` — roadmap Phase 2 (G2).

### Changed / fixed (live-GUI verification, 2026-06-01)
- Node resolver tries `import metaldock_modules` first, then `METALDOCK_SRC` →
  bundled `scripts/` → vendored `_vendor/` → repo `src/`. The package now vendors
  `metaldock_modules` under `_vendor/` so the *installed* copy is self-contained.
- `pixi.toml` platforms narrowed to `["linux-64", "osx-64"]` (`autodock`/`mgltools`
  have no `osx-arm64` builds; Apple Silicon runs under Rosetta/osx-64).
- Added `registry.json` (shelf-source manifest) so the package installs into a
  BoCoFlow dev stack, and `operations/fix_metaldock_env.sh` to repair the
  mgltools/python3 conflict in the built env (see dev-note).
- **Verified live** (5 Playwright tests, all PASS in the GUI): build/wire/configure;
  node-by-node prep execution; full node-by-node chain (ORCA + AutoDock, ΔG −5.55);
  orchestrated prep workflow (Run Workflow → execute_async job); and the full
  orchestrated 6-node workflow via a single Run Workflow click (ORCA + AutoDock,
  ΔG −5.55, ≈7 min via `orchestrate_workflow_parallel`).
- Note: orchestrated runs resolve node `Output Directory` under the workflow working
  path — use a `rel:` prefix (bare-relative resolves against the worker CWD).

### Known limitations
- ORCA must be supplied externally (not a conda package); pass `orca_path` to `mdock_qm_charges`.
- `mdock_autodock_run` requires `box_center` for targeted docking (no blind-docking auto-center yet).
- AutoDock4's 32-torsion hard limit applies; `mdock_ligand_pdbqt` freezes metal-proximal bonds to comply.
- The `metaldock_vwf` env needs `operations/fix_metaldock_env.sh` run once after install
  (mgltools claims `bin/python`→py2.7). Re-run after any env rebuild.
- Full QM + docking chain (nodes 3–6) needs ORCA + AutoDock4; not yet run end-to-end in the GUI.
