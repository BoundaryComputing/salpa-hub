# Changelog

All notable changes to the `metaldock-vwf` package are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/); versioning is
[SemVer](https://semver.org/).

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
