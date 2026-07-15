# Changelog

All notable changes to the `metaldock-vwf` package are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/); versioning is
[SemVer](https://semver.org/).

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
