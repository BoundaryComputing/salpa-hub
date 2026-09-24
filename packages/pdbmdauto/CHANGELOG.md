# Changelog — pdbmdauto

All notable changes to this package. The format follows [Keep a Changelog](https://keepachangelog.com/);
versions are the `[package].version` in `package.toml`, which is what the Marketplace's Updates tab
compares against.

## [1.2.3] — 2026-09-24

The package's own tests. No node changed. `pixi.toml` changed, so an installed environment is
rebuilt once.

### Fixed
- **`pixi run test` runs every test in the package.** It ran only the 10 tests in the package's
  own `tests/`, in an environment without the node-authoring API, so none of the 214 tests in the
  nodes' own `tests/` ran. It now runs all of them, from the package root, in the `test`
  environment; `pytest.ini` lists them. The `dev` environment, which only that task used, is gone.
- **Three node test files load again.** pdb2pqr's and pdb_fasta_biopython's imported their nodes
  through a module layout the app no longer installs. gmx_solv_ion's worked only when run from the
  package root.
- **Thirteen tests that had drifted from their nodes are updated.** They now expect the options
  and output locations the nodes have: the inherited `force_to_run`, gmx_mdrun_local's
  `output_folder`, pdb_fasta_biopython's case folder, and its `chain_<id>.fasta` file names.

## [1.2.2] — 2026-09-24

Two node fixes and two corrections to what the package says. The environment is unchanged, so
nothing is rebuilt.

### Fixed
- **PDB to PQR Converter runs with Optimize Hydrogens off.** It passed `--no-optimize`, an option
  PDB2PQR does not have: PDB2PQR's is `--noopt`. PDB2PQR's parser is strict, so every run with the
  option off stopped at once with "unrecognized arguments". The node's tests now check the command
  against PDB2PQR's own parser, and run it once with the option off.
- **pKa + GROMACS EM fails when PDB2PQR fails.** It used to carry on with the unprotonated structure
  and report success, with "0 protonation change(s)". pdb2gmx then chose the protonation states,
  not PDB2PQR at the requested pH, and the only record went to a log file the app does not show.
  Now the step stops. The Log Center shows why, from the end of PDB2PQR's output, where the reason
  is, and so does the node's error. To let pdb2gmx choose the states deliberately, turn off
  *Run PDB2PQR*.
- **Fix Missing Residues is described as it runs.** The README, the walkthrough and the node's own
  description said the step minimises the model and runs `pm build-model`. It does neither: it
  drives ProMod3's Python API, and the model is first minimised in the next step. The node's error
  messages no longer name `build-model`. Its description no longer says DNA and RNA chains are kept:
  they are left out of `fixed.pdb`, as its warning says.
- **PDB2PQR's citation.** The PDB to PQR Converter gave the 2004 paper's title with the 2007 paper's
  DOI. It now gives the 2004 paper's own, 10.1093/nar/gkh381, as the README's References do.

## [1.2.1] — 2026-09-14

The bundled template only. No node code changed, and nothing that runs changed.

### Fixed
- **Every parameter in `workflows/pdbmdauto-pipeline.json` now holds one value.** A template
  records each parameter in four places. Loading reads one of them, and that one was right; two of
  the others held values the pipeline does not use. There, `case_name` was empty or null on all 11
  nodes, against `e2e_4z8j`; PDB FASTA Parser's `pdb_id` and `output_dir` were null, against `4Z8J`
  and `rel:pdbmdauto-e2e-full`; `check_pdb_header` was false, against true; and the second
  relaxation's `protocol` was `full_4step`, the node's default, against `em_only`. A reader of
  those copies saw a different pipeline from the one that runs. The file is now what the app
  writes when the template is loaded and exported again.

## [1.2.0] — 2026-09-03

Documentation and metadata. No node code changed.

### Added
- `TUTORIAL.md` — the plain-language walk for someone who has never run a simulation. Kept apart
  from the README on purpose: one document is for citing, the other for reading.
- `workflows/pdbmdauto-pipeline.md` — the walkthrough for the bundled template, with three figures
  from a real run (`workflows/figures/`). The app shows it as *Read the walkthrough* on the template
  picker; salpa.app renders the same file as a docs page.
- `NOTICE` — what is ours (MIT), what is bundled as demonstration data and where it came from, and
  which third-party tools the environment installs but this package does not redistribute.
- `CHANGELOG.md` (this file).
- `[package.documentation] tutorial = "TUTORIAL.md"`.

### Changed
- **Platform statement corrected.** Apple Silicon is native (GROMACS and ProMod3 ship arm64 builds;
  no Rosetta), and Windows is served through WSL2, which Salpa sets up itself. `platform_note`,
  the README and every node's `[node.platforms]` block said otherwise — "Docker/WSL2 support
  coming soon", "Apple Silicon under Rosetta" — long after both had stopped being true.
- **`linux-aarch64` declared.** Every dependency has an ARM Linux build and the solve is
  identical to `linux-64`; the platform list was hand-maintained and wrong in this direction
  (bocoflow#105).
- **`openstructure` declared.** `fix_residues_promod3` imports it directly; until now that import
  was satisfied only transitively through `promod3`.
- README: thirteen nodes (not fourteen — `pdb_tools_clean` left in 1.0.x); the deleted "PDB Clean"
  node no longer listed; 4Z8J has **six** unresolved residues, all N-terminal, not five; measured
  runtimes with the machine named; the licensing footer no longer says "being finalized".
- Template `template_info`: the description says which stages run, that the six residues are
  terminal, and that the structure is fetched from RCSB (network needed); `estimated_time` is a
  measured figure.
- `pixi.toml [project].version` now tracks the package version.
- Authors unified to `BoundaryComputing` on the two nodes that still said "BoCoFlow Community".
- Two stale `pixi.lock` files are no longer tracked.

### Environment rebuild
`pixi.toml` changed (platforms, `openstructure`, version). The app hashes the installed manifest
and rebuilds the shared environment when it differs, so **an existing installation rebuilds
`pdbmdauto` once on the next run** — with a warm package cache that is under a minute; a cold one
downloads the packages again. The solved package set is the same as before apart from build-number
bumps of unchanged versions.

### Known
- The *Model Terminal Extensions* option on Fix Missing Residues has no effect: the node always
  models the termini. Documented rather than changed, since this release touches no node code.
- The production run is 2 ps — a pipeline exercise, not sampling.

## [1.1.2] — 2026-09-01
- Pin `openmm >=8.3.1,<8.6`: openmm 8.6.0 removed a symbol promod3's compiled extension links
  against, so every environment solved after 2026-08-19 died at Fix Missing Residues (bocoflow#130).

## [1.1.1] — 2026-08-31
- Test fixture builds its awkward paths instead of hardcoding a home directory.

## [1.1.0] — 2026-08-31
- Every subprocess call is an argv list; no shell is involved anywhere in the package. The
  quoting layer added in 1.0.6–1.0.8 was removed rather than patched (`tests/test_shell_safety.py`
  guards the invariant).

## [1.0.6] – [1.0.8] — 2026-08-31
- Quote every path that reaches a shell (bocoflow#104: the pipeline died at Solvate & Ionize on a
  packaged macOS install, whose path contains a space); the guard was interpreter-dependent and
  hid six more sites in `gmx_md_relax`; its own explanation corrected.

## [1.0.2] – [1.0.5] — 2026-08-01 … 2026-08-03
- Repository points at the public Salpa Hub; authors read `BoundaryComputing`; the template
  carries its own author and a name that is not a save timestamp.

## [1.0.1] — 2026-07-21
- `pdb-tools` and `mdanalysis` restored to the shared environment (pdb2pqr needs them at run
  time; trimming them broke the node, bocoflow#64). `pdb_tools_clean` removed on 2026-07-28
  (bocoflow#73); the package has thirteen nodes since.

## [1.0.0] — 2026-03-31
- First release: fourteen nodes replacing `gromacs-suite` + `pdb-toolkit`, and the
  `pdbmdauto-pipeline` template.
