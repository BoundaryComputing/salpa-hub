# Changelog — pdbmdauto

All notable changes to this package. The format follows [Keep a Changelog](https://keepachangelog.com/);
versions are the `[package].version` in `package.toml`, which is what the Marketplace's Updates tab
compares against.

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
