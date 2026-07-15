# Salpa Hub

Official, first-party scientific node packages for **Salpa** — built, tested, and
maintained by BoundaryComputing. Download them into the Salpa desktop app and run.
**100% free and open.**

Salpa Hub is a public, additive distribution channel: a static repository of
installable packages. Add it as a marketplace source in Salpa and install what you
need. Heavy execution routes to **Salpa Compute** — never gated here.

> **Status:** live. Packages are mirrored from their canonical dev repositories once
> they pass end-to-end verification and ship experimentalist docs. This repo is the
> generated/published mirror — packages are never edited Hub-first.

## Packages

| Package | What it does | License | Platforms |
|---|---|---|---|
| [`hello-world-pipeline`](packages/hello-world-pipeline) | Getting-started starter — a 3-node encrypt/decrypt/reveal pipeline (Caesar cipher) to try the Hub in seconds; pure-Python, no environment — 3 visual nodes | MIT | all |
| [`metalparm-vwf`](packages/metalparm-vwf) | Force-field parameterization for metal-containing systems (built on EasyParm) + metallopeptide fragment fusion + MD preparation — 20 visual nodes | LGPL-2.1 | linux · macOS |
| [`pdbmdauto`](packages/pdbmdauto) | Automated protein structure preparation for MD — homology modeling (ProMod3), protonation, solvation/ionization, staged GROMACS minimization — 14 visual nodes | MIT | linux · macOS |
| [`metaldock-vwf`](packages/metaldock-vwf) | Metal-protein docking — refactored MetalDock pipeline (protein/ligand prep, QM CM5 charges, metal-aware PDBQT, AutoDock4) — 6 visual nodes | ⚠ Academic / non-commercial † | linux · macOS (x86) |

> **†** "Free and open" here means **never paywalled** — not a commercial-use grant. `metaldock-vwf`
> is **academic / non-commercial only** while it depends on MGLTools (Scripps license); the planned
> Meeko swap will make it commercially clean. See [`LICENSING.md`](LICENSING.md).

## Install (in Salpa)

**Marketplace → Sources → Add Source** → point at this repository → **sync** →
**install** the package. Restart the server and worker so the registry discovers the
new nodes.

## Repository layout

- `packages/<pkg>/` — one folder per package, mirrored from its canonical dev
  repository (never authored here).
- `registry.json` — the generated catalog that Salpa reads.
- `LICENSING.md` — per-package license posture and attribution.

## Licensing

See [`LICENSING.md`](LICENSING.md). Every package ships free and open; the license is
finalized per package as its dependency chain settles.
