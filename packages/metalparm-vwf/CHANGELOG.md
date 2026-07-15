# Changelog

All notable changes to **metalparm-vwf** are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/) and the version
matches `package.toml`. Dates are ISO-8601.

## [1.29.0] — 2026-07-05

**Removed: the deprecated monolithic `ep_fragment_fuse` node.** It was
deprecated in v1.12.0 when it was split into `ep_fragment_fuse_topology`
(topology → `complex.prmtop`) + `ep_apply_coords` (coords →
`complex.rst7`/`complex.pdb`), and scheduled for removal in v1.13.0. The
removal is now done — 15 releases late. The two-node split is a strict
superset (reusable topology, decoupled re-runs, and the v1.19.0
interface-charge rebalance the monolith never gained), and the full
end-to-end `snp-complete-pipeline.spec.ts` has driven the split rather
than the monolith since v1.12.0.

- Deleted `ep_fragment_fuse/` and dropped it from `package.toml`
  (`[package.nodes]` now lists **20** nodes) and from `registry.json`.
- `fragment_align/node.py` now imports `DEFAULT_INTERFACE_BONDS` /
  `parse_interface_bonds` from `ep_fragment_fuse_topology.fuse_helpers`
  (was `ep_fragment_fuse.fuse_helpers`; the local fallback is unchanged).
- Descriptions/comments in `peptide_builder`, `fragment_align`, and
  `snp_builder` that named `ep_fragment_fuse` now point at
  `ep_fragment_fuse_topology`.
- Tests: removed `ep_fragment_fuse` from `test_node_packages.py`
  (`ALL_NODES`, `PURE_PYTHON_NODES`, `EXPECTED_DEMO_DATA`) and from
  `test_pipeline_e2e.py`'s node list (added the two split nodes there).
  `test_fragment_fuse_helpers.py` already loaded the topology node's
  helpers, so it is unaffected.
- Specs: deleted the standalone `snp-fragment-fuse.spec.ts` (it drove the
  monolith); `snp-complete-pipeline.spec.ts` retains coverage of the split
  fuser. Its dead `nodeIds['EpFragmentFuse']` Phase-1 deferral check now
  targets `EpFragmentFuseTopology` so it actually fires.

**Migration:** any saved workflow that wired `ep_fragment_fuse` must
replace it with `ep_fragment_fuse_topology → ep_apply_coords`. Wire
`ep_library_generation` + `fragment_align` into Fuse Topology, then Fuse
Topology + `fragment_align` into Apply Coords. Restart the BoCoFlow
server + worker after upgrading so the registry drops the removed node.

### Fixed — release-hardening from the node robustness audit

Three "silent success on failure" blockers where a node reported success
while emitting empty/garbage output, plus the cross-cutting error-message
gap behind them:

- **`ep_forcefield_assembly`**: `parmchk2` ran with no exit-code/output
  check, and the internal `run_script` helper swallowed each step's
  stderr. A parmchk2 (or step 03-13) failure now aborts with the tool's
  stderr tail, and the node refuses to report success unless
  `COMPLEX.frcmod` was actually produced (previously it emitted a phantom
  `output_frcmod` path pointing at a nonexistent file).
- **`ep_library_generation`**: `tleap` ran unchecked and `result.success`
  was hard-set `True` regardless of outcome. Success is now judged by a
  non-empty `COMPLEX.lib` (tleap exits nonzero on benign warnings, so
  exit code alone is unreliable); if the lib is missing the node aborts
  with tleap's log. `12_generate_lib.py` failures now surface stderr too.
- **`ep_mol2_generation`**: the antechamber / `xyz_to_pdb` / atom-typing
  subprocess calls funneled every failure through a bare exit code; they
  now route through a `run_step` helper that surfaces the tool's stderr.
- **`md_membrane_build`**: the `lipid` option advertised arbitrary
  packmol-memgen lipids but the residue classifier + tleap step are
  hardcoded to Lipid21 DPPC — a non-DPPC value (e.g. POPC) silently
  misrouted its tail residues into the solute and crashed on the
  atom-count guard. The node now rejects any non-DPPC `lipid` up front
  with an actionable message, and the option docstring says so.

## [1.28.0] — 2026-06-04

**Added: `md_membrane_build` exposes packmol-memgen box/iteration knobs.**
Two new OPTIONS surface the packmol-memgen flags needed when its auto-sized
XY box is too small for the lipid count, or its all-together packing
loop converges slowly:

| Option | packmol-memgen flag | Default | When to set |
|---|---|---|---|
| `xy_box_A` (Float) | `--distxy_fix` | `0` (auto) | Force XY box to fixed Å (bypass auto-size from `dist` + solute extent) |
| `nloop_all` (Int) | `--nloop_all` | `0` (default) | Bump all-together iterations when packing converges slowly |

Symptom this addresses: with a contracted-helix solute (longest axis
~33–34 Å), packmol-memgen auto-sized the XY box to ~54 Å — too tight
for ~150 DPPC molecules — and the all-together packing loop ran without
lowering its function value. Past Case-2 long-peptide runs with a
38.6 Å helix auto-sized to ~80 Å and converged with the same `--dist`.
Setting `xy_box_A=80` reproduces that geometry regardless of solute
extent.

`run_packmol_memgen(...)` in core.py grows two optional kwargs
(`xy_box_A`, `nloop_all`) that emit the corresponding flags only when
positive; passing them as `None`/`0` preserves the existing behaviour.
No change to the default packing for callers that don't opt in.

## [1.27.0] — 2026-06-02

**Changed: aligned `ep_orca_run`'s SLURM template handling with
`pdbmdauto/gmx_mdrun`.** Walked back the v1.26.0 backend prefill
override on `slurm_script` — the BoCoFlow front-end's
`NodeInstanceConfig.tsx` initialises every HPC textarea from
`node.config?._hpc_slurm_script ?? ''` and never consults the
schema's `default_value`, so v1.26.0's override was invisible to
users despite reaching `/api/node-schemas/search` and `/api/getnodes`
correctly. See `dev-notes/slurm-script-prefill-not-yet-supported.md`
for the GUI-layer gap, the backend-pipeline-is-fine evidence, and the
two-line upstream patch we are *not* taking right now.

What this release actually does:

- Removed the `_DEFAULT_SLURM_SCRIPT` read at import time and the
  `OPTIONS["slurm_script"]` override in `ep_orca_run/node.py`; the
  class now falls back to the inherited empty default from
  `HPCNodeBase.HPC_OPTIONS["slurm_script"]` — identical to
  `gmx_mdrun`'s practice.
- Deleted `ep_orca_run/demo_data/default-slurm.sh` (was only there
  to feed the override). `demo_data/OPTIMIZED.xyz` stays — it's the
  real demo XYZ.
- Kept `ep_orca_run/templates/default-slurm.sh` unchanged: the formal,
  Snellius-validated 94-atom Sn(IV)-porphyrin opt+freq+CHELPG SLURM
  script with `module load 2025` + `ORCA/6.1.0-gompi-2025a-avx2`,
  `genoa` partition, 128 GB memory, `${EBROOTORCA}/bin/orca`, and
  `$TMPDIR` node-local scratch. This is the file users now copy-paste
  into the GUI's SLURM Job Script field.
- Reverted `node.py` module docstring + `core.py` doc reference +
  `README.md` HPC paragraph to describe the copy-paste workflow
  explicitly. Each location notes the prefill gap and points at the
  dev-note.

No change to the SLURM template *content* — the substantive artifact
from v1.26.0 is preserved; only the prefill *mechanism* is reverted.

## [1.26.0] — 2026-06-02

**Changed: `ep_orca_run` ships a prefilled, formal SLURM template.** The
**SLURM Job Script** field on the ORCA Run node is no longer empty by
default — it is prefilled at module-import time from
`ep_orca_run/demo_data/default-slurm.sh`, a formal template validated on
Snellius (SURF) for a 94-atom Sn(IV)-porphyrin DFT opt+freq+CHELPG job.
The template:

- uses `{{VARIABLE}}` placeholders (`{{JOB_NAME}}`, `{{WORKING_DIR}}`,
  `{{ORCA_INPUT_FILE}}`, `{{ORCA_OUTPUT_FILE}}`, `{{RUN_LABEL}}`,
  `{{NPROCS}}`, ...) substituted by `HPCNodeBase` at submit time;
- invokes ORCA via `${EBROOTORCA}/bin/orca` so a silently-failed
  `module load` cannot turn into a `$(which orca)` empty-expansion that
  runs the `.inp` as a shell command;
- runs ORCA in node-local scratch (`$TMPDIR`) and copies only the small
  artifacts (`.out`, `.hess`, `.xyz`, CHELPG) back to the submit dir, so
  $HOME quota cannot trigger `Unable to write data in TBasis::WriteElement!`;
- requests 128 GB RAM and `--ntasks={{NPROCS}}` on the `genoa` partition
  by default (Snellius-specific lines flagged in comments for adjustment).

`templates/default-slurm.sh` is kept in sync with `demo_data/` as a
back-compat reference for documentation. Users can still paste their own
script — the prefilled default is a starting point, not a constraint.

Doc references in `node.py`, `core.py`, and `README.md` updated to point
to the new canonical location. Package `slurm_script` parameter docstring
explains the prefill behaviour.

## [1.25.0] — 2026-06-01

**Changed: short, slide-style node display labels.** Renamed 13 of the 15
demo-canvas nodes to drop the `EasyParm:` / `MD:` category prefix and
parenthetical qualifiers. The dropped qualifiers stay in each node's
tooltip (no information loss), and the stable `node_key` (e.g.
`EpBondDetection`) is unchanged — registry lookup, port ids, and edge
wiring are unaffected.

| `node_key` | old `display_name` | new `display_name` |
|---|---|---|
| `SnpBuilder` | SnP Fragment Builder | **SnP Builder** |
| `XtbOpt` | xTB Geometry Opt | **xTB Opt** |
| `EpBondDetection` | EasyParm: Bond Detection | **Bond Detection** |
| `EpMol2Generation` | EasyParm: MOL2 Generation | **MOL2 Generation** |
| `EpOrcaRun` | EasyParm: ORCA Run (local/HPC) | **ORCA Run** |
| `EpSeminarioOrca` | EasyParm: Seminario (ORCA) | **Seminario** |
| `EpForcefieldAssembly` | EasyParm: Force Field Assembly | **FF Assembly** |
| `EpLibraryGeneration` | EasyParm: Library Generation | **Library Generation** |
| `EpCharges` | EasyParm: Charges (RESP / CHELPG) | **Charges** |
| `EpFragmentFuseTopology` | Fragment Fuse — Topology (tleap) | **Fuse Topology** |
| `EpApplyCoords` | Apply Coordinates (ParmEd) | **Apply Coords** |
| `EpAmberToGromacs` | EasyParm: AMBER → GROMACS | **AMBER → GROMACS** |
| `MdMembraneBuild` | MD: Membrane Build (DPPC) | **Membrane Build** |

`PeptideBuilder` and `FragmentAlign` were already short — unchanged.

Spec files (`addNodeBySearch` + `clickNode` call sites), demo `.bcflow`
files (`pyworkflow.graph.nodes[i].name`), and the project node table in
`CLAUDE.md` are updated to match. Restart the BoCoFlow server + worker
after upgrading so the registry refreshes.

## [1.24.0] — 2026-06-01

**Added: `md_traj_center` `extract_first` option.** Emit the first
PBC-corrected frame as a standalone PDB alongside the centered
trajectory — useful as a topology/reference structure for downstream
analysis nodes and for visual sanity checks of the unwrap+center result.
(Backfilled changelog entry — this release was shipped without one.)

## [1.23.0] — 2026-05-22

**Added: `ep_charges` node** — inject QM-derived partial charges into the
zero-charge MOL2 from `ep_mol2_generation`, restoring the EasyParm charge
step the node refactor originally dropped (antechamber was wired but the
RESP/CHELPG charge assignment was not, so metal-fragment MD ran on
zero/near-zero partial charges). The node sits between MOL2 Generation and
FF Assembly: it awks CHELPG out of the ORCA `.out` (auto-discovered from
`ep_orca_run` via `output_out`) or runs `!RESP` / classic `.vpot`+`resp`,
then overwrites the MOL2 charges that `ep_forcefield_assembly` consumes.
Bundles `RESP_ORCA.py` + `Retrieve_RESP_Charges.py`; ships `core.py` and 7
unit tests (`tests/test_charges_core.py`). `package.toml` → v1.23.0,
`ALL_NODES` + README node table/topology synced. (Backfilled changelog
entry — this release was shipped without one.)

## [1.22.1] — 2026-05-21

**Fixed: `md_solvate_gmx` (and `md_membrane_build`) in single-canvas
workflows.** Two latent bugs surfaced only when one BoCoFlow workflow
owns the whole pipeline, so `ep_amber_to_gromacs` and the solvate/membrane
node share a single `working_path`. (Standalone use — separate workflows
with separate working paths, as in the per-step operation specs — was
unaffected, which is why this was not seen earlier.)

1. **Cross-node ITP-split import.** `save_gromacs_outputs` borrows
   `_split_top_into_itp` / `normalize_itp_basename` from the sibling
   `ep_amber_to_gromacs` node. Its script-mode fallback did
   `from core import …`, which resolves to *this* node's already-imported
   `core` (which lacks the helper) → `ImportError`. Now the sibling's
   `core.py` is loaded by file path under a unique module name. Applied to
   both `md_solvate_gmx/core.py` and `md_membrane_build/core.py`.

2. **Staging same-file copy** (`md_solvate_gmx/node.py`). `execute()`
   stages the input `.itp` to `output_dir/<name>.itp`; when the input is
   already in `output_dir` (shared working path) this is
   `shutil.copyfile(x, x)` → `shutil.SameFileError`. Staging now skips
   when source and destination resolve to the same file.

Regression coverage in `tests/test_md_solvate_gmx_core.py` exercises the
sibling ITP-split import in script mode (the failing context) and the
same-file staging guard.

## [1.22.0] — 2026-05-19

**Added: `md_traj_center`** — PBC-correct an MD trajectory before
analysis.

```
gmx_mdrun(_local) → md_traj_center → md_analysis_helix / md_analysis_distance
```

GROMACS writes coordinates wrapped into the periodic box, so a solute
that drifts across a box face is split across the boundary. DSSP (the
`md_analysis_helix` node) then mis-assigns secondary structure —
backbone H-bond geometry is wrong on a split molecule. (A
minimum-image *distance*, as in `md_analysis_distance`, is immune; the
helix content is not.) This node:

1. **Unwrap** — makes the solute whole again across PBC using the bond
   graph (so the topology must be a `.tpr` — a bare `.gro` carries no
   bonds).
2. **Center** — translates the solute centre of mass to the box
   centre.
3. **Wrap** — wraps the solvent / membrane back into the box.

It is the MDAnalysis equivalent of `gmx trjconv -pbc whole -center
-pbc mol`, done with MDAnalysis transformations so no GROMACS binary
is needed and the corrected trajectory is exactly what the (also
MDAnalysis-backed) analysis nodes see. Emits a PBC-corrected `.xtc` +
a first-frame `.gro`; passes the input `.tpr` through (atom order
unchanged) so the analysis nodes auto-discover a topology.

Recommended placement: insert between an MD-run node and the analysis
nodes so DSSP / distance always run on an unbroken metallopeptide.

**Fixed: `md_membrane_build` periodic box caused the first EM step to
diverge.** The box was taken from packmol's `inside box` constraint,
but packmol packs *non-periodically*: it never checks an atom against
a periodic image, and `inside box` is a *soft* (penalty) constraint
atoms leak 1–2 Å past. So atoms near opposite box faces became
periodic neighbours ~0.1 Å apart — the r⁻¹² wall gave an infinite
force and `gmx mdrun` energy minimisation diverged at step 0
(`LJ (SR) = 1e17`).

`assemble_membrane_system` now derives the periodic box from the
*actual packed-coordinate bounding box* + the packmol `tolerance` as a
margin (new `read_packmol_tolerance` helper), and shifts the
coordinates to sit centred in it. Every atom then lies within an
extent of `box − margin`, so no atom can approach a periodic image
closer than the tolerance — the same minimum packmol enforced for
in-box pairs. The box ends a few % larger than packmol-memgen
intended; the thin vacuum seam closes within the first ps of NPT.
The now-unused `read_packmol_box` (which read the wrong, too-small
constraint box) is removed.

Verified: a Case 2 DPPC system that diverged at EM step 0 now
minimises cleanly and runs MD.

## [1.21.1] — 2026-05-19

**Fixed: `md_membrane_build` reported atom counts as residue counts.**
`split_packed_pdb`'s `counts` are *atom* counts (per its docstring),
but `node.py` surfaced them as `n_lipid_residues` / `n_water` /
`n_ion` and logged "lipid residues" — so a 152-lipid DPPC bilayer was
reported as "19760 lipid residues" (its atom count). Added
`count_membrane_residues(system)` to `core.py`, which counts molecules
from the assembled ParmEd topology (Lipid21 DPPC is a split residue —
the PC headgroup count is the lipid count). The node result now
carries `n_lipid` / `n_water` / `n_ion` as true molecule counts.
Verified on the Case 2 long peptide: 152 DPPC, 6214 water, 28 ions.

## [1.21.0] — 2026-05-18

**Added: `md_membrane_build`** — embed a dry metallopeptide GROMACS
topology in a DPPC bilayer (transmembrane) and emit a solvated GROMACS
topology for membrane MD. The Case 2 counterpart of `md_solvate_gmx`.

### Pipeline

```
ep_amber_to_gromacs → md_membrane_build → (gmx grompp / mdrun)
```

1. **Orient** — align the peptide helix axis to z (PCA on the Cα
   atoms). MEMEMBED (packmol-memgen's default orienter) fails on a
   metallopeptide — the non-standard fragment residue breaks it — so
   the node pre-orients geometrically and runs packmol-memgen
   `--preoriented`.
2. **Pack** — packmol-memgen builds the DPPC bilayer + water + ions
   around the solute, geometry only (no `--parametrize`). The solute
   passes through with `--notprotonate --nottrim`.
3. **Split** — partition the packed PDB into the solute part and the
   membrane part (DPPC + water + ions — all *standard* residues).
   `TER` records are preserved so tleap does not bond adjacent lipids
   into one chain.
4. **Parametrise** — tleap on the membrane part only
   (`leaprc.lipid21` + `leaprc.water.opc`). Standard residues, so no
   re-derivation hazard; the non-standard solute is never tleap'd.
5. **Assemble** — ParmEd-concatenate the *preserved* solute topology
   (from the GROMACS `.top`) with the membrane topology; box read from
   packmol-memgen's `packmol.inp` `inside box` constraints.
6. **Export** — GROMACS `.top` + `.gro` (+ optional `.itp` split).

Like `md_solvate_gmx`, the solute topology is the source of truth and
is never round-tripped through tleap — fragment-fused residues survive.

### Verified

Live on `03-fusion-long-peptide`'s charge-neutral topology: 448-atom
solute → DPPC bilayer (152 lipids) + 6214 OPC water + 28 ions =
**45 092 atoms**, box **7.30 × 7.30 × 8.10 nm**, system charge
+0.00004 e (neutral). Regression test for `read_packmol_box` in
`tests/test_md_membrane.py`.

## [1.20.0] — 2026-05-18

**Added: two MD trajectory-analysis nodes** — `md_analysis_helix` and
`md_analysis_distance`. Together they turn a finished MD trajectory
into the SnP-peptide case's scientific deliverable.

### `md_analysis_helix` — α-helix content

Per-frame DSSP secondary structure of the peptide along the
trajectory. Reads `.xtc`/`.trr` + `.tpr`/`.gro` (set explicitly or
auto-discovered from an upstream MD-run node), runs MDAnalysis' DSSP,
and writes `<case>_helix.csv` (frame, time, helix-residue count,
helix fraction) plus a per-residue helix propensity and a
mean/std/min/max summary. Backbone-incomplete residues (ACE/NME caps,
ions) are auto-dropped — DSSP requires a complete N/CA/C/O backbone
per residue, so the user need not hand-tune the selection.

### `md_analysis_distance` — residue ↔ metal distance

Per-frame minimum-image (PBC-aware) distance from a set of probe
atoms — by default every tyrosine hydroxyl O and histidine ring N,
the PCET quenchers — to the metal centre (Sn). Writes
`<case>_distance.csv` (one Å column per probe + the per-frame
minimum) and a closest-approach summary (value / probe / frame).

### Dependency

`mdanalysis >=2.8` added to the `metalparm_vwf` pixi env (the DSSP
module needs ≥2.8). Run `pixi install` on the env after pulling.

### Verified

Pure helpers covered by 9 unit tests (`tests/test_md_analysis.py`).
Both `run_*_analysis` passes verified live on
`09-md-case1-relax-and-run/md.{tpr,xtc}` (the heptapeptide 100 ps
relaxation): helix mean 67.5 %, distance probes HID2/TYR3/TYR4,
closest approach 7.6 Å.

## [1.19.0] — 2026-05-18

**Added: interface charge rebalancing in `ep_fragment_fuse_topology`** —
a new `charge_rebalance` parameter (default **on**) that makes every
residue's net charge — and the complex total — an exact integer after
the fuse.

### Why

`ep_fragment_fuse_topology` forms the linkage bond by deleting interface
cap atoms (an interface bond's `pep_remove` / `frag_remove`). Each
deleted atom carries a partial charge, so the residue it was removed
from is left non-integer — and `tleap`'s `saveamberparm` does not
repair it. The SnP metallopeptide lands at a net **−0.1812 e**, all of
it on the linkage GLU: deleting `OE2` (the −0.8188 e carboxylate
oxygen) from a −1 glutamate template leaves `−1.000 + 0.8188 = −0.1812`.
`grompp` then warns *"System has non-zero total charge"* and PME runs
with a neutralising background — a real artefact right at the
chemically interesting GLU–porphyrin junction.

Full analysis: `dev-notes/fragment-fusion-interface-charge.md`,
metal-md issue #4.

### What it does

After `tleap`, the node loads `complex.prmtop` with ParmEd and, for
every residue, redistributes the remainder `round(q) − q` equally over
the residue's atoms (`fuse_helpers.redistribute_to_integer` /
`rebalance_residue_charges`). The SnP GLU goes `−0.1812 → 0`; the
complex total goes to an exact `0`. A residue more than `0.4 e` from
any integer aborts the node — that is a genuine parameterisation error
(wrong protonation state), not a small deletion remainder, and must not
be silently snapped to the wrong integer.

This is a charge **rebalance**, not a re-derivation: the per-atom
correction is tiny (the SnP GLU's −0.1812 e spread over 14 atoms =
−0.013 e/atom). A proper RESP re-derivation of the linked residue, and
re-typing the now-amide `CD`/`OE1` from carboxylate types, remain
possible future refinements (tracked in the dev-note / issue).

### Tests

6 new unit tests for `redistribute_to_integer` in
`tests/test_fragment_fuse_helpers.py` (no-op when already integer,
remainder spread, nearest-integer rounding, large-imbalance abort,
empty, input not mutated). The test module now targets the live
`ep_fragment_fuse_topology` helpers. 22 fuse-helper tests pass; live
`rebalance_residue_charges` verified against the SnP `complex.prmtop`
(−0.18120 → 0.000000, GLU6 → 0).

## [1.18.0] — 2026-05-18

**Added: `cap_peptide_termini` parameter on `snp_builder`** — when
"Extract Peptide" is on, geometrically adds an **ACE** N-cap and an
**NME** C-cap to the extracted chain (`Ac-…-NH-CH₃`) and drops the
charged-terminus atoms. Default off (verbatim termini).

### Why

PDB-mode `peptide_builder` runs only `loadpdb` — it never applies
ACE/NME; its `n_term`/`c_term` cap options exist only for sequence
mode. So the genuine Case 1 heptapeptide carved by `peptide_residue_range`
(v1.17.0) came out with **charged zwitterionic termini** (HID1 = NH₃⁺,
ALA7 = COO⁻), not the case spec's `Ac-…-NH₂`. For a 7-residue peptide
whose folding is the scientific question, terminal charges oppose the
helix macrodipole and suppress helix propensity — a confound. Capping
upstream, in `snp_builder`'s extraction, fixes it while keeping the
helical truncation coordinates.

### What it does

- `extract_peptide` (`metal_fragment.py`) gained a `cap_termini`
  keyword. When True it places an ACE before the first written residue
  and an NME after the last, via NeRF internal-coordinate placement
  (`_place_from_internal`) off the backbone N/CA/C atoms — bond lengths
  exact, helical-ish backbone dihedrals (φ ≈ -60°, ψ ≈ -47°). Cap atoms
  use the exact ff19SB template names (ACE: `C O CH3 H1 H2 H3`; NME:
  `N H C H1 H2 H3` — the NME methyl carbon is `C`, not `CH3`).
- The redundant charged-terminus atoms are dropped: N-terminal
  `H1/H2/H3` (tleap rebuilds the single amide H), C-terminal `OXT`.
- All serials are renumbered sequentially in the capped output.
- The returned residue count includes the 2 caps.

Verified: a capped extraction loads cleanly through `tleap`
(`leaprc.protein.ff19SB` → `saveamberparm`, 0 errors).

### Tests

6 new unit tests in `tests/test_snp_builder.py` (ACE/NME added,
template atom names for each cap, charged N-terminal protons dropped,
amide-bond lengths, default-off). All 35 `snp_builder` tests pass.

## [1.17.0] — 2026-05-18

**Added: `peptide_residue_range` parameter on `snp_builder`** — when
"Extract Peptide" is on, the extracted chain can be restricted to an
inclusive residue-number span written `lo-hi` (e.g. `1-7`). Empty =
whole chain (unchanged default behaviour).

### Why

`collect/simulation-case-request/snpp.pdb` is the **23-residue long
peptide** (Case 2's solute) — confirmed by 23 `SEQRES` residues, `TER`
at `ALA A 23`, and a folded α-helix. The verbatim request note states
*"the short peptide-SnPP is a part of the long peptide-SnPP one"*: the
user supplied only the long structure, and the genuine Case 1
heptapeptide (`Ac-HYYLA-E[Sn(OCH₃)₂P]-A-NH₂`) is **residues 1–7** of
that chain. It has no structure file of its own.

`peptide_residue_range` carves it out directly — `snp_builder` with
`peptide_residue_range = "1-7"` emits the residues-1–7 sub-peptide,
keeping the YASARA helical coordinates so the Case 1 MD starts from a
realistic (helical) conformation. This is a general node enhancement,
not a hand-edited `snpp_short.pdb`: any "model a sub-segment of a
longer chain" case is now covered.

### What it does

- `extract_peptide` (`metal_fragment.py`) gained a `residue_range`
  keyword. When given `(lo, hi)`, only ATOM/HETATM records whose
  `resSeq` parses to an int in `[lo, hi]` are written.
- The chain's original `TER` (at the true chain end) is dropped and a
  single **synthetic `TER`** is emitted after the kept atoms, so the
  truncated chain terminates cleanly for `tleap loadpdb`.
- A range matching no residues raises `ValueError` (not a silent empty
  file); `lo > hi` is rejected at parse time with a clear setup error.

### Tests

5 new unit tests in `tests/test_snp_builder.py` (span-only kept,
coordinates preserved, synthetic-TER placement, no-match raises,
inverted-range raises). All 29 `snp_builder` tests pass.

## [1.16.0] — 2026-05-12

**Added: `md_solvate_gmx` node** — GROMACS-side solvator that bypasses
the tleap re-derivation bug in `md_solvate_packmol`. Consumes the dry
GROMACS topology from `ep_amber_to_gromacs` (`complex.top` +
`complex.gro` + optional `metallopeptide.itp`) and produces a fully
solvated `complex.{top,gro}` + `metallopeptide_solv.itp` ready for
`gmx grompp`.

### Why

`md_solvate_packmol` (v1.15.0) wraps `packmol-memgen`, which
internally invokes `tleap` to (re-)build the protein topology before
packing solvent. tleap re-types every residue against the standard
AMBER library; for non-standard residues (e.g. a fragment-fused GLU
that has lost OE2 in step 03), tleap auto-completes the "missing"
atom from the standard library at the standard geometry — landing
0.13 Å from the fragment N1 with no covalent bond. Step 0 of MD
explodes: LJ-SR = 1.07 × 10¹⁷ kJ/mol, EM diverges to inf force on the
re-introduced atom.

The user rejected "rename the fused residue to GL5" as a layering
workaround that doesn't generalise. The proper fix is to **never
re-derive solute chemistry downstream of step 03**. `md_solvate_gmx`
implements that: solute topology from `ep_amber_to_gromacs` is the
source of truth, raw `packmol` does coordinate placement only, and
ParmEd assembles the final Structure via `Structure + Structure` /
`Structure * N` — no tleap round-trip on the solute, no re-typing.

### Pipeline slot

```
Old (md_solvate_packmol):
ep_apply_coords ─▶ md_solvate_packmol ─▶ ep_amber_to_gromacs ─▶ (gmx)
                    [tleap re-types the solute — bug for non-standard residues]

New (md_solvate_gmx):
ep_apply_coords ─▶ ep_amber_to_gromacs ─▶ md_solvate_gmx ─▶ (gmx)
                    [solute → GROMACS once,         [pack solvents
                     never re-derived again]         around it]
```

`md_solvate_packmol` is kept indefinitely for standard-residue protein
systems where the tleap re-typing is harmless and the
packmol-memgen-bundled cosolvents are convenient.

### What it does

1. Loads the dry solute as a `parmed.Structure` via
   `parmed.load_file(top, xyz=gro)` (preserves topology bytewise).
2. Builds single-molecule prmtops for the requested solvents
   (OPC water from `solvents.lib` via `leaprc.water.opc`; methanol
   from `MEOHBOX` in the same lib, extracted via ParmEd `:1`; ions
   from `atomic_ions.lib`) in one `tleap` call (~3-5 s).
3. Computes solvent counts from a box-volume target + a
   user-supplied **molar** ratio (bulk densities: MeOH 14.8 mol/nm³,
   H₂O 33.4 mol/nm³), and ion counts from the solute net charge
   plus `saltcon_M` salt concentration.
4. Writes a `pack.inp` with the same packmol options as
   `packmol-memgen` uses (tolerance 2.0, `nloop 20`, `radius 1.5` on
   solute, `add_amber_ter`), centred at the origin.
5. Runs raw `packmol` to place 1 solute + N solvents + ions in the
   box.
6. Replicates each per-unit Structure via `Structure * N`,
   concatenates via `Structure + Structure` (preserves all atom-types
   + bonded params; no re-derivation), then overwrites the
   coordinates from packmol's PDB by atom index.
7. Saves as GROMACS via `Structure.save(..., format="gromacs")` and
   splits the moleculetype block into a sibling `.itp` (reusing
   `ep_amber_to_gromacs.core._split_top_into_itp`).

### OPTIONS schema

| Field | Type | Default | Notes |
|---|---|---|---|
| `case_name` | StringParameter | `"complex"` | |
| `output_dir` | FolderParameter | — | required |
| `input_top` | FileParameterEdit | `""` | auto-discovers `output_top` from predecessor |
| `input_gro` | FileParameterEdit | `""` | auto-discovers `output_gro` |
| `input_itp` | FileParameterEdit | `""` | optional; auto-discovers `output_itp` |
| `solvents` | StringParameter | `"MOH:WAT"` | currently 'WAT' or 'MOH:WAT' |
| `solvent_ratio` | StringParameter | `"2:1"` | **molar** ratio (explicit, not basis-unspecified) |
| `padding_A` | FloatParameter | `12.0` | |
| `water_model` | SelectParameter | `"opc"` | opc / opc3 / tip3p / tip4pew / spce |
| `cation` | SelectParameter | `"K+"` | K+ / Na+ |
| `anion` | SelectParameter | `"Cl-"` | |
| `saltcon_M` | FloatParameter | `0.0` | above neutralisation |
| `random_seed` | IntegerParameter | `-1` | -1 = packmol picks |

### Forwarded data (matches `ep_amber_to_gromacs`'s keys)

`output_top`, `output_gro`, `output_itp`, `case_name`, `working_path`,
`solvent_counts`, `box_dimensions_nm`, `solute_charge`.

### Tests

24 new unit tests in `tests/test_md_solvate_gmx_core.py` cover the
pure-Python math (parse_solvent_ratio, compute_box_dimensions_nm,
compute_solvent_counts, compute_ion_counts, write_packmol_input).
The tleap/ParmEd-dependent functions are exercised by the live
BoCoFlow spec.

### Live verification

`simulations/snp-peptide-md/operations/05b-md-case1-meoh-water-gmx.spec.ts`
runs the node against `03-fusion-long-peptide`'s outputs. The
follow-up `06-md-case1-relax-and-run` spec then runs EM + 100 ps NPT
MD against the solvated topology — **EM must converge to finite max
force** (no inf force on atom 103, in contrast to the
`md_solvate_packmol` output where it does).

### Changed — md_solvate_packmol: drop input staging, pass abs paths

Refactor follow-up after the v1.15.0 first-cut: the node was
defensively copying `input_pdb` / `fragment_lib` / `fragment_frcmod` /
`linkage_frcmod` into its `output_dir` and then passing basenames to
packmol-memgen. Empirical test (2026-05-11) confirmed packmol-memgen
accepts **absolute paths** for `--pdb`, `--ligand_param FRCMOD:LIB`,
and abs-path references inside `--leapline`. So we now resolve the
prefix-tagged BoCoFlow paths via `self.resolve_path()` and pass the
resulting filesystem paths straight to packmol-memgen — no copy.

Why: the BoCoFlow `abs:` / `rel:` / `node:` model already gives nodes
a real filesystem path after `resolve_path()`. Any tool that accepts
abs paths (which is almost all of them, including the entire AMBER
toolchain) needs no staging. See
[`dev-notes/node-io-helpers-duplication.md`](../../dev-notes/node-io-helpers-duplication.md)
for the full discussion + the duplication audit across 13 of 15
metalparm-vwf nodes.

User-visible impact: bytes-identical solvated topology, but the
node's output directory is now leaner — no `fragment.lib` /
`fragment.frcmod` / `linkage.frcmod` copies polluting the run dir.
Upstream artifacts read in place. Re-verified live in
`05-md-case1-meoh-water.spec.ts` (workflow `b4a5cf84-…`, 1.4 min).

### Removed

- `md_solvate_packmol/helpers.py::ensure_in_workdir()` — unused after
  the abs-path refactor; replaced by an explanatory comment pointing
  at the dev-note. `get_from_predecessors()` stays; it's a separable
  predecessor-data convenience that's broadly useful.

## [1.15.0] — 2026-05-11

**Added: `md_solvate_packmol` node** — first member of the MD-prep
sub-area of the package. Wraps AmberTools' `packmol-memgen` to
solvate a metallopeptide AMBER topology in a single- or mixed-solvent
box with proper ions, all in one tool call.

### Why

The case-folder for SnP-peptide Case 1 ends at
`ep_amber_to_gromacs` with the dry metallopeptide topology
(`complex.top` + `complex.gro` + `metallopeptide.itp`). To run real
production MD in a MeOH:H₂O 2:1 mixed solvent, we need a box-build +
solvent-pack + ion-neutralize step that:

1. Accepts ff19SB-compatible inputs (the existing
   `bocoflow-nodes/installed/gmx_solv_ion` hardcodes
   `[amber99sb, charmm27, oplsaa, gromos53a6]` and rebuilds via
   pdb2gmx, discarding our SnP fragment FF).
2. Supports **mixed** solvents at arbitrary ratios.
3. Pairs ff19SB with **OPC water** (the FF's design pairing per
   Tian et al. 2020) by default.

### How

Research into the AMBER ecosystem (recorded in
[`dev-notes/solvent-ff-amber-ff19sb-research.md`](../../dev-notes/solvent-ff-amber-ff19sb-research.md))
showed `packmol-memgen` (already bundled in AmberTools, available
in this package's `metalparm_vwf` pixi env) ships an AMBER-
compatible cosolvent library covering MOH (Cieplak et al. 2001),
CL3, DMS, NMA, ACN, ACT, BNZ, IPH, TFE — and natively supports
`--solvents A:B --solvent_ratio X:Y` for mixed boxes. So the node
is a thin (~280 line) wrapper that translates BoCoFlow OPTIONS
into the right `packmol-memgen` CLI flags, runs the tool, parses
the log for molecule counts, and forwards `output_prmtop` +
`output_rst7` + `output_pdb` to the next node.

### Node placement in the workflow

```
ep_apply_coords (dry metallopeptide AMBER topology)
    → md_solvate_packmol (solvated AMBER topology)
        → ep_amber_to_gromacs (GROMACS conversion)
            → gmx_md_relax (EM + equilibration)
                → gmx_mdrun_local (production MD)
```

The node sits **before** `ep_amber_to_gromacs`, not after — by
keeping the solvation in AMBER, tleap handles the topology merge
(no Python-side topology stitching needed), and the existing
`ep_amber_to_gromacs` cleanly converts the unified solvated
topology to GROMACS for downstream MD nodes.

### Added

- `md_solvate_packmol/` (node.py + meta.toml + demo_data/)
- `package.toml` nodes list grew to 15 entries.

### OPTIONS

`case_name`, `output_dir`, `input_pdb`, `fragment_lib`,
`fragment_frcmod`, `linkage_frcmod`, `solvents`, `solvent_ratio`,
`padding_A`, `ff_protein` (default ff19SB), `water_model` (default
**opc**), `saltcon_M`, `extra_leaplines`. Auto-discovers all input
files from upstream predecessor data when blank.

### Outputs

`output_prmtop`, `output_rst7`, `output_pdb`, `solvent_counts`
(dict of per-solvent + per-ion molecule counts parsed from the
packmol-memgen log), `box_dimensions_A` (lx, ly, lz from rst7),
`packmol_memgen_log` path.

### Verified during development

Smoke-tested 2026-05-11 against the existing Case 1 outputs at
`simulations/snp-peptide-md/03-fusion-long-peptide/`:
`packmol-memgen --solvents MOH:WAT --solvent_ratio 2:1 --dist 12
--ffprot ff19SB --ffwat opc` produced a solvated topology with
3036 MeOH + 4540 H₂O + 1 K⁺ neutralizer (the existing K⁺ default
of packmol-memgen; switch to NaCl via `--salt_c Na+ --salt_a Cl-`
in a future extension). End-to-end pipeline live verification
pending (paired narrative MD + Playwright spec for the case
folder lands separately).

## [1.14.0] — 2026-05-09

**Added: bundled `amide_glh_gaff2_n.frcmod` linkage** for the
GLH (protonated/neutral glutamic acid) form of the cross-FF amide
patch.

### Why

`amide_glu_gaff2.frcmod` (v1.9.0) and `amide_glu_gaff2_n.frcmod`
(v1.13.0) target the GLU residue, where ff19SB types the carboxyl C
as `CO` and the carboxylate O as `O2`. But sequence-mode
`peptide_builder` building a peptide with neutral side chains (the
typical case at neutral pH for downstream MD with explicit
protonation) uses the **GLH** residue template — and ff19SB types
GLH's carboxyl C as `C` (not `CO`) and OE1 as `O` (not `O2`). The
GLU-keyed linkages don't apply, and tleap aborts with `No torsion
terms for atom types: 2C-C-n-hn` (or the angle/improper analogues).

The new file is the GLH analogue of v1.13.0's `_n` linkage: same
numerical values, only the carboxyl-side atom types change
(`CO → C`, `O2 → O`). N-side stays GAFF2 type `n` (same antechamber-
on-real-DFT typing situation).

### How to choose

Pick by your peptide.lib's residue template AND your fragment.lib's
linkage-N type:

| Peptide residue | Glu's CD type | Fragment linkage-N type | Use linkage |
|---|---|---|---|
| GLU (charged) | CO | `ns` (sp2 amide) | `amide_glu_gaff2.frcmod` |
| GLU (charged) | CO | `n`  (sp3 amine) | `amide_glu_gaff2_n.frcmod` |
| GLH (protonated) | C  | `n`  (sp3 amine) | `amide_glh_gaff2_n.frcmod` (NEW) |

The `peptide_builder` node's `peptide_mode = sequence` with `GLH`
written in the sequence string produces the third row above — the
typical sequence-built case.

### Added

- `ep_fragment_fuse_topology/demo_data/linkages/amide_glh_gaff2_n.frcmod`
- README + `dev-notes/cross-ff-linkage-frcmods.md` updated with the
  expanded type-selection table.

### Migration

No code change required for existing workflows already using
`amide_glu_gaff2{,_n}.frcmod`. Sequence-mode peptides written with
GLH at the SnP-bearing residue (the natural choice at neutral pH)
should now use `amide_glh_gaff2_n.frcmod` instead. First production
verification:
`metal-md/simulations/snp-peptide-md/operations/04-fusion-case2-long.spec.ts`
(2026-05-09).

## [1.13.0] — 2026-05-09

**Added: bundled `amide_glu_gaff2_n.frcmod` linkage** for the
`n`-typed (sp3 amine) variant of the GLU-amide cross-FF patch.

### Why

`amide_glu_gaff2.frcmod` (added in v1.9.0) targets GAFF2 atom type
`ns` — the demo SnP fragment lib uses that type for its linkage
nitrogen. But when the easyPARM FF leg is run against a real-DFT
QM Hessian (the standard production path for new metals), antechamber
routinely types the linkage N as `n` (sp3 amine) instead of `ns`
(sp2 amide). Without an `n`-typed parallel linkage, tleap aborts at
`ep_fragment_fuse_topology` with `No torsion terms for atom types:
2C-CO-n-hn` and the workflow halts.

The new file is the same numerical patch keyed to type `n` instead of
`ns` — same bond/angle/dihedral/improper values, only the atom type
strings change. Use whichever linkage matches your fragment lib's
typing.

### How to choose

Open your fragment `.lib` (or `.mol2`) and look at the linkage
nitrogen's GAFF2 type:

| Linkage-N type in fragment lib | Use this linkage |
|---|---|
| `ns` (sp2 amide N — typical of demo SnP / hand-edited libs) | `amide_glu_gaff2.frcmod` |
| `n`  (sp3 amine N — typical of antechamber-on-real-DFT)     | `amide_glu_gaff2_n.frcmod` |

Both ship under `ep_fragment_fuse_topology/demo_data/linkages/` and
are accessible via the `Linkage Frcmod` parameter as
`node:demo_data/linkages/<filename>`.

### Added

- `ep_fragment_fuse_topology/demo_data/linkages/amide_glu_gaff2_n.frcmod`
- README + `dev-notes/cross-ff-linkage-frcmods.md` updated with the
  type-selection rule and the new file's role.

### Changed

- (none — additive only)

### Migration

No code change required for existing workflows already using
`amide_glu_gaff2.frcmod`. Users hitting the
`No torsion terms for … 2C-CO-n-…` tleap error on a real-DFT FF leg
should switch their `Linkage Frcmod` to
`node:demo_data/linkages/amide_glu_gaff2_n.frcmod`. First production
verification:
`metal-md/simulations/snp-peptide-md/operations/03-fusion-long-peptide.spec.ts`
(workflow `399211cc-…`, 2026-05-09).

## [1.12.0] — 2026-05-07

**Split `ep_fragment_fuse` into topology-only + coordinate-applier nodes**
(closes metal-md#2).

The old monolithic node conflated two genuinely separable concerns:
1. *Topology fusion* — merge peptide.lib + fragment.lib into
   complex.prmtop via tleap.
2. *Coordinate assembly* — write the merged coords into complex.rst7 /
   complex.pdb.

Splitting them gains:

- **DAG honesty.** The old `fragment_align → ep_fragment_fuse` edge
  *looked* like a topology dependency but was actually carrying both
  peptide passthrough and aligned-coords through the lib's positions
  table. The split version makes each edge mean exactly one thing.
- **Decoupled re-runs.** Tweak interface_bonds → re-run topology only
  (cheap). Tweak alignment → re-run coords only (cheaper still). Today
  any change forced a full re-fuse.
- **Reusable topology.** One complex.prmtop can now be applied to
  multiple starting coordinate sets (50 docking poses, replica-exchange
  seeds, free-energy decoys) by re-running `ep_apply_coords` alone.

### Added

- **`ep_fragment_fuse_topology`** (`v1.0.0`) — same `OPTIONS` and tleap
  script as the old fuse, but only forwards `output_prmtop` (plus
  `peptide_residues` / `interface_bonds` / `case_name` / `working_path`
  for downstream chaining). The auto-generated rst7+pdb that tleap's
  `saveamberparm` always writes are renamed to `<case>_initial.rst7`
  and `<case>_initial.pdb` to make clear they're not the canonical
  coords (which come from `ep_apply_coords`).

- **`ep_apply_coords`** (`v1.0.0`) — ParmEd-based coord-applier. Loads
  a topology-only prmtop, walks its atom order, looks up each atom in
  the aligned source PDBs by `(residue_idx, atom_name)` (peptide
  source for residues 1..N, fragment source for residues N+1..),
  assigns coords, writes the final `complex.rst7` + `complex.pdb`.
  Atoms that tleap's `remove cpx.X.Y` operations dropped (caps OE2/HE2
  on peptide GLU, CM/HM1-3/CAP/OAP on the SnP fragment) are simply
  skipped — they aren't in the prmtop's atom list.

- 15 unit tests in `tests/test_apply_coords_core.py`: PDB parsing,
  atom-index dedup, residue-bounded source-PDB switch, removed-atom
  tolerance, missing-atom error reporting (peptide-side and
  fragment-side), translation round-trip, fragment-resid offset.

### Changed

- **`peptide_builder`** now also emits the role-specific
  `output_peptide_pdb` key alongside the generic `output_pdb`. Allows
  `ep_apply_coords` to disambiguate when both `peptide_builder` and
  `fragment_align` are upstream of the same consumer (the latter also
  emits `output_pdb` for the aligned fragment).

- **`fragment_align`** passthrough list grows to include
  `output_peptide_pdb`. With this, `ep_apply_coords` only needs two
  direct predecessors (`ep_fragment_fuse_topology` + `fragment_align`)
  rather than three (which would otherwise need to wire
  `peptide_builder → ep_apply_coords` separately).

- **`operations/snp-complete-pipeline.spec.ts`** updated from 12 nodes
  / 13 edges to 13 nodes / 14 edges. Old `[10, 11]` (fuse → amber) is
  now `[10, 11]` (fuse_topology → apply_coords) + `[9, 11]`
  (fragment_align → apply_coords) + `[11, 12]` (apply_coords → amber).
  The on-disk artifact set is unchanged
  (`complex.{prmtop,rst7,pdb,top,gro,itp}`).

- **`CLAUDE.md`** package descriptions and node table updated to
  reflect the 14-node count and the deprecation of the old fuse.

### Deprecated

- **`ep_fragment_fuse`** is marked deprecated (`v0.2.0`, display name
  "Fragment Fuse (deprecated)") and will be removed in `v1.13.0`.
  Existing saved workflows that contain this node continue to work
  through `v1.12.0`. New workflows should use the 2-node split.

### Migration

For users with workflows that wire today's `ep_fragment_fuse`:

```
old:  ep_library_generation ─┐                          ┌→ ep_amber_to_gromacs
                              ▼                          │
                         ep_fragment_fuse ────────────────┘
                              ▲
              fragment_align ─┘

new:  ep_library_generation ─┐
                              ▼
                         ep_fragment_fuse_topology ─┐
                              ▲                     ▼
              fragment_align ─┴───────────────► ep_apply_coords ─→ ep_amber_to_gromacs
```

Replace the single `ep_fragment_fuse` node with `ep_fragment_fuse_topology`
+ `ep_apply_coords` in series. Wire `fragment_align` to BOTH (the
topology side keeps it as a predecessor for peptide-side passthrough,
the coord side reads aligned PDBs from it). Configuration carries
over — case_name / output_dir / linkage_frcmod / interface_bonds stay
on the topology node; apply_coords only needs case_name + output_dir.

## [1.11.0] — 2026-05-07

**Added: optional `.itp` export from `ep_amber_to_gromacs`** for mixing
the metallopeptide with other molecules in downstream GROMACS MD setups.

### Added

- New `ITP Filename` (`itp_filename`) `StringParameter` on
  `ep_amber_to_gromacs`. When set, the node post-processes the
  ParmEd-generated `.top` and emits a self-contained
  `<itp_filename>.itp` containing the metallopeptide's `[ moleculetype ]`
  block (atoms, bonds, pairs, angles, dihedrals, exclusions, position
  restraints — whatever ParmEd writes between `[ moleculetype ]` and
  `[ system ]`). The master `.top` then `#include`s the `.itp`
  alongside the unchanged `[ defaults ]`, `[ atomtypes ]`, `[ system ]`,
  and `[ molecules ]` sections — so it can be merged into a multi-
  molecule master topology (water, ions, lipids, second peptide copy,
  etc.) without re-engineering by hand. Both `'complex'` and
  `'complex.itp'` resolve to `complex.itp` (trailing extension is
  stripped before re-appending). Path separators are rejected — the
  param is a basename only; the .itp lives next to the .top.

- `output_itp` is forwarded on `result.data` and `result.files["output"]`
  (`None` / absent when the split is skipped).

- 6 new tests in `tests/test_amber_to_gromacs_core.py`
  (`test_normalize_itp_basename_strips_extension`,
  `test_normalize_itp_basename_rejects_path_separator`,
  `test_split_top_into_itp_unit`,
  `test_split_top_into_itp_missing_block_raises`,
  `test_round_trip_with_itp_split`,
  `test_itp_filename_none_skips_split`,
  `test_itp_filename_empty_string_skips_split`)
  + 1 new E2E test `TestAmberToGromacsE2E.test_tleap_to_gromacs_with_itp_split`.

### Changed

- Default behavior unchanged when `itp_filename` is empty: identical
  `.top` + `.gro` to v1.10.0. No DAG / config-hash impact for existing
  workflows.

### Why

User requested mixing the metallopeptide with solvent / membrane / ions,
which requires `#include`-able topology rather than the monolithic
`.top`. This is also how AMBER → GROMACS workflows are conventionally
consumed in published MD protocols (CHARMM-GUI, `gmx pdb2gmx`, OPLS-AA
tutorials all expect per-molecule `.itp`).

### Why `[ atomtypes ]` stays in the master `.top`

GROMACS requires atomtypes at top-level scope **before any
`[ moleculetype ]` directive**. Putting them inside the `.itp` makes the
`.itp` non-`#include`-able after another moleculetype has already
opened in the master topology. Keeping them in the `.top` makes the
`.itp` a portable fragment that the user can include from any master
`.top` without ordering surprises. (A future split that emits a separate
`<itp>_atomtypes.itp` is on the table for multi-custom-molecule
topologies but out of scope here.)

## [1.10.0] — 2026-05-05

**New node `ep_amber_to_gromacs`** — converts the AMBER topology pair
produced by `ep_fragment_fuse` (`complex.prmtop` + `complex.rst7`) into
GROMACS format (`complex.top` + `complex.gro`) for downstream GROMACS
MD simulation. Implementation mirrors easyPARM's
`collect/easyPARM/scripts/amber_converter.py:71-83`: load the AMBER
pair with [ParmEd](https://github.com/ParmEd/ParmEd) (already an
AmberTools dependency), call `Structure.save(format='gromacs')` and
`Structure.save(format='gro')`. Node count: 11 → 12.

### Added

- **`packages/metalparm-vwf/ep_amber_to_gromacs/`** — new node:
  - `core.py` (~85 lines): pure helper
    `convert_amber_to_gromacs(prmtop, rst7, output_prefix, *,
    add_box_if_absent=True, box_padding=10.0)`. Returns a stats dict
    with paths + atom count + box-handling info.
  - `node.py` (~165 lines): BoCoFlow wrapper. 3-tier predecessor
    resolution (explicit option → `output_prmtop` / `output_rst7` from
    a `ep_fragment_fuse` predecessor → setup error). Forwards
    `output_top` / `output_gro` for downstream consumers.
  - `meta.toml` with the standard node metadata (display name "EasyParm:
    AMBER → GROMACS", shared environment `metalparm_vwf`).

- **Cubic-box auto-add** when the prmtop is non-periodic (the typical
  fuse output). Sized to the molecule's extent + 2 × `box_padding` (Å).
  Default padding 10 Å matches AMBER tutorial conventions and easyPARM.
  Disabled via the `Add Cubic Box if Non-Periodic` BooleanParameter
  for users who plan to solvate via `gmx editconf` / `gmx solvate`
  themselves.

- **`tests/test_amber_to_gromacs_core.py`** — 6 tests:
  - missing-input raises (no AmberTools dependency)
  - round-trip on a tleap-built ALA peptide (default box added,
    non-periodic input)
  - `add_box_if_absent=False` keeps non-periodic
  - `box_padding=15` produces a box ~10 Å larger than `box_padding=10`
    (verified via `.gro` last-line box-vector parse)
  - re-running on same prefix overwrites cleanly

- **`tests/test_pipeline_e2e.py::TestAmberToGromacsE2E`** — 1 E2E
  test: tleap-built ACE-ALA-NME → ParmEd → assert the produced `.top`
  has the standard GROMACS sections (`[ atomtypes ]`,
  `[ moleculetype ]`, `[ atoms ]`, `[ bonds ]`) and the `[ atoms ]`
  block atom count matches the prmtop's.

### Changed

- **`package.toml`**: `nodes` list extended (11 → 12);
  `description` updated; new hashtags `gromacs`, `parmed`,
  `format-conversion`.
- **`tests/test_node_packages.py:ALL_NODES`** and
  **`tests/test_pipeline_e2e.py:TestNodePackages.NODE_NAMES`** updated
  to include `ep_amber_to_gromacs` (structural validation passes).
- **CLAUDE.md** + **README.md** node count and pipeline diagram
  updated; node table gains a row for the converter.

### Verified

- Real SnP fuse output (`complex.prmtop` 219 KB, `complex.rst7` 16 KB,
  448 atoms) → `complex.top` 247 KB + `complex.gro` 20 KB. All
  expected GROMACS sections present. Custom types (`n2-Sn`, `Sn-os`,
  `ns`, etc.) and the v1.5.1 metal MASS+NONBON entries round-trip
  faithfully via ParmEd's `[ atomtypes ]` translation. Cubic box
  auto-added with default 10 Å padding.
- Test count: 366 → 376 (+10: 6 unit + 1 E2E + 3 from extending the
  TestNodePackages NODE_NAMES list).

### Why a separate node

Three integration paths were considered (CHANGELOG-internal note —
documented in the dev discussion that prompted this work):

1. *Run ParmEd ad-hoc as a one-liner snippet in the README.* Lowest
   effort but invisible in the visual graph.
2. *Inline as a final step in `ep_fragment_fuse` gated by a flag.*
   Keeps node count at 11 but couples the format conversion to the
   fuse step (less reusable for third-party prmtops).
3. **Separate node (this release).** Surfaces the conversion as an
   explicit pipeline step; users can also feed in third-party
   `complex.prmtop` files (e.g. from MCPB.py or MetalDock outputs)
   without going through fuse. Matches easyPARM's "separate
   converter script" pattern.

## [1.9.1] — 2026-05-05

**fragment_align primary-rotation sign bug — discovered via VMD
visualization of v1.9.0's complex.pdb.** The user loaded the produced
prmtop in VMD, which flagged "Unusual bond between residues 6/10/14
and 24" — atoms within 0.8-1.6 Å, indicating severe steric overlap
between the SnP porphyrin core and peptide LEU residues 10 and 14.

### Root cause

`compute_rigid_transformation()` calls
`rotation_aligning(frag_outward, -pep_outward)`, which returns the
unique minimum-angle rotation. But for sp2/sp3 anchors, the
`compute_outward_direction()` sign is arbitrary (cross-product order
of the bonded atoms), so the primary rotation can place the bulk of
the fragment on the **wrong side** of the bond axis.

For the SnP demo: snpp.pdb's YASARA-designed input has SnP correctly
oriented (NH2 already 1.51 Å from GLU.CD, porphyrin facing away from
the helix). The primary rotation flipped the porphyrin 150° around
the bond axis, sending it INTO the helix — producing 0.8-1.2 Å
overlaps with LEU.10 (HD23/HD11/HG/CG) and LEU.14 (HD23/HD12/HD22).

### Fix

`compute_rigid_transformation()` now tries both candidate rotations —
the primary one and a 180°-flip around the bond axis — and picks
whichever has fewer peptide-fragment **clash pairs** (atom pairs at
< 2.0 Å). The clash count *excludes* atoms within 2.0 Å of the
peptide anchor (which include the new amide-bond partner and existing
GLU side-chain atoms — chemically expected close contacts). This
metric is robust where simple min-distance fails: both candidates
have the new amide bond at ~1.3 Å in their min-distance, so a min-
distance metric can't distinguish them; counting clashes elsewhere
on the fragment surface does. The existing `find_clash_free_rotation`
secondary scan still runs after, providing fine-grained adjustment
around the bond axis. Net additional cost: one matrix multiply + one
O(|pep|×|frag|) distance scan, negligible.

### Verified

- 1 new regression test
  (`test_rigid_transformation_picks_clash_free_side`) constructs a
  synthetic peptide with a LEU bulk at +z and a fragment with phenyl
  bulk; asserts the alignment places the phenyl on the clash-free
  side (min-dist > 2 Å, was ~1 Å before fix).
- The full snp-fragment-fuse Playwright spec re-run still produces
  `complex.{prmtop,rst7,pdb}` cleanly. Specific clash distances
  before/after the fix:

  | Residue | Atom pair | Before  | After   |
  |---------|-----------|---------|---------|
  | res 10  | LEU.HD23 ↔ SnP.H8  | 0.809 Å ❌ | 4.014 Å ✓ |
  | res 14  | LEU.HD23 ↔ SnP.C15 | 1.248 Å ❌ | 5.979 Å ✓ |
  | res 6   | GLU.CD ↔ SnP.NH2   | 1.512 Å (the expected new amide bond) | 1.512 Å |

  VMD now reports only the expected res 6 → res 24 amide bond; the
  spurious res 10/14 → 24 contacts are gone.
- All 28 fragment_align unit tests pass; full suite 365 → 366.

### Why VMD's distance-based bond detection caught it

VMD's PSF builder uses 0.6 × sum-of-covalent-radii as bond cutoff
(~1.6-1.8 Å for typical pairs). Atoms at 0.8-1.5 Å register as
"unusual bonds" — i.e. severe overlap that no MD step would relax
gracefully. tleap's `saveamberparm` doesn't check geometry (it only
checks parameter coverage), so v1.9.0 produced a parametrically-valid
prmtop with broken structure. v1.9.1 closes that gap by improving
the geometric placement at fragment_align time.

## [1.9.0] — 2026-05-05

**Closes the cross-FF parameter gap at fuse-time amide attachments.**
The v1.8.1 spec run found that fusing a peptide (ff19SB types: `2C`,
`CO`, `O2`, `HC`) to a metal-fragment ligand (GAFF2 types: `ns`, `ca`,
`hn`) via a covalent amide bond fails `saveamberparm` with ~17 missing
BOND/ANGLE/DIHE/IMPROPER errors. The cross-types (`CO-ns`, `2C-CO-ns`,
…) aren't pre-tabulated in either FF, even though chemically they're
just standard amide bonds.

The published practice (Glycam, MCPB.py, parmed-based pipelines) is to
provide a small **bond-type-specific "linkage frcmod"** with
canonical AMBER amide values. v1.9.0 ships that linkage as a bundled
demo asset and wires it through `ep_fragment_fuse`.

### Added

- **`packages/metalparm-vwf/ep_fragment_fuse/demo_data/linkages/amide_glu_gaff2.frcmod`**
  — cross-FF amide linkage frcmod with values verified against parm10:

  | Term | Type pair | Value | Source |
  |------|-----------|-------|--------|
  | BOND | `CO-ns` | 490.0 / 1.335 | parm10 `C-N` AA general |
  | ANGLE | `2C-CO-ns` | 70.0 / 116.6 | parm10 `CT-C-N` |
  | ANGLE | `O2-CO-ns` | 80.0 / 122.9 | parm10 `N-C-O` |
  | ANGLE | `CO-ns-hn` | 50.0 / 120.0 | parm10 `C-N-H` (NMA) |
  | ANGLE | `CO-ns-ca` | 50.0 / 121.9 | parm10 `C-N-CT` |
  | DIHE  | 7 cross-FF torsions | 0.000 amplitude | parm10 `X-CT-N-X = 0.00` convention |
  | IMPRO | `O2-2C-CO-ns` | 10.5 / 180 / 2 | parm10 `X-X-C-O` (line 893) |

  The torsion barriers are zero by convention — the planar restraint is
  carried by the `X-X-C-O` improper, exactly mirroring how parm10 sets
  `X-CT-N-X = 0.00` and lets `X-X-C-O = 10.5` lock amide planarity.

- **`linkage_frcmod` parameter** on `EpFragmentFuse` (FileParameterEdit,
  optional) — points at a cross-FF linkage frcmod. Auto-discovers
  `output_linkage_frcmod` from a predecessor when set. Staged into
  the work directory so tleap finds it via basename.

- **`linkage_frcmods` argument** on `fuse_helpers.build_tleap_script()`
  — list of additional frcmods that get `loadamberparams`'d after
  `fragment.frcmod` (so linkage entries take precedence over any that
  fragment.frcmod might have defined).

- **2 new unit tests** in `tests/test_fragment_fuse_helpers.py`:
  - `test_script_with_linkage_frcmod` — asserts the linkage is loaded
    after fragment.frcmod via basename only.
  - `test_script_no_linkage_when_not_supplied` — asserts default path
    (no linkage) emits exactly one `loadamberparams` line.

### Changed

- **Bundled `snp.{lib,frcmod}`** in both
  `ep_fragment_fuse/demo_data/` and `fragment_align/demo_data/`
  regenerated with `--at-type gaff2` (was `gaff` in v1.6.0). The
  amide-N type is now `ns` (gaff2) instead of `n` (gaff), matching
  the linkage frcmod's type vocabulary. End-to-end verified: fresh
  tleap session producing `complex.prmtop` (~219 KB, 0 errors,
  6 benign warnings about non-integer charges from antechamber BCC
  fitting).
- **`operations/snp-fragment-fuse.spec.ts`** wires the linkage frcmod
  into the fuse node's `Linkage Frcmod` field. The previous v1.8.1
  spec failure (cross-FF torsion errors) is resolved.

### Test count

365 unit + E2E tests passing (363 → 365 with the 2 new linkage tests).

### Known limitations

- **Linkage frcmod is bond-type-specific.** The bundled
  `amide_glu_gaff2.frcmod` covers GLU/ASP-amide attachments to a
  GAFF2-typed amine. Other linkage types (His-coordination to Zn,
  Cys-disulfide to a cofactor, Pt-N4 square-planar) need separate
  files. The directory `ep_fragment_fuse/demo_data/linkages/` is
  the canonical home; expand as new cofactor cases land.
- **GAFF (not GAFF2) typing requires a parallel `amide_glu_gaff.frcmod`** —
  same body but with `n` instead of `ns`. Not bundled yet (gaff2 is
  the recommended path for ff19SB compatibility per the upstream
  AMBER documentation).

## [1.8.1] — 2026-05-05

**Bugfix found via the full SnP Playwright spec.** The v1.8.0
preprocessor's `drop_heteroatoms=True` filter was too aggressive — it
dropped any record with `HETATM` record type, regardless of resname.
YASARA-style PDBs use `HETATM` for **standard** residues too (the SnP
demo's `snpp.pdb` has `HETATM ... GLU A 6 ...`), so the GLU-residue at
the SnP-fuse interface bond got silently dropped during preprocessing,
then `fragment_align` failed with `"atom not found: resseq=6 atom=CD"`.

### Fixed

- **`peptide_pdb_preprocess`** now gates the heteroatom-drop on the
  **residue name** (against `STANDARD_RESIDUES`), not the **record
  type**. Drops only records whose resname isn't in the ff19SB template
  set. So `HETATM ... GLU A 6 ...` is kept (standard resname) while
  `HETATM ... HOH A 100 ...` is still dropped (non-standard).
- **`tests/test_peptide_pdb_preprocess.py:test_preprocess_keeps_hetatm_with_standard_resname`** —
  regression test pinning the new behavior. Total preprocess tests:
  22 → 23.

### Verified end-to-end

Re-ran `operations/snp-fragment-fuse.spec.ts` after the fix:
- ✅ `snp_builder` extracts the SnP metal site + emits a peptide PDB.
- ✅ `peptide_builder` (PDB mode, v1.8.0 preprocessing) preserves all
  22 residues including GLU at resseq 6. Without the v1.8.1 fix this
  step silently dropped GLU.
- ✅ `fragment_align` finds GLU.CD as the peptide anchor and computes
  the rigid-body transform. Emits `<case>_aligned.{lib,pdb}`.
- ❌ `ep_fragment_fuse` `saveamberparm` fails with 17 missing-torsion
  errors at the GLU.CD-SnP.NH2 amide boundary: `HC-2C-CO-n`,
  `2C-CO-n-hn`, `2C-CO-n-ca`, `2C-2C-CO-n`, `O2-2C-CO-n` improper.
  This is a **cross-FF parameter gap**: peptide.lib uses ff19SB types
  (`2C`, `CO`, `O2`), fragment.lib uses GAFF types (`n`). Tleap has no
  fallback torsions for the cross-FF type pairs. Out of scope for
  v1.8.1 — see Known limitations.

### Known limitations (added)

- **Cross-FF torsion parameters at fuse interface bonds.** When fuse
  creates a bond between a peptide atom (ff19SB types: `2C`/`CO`/`O2`)
  and a fragment atom (GAFF types: `n`/`c2`), tleap needs torsion
  parameters spanning the cross-FF boundary. Neither force field
  defines these out of the box. Workarounds:
  - **Hand-patch `peptide.frcmod`** with the boundary torsions —
    add lines like `2C-CO-n-hn 0.0 0.0 0.0` to silence tleap (works
    but the dynamics through that bond are then unrestrained, which
    is usually fine for a covalent attachment).
  - **Build the peptide with GAFF types** instead of ff19SB — change
    `peptide_builder` to use a GAFF `loadpdb` path. Loses ff19SB's
    backbone parameters but unifies the type set.
  - **Use ff19SB types throughout** — re-parameterize the metal
    fragment with ff19SB types instead of GAFF. Requires antechamber
    `-at amber` and corresponding ff19SB-style frcmod.
  Pending evaluation; for now the metalloprotein workflow is usable
  end-to-end through the lib-generation step but the fuse-step
  `saveamberparm` is a known-broken artifact for hybrid-FF complexes.

## [1.8.0] — 2026-05-05

**Phase 2 — peptide_builder PDB-mode handles AlphaFold / ProteinMPNN /
experimental PDBs without external preprocessing.** Adds three new
options + a single PDB-cleanup pass that fronts the existing tleap
load.

### Added

- **`peptide_builder/core.py:peptide_pdb_preprocess()`** — single-pass
  PDB cleanup that runs in PDB mode before tleap. Six steps in order:
  1. **MODEL 1 selection** — multi-model NMR / AlphaFold-multimer /
     ensemble PDBs are common; tleap doesn't handle them. Take first.
  2. **Chain filter** (e.g. `"A"` or `"A,B"`; empty / `"all"` keeps
     everything).
  3. **Residue range** (e.g. `"5-30"`; empty / `"all"` keeps
     everything). Useful for trimming AlphaFold predictions to a
     binding loop.
  4. **Heteroatom drop** (default ON). Drops `HETATM` records and
     `ATOM` records with non-standard residue names — waters, ligands,
     post-translationally modified residues that ff19SB doesn't have
     templates for. The standard residue set covers the 20 AAs +
     ff19SB protonation / tautomer variants (HID/HIE/HIP, ASH/GLH,
     LYN, CYM/CYX) + caps (ACE/NME/NHE) + N-/C-terminal forms.
  5. **YASARA atom rename** — `GLU.COOH → HE2`, `SER.HO → HG`, etc.
     Mirror of the table in `snp_builder/metal_fragment.py` (per the
     v1.7.0 cross-node-import policy: each node bundles what it
     needs, kept in sync).
  6. **HIS tautomer inference** — rewrite `HIS` residue names to
     `HID` (HD1 only) / `HIE` (HE2 only — leave as `HIS`, ff19SB
     default) / `HIP` (both) based on which protonation hydrogens
     are present.

  Returns a stats dict (`models_dropped`, `chains_dropped`,
  `residues_out_of_range`, `het_dropped`, `his_renamed`, `residues`,
  `atoms`) for logging / testing. Raises `ValueError` on unparseable
  range / chain spec or missing input.

- **3 new BoCoFlow node OPTIONS** on `PeptideBuilder`:
  - `chain_filter: StringParameter` (default `""` = all)
  - `residue_range: StringParameter` (default `""` = all)
  - `drop_heteroatoms: BooleanParameter` (default `True`)

  PDB mode now stages the user-supplied PDB through
  `peptide_pdb_preprocess()` before the tleap loadpdb. Sequence mode is
  unchanged.

- **`peptide_builder/demo_data/alphafold_peptide.pdb`** — synthesized
  AlphaFold-style demo derived from the existing `peptide.pdb`. Wraps
  the 9-residue peptide in a `MODEL 1` block, adds a 4-atom chain B
  stub, a 3-atom HOH water, and a `MODEL 2` stub atom. Exercises every
  preprocessing step.

- **`tests/test_peptide_pdb_preprocess.py`** — 22 unit tests covering:
  - Range / chain spec parsers (3 + 3, incl. empty / "all" /
    invalid / reversed-range raises).
  - YASARA atom rename (1).
  - HIS tautomer inference (3 covering HD1-only / HE2-only / both).
  - Full preprocessor pipeline (12, covering: default flags on the
    bundled demo, chain-only filter, residue-range trim,
    drop-heteroatoms-off keeps water, non-standard residue drop,
    YASARA rename + disable, HIS rename, combined filters,
    missing-input raise, invalid-range raise).
- **`tests/test_pipeline_e2e.py::TestPeptideBuilderPdbPreprocess::test_alphafold_demo_passes_through_tleap`** —
  drives the bundled `alphafold_peptide.pdb` through the
  preprocessor + tleap; asserts a clean `peptide.pdb` and `peptide.lib`
  are produced and the HIS-with-HE2 has been auto-resolved to `HIE`.

### Changed

- **No behavioral change** to sequence-mode peptide builds. PDB mode's
  output for already-clean inputs (e.g. tleap-emitted `peptide.pdb` from
  a prior run) is byte-equivalent except for an explicit `END` record
  at the bottom.

### Test count

361 → 362 unit + E2E tests passing (361 unit + 23 new in
`test_peptide_pdb_preprocess.py` and `test_pipeline_e2e.py`).
Effective totals: 362 with the new `TestPeptideBuilderPdbPreprocess`
class.

### Known limitations (unchanged)

- AlphaFold / ProteinMPNN peptides with **non-standard residues**
  (phosphoryl, methyl, hydroxyproline, …) are dropped by the
  heteroatom filter. ff19SB doesn't have templates for these; users
  with PTMs need to either pre-resolve the residues to standard +
  custom ff fragments or build via a sequence with `XAA`-style
  custom-residue codes. Out of scope for v1.8.0.
- **`metal_fragment_builder` refactor** (Phase 3b — resolved in v1.7.0
  with the metal_fragment.py shared library).

## [1.7.1] — 2026-05-05

**Vanilla easyPARM force-constant scaling now applied by default.**
Closes the last documented divergence from upstream easyPARM:
01_easyPARM.sh lines 1269-1299 has an awk block that boosts weak
BOND/ANGLE force constants by empirical multipliers (Seminario
underestimates dative metal-ligand bonds). Our pipeline previously
emitted raw Seminario constants — viable for Sn / Ru / strong-bond
cases but suboptimal for weakly-bound metals (Zn²⁺, Cu²⁺ centers in
the published validation set).

### Added

- **`packages/metalparm-vwf/ep_forcefield_assembly/scripts/apply_fc_scaling.py`** —
  pure-Python port of the awk scaling block. Reads a frcmod, scales
  BOND and ANGLE force constants in place using:
  ```
  ANGLE k:    k<5 → ×11.599   k<10 → ×7.799
              k<20 → ×3.599   k<29 → ×2.699   k≥29 → unchanged
  BOND k:     k<20 → ×4.599   k≥20 → unchanged
  ```
  Strict-less-than boundaries match awk's `<` operator. Idempotent for
  k ≥ cutoffs (the common case after Seminario averaging on strong
  bonds, e.g. SnP's na-Sn k≈115 stays 115). DIHE / IMPROPER / NONBON
  sections are untouched.
- **`src/steps/step_fc_scaling.py`** — shallow Python wrapper
  (subprocess) for the script.
- **`scale_fc` toggle** plumbed through:
  - `BooleanParameter("Apply easyPARM Force-Constant Scaling",
    default=True)` on `EpForcefieldAssembly` node — runs after
    `metal_nonbon_fill` and before tleap.
  - `--scale-fc` / `--no-scale-fc` mutually-exclusive CLI flags on
    `run_pipeline.py`, `run_pipeline_metalloprotein.py`, and
    `run_pipeline_metalloprotein_gaussian.py` (default ON).
  - `scale_fc=True` parameter on `run_pipeline()`,
    `step_fc_scaling.run()`.
- **`tests/test_apply_fc_scaling.py`** — 16 unit tests covering all
  cutoffs + multipliers + boundary conditions (strict-less-than),
  column preservation, section-header / blank-line / DIHE / IMPROPER
  / NONBON immunity, full-file pass with mixed-magnitude entries,
  idempotency for above-cutoff values.

### Verified

- **SnP regen** (default `--scale-fc`): bundled
  `snp.{lib,frcmod}` are byte-equal to v1.6.0 — all SnP constants
  (na-Sn=115, Sn-os=168, c2-na-Sn=78, na-Sn-na=49, na-Sn-os=31,
  Sn-os-c3=47, etc.) are above the scaling cutoffs, so scaling is a
  no-op. Fresh-tleap saveamberparm: 0 errors, 0 warnings, prmtop OK.
- **All 339 tests pass** (323 + 16 new). E2E pipeline run on the SnP
  fixtures (`collect/snp_ff_v1/orca/`) produces a clean
  `COMPLEX.{mol2,frcmod,lib}`.

### Affects which cases?

After Seminario averaging:
- **Strong covalent bonds** (k_bond > 200, k_angle > 60) → unchanged.
- **Sn / Ru / similarly-bound metals**, post-averaging — typically
  unchanged (SnP is the worked example).
- **Zn²⁺, Cu²⁺, weak-dative-bond cases** — scaling kicks in. If you
  were comparing our output to vanilla easyPARM previously and saw
  k_bond ~10× smaller, that's the gap this closes.

To replicate v1.7.0 behavior (raw Seminario): pass `--no-scale-fc` on
the CLI or untick the assembly node's "Apply easyPARM Force-Constant
Scaling" checkbox.

## [1.7.0] — 2026-05-05

**Phase 3b refactor — `snp_builder` split into a thin SnP-specific
orchestrator + a generic `metal_fragment` shared library.** Lays the
groundwork for sibling builders (Zn-finger, heme, Pt-N4, …) without
forking porphyrin-specific code.

### Added

- **`packages/metalparm-vwf/snp_builder/metal_fragment.py`** — generic
  shared library (~330 lines). The canonical home for:
  - **Data records**: `AtomRec`, `BuildResult` with `.to_xyz()` /
    `.to_pdb()` writers.
  - **Geometry**: `ring_plane` (planar SVD), `perpendicular_unit`,
    `tetrahedral_h` (methyl-H placement), `pdb_atom_name` (column
    formatting helper).
  - **Residue extraction**: `extract_residue_atoms(pdb, resname)` —
    generalizes the previous "find UNK in PDB" pass.
  - **Metal swap**: `swap_metal_at_centroid(atoms, centroid, in, out)`
    — drop placeholder, insert new metal at centroid; preserves the
    "PDB-name uppercase, element-symbol mixed-case" convention.
  - **Axial ligand placement**: `place_axial_ligand(centroid, normal,
    ligand, bond_len, sign, tag)` — supports `OMe`, `OH`, `Cl`,
    `none`. Per-axial tagging so multiple axials don't name-collide.
  - **Peptide cleanups**: `extract_peptide`, `infer_his_tautomer_renames`,
    `YASARA_ATOM_RENAMES`, `_rename_atom_in_pdb_line` — all moved
    here unchanged from the previous `core.py`.
  - **File writers**: `write_outputs(result, out_dir, basename, resname)`.
- **`tests/test_metal_fragment.py`** — 22 unit tests covering the
  generic helpers in isolation: PDB column formatting (3 tests),
  AtomRec data class (1), ring_plane SVD (2), perpendicular_unit (1),
  tetrahedral_h (2), extract_residue_atoms (3 incl. missing-resname
  raise), swap_metal_at_centroid (2), place_axial_ligand (5 covering
  OMe/OH/Cl/none/unknown), YASARA renames (3).
- **`dev-notes/metal-fragment-builders.md`** — guide for adding sibling
  builders. Includes:
  - What's shared in `metal_fragment.py` vs what's porphyrin-specific
    in `core.py`.
  - Sketches for `zn_finger_builder` (no axials, tetrahedral
    Cys/His coordination) and `heme_builder` (porphyrin-like with FE
    + axial His from peptide).
  - 8-step recipe for shipping a new builder (copy
    `metal_fragment.py`, write `core.py` orchestrator, mirror
    `node.py`, demo data, meta.toml, package.toml entry, tests,
    interface_bonds template).

### Changed

- **`snp_builder/core.py`** is now ~270 lines (was 511) and contains
  only the porphyrin-specific orchestrator: `PYRROLE_N_NAMES`,
  `SN_O_AXIAL`, `build_snp_fragment`, `_apply_cap` (the aniline-NH₂
  amide cap, keyed on the `NH2` atom name + porphyrin ring normal
  fallback). Generic helpers are imported from
  `metal_fragment` via a 3-tier fallback that mirrors the
  `node.py` import pattern (package-relative → direct-path →
  `importlib.util` for the test suite's standalone load mode).
- **Back-compat re-exports**: `core.py` still exposes
  `AtomRec`, `BuildResult`, `extract_peptide`, `infer_his_tautomer_renames`,
  `write_outputs`, `_pdb_atom_name`, `_ring_plane`, all chemistry
  constants (`O_C_METHOXY`, `C_H_METHYL`, etc.) and the
  `YASARA_ATOM_RENAMES` table — anyone still importing them off
  `snp_builder.core` keeps working. `write_outputs` is a thin SnP-
  flavoured wrapper (default `basename="snp_frag"`, `resname="SNP"`)
  around the generic version.
- **`tests/test_snp_builder.py`** unchanged — all 24 tests still pass.
  The test loads `core.py` via `importlib.util.spec_from_file_location`,
  which exercises the standalone-load fallback in `core.py`'s import
  block (registering the dynamically-loaded `metal_fragment` module
  in `sys.modules` before `exec_module` so the dataclass decorator's
  `cls.__module__` lookup works).

### Test count

323 unit + E2E tests passing (was 301; +22 new in
`test_metal_fragment.py`).

### What's NOT in this release

- **No new builder nodes**. `zn_finger_builder` / `heme_builder` /
  `pt_n4_builder` are still gated on real workflow triggers per the
  v1.4.2 strategy doc — when one of those cases lands, the recipe in
  `dev-notes/metal-fragment-builders.md` documents the path. The
  shared module is the unblocker, not the deliverable.
- **No marketplace-install changes**. The shared `metal_fragment.py`
  ships inside `snp_builder/`. Sibling builders that want to reuse it
  should copy the file (~330 lines, single dependency on
  numpy + Bio.PDB) rather than rely on cross-node imports — the
  bocoflow node-install process doesn't yet handle shared libraries
  cleanly. This is a known temporary cost; revisit when (a) two
  sibling builders coexist and the duplication is felt, or (b) the
  marketplace install supports shared per-package libs.

## [1.6.0] — 2026-05-05

**Atom names from `snp_builder` (and any fragment builder that emits a
PDB) now survive QM optimization → mol2 generation → final lib.** The
"two `interface_bonds` JSON variants" footgun is gone — `fragment_align`
and `ep_fragment_fuse` consume the same JSON.

### Root cause

`xyz_to_pdb.py` (the first step of `ep_mol2_generation` and
`step_antechamber`) rebuilds the PDB from a post-QM XYZ. XYZ files
contain only element symbols + coordinates, so the script generated
atom names as `<element><1-based-counter>` — losing the meaningful
names (`NH2`, `CAP`, `OAP`, `CM`, `HM1-3`) that snp_builder writes into
its emitted PDB. The earlier suspect (antechamber's `-j 5` flag) turned
out to be a red herring: antechamber preserves names that start with a
valid element symbol.

### Added

- **`--template <pdb>` flag on the bundled
  `packages/metalparm-vwf/ep_mol2_generation/scripts/xyz_to_pdb.py`**.
  When supplied, atom names from the template overlay onto the XYZ
  coords (atom-by-atom, in order). Validates atom count + per-position
  element match; falls back to legacy `<element><counter>` naming with
  a stderr warning on mismatch.
- **`pdb_template` parameter on `EpMol2Generation`** node — optional
  `FileParameterEdit`. Auto-discovers from a predecessor's `output_pdb`
  (which `snp_builder` already exposes), so wiring `snp_builder →
  ep_mol2_generation` requires no extra user action.
- **`--pdb-template` CLI flag on `run_pipeline.py`,
  `run_pipeline_metalloprotein.py`, and
  `run_pipeline_metalloprotein_gaussian.py`**, plumbed via
  `step_antechamber.run(..., pdb_template=...)` and
  `step_antechamber_gaussian.run(..., pdb_template=...)`.
- **`tests/test_xyz_to_pdb_template.py`** — 7 unit tests for the
  template-overlay logic: legacy fallback, successful preservation,
  count mismatch (warns + falls back), element mismatch (warns + falls
  back), 4-char atom-name column alignment, PDB element-column
  fallback to atom-name prefix, end-to-end CLI invocation.
- **`tests/test_pipeline_e2e.py::test_atom_names_preserved_through_antechamber`** —
  drives `step_antechamber.run()` with a template PDB through real
  antechamber; asserts the template's atom names appear in
  `COMPLEX.mol2`.

### Changed

- **`step_antechamber.py` and `step_antechamber_gaussian.py`** now
  prefer the bundled `xyz_to_pdb.py` (which supports `--template`)
  over the unchanged copy in the easyPARM submodule — same pattern
  used by v1.5.0 (`05_prepare_mol2_frcmod_passthrough.py`) and v1.5.1
  (`metal_nonbon_fill.py`). The submodule itself is read-only.

### Constraint (worth flagging)

Antechamber re-normalizes atom names whose leading character isn't a
valid element symbol (e.g., it would rewrite `Z..` to `C..` for a
carbon atom). For real-world fragment builders this is a non-issue —
snp_builder's names (`NH2`, `CAP`, `CM`, `HM1-3`) all start with
valid element symbols and are preserved verbatim.

### Known limitations (no longer including atom-name preservation)

- **Force-constant scaling** (vanilla easyPARM CLI lines 1269-1299) —
  resolved in v1.7.1 (now applied by default).
- **AlphaFold/ProteinMPNN peptide handling** (Phase 2 — wait for trigger).
- **`metal_fragment_builder` refactor** (Phase 3b — resolved in v1.7.0
  with the metal_fragment.py shared library).

### Test count

276 unit tests (269 + 7 new). 25 pipeline E2E tests (24 + 1 new). All
passing.

## [1.5.1] — 2026-05-05

**Closes the metal-element vdW gap.** A fresh tleap session loading the
produced `lib + frcmod` now succeeds with `saveamberparm` — no more
`could not find vdW (or other) parameters for type (Sn|Ru|Zn|...)`. This
was the last known blocker for the SnP demo end-to-end.

### Added

- **`metal_nonbon_fill.py`** (bundled in
  `packages/metalparm-vwf/ep_forcefield_assembly/scripts/`) — appends
  MASS + NONBON entries to the final frcmod for any metal element type
  (atomic number > 10) present in the mol2 but missing from the frcmod.
  - **MASS** entries (`<type>  <atomic_mass>`) come from a built-in
    table of standard atomic weights (IUPAC 2021).
  - **NONBON** entries come from `uff_data.txt` via the existing
    UFF-coordination lookup logic. Format: 2-space indent + R_min/2 + ε,
    matching `gaff.dat` and `frcmod.tip3p` convention.
  - Inserts NONBON entries **immediately after the section header** —
    blank lines terminate the section in tleap's parser, so trailing
    appends were silently dropped (silent only because tleap warned
    "Unknown keyword" but parsed on; saveamberparm then failed).
- **`src/steps/step_metal_vdw.py`** — shallow Python wrapper
  (subprocess) for the script.
- **Wiring**: the new step runs in all three CLI pipelines
  (`run_pipeline.py`, `run_pipeline_metalloprotein.py`,
  `run_pipeline_metalloprotein_gaussian.py`) and in
  `ep_forcefield_assembly/node.py`, immediately after the COMPLEX.frcmod
  rename and before tleap is invoked.
- **`tests/test_metal_nonbon_fill.py`** — 13 unit tests covering: MASS +
  NONBON insertion, dedup against existing entries, multi-metal,
  organic-atom skip, missing-input guards, the 2-space indent format
  (regression), and post-header insertion (regression).
- **`tests/test_pipeline_e2e.py::test_fresh_tleap_saveamberparm_succeeds`** —
  acceptance test: full Ru pipeline → fresh tleap session →
  `saveamberparm` produces a non-empty prmtop with no
  "could not find" errors. Locks in the v1.5.x behavior end-to-end.

### Verified

- `pixi run -e pipeline pipeline ... --no-uls` on
  `examples/Ru_orca/` produces a `COMPLEX.frcmod` with
  ```
  MASS
  Ru     101.070

  ...

  NONBON
    Ru          1.4815  0.0560
  ```
  Fresh tleap session loads + `saveamberparm` succeeds (exit 0,
  no errors, prmtop ~34kB).
- 269 unit tests pass (256 + 13 new). 24 pipeline E2E tests pass
  (23 + 1 new acceptance test).

### Known limitations (no longer including the metal vdW gap)

- **Force-constant scaling** (vanilla easyPARM CLI lines 1269–1299) is
  still not applied. Pre-existing, not a regression.
- **Atom-name preservation through `ep_mol2_generation`** — root cause
  identified (`xyz_to_pdb.py` losing PDB names) and resolved in
  v1.6.0; see that section for details.
- **MASS values** in v1.5.1 use a built-in IUPAC 2021 table. UFF defines
  metal masses too but with slightly different precision; the
  difference (≤ 0.1 amu) is negligible for force-field MD but worth
  noting for high-precision applications.

## [1.5.0] — 2026-05-05

**ULS becomes opt-in.** Step 05 of the easyPARM pipeline used to rename
metal-coordinating atom types to unique labels (`na → n1/n2/n3/n4`,
`os → o1/o2`) so multiple distinct metal complexes wouldn't collide. For
single-metal-cofactor cases (the entire metalparm-vwf roadmap — SnP,
Zn-finger, heme, Pt-N4, …) ULS isolation is unused, and the rename created
an atom-type vocabulary divide between the lib and the frcmod that v1.4.3
patched at the lib-generation step. v1.5.0 makes ULS opt-in: by default,
step 05 is now a pass-through and the lib + frcmod use GAFF types
throughout, matching by construction.

### Added

- **`05_prepare_mol2_frcmod_passthrough.py`** — pass-through variant of
  step 05. Bundled at
  `packages/metalparm-vwf/ep_forcefield_assembly/scripts/` and mirrored
  into `collect/easyPARM/scripts/`. Copies `COMPLEX_modified.mol2` to
  `NEW_COMPLEX.mol2` unchanged; writes identity-mapping
  `new_atomtype.dat` and `metalloprotein_atomtype.dat`; writes an empty
  `addAtomTypes { }` block in `Hybridization_Info.dat`. Steps 06–13
  operate on GAFF types throughout — no downstream changes needed.
- **`use_uls` toggle**:
  - `BooleanParameter("Use Unique Labeling Strategy (ULS)", default=False)`
    on the `ep_forcefield_assembly` node.
  - `--use-uls` / `--no-uls` mutually exclusive flags on
    `src/steps/run_pipeline.py` (default `--no-uls`).
  - `use_uls=False` parameter on `src/steps/step_05.run()`.
- **`tests/test_step_05_passthrough.py`** — 8 unit tests for the
  pass-through script: byte-equal mol2 copy, identity metadata files,
  empty `addAtomTypes` block, missing-input handling, mol2 parser, metal-
  bonded atom finder.

### Changed

- **`step_tleap.build_tleap_input()` and
  `ep_library_generation/node.py`** — only inline the `addAtomTypes`
  block from `Hybridization_Info.dat` when it has actual entries
  (between-braces non-empty). An empty block makes tleap fail with
  *"Argument #1 is of type ?? Unknown type ?? must be of type: [list]"*.
  In pass-through mode the block is empty by design, so the tleap input
  no longer sources it. Both consumers prefer `NEW_COMPLEX.mol2` if it
  exists (it equals `COMPLEX.mol2` in pass-through mode, ULS-renamed in
  ULS mode) — convergent fall-through.
- **Test count**: 247 → 256 (+9 — passthrough tests + 1 new step_tleap
  test for the empty-addAtomTypes case).

### Verified end-to-end

- **No-ULS pipeline on Ru example**:
  ```
  pixi run -e pipeline pipeline --xyz examples/Ru_orca/OPTIMIZED.xyz \
    --hess examples/Ru_orca/freq_chelpg.hess --log examples/Ru_orca/freq_chelpg.out \
    --charge 2 --mult 1 --at-type gaff --no-uls --work-dir /tmp/ru_no_uls
  ```
  Produces `COMPLEX.lib` with `"N1" "n2"` for all coordinating N's and
  `COMPLEX.frcmod` with `n2-Ru` BOND / `n2-Ru-n2` ANGLE — vocabularies
  match by construction.
- **Fresh-tleap saveamberparm**: `loadoff COMPLEX.lib` + `loadamberparams
  COMPLEX.frcmod` + `saveamberparm` succeeds with no atom-type vocabulary
  errors. The single remaining `"could not find vdW for type (Ru)"` error
  is the pre-existing metal-element-vdW gap (see Known limitations);
  identical effect in ULS mode.
- **Regression on ULS-on path** (`--use-uls`): same pipeline produces the
  legacy ULS-renamed artifacts (`n1-Ru, n3-Ru, …`); the artifact shape
  matches the v1.4.x output. ULS code path untouched.

### Known limitations

- **Metal-element vdW gap**: `parmchk2` does not emit NONBON entries for
  metals (`Ru`, `Sn`, `Zn`, `Fe`, `Pt`, …). easyPARM step 11 only patches
  *existing* NONBON entries via UFF; it never adds new ones. So a fresh
  tleap session loading just `lib + frcmod` always fails with `"could
  not find vdW for type (<metal>)"`. Both ULS and no-ULS modes are
  affected equally — this is a pre-existing pipeline gap, independent of
  the ULS choice. Workaround: add metal NONBON entries to the frcmod by
  hand (or manually source UFF Lennard-Jones data). Future task: have
  step 11 (or a new step) emit a NONBON entry for each metal type.
- **Force-constant scaling (vanilla easyPARM CLI lines 1269–1299)** is
  *not* applied in `run_pipeline.py` or the BoCoFlow chain. Vanilla
  easyPARM scales BOND/ANGLE force constants by 4.6× / 11.6× / 7.8×
  /3.6× / 2.7× depending on magnitude after step 13. We emit unscaled
  Seminario constants. Not a regression from v1.4.3 — same behavior as
  earlier versions; flagged here since users comparing our output to
  vanilla easyPARM will see the difference.

## [1.4.3] — 2026-05-05

Phase 4 from the v1.4.2 "Known limitations" — fixes the lib/frcmod
atom-type vocabulary mismatch that broke `tleap saveamberparm` on the
SnP demo. The lib generation step now emits a lib whose atom types match
the frcmod's ULS labels (`n1-n4`, `o1-o2`, …), so the bundled
`snp.lib + snp.frcmod` can be loaded together for a complete prmtop.

### Fixed

- **`step_tleap.run()` and `ep_library_generation` now load
  `NEW_COMPLEX.mol2` + source `Hybridization_Info.dat`** when both are
  present in the work directory, instead of always loading
  `COMPLEX.mol2`. `NEW_COMPLEX.mol2` is the ULS-renamed mol2 produced by
  step 05 (`05_prepare_mol2_frcmod.py`); `Hybridization_Info.dat`
  contains the `addAtomTypes` block declaring the new types to tleap.
  Without sourcing it, tleap saved a lib with the antechamber-default
  GAFF types (`na`, `os`) which didn't match the frcmod's ULS labels
  (`n1-n4`, `o1-o2`), so `loadamberparams + saveamberparm` failed with
  `Could not find angle parameter for atom types: c3 - os - Sn`. Falls
  back to the legacy `COMPLEX.mol2` path when the ULS files are absent
  (non-metal cases).
- **`ep_forcefield_assembly` no longer overwrites `COMPLEX.mol2` with
  `NEW_COMPLEX.mol2`.** Both files now coexist in the work directory,
  with `output_mol2` pointing at the antechamber-typed `COMPLEX.mol2`
  and new `output_new_mol2` / `output_hybridization_info` keys pointing
  at the ULS-renamed mol2 and addAtomTypes block respectively. This
  matches the CLI `run_pipeline.py` behavior and makes the data flow
  explicit.
- **Bundled `snp.lib` regenerated with ULS-aligned atom types** in both
  `ep_fragment_fuse/demo_data/` and `fragment_align/demo_data/`. The
  bundled lib previously had GAFF defaults (`N2-N5: na`, `O1-O2: os`)
  even though the bundled frcmod was ULS-renamed (`n1-n4`, `o1-o2`).
  Patched mapping (deterministic from step 05's atom-id-order rename
  logic): `N2 → n1, N3 → n2, N4 → n3, N5 → n4, O1 → o1, O2 → o2`.
- **Bundled `snp.frcmod` gap-fill** added the four GAFF-equivalent
  bonded entries that were dropped during ULS rename:
  - `BOND: o1-c3 314.0 1.4170` and `o2-c3 314.0 1.4170`
  - `ANGLE: o1-c3-h1 50.840 109.470` and `o2-c3-h1 50.840 109.470`
  Values are GAFF2 defaults for `os-c3` / `os-c3-h1` (the parent types
  before the ULS rename). Without these, `saveamberparm` fails on the
  methoxy O-CH3 bonds and the methoxy O-C-H angles.

### Added

- **`tests/test_step_tleap.py`** — 5 tests covering the smart tleap-input
  builder: ULS-on path (uses `NEW_COMPLEX.mol2` + `addAtomTypes`), legacy
  path (uses `COMPLEX.mol2`, no `addAtomTypes`), partial-ULS fallback
  (only one of the two files present → legacy path), `at_type=amber`
  (`leaprc.ff19SB`), and unknown `at_type` defaulting to GAFF.

### Known limitations

- **The deeper pipeline gap-fill remains unfixed.** EasyParm step 11
  (`11_retrieve_uffdata.py`) only fills NONBON for ULS-renamed types
  via UFF; it does not propagate GAFF's bonded entries (`os-c3`,
  `os-c3-h1`, etc.) to the renamed types (`o1-c3`, `o2-c3`,
  `o1-c3-h1`, `o2-c3-h1`). For now, the bundled SnP demo frcmod has
  these four entries patched in by hand. Future runs through the full
  pipeline on a new metal cofactor will surface analogous gaps for
  that cofactor's local connectivity — fix at the pipeline level when
  this becomes a recurring problem.

## [1.4.2] — 2026-05-04

Phases 1b + 3a from the strategy doc. Phase 1a (regenerate `snp.frcmod`
with complete coverage) deferred — investigation revealed it requires a
deeper pipeline-level fix (lib/frcmod atom-type vocabulary mismatch + missing
GAFF entries for ULS-renamed types; covered in detail under "Known limitations"
below).

### Added

- **HIS tautomer rename in `snp_builder.extract_peptide`.** New helper
  `infer_his_tautomer_renames()` inspects each `HIS` residue and rewrites
  its name to `HID` / `HIE` / `HIP` based on which protonation hydrogens
  (`HD1`, `HE2`) are present. Solves the
  `Atom .R<NHIE 1>.A<HD1 20> does not have a type` failure when YASARA-
  extracted PDBs carry an `HD1` atom but tleap defaults to the HIE/NHIE
  template. New `rename_his_tautomers` flag (default True) toggles the
  behavior. Added 7 unit tests covering all four protonation cases plus
  the rename-disabled path. `infer_his_tautomer_renames()` only inspects
  `HIS`-named residues; `HID`/`HIE`/`HIP` are left alone.

- **General-metal support in `fragment_align/core.py`** (Phase 3a):
  - `STANDARD_BOND_LENGTHS` extended with Fe-N (2.0 Å), Fe-S (2.3),
    Fe-O (2.0), Cu-N/Cu-S/Cu-O, Ni-N/Ni-S, Mg-O/Mg-N, Ca-O, Mn-O/Mn-N,
    Co-N/Co-O, Pt-N (2.0)/Pt-S (2.3)/Pt-Cl, Ru-N/Ru-O.
  - `FRAGMENT_HYBRIDIZATION_TABLE` adds entries for the common metal
    anchor names (`("*","ZN")` → `metal_tetrahedral`, `("*","FE")` →
    `metal_axial`, `("*","CU")`, `("*","NI")`, `("*","MG")`, `("*","CA")`,
    `("*","MN")`, `("*","CO")`, `("*","PT")` → `metal_square_planar`,
    `("*","RU")` → `metal_octahedral`, `("*","SN")`).
  - `compute_outward_direction()` learns three new cases:
    `metal_tetrahedral` (4-coord; outward completes the tetrahedron),
    `metal_octahedral` (6-coord; outward to the 6th vertex),
    `metal_square_planar` (4-coord; in-plane outward with 3 ligands,
    plane normal with 4). 7 new unit tests.

- **`dev-notes/interface-bonds-templates.md`** — recipes for SnP-amide
  (current default), Zn-finger 2C2H, heme b axial His, disulfide-linked
  cofactor, Pt-N4 square-planar. Each recipe is ready-to-paste JSON for
  `fragment_align` + `ep_fragment_fuse` with chemistry notes and caveats.
  Linked from the package README's new "Cookbook" section.

### Known limitations

- **Bundled `snp.frcmod` parameter coverage** (was already a known
  limitation; clarified now): The lib and frcmod ship with mismatched
  atom-type vocabularies — the lib uses GAFF defaults (`na`, `os`)
  while the frcmod uses ULS-relabeled types (`n1`/`n2`/`n3`/`n4`,
  `o1`/`o2`). Even regenerating both via `src/steps/run_pipeline.py`
  produces the same mismatch (verified — md5 identical to the bundled
  artifact). Fixing requires either (a) re-saving the lib via tleap
  from `NEW_COMPLEX.mol2` (which has correct ULS types but exposes
  a second gap: missing entries for ULS-renamed-type interactions
  with non-metal atoms, e.g. `o2-c3-h1`); or (b) augmenting the
  pipeline's `13_final_clean.py` / `11_retrieve_uffdata.py` step to
  emit a frcmod with both metal-coordination params AND the GAFF-
  equivalent entries for every ULS-renamed type's interactions. This
  is Phase 4 deeper-pipeline work; the architectural pipeline succeeds
  through interface-bond creation, only the final saveamberparm step
  is blocked on the demo path.

## [1.4.1] — 2026-05-04

End-to-end Playwright run of `snp-fragment-fuse.spec.ts` exposed five
real bugs in the v1.4.0 architecture. Each fix is small and independent;
together they make the topology-coupled fuse pipeline run all the way
through the tleap interface-bond step.

### Fixed

- **Predecessor key collision (`peptide_lib` vs `fragment_lib`).** Both
  `peptide_builder` and `fragment_align` were forwarding the generic
  `output_lib`, so `ep_fragment_fuse._resolve` grabbed the same source
  for both peptide and fragment. Added role-specific keys —
  `output_peptide_lib`/`output_peptide_frcmod` from `peptide_builder`,
  `output_fragment_lib`/`output_fragment_frcmod` from `fragment_align`.
  `fragment_align` also passes through the peptide-side keys so the
  single-predecessor chain (peptide_builder → fragment_align → fuse)
  preserves both libs. `ep_fragment_fuse` prefers role-specific keys and
  falls back to generic `output_lib` for back-compat with consumers that
  don't yet emit role-specific names.

- **Module-introspection import failure (`No module named 'core'`).**
  `peptide_builder`, `fragment_align`, and `ep_fragment_fuse` used a
  2-tier `try/except` for their core/helpers imports. At server-side
  introspection (no `sys.path` tweak yet) BOTH tiers failed; the module
  didn't load, and `OPTIONS` came back empty — disabling the panel's
  Save button. Added a third-tier `name = None` fallback (matching the
  `snp_builder` pattern) with a clean `NodeException` at execute() if
  the helpers really aren't available.

- **`FileParameterEdit` options not optional.** `peptide_pdb`,
  `peptide_lib`, `fragment_lib`, `fragment_frcmod`, `fragment_pdb`
  defaulted to required; the UI disabled Save when they were empty,
  even though the predecessor-discovery path would fill them at execute
  time. Marked all predecessor-discoverable file inputs `optional=True`.

- **tleap unit-name collision + saveoff-append behavior.** `tleap`'s
  `saveoff varname filename` writes the unit under tleap's *internal*
  name (default "mol"), NOT the variable name. With both libs named
  "mol", the second `loadoff` clobbered the first.  AND `saveoff`
  *appends* to existing files rather than overwriting, so prior runs in
  the shared work_dir produced libs with two units.  Added
  `rename_lib_unit()` in `peptide_builder/core.py` to rewrite the
  peptide lib's unit name from "mol" to "pep" after tleap. All three
  nodes (peptide_builder, fragment_align, ep_fragment_fuse) now wipe
  stale output files in their work_dir before invoking tleap.

- **tleap `remove` syntax + YASARA atom names.** Two related fixes:

  - `remove cpx.6.OE2` failed — tleap requires `remove <unit> <atom>`
    with the *unit* as the first arg. `fuse_helpers.build_tleap_script`
    now emits `remove cpx cpx.<resid>.<atom>`.

  - The YASARA-extracted GLU residue uses `COOH` for the acid proton
    (not ff19SB's `HE2`), so `remove cpx cpx.6.HE2` failed too. Added
    `YASARA_ATOM_RENAMES` in `snp_builder/extract_peptide()` that
    rewrites `GLU.COOH → GLU.HE2` (plus a handful of other common
    YASARA hydroxyl/thiol naming quirks for SER/THR/TYR/CYS/ASP). New
    `rename_yasara_atoms` flag (default True) toggles the behavior.

### Known limitations

- **Atom-name preservation through `ep_mol2_generation`.** `snp_builder`
  writes atoms with chemical names (`NH2`, `CAP`, `OAP`, `CM`, `HM*`)
  but `ep_library_generation` runs `antechamber` which renames them
  generically (`N1`, `C47`, `O3`, `C48`, `H35-37`). The default
  `DEFAULT_INTERFACE_BONDS` uses snp_builder names — they work for
  `fragment_align` (which reads `snp_frag.pdb`) but not for fuse (which
  loads the antechamber-renamed lib). Workflows that target a
  `fragment_align`-output lib must override `interface_bonds` in fuse's
  panel with the actual antechamber names. The `snp-fragment-fuse`
  Playwright spec demonstrates this pattern.  Permanent fix is to
  preserve atom names through `ep_mol2_generation` — separate task.

- **Bundled `snp.frcmod` parameter coverage.** The fragment frcmod
  shipped in the `fragment_align`/`ep_fragment_fuse` `demo_data/` (from
  `collect/snp_ff_v1/`) lacks some angle parameters (e.g. `c3-os-Sn`)
  and atom-typing for HD1 on the N-terminal HIS, so end-to-end fuse via
  `saveamberparm` doesn't produce a valid prmtop on the demo path.
  Architectural pipeline succeeds through interface-bond creation;
  `saveamberparm` is the failure point.  Out of scope: regenerate
  `snp.frcmod` with complete coverage — separate task.

## [1.4.0] — 2026-05-04

### Added
- **`fragment_align` node** (11th node): pure-geometry rigid-body placement of the parameterized fragment onto the peptide anchor. Avogadro-style algorithm — outward direction from hybridization rules + standard bond lengths + Rodrigues rotation + optional clash-free secondary rotation, applied via tleap `transform` + `translate` + `saveoff`. The output `<case>_aligned.lib` carries pre-aligned coordinates so `ep_fragment_fuse` produces a chemically sensible interface bond geometry without any minimization. Bundled `demo_data/` ships a peptide.pdb + snp_frag.pdb (with SnP-style atom names) + snp.lib + snp.frcmod.
- `peptide_builder` now ALSO emits `peptide.lib` (via `saveoff pep peptide.lib`) and a `peptide.frcmod` placeholder. Both modes (sequence + pdb) produce a parameterized peptide topology for the fuse boundary. Forwarded data: `output_lib` + `output_frcmod`.
- `snp_builder` gains an `extract_peptide` option (default `true`) — emits the peptide chain (everything except the UNK residue) as `peptide_from_pdb.pdb`. Forwarded as `output_peptide_pdb`. peptide_builder PDB-mode auto-discovers it via the standard 3-tier predecessor pattern. Result: SnP-via-snpp.pdb path now produces consistent fragment + peptide coordinates from a single source.
- `tests/test_fragment_align_core.py`: unit tests for outward-direction computation (sp2 carbonyl bisector, sp3 tetrahedral completion), rigid-body transformation (orthonormal rotation, correct translation), clash-free rotation scan, and tleap-script syntax.

### Changed
- **`ep_fragment_fuse` is now strictly a topology fuser.** It loads two parameterized OFF libraries (`peptide.lib` + `fragment.lib`), combines them, applies cap removals + interface bonds, and saves `complex.prmtop / .rst7 / .pdb`. The peptide PDB is no longer an input — peptide topology is consumed exclusively as a `.lib`. Auto-detects unit names from each lib's index header.
- `fuse_helpers.build_tleap_script`: signature changed — replaced `peptide_pdb` with `peptide_lib` + `peptide_resname`; emits `loadoff peptide.lib` + `pep = copy <peptide_resname>` instead of `loadpdb`. Existing workflows that wired peptide_builder → ep_fragment_fuse must update — peptide_builder's output_lib now flows in instead of output_pdb.

### Architecture
- The package now spans **11 nodes**. The two metallopeptide paths cross the fuse boundary via parameterized topology on both sides:
  1. **SnP-from-PDB**: `snp_builder` (with `extract_peptide=true`) → peptide_builder (PDB mode) + parallel metal-FF chain → fragment_align (no-op since coords aligned) → ep_fragment_fuse
  2. **Fresh peptide**: peptide_builder (sequence mode) + metal-FF chain → fragment_align (computes rigid-body transform) → ep_fragment_fuse

## [1.3.0] — 2026-05-04

### Added
- **`peptide_builder` node** (10th node): builds a standalone peptide PDB from a sequence (via tleap, ff19SB or ff14SB, with optional ACE/NME caps) OR accepts a user-supplied PDB (mutations, multiple chains, pre-relaxed peptide). Forwards `output_pdb` + `forcefield` so `ep_fragment_fuse` consumes them via the standard 3-tier predecessor pattern. Bundled `demo_data/peptide.pdb` (9-residue ACE-HIS-TYR-TYR-LEU-ALA-GLU-ALA-NME built via ff19SB).
- `tests/test_peptide_builder_core.py`: unit tests for sequence normalization, residue counting, tleap-script emission, and user-PDB validation.

### Changed
- **BREAKING — `ep_fragment_fuse` is now a pure fuser.** Removed options: `peptide_mode`, `peptide_sequence`, `peptide_pdb` (in build mode), `n_term`, `c_term`. `peptide_pdb` is now an input only — it's auto-discovered from a `peptide_builder` predecessor (or set explicitly). Sequence-handling logic moved into `peptide_builder/core.py`. Existing workflows that drove fuse with `peptide_mode=sequence` must add a `peptide_builder` upstream.
- `fuse_helpers.build_tleap_script`: dropped `peptide_mode`, `peptide_sequence`, `n_term`, `c_term` arguments. Always uses `pep = loadpdb peptide.pdb`. `peptide_residue_count` moved to `peptide_builder/core.py`.
- `pixi.toml` (in `packages/metalparm-vwf/`): `[project]` → `[workspace]` (pixi 0.50+ rename).

## [1.2.0] — 2026-05-01

### Added
- xtb_opt: bundled `demo_data/OPTIMIZED.xyz` (Ru(bpy)₃ fixture, same as ep_bond_detection) so the node is testable standalone in the marketplace UI.
- ep_orca_run: bundled `demo_data/OPTIMIZED.xyz` + `demo_data/default-slurm.sh`.
- ep_fragment_fuse: bundled `demo_data/snp.lib` + `demo_data/snp.frcmod` (from `collect/snp_ff_v1/ff/`) so the fuse node has a self-contained demo path.
- `tests/test_node_packages.py`: structural smoke tests covering all 9 nodes — meta.toml validity, node.py compileability, bundled-scripts presence, demo_data presence.
- `operations/snp-fragment-fuse.spec.ts`: single-node Playwright operation symmetric with `snp-builder` and `snp-preopt`.

### Removed
- Empty `demo_data/` directories from `ep_forcefield_assembly` and `ep_library_generation` (downstream consumer nodes — inputs come from predecessors via auto-discovery, so a standalone demo_data dir was misleading).

### Fixed
- Node-count and path drift in `CLAUDE.md` and `README.md` (5/8 → 9 nodes; `metalparm_vwf/` paths → `packages/metalparm-vwf/`).
- `.gitignore`: ignore root-level pipeline scratch (`angle.dat`, `dihedral.dat`, `distance.dat`, `distance_type.dat`, `metal_number.dat`), `examples/*/test_output/`, ORCA `*.property.txt` and `*.chelpg.xyz`.

### Changed
- **Restructure**: package moved from repo-root `metalparm_vwf/` to `packages/metalparm-vwf/`. Aligns with the bocoflow-marketplace package layout. Operations moved from `metalparm_vwf/operations/` to repo-root `operations/`.
- Added lifecycle install spec (`operations/install-metalparm-vwf.spec.ts`) — idempotent end-to-end install verification: add metal-md as a marketplace source → sync → install metalparm-vwf → verify all 9 nodes register.

## [1.1.0] — 2026-04-23

### Added
- **xtb_opt** node: GFN2-xTB pre-optimizer bridging `snp_builder` and `ep_orca_run`. Cuts DFT convergence time substantially by seeding with a near-minimum geometry.
- **Auto-ECP** in `ep_orca_run`: leaving `ecp_block` empty triggers a `def2-ECP` block for atoms with Z ≥ 37 (Sn, Pb, Ru, Pt, …).
- `operations/` subdirectory with the first Playwright operation script (`snp-builder.spec.ts`).

## [1.0.0] — 2026-04-23

### Added
- **Consolidation**: previous `easyparm_vwf` (6 EasyParm nodes) and `metalpep_vwf` (2 metallopeptide nodes) merged into a single `metalparm-vwf` package. The shared environment is `metalparm_vwf` (AmberTools + xtb + numpy + biopython + redis-py).
- 8 nodes ship: `ep_bond_detection`, `ep_mol2_generation`, `ep_orca_run`, `ep_seminario_orca`, `ep_forcefield_assembly`, `ep_library_generation`, `snp_builder`, `ep_fragment_fuse`.

<!-- Version-compare/release links removed: they pointed at the private metal-md dev repo. -->

