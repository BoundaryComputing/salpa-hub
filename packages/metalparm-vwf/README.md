# metalparm-vwf

Visual workflow nodes for **metal force-field parameterization** in AMBER, plus **metallopeptide fragment-fusion** — one installable package bundling the EasyParm QM/Seminario pipeline and the tleap-based fragment builder + fuser.

## What it does

**Core EasyParm pipeline** — parameterizes metal-containing molecular systems (transition metal complexes, metalloproteins, metallonucleic acids) for classical MD in AMBER. Derives bonded parameters from quantum-mechanical (QM) Hessian data via the Modified Seminario method, fills gaps with UFF, and produces ready-to-use `.frcmod` and `.lib` files.

**Metallopeptide fusion** — for peptide + custom-cofactor systems (Sn-porphyrin labels, metalloenzymes, covalent inhibitors), running QM + Seminario on the *entire* complex every time the peptide sequence changes is wasteful (each DFT Hessian can take hours). The fragment+fuse pattern splits the work:

- **Parameterize the cofactor ONCE** → reusable AMBER library (`.lib` + `.frcmod`)
- **Build the peptide** with standard `ff19SB` (no QM)
- **Fuse** in tleap at user-specified interface bonds (amide, thioether, coordination)

Swapping peptide sequences then re-runs only the cheap fuse step.

## Pipelines

### Core EasyParm (metal complex → reusable FF)

```
XYZ geometry ─┬─> [ep_bond_detection]
              │        │
              │        └─> distance.dat, angle.dat, dihedral.dat, distance_type.dat, metal_number.dat
              │
              └─> [ep_mol2_generation] <── bonds/angles
                        │
                        └─> COMPLEX.mol2 (metal-aware atom typing)

Optimized XYZ ─> [ep_orca_run] ─> .out, .hess, optimized .xyz (local or SLURM)
                        │
                        ▼
[ep_seminario_orca] <── bonds/angles/dihedrals
                        │
                        └─> bond_angle_dihedral_data.dat (force constants)

mol2 + force_constants + bonds ─> [ep_forcefield_assembly]
                        │
                        └─> COMPLEX.frcmod, COMPLEX.mol2 (labeled, ULS)

mol2 + frcmod + distance_type ─> [ep_library_generation]
                        │
                        └─> COMPLEX.lib (AMBER library with correct atomic numbers + bonds)
```

### Metallopeptide fragment + fuse (for the SnP case) — recommended DAG

```
snp_builder ────► output_pdb (fragment, capped) ──► xtb_opt ──► [bond_det, mol2_gen, orca_run]
       │                                                                       │
       │ output_peptide_pdb (extracted peptide chain                            ▼
       │  at same YASARA coords)                                       ep_seminario_orca
       │                                                                       │
       └──► peptide_builder (PDB mode)                                          ▼
                  │                                                  ep_forcefield_assembly
                  │  peptide.{pdb,lib,frcmod}                                   │
                  │                                                             ▼
                  │                                                  ep_library_generation
                  │                                                             │
                  ▼                                                             │
            fragment_align ◄──── fragment_pdb + fragment.lib + fragment.frcmod ◄┘
                  │
                  │  <case>_aligned.lib (coords transformed)
                  ▼
            ep_fragment_fuse_topology  ──► complex.prmtop
                  │
                  ▼
            ep_apply_coords  ──► complex.{rst7, pdb}
                  │
                  │  (optional)
                  ▼
            ep_amber_to_gromacs ──► complex.{top, gro}  for GROMACS MD
```

Key design points:
- **Two crossings of the parameterized topology boundary**: peptide side and fragment side both ship `.lib + .frcmod`. `ep_fragment_fuse_topology` consumes parameterized topology on both sides (no implicit `loadpdb` of the peptide).
- `snp_builder`'s `extract_peptide=true` produces the peptide PDB at the same YASARA coordinates as the metal fragment — `peptide_builder` PDB-mode then auto-discovers it. **Coordinates are aligned by construction** along the SnP path; `fragment_align` becomes a near-no-op.
- For the **fresh-peptide case** (sequence-built peptide_builder, no aligned input PDB): `fragment_align` does pure-geometry rigid-body placement (Avogadro/GaussView style — hybridization rules + standard bond lengths + Rodrigues rotation) so the fuse produces a chemically sensible interface bond geometry without any minimization step.
- `xtb_opt` sits right after `snp_builder` so every downstream XYZ consumer sees a relaxed geometry.
- `ep_bond_detection`, `ep_mol2_generation`, and `ep_orca_run` run in **parallel** — all three take the relaxed XYZ but have no dependency on each other. `ep_orca_run` is the only expensive branch; the two topology branches finish in seconds.
- `ep_orca_run` **auto-generates the def2-ECP block** for atoms with Z ≥ 37 (Rb and heavier — e.g. Sn, Pb, Ru, Pt).
- `ep_seminario_orca` is the only node that blocks on the DFT Hessian; the force-field assembly step converges everything back together.

### Geometry placement (`fragment_align`)

For the fresh-peptide case, the peptide is built by tleap from a sequence in *whatever default coordinates tleap chooses*. The fragment library carries the QM-parameterized geometry. Without intervention, after `combine + bond` in tleap the interface bond can span an arbitrary distance — physically meaningless until a minimization step relaxes it.

`fragment_align` solves this with pure linear algebra (no FF, no QM):

1. Read peptide.pdb. Look up the peptide anchor's hybridization (e.g. `("GLU","CD") → "sp2_open"` after OE2/HE2 are removed). For sp2 with two existing bonds: outward = `-unit(unit(CD→CG) + unit(CD→OE1))` — the bisector pointing away from the existing bonds. For sp3: complete-the-tetrahedron formula.
2. Read fragment.pdb. Compute the same outward direction at the fragment anchor (after `frag_remove` atoms are dropped).
3. Rotation: align fragment outward anti-parallel to peptide outward, via Rodrigues' formula on the cross-product axis.
4. Translation: place fragment anchor at `pep_anchor + L · pep_outward` where `L` is from `STANDARD_BOND_LENGTHS` (1.33 Å for C–N amide, 2.05 Å for S–S, etc.) or auto-inferred from the anchor element pair.
5. Optional secondary rotation around the new bond axis (12-step scan) to maximize the minimum peptide↔fragment atom distance.
6. Apply via tleap: `transform mol <rotation>; translate mol <translation>; saveoff mol <case>_aligned.lib`.

The aligned lib feeds directly into `ep_fragment_fuse_topology`, which sees the fragment already at correct geometry and just needs to combine + bond; `ep_apply_coords` then transfers the aligned source-PDB coordinates onto the resulting prmtop.

## Nodes

### Core EasyParm

| Node | Display name | Purpose |
|------|--------------|---------|
| `ep_bond_detection` | Bond Detection | Detect bonds/angles/dihedrals from XYZ via covalent radii (metal-aware tolerance) |
| `ep_mol2_generation` | MOL2 Generation | antechamber + metal-aware atom type correction (GAFF/GAFF2/AMBER). Emits a **zero-charge** MOL2 — charges come from `ep_charges` |
| `ep_charges` | Charges | Inject QM-derived partial charges (ORCA 6 native `!RESP`, CHELPG, or classic `.vpot`+`resp`) into the MOL2. Sits between MOL2 Generation and Force Field Assembly; consumes the ORCA `.out` |
| `ep_orca_run` | ORCA Run | Run ORCA QM (opt + freq + CHELPG, optional `!RESP`) locally or via SLURM — produces the .hess + .out the Seminario and Charges nodes need |
| `ep_seminario_orca` | Seminario | Force constants from ORCA Hessian via Modified Seminario method |
| `ep_forcefield_assembly` | FF Assembly | Unique Labeling Strategy + Seminario merge + UFF gap-fill |
| `ep_library_generation` | Library Generation | tleap → `.lib` with metal atomic numbers and bonds fixed |

### Pre-optimization

| Node | Display name | Purpose |
|------|--------------|---------|
| `xtb_opt` | xTB Opt | Fast semi-empirical GFN2-xTB relaxation. Seeds the expensive DFT step with a near-minimum structure and gives `ep_bond_detection`/`ep_mol2_generation` cleaner distances to work with. |

### Metallopeptide fusion

| Node | Display name | Purpose |
|------|--------------|---------|
| `snp_builder` | SnP Builder | From a YASARA-style PDB (UNK = ZnPP), produce a capped Sn(IV)(OMe)₂-porphyrin fragment XYZ + PDB ready for the pipeline. Default: ACE cap on the aniline-NH side. With `extract_peptide=true`, also emits the peptide chain (everything except UNK) at the same coordinates so the downstream fuse starts with consistent geometry. `peptide_residue_range` (e.g. `1-7`, v1.17.0) carves a residue-number sub-span out of that chain — used to model a shorter peptide that is part of a longer one. `cap_peptide_termini` (v1.18.0) adds ACE/NME caps to the extracted chain (`Ac-…-NH-CH₃`). |
| `peptide_builder` | Peptide Builder | Build a standalone peptide topology — peptide.{pdb,lib,frcmod} — from a sequence (via tleap) under the chosen forcefield (ff19SB / ff14SB) with optional ACE/NME caps, OR load a user-supplied PDB (mutations, multiple chains, pre-relaxed peptide). The lib is the parameterized boundary `ep_fragment_fuse_topology` consumes. |
| `fragment_align` | Fragment Align | Pure-geometry rigid-body placement of the parameterized fragment onto the peptide anchor. Avogadro/GaussView style — hybridization rules + standard bond lengths + Rodrigues rotation + optional clash-free secondary rotation, applied via tleap `transform` + `translate` + `saveoff`. Output: aligned `<case>_aligned.lib` with pre-aligned coords. |
| `ep_fragment_fuse_topology` | Fuse Topology | Topology side of the fuser split: combine a peptide.lib (from `peptide_builder`) with a fragment.lib (raw from `ep_library_generation` or aligned from `fragment_align`), apply one or more user-specified interface bonds, and write `complex.prmtop` only. Rebalances the non-integer net charge left by interface-atom deletion so the complex total is an exact integer. |
| `ep_apply_coords` | Apply Coords | Coord side of the fuser split: takes the `complex.prmtop` from `ep_fragment_fuse_topology` plus the aligned peptide/fragment PDBs from `fragment_align`, and transfers those source-PDB coordinates onto the prmtop via ParmEd → final `complex.rst7` + `complex.pdb`. |
| `ep_amber_to_gromacs` | AMBER → GROMACS | Convert the produced `complex.prmtop`+`complex.rst7` to GROMACS `.top`+`.gro` via [ParmEd](https://github.com/ParmEd/ParmEd). All custom atom types, cross-FF linkage parameters, and metal MASS+NONBON entries round-trip faithfully. Auto-adds a cubic box for non-periodic input (matches easyPARM's `amber_converter.py` behavior). |

### MD preparation

| Node | Display name | Purpose |
|------|--------------|---------|
| `md_solvate_packmol` | MD: Solvate (packmol-memgen) | Wrap AmberTools' `packmol-memgen` — single command line for box + (mixed) solvent + ions + tleap-built AMBER topology. **Warning**: packmol-memgen's internal tleap re-types residues from the input PDB; for non-standard residues (fragment-fused GLU, etc.) it re-introduces auto-completed atoms that collide with the fragment. Use `md_solvate_gmx` instead for non-standard residue systems. |
| `md_solvate_gmx` | MD: Solvate (GROMACS-side) | Solvate a dry GROMACS topology (output of `ep_amber_to_gromacs`) via **raw packmol + ParmEd** — never re-derives the solute. Build solvent moleculetypes (OPC water from `solvents.lib`, Cieplak methanol from `MEOHBOX`, ions from `atomic_ions.lib`) at runtime via tleap, then assemble the final Structure entirely in ParmEd's `Structure + Structure` / `Structure * N` algebra. Safe for fragment-fused residues. |
| `md_traj_center` | MD Analysis: Trajectory Center | PBC-correct a trajectory before analysis: unwrap the solute (metallopeptide) so it is whole across periodic boundaries, centre it in the box, wrap solvent/membrane back in. MDAnalysis equivalent of `gmx trjconv -pbc whole -center -pbc mol`. Insert between an MD-run node and the analysis nodes — DSSP mis-assigns secondary structure on a solute split across a box face. Needs a `.tpr` (unwrapping uses the bond graph). |
| `md_analysis_helix` | MD Analysis: α-Helix Content | Per-frame DSSP secondary structure of the peptide (MDAnalysis). Reads `.xtc`/`.trr` + `.tpr`/`.gro`, writes `<case>_helix.csv` (helix fraction over time) + a per-residue helix propensity. Backbone-incomplete residues (ACE/NME caps) are auto-dropped. Run on a `md_traj_center` output for correct DSSP. |
| `md_analysis_distance` | MD Analysis: Residue–Metal Distance | Per-frame minimum-image distance from probe atoms (default: every Tyr hydroxyl O + His ring N — the PCET quenchers) to the metal centre (Sn), via MDAnalysis. Writes `<case>_distance.csv` (Å) + a closest-approach summary. |
| `md_membrane_build` | Membrane Build | Embed a dry metallopeptide GROMACS topology in a DPPC bilayer (transmembrane) via packmol-memgen + ParmEd. Pre-orients the helix axis along z (MEMEMBED can't orient a metallopeptide), packs the bilayer + water + ions geometry-only, tleap-parametrises just the standard membrane part, then ParmEd-concatenates the **preserved** solute topology — never re-derives the solute. Case 2's counterpart of `md_solvate_gmx`. v1.28.0 surfaces two packmol-memgen tuning knobs: `xy_box_A` (force `--distxy_fix` Å) and `nloop_all` (`--nloop_all` iterations) — set when the auto-sized XY box is too small for the lipid count and packmol's all-together loop won't converge. |

## Capping strategy (SnP case)

The aniline-N of the SnP dye forms an amide bond with the GLU sidechain Cδ carbonyl at runtime. During QM it must see the **amide** electronic environment, not a free aniline (-NH₂) or an N-methyl (-NHMe):

- **ACE cap** (default): prepend `CH3-C(=O)-` in front of the aniline N. Adds 6 atoms (CM + 3 HM + CAP + OAP). The methyl stands in for the GLU-Cβ during QM; electronics at the aniline-N match the real amide.
- **H cap**: keep the aniline as -NH₂. Simplest geometry; electronically wrong (amine donor, not amide). Only use if you want a quick pass without careful parametrization.
- **NHMe cap**: small methyl on the N side. Electronically between H-cap and ACE-cap. Uncommon.

At fuse time, `ep_fragment_fuse_topology` removes the cap atoms (default: `CM`, `HM1`, `HM2`, `HM3`, `CAP`, `OAP`) on the fragment side and the `OE2`/`HE2` on the GLU side, then creates `bond GLU.CD SNP.NH2` to close the amide.

## `interface_bonds` schema

`ep_fragment_fuse_topology` takes a list of interface-bond specs — one entry per covalent link between the peptide and the fragment. Default for the SnP case:

```json
[
  {
    "pep_resid": 6,
    "pep_atom": "CD",
    "frag_resid": 1,
    "frag_atom": "NH2",
    "pep_remove": ["OE2", "HE2"],
    "frag_remove": ["CM", "HM1", "HM2", "HM3", "CAP", "OAP"]
  }
]
```

Generalizes to multi-bond cofactors — e.g., a Zn-finger with 2 Cys + 2 His coordination passes four entries, one per coordinating atom.

## Requirements

- **AmberTools ≥ 24** (`antechamber`, `parmchk2`, `tleap`) — bundled via `pixi.toml`
- **ORCA** — user-installed. Local: set `ORCA_BIN` env var or put `orca` on `$PATH`. HPC: `module load orca/...` inside the SLURM script passed to `ep_orca_run`.
- **Python ≥ 3.11** with `numpy`, `scipy`, `periodictable`, `biopython`

Windows is not a supported conda platform for AmberTools; use Docker or WSL2.

## HPC execution (`ep_orca_run`)

`ep_orca_run` inherits from `HPCNodeBase` — the same base class used by pdbmdauto's `gmx_mdrun`. Toggle between `local` and `remote` via the **Execution Mode** option. Remote runs use an HPC profile (configured in Salpa **Settings → HPC Profiles**) plus a user-provided SLURM script. A formal, Snellius-validated reference template ships at `ep_orca_run/templates/default-slurm.sh` — open it, copy the contents into the GUI's **SLURM Job Script** field, and edit the cluster-specific lines (partition, module tree, memory, scratch). `{{VARIABLE}}` placeholders are substituted at submit time. Note: the textarea does not auto-prefill yet — the front-end currently only reads from saved-config, not from schema defaults.

## Installation

### In Salpa/BoCoFlow

1. Open Marketplace → Add Source → point to this repository or the published marketplace
2. Install `metalparm-vwf` — all 20 nodes appear under "Force Field Parameterization"

### Standalone (development)

```bash
cd packages/metalparm-vwf
pixi install
```

## Data flow (predecessor pattern)

Nodes use the 3-tier priority for file inputs:

1. **Explicit user config** — user sets the file path in the node panel
2. **Predecessor auto-discovery** — read upstream node's declared output (`output_mol2`, `output_frcmod`, etc.)
3. **Default filename** — fallback to a fixed name in the working directory

Connecting `ep_bond_detection → ep_mol2_generation → ep_charges → ep_forcefield_assembly → ep_library_generation` wires everything automatically. `ep_charges` also takes the ORCA `.out` (auto-discovered via `output_out`) and overrides the zero-charge `output_mol2` with the charged copy, which `ep_forcefield_assembly` then picks up. `ep_seminario_orca` branches off `ep_bond_detection` and feeds into `ep_forcefield_assembly`. `ep_fragment_fuse_topology` accepts `fragment_lib` / `fragment_frcmod` from `ep_library_generation` via `output_lib` / `output_frcmod`; `ep_apply_coords` then consumes its `complex.prmtop` plus the aligned PDBs from `fragment_align`.

## Bundled scripts

Each EasyParm node bundles only the upstream easyPARM Python scripts it actually calls in its own `scripts/` subdirectory (`ep_bond_detection/scripts/02_get_bond_angle.py`, etc.). This keeps each node self-contained at install time — no shared package directory to copy.

Nodes resolve their scripts via `_find_scripts_dir()`, which checks:

1. `EASYPARM_SCRIPTS` env var (dev override — single dir with all scripts)
2. `<node_dir>/scripts/` (bundled with node, used at install)
3. `<node_dir>/../../collect/easyPARM/scripts/` (metal-md source-tree fallback)

Note: `uff_data.txt` is bundled alongside `11_retrieve_uffdata.py` in `ep_forcefield_assembly/scripts/` because the script reads it relative to its own location.

`xtb_opt`, `snp_builder`, `peptide_builder`, `fragment_align`, `ep_fragment_fuse_topology`, and `ep_apply_coords` are pure-Python / pure-tleap wrappers and ship no external scripts.

## Shared environment

Declared in `pixi.toml` as `metalparm_vwf`. All 20 nodes share this env (AmberTools 24, xtb, numpy, scipy, biopython, periodictable, redis-py).

## Cookbook

Three companion documents in the metal-md repository (`dev-notes/`) cover the cofactor-fuse extension surface:

- **`interface-bonds-templates.md`** — `interface_bonds` JSON recipes for common topologies (SnP-amide, Zn-finger 2C2H, heme b axial His, disulfide-linked cofactor, Pt-N4 square-planar). Ready to paste into `fragment_align` and `ep_fragment_fuse_topology` panels; lists the matching `STANDARD_BOND_LENGTHS` / `FRAGMENT_HYBRIDIZATION_TABLE` entries (v1.4.2+).
- **`metal-fragment-builders.md`** — recipe for adding a new metal-cofactor builder node (extending the v1.7.0 shared-core split). Sketches for `zn_finger_builder`, `heme_builder`.
- **`cross-ff-linkage-frcmods.md`** — design notes for cross-FF linkage frcmods (v1.9.0). Explains why fusing GAFF2-typed cofactors to ff19SB peptides via covalent bonds needs a small parameter patch, the parm10-verified value table, and the recipe for adding new linkage types when new cofactor topologies land.

## Known limitations

### YASARA-style PDB atom names

`snp_builder.extract_peptide()` and `peptide_builder` (PDB-mode preprocessing, v1.8.0+) rewrite a small handful of YASARA-style atom names to ff19SB equivalents — `GLU.COOH → GLU.HE2`, `ASP.COOH → ASP.HD2`, `SER.HO → SER.HG`, `THR.HO → THR.HG1`, `TYR.HO → TYR.HH`, `CYS.HS → CYS.HG`. Without this, `tleap` keeps both the YASARA and the standard atom in the residue (close-contact warning), and downstream removals targeting standard ff19SB names (e.g. `HE2`) silently miss. Pass `rename_yasara_atoms=False` to opt out.

### Linkage frcmod is bond-type-specific

`ep_fragment_fuse_topology`'s `linkage_frcmod` parameter (v1.9.0+) closes the cross-FF parameter gap at covalent attachments between a GAFF2-typed cofactor and an ff19SB peptide. Two bundled patches cover GLU/ASP-amide bonds — pick by the GAFF2 atom type of the linkage nitrogen in your fragment `.lib`:

| Glu protonation | Glu CD type | Linkage-N type | Bundled file |
|---|---|---|---|
| GLU (charged) | `CO` | `ns` (sp2 amide N) | `amide_glu_gaff2.frcmod` |
| GLU (charged) | `CO` | `n`  (sp3 amine N) | `amide_glu_gaff2_n.frcmod` (v1.13.0+) |
| GLH (protonated) | `C`  | `n`  (sp3 amine N) | `amide_glh_gaff2_n.frcmod` (v1.14.0+) |

`peptide_builder` in PDB mode loading a PDB whose Glu has the HE2
hydrogen (e.g. YASARA-built structures) typically yields row 1/2.
`peptide_builder` in sequence mode with `GLH` written in the
sequence string yields row 3.

Both ship under `ep_fragment_fuse_topology/demo_data/linkages/` and are addressable as `node:demo_data/linkages/<filename>`. Other linkage types — Zn-finger Cys/His coordination, heme axial-His, disulfide-linked cofactors, Pt-N4 — need separate frcmods following the same convention. The recipe is in `dev-notes/cross-ff-linkage-frcmods.md` in the metal-md repository.

### Resolved in earlier versions (kept for traceability)

- **Atom-name preservation through `ep_mol2_generation`** — resolved in **v1.6.0** via the `xyz_to_pdb.py --template` flag; meaningful atom names from `snp_builder` (`NH2`, `CAP`, `OAP`, `CM`, `HM1-3`) now survive into `COMPLEX.lib` automatically. `DEFAULT_INTERFACE_BONDS` works directly without overrides.
- **`snp.frcmod` parameter coverage** — resolved in **v1.5.1** (metal-element vdW NONBON) + **v1.5.0** (lib/frcmod vocabulary alignment via no-ULS pass-through) + **v1.9.0** (cross-FF amide linkage at fuse time). The full SnP fragment-fuse Playwright spec passes end-to-end producing `complex.{prmtop,rst7,pdb}` (verified 1.9 min total, 0 errors).

## License

LGPL-2.1 — inherited from [easyPARM](https://github.com/abenmb/easyPARM). The bundled scripts in each node's `scripts/` directory are redistributed under the same license. Node wrappers (`*/node.py`), `snp_builder`, `ep_fragment_fuse_topology`, and `ep_apply_coords` are original work by the metal-md authors.

## Citations

- Abdelazim M. A. Abdelgawwad, Antonio Francés-Monerris. *easyPARM: Automated, Versatile, and Reliable Force Field Parameterization Workflow for Metal-Containing Molecules with Unique Labeling of Coordinating Atoms.* J. Chem. Theory Comput. 2025, 21, 4, 1817–1830.
- Allen, A. E. A.; Payne, M. C.; Cole, D. J. *Harmonic Force Constants for Molecular Mechanics Force Fields via Hessian Matrix Projection.* J. Chem. Theory Comput. 2018, 14, 1, 274–281.
- Case, D.A. et al. *AMBER 2024* — `tleap`, `ff19SB`, `gaff2`.
