# Licensing — Salpa Hub

Salpa Hub ships **100% free and open**. There is no single repo-wide license: the
license is finalized **per package** as each package's dependency chain settles.
This file records the current posture; each package also carries its own `NOTICE`
and/or `LICENSE`.

Premium execution lives entirely on **Salpa Compute**. Nothing in this Hub is ever
paywalled.

## hello-world-pipeline

- **License:** MIT.
- **Attribution:** original BoCoFlow work; pure-Python stdlib, no bundled third-party code.
- **Platforms:** all (stdlib-only, runs IN_PROCESS — no conda environment).
- **Third-party tools invoked (not redistributed):** none.

## metalparm-vwf

- **License:** LGPL-2.1 — inherited from [easyPARM](https://github.com/abenmb/easyPARM).
- **Attribution:** see [`packages/metalparm-vwf/NOTICE`](packages/metalparm-vwf/NOTICE).
  The bundled easyPARM scripts retain their original copyright/license headers; the
  node wrappers and metallopeptide fusion nodes are original metal-md work, also
  under LGPL-2.1.
- **Upstream citation:** Abdelgawwad & Francés-Monerris, *easyPARM*, J. Chem. Theory
  Comput. 2025, 21, 4, 1817–1830.
- **Platforms:** linux-64, osx-64, osx-arm64. **No Windows** — AmberTools has no
  Windows conda package (use Docker or WSL2).
- **Third-party tools invoked (not redistributed):** AmberTools, ORCA (free for
  academia, user-supplied), xtb, ParmEd, MDAnalysis, packmol/packmol-memgen.

## pdbmdauto

- **License:** MIT.
- **Attribution:** original BoCoFlow work; no bundled third-party source.
- **Platforms:** linux-64, osx-64, osx-arm64. **No Windows** — GROMACS and ProMod3
  have no Windows conda packages (use Docker or WSL2).
- **Third-party tools invoked (not redistributed):** GROMACS, ProMod3 (Apache-2.0),
  PDB2PQR, PROPKA, pdb-tools, Biopython; SLURM (`sbatch`) for the HPC `gmx_mdrun` node.

## metaldock-vwf

- **License:** ⚠ **Academic / non-commercial use only.** The AutoDock file-prep steps
  depend on **MGLTools / AutoDockTools**, distributed under the Scripps "MGLTOOLS
  SOFTWARE LICENSE AGREEMENT" (**non-commercial use only**). Until the planned **Meeko**
  (Apache-2.0) swap removes that dependency, this package must **not** be used
  commercially and carries no permissive/commercial label — regardless of the LGPL-2.1
  license of its own wrapper code. _Not legal advice._
- **Attribution:** see [`packages/metaldock-vwf/NOTICE`](packages/metaldock-vwf/NOTICE).
  Refactored from **MetalDock** (MIT); node wrappers are original metal-docking work.
- **Upstream citation:** Hakkennes et al., *MetalDock: An Open-Source Docking Tool for
  Metal-Organic Compounds*, J. Chem. Inf. Model. 2023, doi:10.1021/acs.jcim.3c01582.
- **Platforms:** linux-64, osx-64 only. **No osx-arm64 / Windows** — MGLTools and
  AutoDock have no builds there (Apple Silicon runs via Rosetta/osx-64).
- **Third-party tools invoked (not redistributed):** MGLTools/AutoDockTools
  (non-commercial), AutoDock4/AutoGrid4 (GPL-2), PDB2PQR, OpenBabel, ASE; ORCA (free for
  academia, user-supplied) for QM/CM5 charges.

---

_More packages are added to this table as they are mirrored in, each with its own
license posture and attribution._
