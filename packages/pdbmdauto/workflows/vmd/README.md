# Regenerating the figures

The three figures in `../pdbmdauto-pipeline.md` come from one run of the template. Two are VMD
renders made headless with the Tachyon software renderer (no display, no window); the third is a
matplotlib plot of the energy files GROMACS wrote.

```bash
VMD=/Applications/VMD*/Contents/vmd/vmd_MACOSXX86_64
GMX=~/.bocoflow/pixi/environments/pdbmdauto/.pixi/envs/default/bin/gmx      # the package's own GROMACS
PY=~/.bocoflow/pixi/environments/pdbmdauto/.pixi/envs/default/bin/python    # has matplotlib

# copy the inputs beside these scripts (paths relative to the run's working directory)
cp "<working_dir>/pdbmdauto-e2e-full/e2e_4z8j/Merge/fixed.pdb" fixed.pdb
cp "<working_dir>/pdbmdauto-e2e-full/e2e_4z8j/gmx/ion.gro" ion.gro
cp "<working_dir>/pdbmdauto-e2e-full/e2e_4z8j/gmx/em_hbonds.edr" .
cp "<working_dir>/pdbmdauto-e2e-full/e2e_4z8j/gmx/em.edr" .

$VMD -dispdev text -e fixed.tcl      # -> fixed.tga
$VMD -dispdev text -e box.tcl        # -> box.tga
printf 'Potential\n0\n' | $GMX energy -f em_hbonds.edr -o em_vacuum.xvg
printf 'Potential\n0\n' | $GMX energy -f em.edr        -o em_solvated.xvg
$PY em_energy.py em_vacuum.xvg em_solvated.xvg ../figures/pdbmdauto-pipeline-energy.jpg
ffmpeg -y -i fixed.tga -q:v 4 ../figures/pdbmdauto-pipeline-rebuilt.jpg
ffmpeg -y -i box.tga   -q:v 4 ../figures/pdbmdauto-pipeline-box.jpg
```

Two things worth knowing:

- **ProMod3 renumbers from 1.** `Merge/merge.pdb` keeps 4Z8J's numbering (chain A 38–133, chain B
  587–593); `Merge/fixed.pdb` comes back as A 1–101 and B 1–8. The six rebuilt residues — the
  entry's REMARK 465 records A 33–37 and B 586 — are therefore `chain A and resid 1 to 5` and
  `chain B and resid 1` in `fixed.tcl`. Downstream, `ori_ndx_builder` still excludes exactly those
  six: its OriBackBone group is 414 atoms = 103 resolved residues × 4 backbone atoms + 2 OXT.
- **Render at the size you want.** `display resize` before `render` gives a sharp image;
  upscaling a small render afterwards makes a larger file that is no sharper.

The copied inputs (`fixed.pdb`, `ion.gro`, `*.edr`, `*.xvg`, `*.tga`) are not tracked.
