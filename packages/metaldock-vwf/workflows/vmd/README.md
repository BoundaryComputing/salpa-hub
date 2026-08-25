# Regenerating the figures

The two renders in `../hsa_ferrocene_walkthrough.html` come from VMD, run headless
with the Tachyon software renderer — no display, no window, no screenshot.

```bash
VMD=/Applications/VMD*/Contents/vmd/vmd_MACOSXX86_64
$VMD -dispdev text -e ligand.tcl     # ferrocene alone      -> ligand.tga
$VMD -dispdev text -e pocket.tcl     # the pose in the site -> pocket.tga
ffmpeg -i ligand.tga ligand.jpg
```

`pocket.tcl` expects two files beside it, both produced by a run of the template:

| copy from | to |
|---|---|
| `<working_dir>/protein/clean_1ao6_A.pdb` | `rec.pdb` |
| `<working_dir>/docking/hsa_fe_ligand_1.pdbqt` | `pose1.pdb` — strip the PDBQT charge/type columns first, or load with `type pdb` |

Two things that cost time when this was first set up:

- **`display resetview` frames the *top* molecule.** Making the 21-atom ligand top and
  then scaling *in* fills the frame with a single iron atom. Scale *out* from it
  (`scale by 0.30`) to bring the pocket into view.
- **Render at the size you want.** `display resize 880 880` before `render` gives a
  sharp image; upscaling a 512 px render afterwards makes a larger file that is no
  sharper.
