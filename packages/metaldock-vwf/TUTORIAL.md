# Metal docking, without the command line

**Who this is for:** you work at a bench, you have a metal-containing compound and a
protein you think it binds, and you would like a testable hypothesis about *where* and
*how* before you spend another month at the bench. No simulation background assumed.
Every acronym is defined the first time it appears.

**Time:** about five minutes of computer time. The example is bundled — nothing to
download, no accounts, no command line.

> Prefer the formal version? `README.md` is the methods-section treatment: scope,
> per-stage methodology, node table, references. This page is the plain-language one.

---

## What it does

You give it a protein structure and a metal-containing molecule. It gives you back a
set of predicted **poses** — specific positions and orientations of that molecule
inside the protein — each with a score saying how good it looks, and a list of which
amino acids it touches.

You get a ranked shortlist of binding hypotheses, not an answer.

## When you'd use it

- You have a metal complex that does something in an assay and you want a plausible
  binding site to design mutants around.
- You want to know which residues to mutate first to test whether a site is real.
- You have two candidate compounds and want a reason to prioritise one.

You would **not** use it to predict how tightly something binds in absolute terms.
More on why in *Reading the results*.

## The idea, in plain terms

Docking is speed dating for molecules. The software tries your compound in the protein
in a very large number of positions and orientations, scores each arrangement with a
quick approximate energy function, and keeps the ones that score best.

That scoring function is a rulebook: how strongly does an oxygen here attract a
hydrogen there, how much does it cost to jam two atoms too close. The rulebook was
written by fitting to thousands of ordinary organic molecules — carbon, nitrogen,
oxygen, the usual cast.

**Metals are not in that cast.** A rhenium or ruthenium atom sitting at the centre of
your compound pulls electrons around in a way the standard rulebook has no entry for.
Dock it naively and the software is guessing about the single most important atom.

This pipeline fixes that by asking a quantum-chemistry program a narrow question first:
*given this exact molecule, how is the electric charge actually distributed across its
atoms?* Those numbers — **partial charges** — are then handed to the docking program in
place of the defaults. The search is conventional; what changes is that the scoring now
has physically sensible numbers near the metal.

One more idea you will meet. A metal atom has a fixed number of "slots" for
neighbouring atoms — its **coordination sphere**. Sometimes a slot is left empty, and
that vacancy is exactly where the protein binds. The bundled example is one of these
cases, and the pipeline marks the empty slot explicitly so the docking search knows to
put something there.

## How it works, step by step

Six nodes, left to right on the canvas. Each does one job and hands its results to the
next.

| # | Step | What it does, plainly | Tool | Why it matters |
|---|------|----------------------|------|----------------|
| 1 | **Protein Prep** | Strips out water and anything that isn't the protein, works out which acidic and basic groups are charged at your chosen pH, adds the hydrogen atoms crystallography can't see, and writes the result in the format the docking program reads. | pdb2pqr, AutoDockTools | A crystal structure is a sketch, not a finished model. Hydrogens and charges are missing and they decide what sticks to what. |
| 2 | **Ligand Prep** | Reads your compound's 3D coordinates and works out which atoms are bonded to which, building a map of the molecule. | OpenBabel | Everything downstream needs to know the molecule's connectivity, not just a cloud of coordinates. |
| 3 | **QM Charges** | Runs a quantum-chemistry calculation on your compound alone and records how much charge sits on each atom. | xtb (default) or ORCA | This is the step that makes the metal behave. See *Fast or careful* below. |
| 4 | **Ligand PDBQT** | Packages the compound with its new charges, decides which bonds are allowed to rotate during the search, and freezes the ones around the metal. Marks any empty coordination slot. | — | Bonds to a metal don't swivel freely. Letting them would waste the search on impossible shapes. |
| 5 | **AutoDock Run** | The search itself. Builds a 3D grid of "how attractive is this spot" around the target region, then tries many positions and orientations, keeping the best. | AutoGrid4, AutoDock4 | This is the docking. Everything before it was preparation. |
| 6 | **Results Analysis** | Ranks the poses, lists which residues each one touches, and — if you gave it a known answer — measures how far off it was. | — | Turns a pile of coordinates into something you can act on. |

### Fast or careful

Step 3 offers a choice, and it is the only one you really need to think about.

| | **xtb** (the default) | **ORCA** |
|---|---|---|
| Speed | Under a second | Minutes to hours |
| Method | Semi-empirical — a fast approximation with experimental shortcuts baked in | Density functional theory (DFT) — a full quantum calculation |
| Setup | Installed with the package | Download separately, free for academic use |
| Use it for | Getting going, screening, exploring | Numbers going into a paper |

On the bundled example, xtb puts the charge on the rhenium atom at **+0.749** where
the full DFT calculation gives **+0.704** — a difference of about 0.05 in units where
a whole electron is 1.0. Close enough to find the same binding site; not what you'd
publish a charge value from.

## Try it: the bundled example

The example is the **1JZI** case from the MetalDock paper: a rhenium complex binding
azurin, a small copper protein from *Pseudomonas aeruginosa*. Both input files ship
inside the package.

1. Open Salpa, and from the workflow templates choose **MetalDock 1JZI Re Pipeline**
   (it's under the *molecular-docking* category).
2. Pick a folder for the results. This is the only thing you have to supply.
3. Press run.

About three minutes later — most of it the docking search in step 5 — you'll have six
folders. The interesting ones:

```
qm/enriched_graph.json         the per-atom charges
pdbqt/1jzi_re_ligand.pdbqt     the compound, ready to dock
docking/…dlg                   the raw search output, ten poses
analysis/1jzi_re_analysis.json the answer
```

## Reading the results

Open `analysis/1jzi_re_analysis.json`. Three things matter.

**`binding_energies`** — one number per pose, in kcal/mol, more negative meaning
better. On this example they land around **−4.6**.

> **Read these as a ranking, never as an affinity.** A docking score is a fast
> approximation, and the error bars on the absolute number are larger than the
> differences you'll typically care about. Pose A scoring better than pose B is
> informative. "−4.6 kcal/mol" as a measured binding strength is not — do not convert
> it to a Kd. Treat the ranking as a hypothesis to test at the bench.

**`interacting_residues`** — which amino acids each pose touches. This is the most
actionable output: it's your mutagenesis shortlist. On this example the top pose
contacts about **12** residues, clustered around one face of the protein.

**`rmsd_values`** — only meaningful if you supplied a known structure to compare
against, as the bundled example does. Root-mean-square deviation is the average
distance between the predicted atoms and the real ones, in ångströms. Here it comes
out **5.5–5.9 Å**.

Is that good? Honestly: it's the same answer the published DFT calculation gives
(5.91 Å), so the pipeline is behaving. But 5–6 Å is a *region*, not a pose — it means
the search found the right neighbourhood and not the exact orientation. That is a
normal and useful outcome for a hard metal-containing case, and it is worth knowing
before you over-read a picture.

### When something looks wrong

| What you see | Likely cause |
|---|---|
| Every pose scores about the same and they're scattered | The search box is too big, or centred on the wrong place. Set the box centre on the site you care about. |
| The compound never goes near the metal site | The empty coordination slot wasn't marked. Check *Vacant Site* is on in step 4. |
| Charges look implausible on the metal | Wrong total charge or spin state on step 3. These come from your chemistry and cannot be guessed. |
| The whole thing is slower than expected | You selected ORCA. That is expected — it's doing real quantum chemistry. |

## Honest limits

- **Scores are rankings.** Said three times because it is the mistake everyone makes.
- **The protein is rigid.** Real proteins move to accommodate what binds them. This
  search does not let them, so a site that only opens up on binding will be missed.
- **No solvent.** Water is modelled crudely or not at all. Water-mediated contacts,
  which are common and often decisive, are invisible here.
- **The starting structure matters.** A poor or wrong crystal structure gives a
  confident and wrong answer. Garbage in, confident garbage out.
- **This is a hypothesis generator.** Its output is a list of things worth testing.

## Going deeper

- Hakkennes, M. et al. *MetalDock: An Open-Source Docking Tool for Metal-Organic
  Compounds.* J. Chem. Inf. Model. 2023. doi:10.1021/acs.jcim.3c01582 — the method this
  package implements, including the 1JZI case above.
- Morris, G. M. et al. *AutoDock4 and AutoDockTools4.* J. Comput. Chem. 2009.
  doi:10.1002/jcc.21256 — the docking engine and its scoring function.
- Bannwarth, C. et al. *Extended tight-binding quantum chemistry methods.* WIREs
  Comput. Mol. Sci. 2021. doi:10.1002/wcms.1493 — what xtb actually does.
- Marenich, A. V. et al. *Charge Model 5.* J. Chem. Theory Comput. 2012.
  doi:10.1021/ct200866d — where the partial charges come from.

For the formal treatment — per-stage methodology, the full node table, platform and
licensing constraints — see `README.md`.
