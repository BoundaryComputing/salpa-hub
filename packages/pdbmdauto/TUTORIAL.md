# From a PDB entry to a simulation-ready protein, without the command line

**Who this is for:** you have a protein structure — a PDB code, or a file someone sent you — and
you want to simulate it, or you want to understand what "preparing a structure for MD" actually
involves before you trust a collaborator's setup. No simulation background assumed. Every acronym
is defined the first time it appears.

**Time:** about five minutes of computer time on a laptop, two on an Apple Silicon Mac, after a
one-time install of a few minutes (about twenty on Windows). The example needs a network
connection, because it fetches the structure from the Protein Data Bank.

> Prefer the formal version? `README.md` is the methods-section treatment: scope, per-stage
> methodology, node table, references. This page is the plain-language one.

---

## What it does

You give it a PDB code. It gives you back a complete protein — missing pieces rebuilt, hydrogens
added, charges decided for your pH — sitting in a box of salt water, relaxed so nothing is
clashing, plus a short test run of molecular dynamics (MD) to prove the whole thing moves. Every
intermediate file stays on disk, so you can look at each stage.

## When you'd use it

- You have a crystal structure and want to run MD on it, and the setup is the part you have never
  done yourself.
- You want a starting point you can hand to a simulation collaborator that is reproducible — same
  inputs, same steps, same result — rather than a folder of files with a story attached.
- You want to see what each preparation step changes, on a small protein, before doing it on yours.

You would **not** use it for a system that needs a bound ligand, a cofactor or a metal: those are
removed early on (step 4), and parameterising them is a separate job. And the test run at the end
is two picoseconds — a smoke test, not a result.

## The idea, in plain terms

A crystal structure is a photograph with pieces missing. Some residues were too mobile to be seen
and are simply absent from the coordinates. Hydrogen atoms are almost never resolved, so the file
has none. Nothing in it says which acidic and basic groups are charged, because that depends on the
pH of the solution, not on the crystal. And there is no water around it. Simulation needs all of
that decided.

- **Rebuilding the missing residues** is homology modelling: the software knows the full sequence,
  sees which residues are absent, and builds them from what similar stretches of protein look like
  in other structures. In this example the missing pieces are at the ends of the chains — a
  five-residue tag left over from expression, and the first residue of the peptide.
- **Protonation** is the pH decision. At pH 7 some histidines are charged and some are not; a tool
  called PROPKA estimates each one's pKa from its surroundings, and PDB2PQR adds the hydrogens
  accordingly.
- **The force field** is the rulebook the simulation uses: how strongly each pair of atoms attracts
  or repels, how stiff each bond is. Building the *topology* means writing that rulebook out for
  your specific protein. This pipeline uses amber99sb, a standard choice for proteins.
- **Energy minimisation** lets the model relax out of bad contacts. Newly built residues and newly
  added hydrogens can start too close to their neighbours; minimisation nudges everything downhill
  until the worst clashes are gone. The number to watch is the potential energy, which should fall
  and flatten.
- **Restraints** hold the atoms you trust while the ones you built settle. The crystal's atoms are
  held in place during the first relaxation; only the rebuilt residues move.
- **Solvation** puts the protein in a box of water with sodium and chloride ions at physiological
  concentration, and enough extra ions to make the whole box electrically neutral.

## How it works, step by step

| # | Step | What it does, plainly | Tool | Why it matters |
|---|---|---|---|---|
| 1 | PDB FASTA Parser | Downloads the entry, reads its sequence and its list of missing residues | Biopython | Everything else starts from this list |
| 2 | Generate Alignment | Lines up what was seen against the full sequence, per chain | — | The gaps in the alignment are what gets built |
| 3 | Multi-Chain Alignment | Merges the per-chain alignments | — | The complex is modelled as one unit |
| 4 | Merge PDB Chains | Writes the chains as one model, dropping water and ligands | Biopython | The model the rest of the pipeline works on |
| 5 | Fix Missing Residues | Builds the missing residues and their sidechains | ProMod3 | A continuous protein, with no holes |
| 6 | pKa + GROMACS EM | Decides charges for pH 7, adds hydrogens, writes the force-field topology, minimises in vacuum | PROPKA, PDB2PQR, GROMACS | The electrostatics of your simulation are decided here |
| 7 | Original Atom Groups | Marks which atoms came from the experiment | GROMACS | So the next step can hold them still |
| 8 | GMX MD Relaxation — in vacuum, restrained | Lets the rebuilt residues settle while the crystal's atoms are held | GROMACS | The built pieces adopt sensible positions without moving what was measured |
| 9 | GMX Solvate & Ionize | Puts everything in a 5 nm box of water with 0.15 M NaCl, neutral | GROMACS | Proteins live in salt water, not vacuum |
| 10 | GMX MD Relaxation — in water | Minimises the whole box | GROMACS | Water settles onto the protein before it moves |
| 11 | GROMACS MD Run (Local) | Runs 2 ps of molecular dynamics | GROMACS | Proves the prepared system runs; the trajectory is the starting point for real work |

## Try it: a worked example

The example is PDB entry **4Z8J**: a small protein domain (a PDZ domain from a sorting protein
called SNX27) holding on to the tail of a hormone receptor. Two chains, about 110 residues, six of
them missing from the crystal.

1. Open Salpa and choose **Start from Template**.
2. Pick **PDBmdAuto Full Pipeline** (under *molecular-dynamics*). Its walkthrough is one click
   away on the same screen — **Read the walkthrough** — with figures from a real run.
3. Choose a working directory. This is the only thing you type. An empty folder is best.
4. Press **Run**. Nothing else needs setting: the PDB code, the pH and the case name are in the
   template.

Watch the eleven steps go green in order. On a 2019 Intel laptop the whole thing took three to
five and a half minutes across three runs; on an Apple Silicon Mac Mini, about two. The first
time, the app also downloads and builds the tools (3.1 GB, about three and a half minutes on a
Mac with a fast connection, about twenty on Windows, where it sets up a Linux environment for you).

Afterwards, open the working directory and go to `pdbmdauto-e2e-full/e2e_4z8j/`. Three files are
worth a look:

- `Merge/fixed.pdb` — the completed protein. Open it in any structure viewer and find the ends of
  the chains: the orange residues in the walkthrough's first figure are the rebuilt ones.
- `gmx/ion.gro` — the protein in its box of water and ions, 12,193 atoms.
- `gmx/md.gro` with `gmx/md.trr` — the final coordinates and the eleven frames of the test run.

## Reading the results

What "good" looks like:

- Six residues rebuilt, and the chains continuous — the step-5 log says chain A went from 96 to
  101 residues and chain B from 7 to 8.
- The minimisation energies fall and flatten. In vacuum (step 6) the run converges in about fifteen
  steps. In water (step 10) the energy drops by a large amount as the water settles onto the
  protein, and is still gently descending when the 500-step limit stops it — that is expected for
  this demonstration; a production study would let it run longer.
- The system is neutral: the topology lists 12 sodium and 11 chloride ions, which cancel the
  protein's charge of −1 at pH 7 and add 0.15 M salt.
- The MD run finishes: `gmx/md.log` ends with a *Performance* line and *Finished mdrun*.

When something looks wrong:

| You see | Likely cause | What to do |
|---|---|---|
| Step 1 fails immediately | No network — 4Z8J could not be fetched | Connect, or point *Input mode* at a local file |
| Step 5 fails | ProMod3 or OpenMM problem in the environment | Check the Log Center; reinstalling the package rebuilds the environment |
| A step stops at Solvate & Ionize | A working directory the tools could not use | Choose a plain, empty folder and rerun |
| The first run is very slow | It is building the environment, or the first run on Windows | Wait; the second run is minutes |
| Numbers differ from the walkthrough's | Normal: rebuilt residues and MD velocities vary between runs | The counts (residues, atoms, ions) should not |

## Honest limits

- The rebuilt residues are predictions. Terminal extensions especially have nothing to pack
  against and should be read with that in mind.
- Ligands, cofactors, metals and crystal waters are removed at step 4. A system that needs them
  needs a different preparation.
- Two picoseconds of dynamics is a smoke test. Real sampling takes nanoseconds to microseconds.
- The force field is an approximation, and protonation is a prediction made once, at the start.
- The *Model Terminal Extensions* switch on Fix Missing Residues has no effect in this version; the
  termini are always built.

## Going deeper

- The walkthrough for this template (`workflows/pdbmdauto-pipeline.md`) names every file each step
  writes and the parameters it uses, with figures from a reference run.
- The tools, each with its own documentation: [GROMACS](https://www.gromacs.org),
  [ProMod3 / SWISS-MODEL](https://swissmodel.expasy.org), [PDB2PQR and PROPKA](https://www.poissonboltzmann.org),
  [Biopython](https://biopython.org). The README lists the papers to cite.
- To prepare your own protein: change the PDB id and the pH on the first node, or switch *Input
  mode* to a local file. Everything else stays the same.
