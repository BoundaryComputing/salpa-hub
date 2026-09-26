# DiffDock Docking (Modal)

Blind protein-ligand docking using [DiffDock](https://github.com/gcorso/DiffDock), run on an A10G GPU by Salpa Compute. Sign in to Salpa to use it; no Modal account is needed.

DiffDock is a generative diffusion model that predicts how small molecules bind to protein targets — a critical step in computational drug discovery. Unlike traditional docking (AutoDock Vina, Glide), DiffDock uses deep learning to simultaneously predict binding pose, position, and orientation without requiring a predefined search box.

**Reference**: Corso et al., "DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking", ICLR 2023 (MIT License)

## Quick Start

1. Drag "DiffDock Docking (Modal)" onto the canvas
2. Select a protein PDB file (or connect from an upstream node like ESMFold)
3. Enter a ligand SMILES string (e.g., `C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1` for Erlotinib)
4. Click Execute — results arrive as ranked SDF pose files

## Timing

| Scenario | Duration | Notes |
|----------|----------|-------|
| **Cold start** | ~10 minutes | ESM-2 model (2.6 GB) loading into GPU memory |
| **Warm start** | ~40 seconds | Container reused within 2-minute window |

**Important**: This node is designed for **interactive single-protein docking**, not batch processing. For docking many protein-ligand pairs, use a local DiffDock installation or HPC pipeline.

## Parameters

### Required (one of)

| Parameter | Description |
|-----------|-------------|
| **Protein PDB** | Protein structure in PDB format. Must have 20+ residues. |
| **Ligand SMILES** | SMILES string for the small molecule (primary input) |
| **Ligand SDF** | Alternative: ligand structure in SDF format (used if SMILES is empty) |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Number of Poses** | 10 | Binding poses to generate (1-40) |
| **Inference Steps** | 20 | Denoising steps (10-40, higher = more accurate but slower) |
| **Samples per Complex** | 10 | Samples per protein-ligand complex (1-40) |
| **Output Folder** | (workflow folder) | Where to save results. A relative folder is placed inside the workflow's folder |
| **Output Prefix** | (auto) | Filename prefix (auto-generates `diffdock_YYYYMMDD_HHMMSS`) |

## Input

- **Protein PDB**: File parameter or predecessor data (`pdb_content`, `protein_pdb`, or `output_file`)
- **Ligand**: Either SMILES string or SDF file

The protein must have at least ~20 residues — very small peptides cause graph construction errors.

## Output

### Files

| File | Description |
|------|-------------|
| `{prefix}.tar.gz` | The full result from DiffDock |
| `{prefix}_rank1_confidence-X.XX.sdf` | Best-ranked docking pose (SDF format) |
| `{prefix}_rank2_confidence-X.XX.sdf` | Second-ranked pose |
| ... | Up to `num_poses` ranked SDF files |

### Result Data

```python
result.data = {
    "output_file": "/path/to/diffdock_20260306_143200.tar.gz",
    "extracted_sdfs": ["/path/to/rank1.sdf", ...],
    "num_poses": 10,
    "top_confidence": -1.18,
    "confidence_scores": [-1.18, -1.46, -1.72, ...],
    "processing_time_seconds": 42.3,
    "protein_pdb": "...",        # Original input (for chaining)
    "ligand_smiles": "...",      # Original input
}
```

## Understanding the Results

### Confidence Scores

DiffDock confidence scores are **negative numbers**. Higher (less negative) = better.

| Score Range | Interpretation |
|-------------|---------------|
| > -1.0 | Excellent — high confidence binding pose |
| -1.0 to -2.0 | Good — likely correct binding site |
| -2.0 to -3.0 | Moderate — plausible but verify |
| < -3.0 | Low confidence — may not represent true binding |

The filename encodes the score: `rank1_confidence-1.46.sdf` means confidence = **-1.46** (the dash after "confidence" is a minus sign).

### Visualizing in PyMOL

```bash
# Load protein and docking poses
pymol

# In PyMOL command line:
load /path/to/your_protein.pdb, protein
show cartoon, protein
color slate, protein

# Load ranked poses
load /path/to/diffdock_rank1_confidence-1.18.sdf, pose1
load /path/to/diffdock_rank2_confidence-1.46.sdf, pose2

# Style ligands
show sticks, pose1
show sticks, pose2
color cyan, pose1
color yellow, pose2

# Zoom to binding site
zoom pose1, 10

# Show binding pocket surface
select pocket, protein within 5 of pose1
show surface, pocket
set transparency, 0.5, pocket
```

### Validating Against Known Structure

If you have a co-crystal structure (e.g., PDB 1M17 for EGFR + Erlotinib):

```bash
# In PyMOL:
fetch 1M17
select crystal_lig, 1M17 and resn AQ4
show sticks, crystal_lig
color green, crystal_lig
# Compare DiffDock prediction (cyan) vs crystal truth (green)
```

## Demo Data

| File | Description |
|------|-------------|
| `demo_data/EGFR_1M17.pdb` | EGFR kinase domain (PDB 1M17, 324 residues, chain A) |
| `demo_data/trp_cage_1L2Y.pdb` | Trp-cage miniprotein (PDB 1L2Y, 20 residues — minimum viable) |

**Recommended demo**: EGFR (1M17) + Erlotinib SMILES `C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1`

This is a classic drug discovery benchmark — Erlotinib (Tarceva) is an FDA-approved EGFR inhibitor. You can validate the predicted pose against the crystal structure in PDB 1M17 (ligand residue AQ4).

## Example Workflows

### Simple Docking
```
[Protein PDB File] → [DiffDock] → [Output SDF Files]
```

### Structure Prediction + Docking Pipeline
```
[Protein Sequence] → [ESMFold / Boltz-2] → [DiffDock] → [Output SDF Files]
```
DiffDock accepts predecessor PDB output automatically via `pdb_content` in predecessor data.

## Large results

A result up to 16 MiB (compressed) comes back directly, which is every DiffDock run so far. A larger one comes back as its poses plus a download link: the node downloads the full archive into the folder, checks its size and SHA-256, and then has the copy on Salpa Compute deleted. An archive that is never downloaded is deleted within 24 hours.

## Timeout

A new node starts with a timeout of 1500 seconds, enough for a cold start and the run. A node saved in a workflow earlier keeps the timeout it was saved with (often 600 seconds): raise it under **Advanced Options > Timeout (seconds)**. The node warns when its timeout is shorter than a run may need.

## Stopping a run

Pressing Stop in Salpa also stops the run on Salpa Compute (Salpa 0.8.3 and later): the GPU work
ends within about half a minute. So does a run the node gives up on when its timeout passes. The
GPU time the run used until then counts toward your monthly GPU hour, and a run whose result you
did not receive is never charged. The same holds for a run that fails because of its input, such
as one that runs out of GPU memory; a failure on our side counts nothing.

You can have two GPU runs in progress at a time. A third is refused until one of them finishes or
is stopped.

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "No edges and no nodes" error | Protein too small (<20 residues). Use a real protein, not a mini peptide |
| Timeout after ~10 minutes | Cold start — try again immediately (container should now be warm) |
| "sign-in has expired" | Sign in to Salpa again |
| "GPU quota exceeded" | Your Salpa Compute quota for this period is used up; it resets at the start of the next period |
| "produced no" poses, or no SDF files | Check that the ligand SMILES is valid (use RDKit to verify) |
| "stopped waiting for DiffDock" | A cold start took too long; run the node again, the container should now be warm |
