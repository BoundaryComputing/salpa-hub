# ProteinMPNN Design (GCP)

Designs amino acid sequences for a protein backbone with
[ProteinMPNN](https://github.com/dauparas/ProteinMPNN), run by Salpa Compute. Sign in to Salpa
to use it; no Google Cloud account is needed.

## Input

- **Input PDB**: a backbone structure, or a PDB from the node before this one.
- **Chains to Design**, **Fixed Positions** and **Omit Amino Acids** narrow the design;
  **Number of Sequences**, **Sampling Temperature**, **Model Variant**, **Backbone Noise
  Level** and **Random Seed** tune it.

## Output

The designed sequences, with their scores and sequence recovery, are in the node's result
data (`sequences`, `native_sequence`). With **Output Folder** set, they are also written to
`{prefix}_designed.fasta` there; a relative folder is placed inside the workflow's folder.

## Timeout

A run takes seconds. A new node starts with a timeout of 600 seconds. The node warns when its
timeout is shorter than that.

A failed run is reported as a failure with its reason, never as an empty success.
