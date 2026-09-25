# Chai-1 Prediction (Modal)

Structure prediction of proteins, nucleic acids, small molecules and their complexes with
[Chai-1](https://github.com/chaidiscovery/chai-lab), run on an H100 GPU by Salpa Compute. Sign
in to Salpa to use it; no Modal account is needed.

## Input

- **FASTA Input** with entity annotations (`>protein|name=…`, `>ligand|name=…` with a SMILES
  line, DNA, RNA), or
- **Protein Sequence**, optionally with **Ligand SMILES**, or a sequence from the node before
  this one.

## Output

Files go to the workflow's folder, or to **Output Folder** (a relative folder is placed inside
the workflow's folder):

| File | What it is |
|---|---|
| `{prefix}_best.cif` | The best-scoring of the five samples |
| `{prefix}_best_scores.npz` | That sample's scores |
| `{prefix}.tar.gz` | The full result: all five samples and their scores |

The node's result data names each file (`structure_file`, `scores_file`, `output_file`) and
carries every sample's aggregate score.

## Large results

A result up to 16 MiB (compressed) comes back directly, which is nearly every Chai-1 run. A
larger one comes back as its structures plus a download link: the node downloads the full
archive into the folder, checks its size and SHA-256, and then has the copy on Salpa Compute
deleted. An archive that is never downloaded is deleted within 24 hours.

## Timeout

A Chai-1 run may take up to 15 minutes on the GPU, plus a cold start. A new node starts with a
timeout of 1500 seconds. A node saved in a workflow earlier keeps the timeout it was saved with
(often 600 seconds): raise it under **Advanced Options > Timeout (seconds)**. The node warns
when its timeout is shorter than a run may need.
