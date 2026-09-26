# Boltz-2 Prediction (Modal)

Protein structure prediction with [Boltz-2](https://github.com/jwohlwend/boltz), run on an
H100 GPU by Salpa Compute. Sign in to Salpa to use it; no Modal account is needed.

## Input

- **Protein Sequence** (one-letter codes), a full Boltz YAML file in **YAML Configuration**, or
  a sequence from the node before this one.
- **MSA Mode**:
  - `empty`: no alignment. Fast (under a minute for a small protein), less accurate.
  - `server`: Boltz-2 builds an alignment with the public ColabFold MSA server
    (api.colabfold.com), which receives the protein sequences to do so. Slower, more accurate.
  - `provided`: your own alignment, in **MSA A3M File**.

## Output

Files go to the workflow's folder, or to **Output Folder** (a relative folder is placed inside
the workflow's folder):

| File | What it is |
|---|---|
| `{prefix}_model_0.cif` | The predicted structure |
| `{prefix}_confidence_model_0.json` | Confidence scores |
| `{prefix}_plddt_model_0.npz` | Per-residue pLDDT |
| `{prefix}_affinity.json` | Binding affinity, when the input asks for it |
| `{prefix}.tar.gz` | The full result: every file Boltz-2 wrote, including the PAE and PDE matrices and the alignment |

The node's result data names each file (`structure_file`, `confidence_file`, `plddt_file`,
`output_file`), and `archive_complete` says whether `{prefix}.tar.gz` holds the full result.

## Large results

A result up to 16 MiB (compressed) comes back directly. A larger one comes back as its main
files plus a download link; MSA-server runs and large complexes are the usual cases, because
the alignment and the PAE/PDE matrices grow quickly. The node downloads the full archive into
the folder, checks its size and SHA-256, and then has the copy on Salpa Compute deleted. An
archive that is never downloaded is deleted within 24 hours.

## Timeout

A Boltz-2 run may take up to 30 minutes on the GPU, plus a cold start and the download. A new
node starts with a timeout of 2400 seconds. A node saved in a workflow earlier keeps the
timeout it was saved with (often 600 seconds): raise it under **Advanced Options > Timeout
(seconds)**. The node warns when its timeout is shorter than a run may need.

## Stopping a run

Pressing Stop in Salpa also stops the run on Salpa Compute (Salpa 0.8.3 and later): the GPU work
ends within about half a minute. So does a run the node gives up on when its timeout passes. The
GPU time the run used until then counts toward your monthly GPU hour, and a run whose result you
did not receive is never charged. The same holds for a run that fails because of its input, such
as one that runs out of GPU memory; a failure on our side counts nothing.

You can have two GPU runs in progress at a time. A third is refused until one of them finishes or
is stopped.

## Limits

On the H100 a complex of 2,340 tokens fits and one of 3,510 does not (measured). A run that
runs out of GPU memory, or produces no structure, is reported as a failure with its reason.
