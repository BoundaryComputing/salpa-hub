# Evo2 DNA Model (GCP)

DNA sequence scoring, embeddings and generation with the Evo2 7B foundation model, run on an L4
GPU by Salpa Compute. Sign in to Salpa to use it; no Google Cloud account is needed.

## Input

- **DNA Sequence** (ACGT), a FASTA or text file in **Sequence File**, or a sequence from the
  node before this one.
- **Mode**: `score` (per-position log-likelihoods), `embed` (a mean-pooled embedding) or
  `generate` (extend the sequence by **Generate Length** bases, with **Temperature** and
  **Top-K**).

## Output

Files go to the workflow's folder, or to **Output Folder** (a relative folder is placed inside
the workflow's folder):

| File | What it is |
|---|---|
| `{prefix}.json` | The full result of the chosen mode |
| `{prefix}.fasta` | The generated sequence (`generate` mode) |

The result data carries the same values, plus `output_file`.

## Timeout

The first call after the service has been idle loads the model, which takes 2-5 minutes; later
calls take seconds. A new node starts with a timeout of 1200 seconds. A node saved in a workflow
earlier keeps the timeout it was saved with: raise it under **Advanced Options > Timeout
(seconds)** if a cold start is cut short. The node warns when its timeout is shorter than a run
may need.

A failed run is reported as a failure with its reason, never as an empty success.

## Stopping a run

Pressing Stop in Salpa also stops the run on Salpa Compute (Salpa 0.8.3 and later): the GPU work
ends within about half a minute. So does a run the node gives up on when its timeout passes. The
GPU time the run used until then counts toward your monthly GPU hour, and a run whose result you
did not receive is never charged. The same holds for a run that fails because of its input, such
as one that runs out of GPU memory; a failure on our side counts nothing.

You can have two GPU runs in progress at a time. A third is refused until one of them finishes or
is stopped.
