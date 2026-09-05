# demo_data

`lac_operator.fasta` — the E. coli *lac* operator O1, 21 bp, `AATTGTGAGCGGATAACAATT`. A
real, textbook sequence a reviewer can verify by eye: it is nearly palindromic, which is why
the Lac repressor binds it as a dimer. Small enough to score in seconds; `DEMO_CONFIG` in
`node.py` points the node at it in *score* mode.

Running this node needs a Salpa account with cloud access: `salpa smoke` reaches the gateway
and stops at authentication. That is expected; the declaration is what `salpa validate` and a
reviewer need.
