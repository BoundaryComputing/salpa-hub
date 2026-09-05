# demo_data

`trp_cage.fasta` — the Trp-cage miniprotein (PDB 1L2Y), 20 residues, `NLYIQWLKDGGPSSGRPPPS`:
the smallest protein with a real fold, so a structure prediction finishes in under a minute
and can be checked against the deposited NMR structure. The same sequence `DEMO_CONFIG` in
`node.py` declares inline, with `msa_mode = "empty"` for the fast path.

Running this node needs a Salpa account with cloud access: `salpa smoke` reaches the gateway
and stops at authentication. That is expected; the declaration is what `salpa validate` and a
reviewer need.
