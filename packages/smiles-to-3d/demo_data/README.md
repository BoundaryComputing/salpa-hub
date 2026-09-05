# demo_data

`ligands.smi` — three small, real molecules (ethanol, aspirin, caffeine) in the `.smi`
convention: SMILES, whitespace, an optional name; `#` lines are comments. It is what
`salpa smoke`, `pixi run test` and a reviewer run the node on.

Hand-checkable expectations from the formulas alone: with hydrogens added, ethanol
C2H6O has 9 atoms, aspirin C9H8O4 has 21, caffeine C8H10N4O2 has 24.
