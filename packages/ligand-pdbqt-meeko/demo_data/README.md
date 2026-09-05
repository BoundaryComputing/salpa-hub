# demo_data

`ligands.sdf` — ethanol, aspirin and caffeine as 3D conformers with explicit hydrogens,
written by the `smiles-to-3d` example node on its own demo data (ETKDGv3, seed 42, MMFF94).
It is what `salpa smoke`, `pixi run test` and a reviewer run this node on.

Hand-checkable expectations, from AutoDock's own rules: Meeko merges non-polar hydrogens, so
the PDBQT keeps heavy atoms plus polar hydrogens — ethanol 3 + 1 = 4 atoms, aspirin 13 + 1 = 14,
caffeine 14 + 0 = 14. Rotatable bonds (TORSDOF): ethanol 1 (C–O; the C–C bond moves only
hydrogens), aspirin 4 (the ester's two single bonds, the ring–carboxyl bond, and C–OH — an
ester is not an amide, so its C(=O)–O bond rotates), caffeine 0 (only methyl rotations, which
move nothing but hydrogens). All three are neutral, so the Gasteiger charges sum to ~0.
