# Figure 1 — the model after Fix Missing Residues: what the experiment resolved (grey), the
# bound peptide (teal), and the six residues ProMod3 rebuilt (orange).
# Input: fixed.pdb, a copy of <working_dir>/pdbmdauto-e2e-full/e2e_4z8j/Merge/fixed.pdb
axes location off
display resize 1000 800
display projection orthographic
display depthcue off
display rendermode GLSL
color Display Background white
mol new fixed.pdb type pdb waitfor all
mol delrep 0 top
# PDZ domain, experimentally resolved part
mol representation NewCartoon 0.30 12.0 4.1 0
mol color ColorID 2
mol selection {chain A and resid > 5}
mol addrep top
# The peptide (chain B)
mol representation NewCartoon 0.30 12.0 4.1 0
mol color ColorID 10
mol selection {chain B}
mol addrep top
# Rebuilt residues: 4Z8J's A 33-37 (GSHGG) and B 586 (Gln), which ProMod3 renumbers to A 1-5 and B 1
mol representation Licorice 0.30 12.0 12.0
mol color ColorID 3
mol selection {(chain A and resid 1 to 5) or (chain B and resid 1)}
mol addrep top
mol representation NewCartoon 0.30 12.0 4.1 0
mol color ColorID 3
mol selection {chain A and resid 1 to 6}
mol addrep top
display resetview
rotate y by 35
rotate x by -20
scale by 1.9
render TachyonInternal fixed.tga
quit
