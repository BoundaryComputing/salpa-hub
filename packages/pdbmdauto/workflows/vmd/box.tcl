# Figure 2 — the solvated, neutralised system: protein cartoon, water as faint lines,
# Na+ / Cl- as spheres, the periodic box outlined.
# Input: ion.gro, a copy of <working_dir>/pdbmdauto-e2e-full/e2e_4z8j/gmx/ion.gro
axes location off
display resize 1000 800
display projection orthographic
display depthcue off
display rendermode GLSL
color Display Background white
mol new ion.gro type gro waitfor all
mol delrep 0 top
# Water: one point per oxygen (a Lines rep of lone oxygens draws nothing — no bonds)
mol representation Points 2.0
mol color ColorID 9
mol selection {resname SOL and name OW}
mol addrep top
mol representation NewCartoon 0.35 12.0 4.1 0
mol color ColorID 2
mol selection {protein}
mol addrep top
# Na+ blue, Cl- green
mol representation VDW 1.0 20.0
mol color ColorID 0
mol selection {resname NA}
mol addrep top
mol representation VDW 1.0 20.0
mol color ColorID 7
mol selection {resname CL}
mol addrep top
# The periodic box, drawn by hand from the .gro box vector (5.0 nm cube = 50 A) — the pbc
# plugin is not loaded in -dispdev text mode, so `pbc box` cannot be relied on here.
set L 50.0
graphics top color black
foreach e {{0 0 0 1 0 0} {0 0 0 0 1 0} {0 0 0 0 0 1} {1 0 0 1 1 0} {1 0 0 1 0 1} {0 1 0 1 1 0} {0 1 0 0 1 1} {0 0 1 1 0 1} {0 0 1 0 1 1} {1 1 0 1 1 1} {1 0 1 1 1 1} {0 1 1 1 1 1}} {
  lassign $e ax ay az bx by bz
  graphics top line [list [expr {$ax*$L}] [expr {$ay*$L}] [expr {$az*$L}]] [list [expr {$bx*$L}] [expr {$by*$L}] [expr {$bz*$L}]] width 2
}
display resetview
rotate x by -65
rotate y by 25
scale by 1.25
render TachyonInternal box.tga
quit
