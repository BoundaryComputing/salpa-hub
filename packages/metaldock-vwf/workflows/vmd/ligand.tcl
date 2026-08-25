# Ferrocene alone — CPK, iron in orange.
axes location off
display resize 880 880
display projection orthographic
display depthcue off
display rendermode GLSL
color Display Background white
mol new ferrocene.xyz type xyz waitfor all
mol delrep 0 top
mol representation CPK 0.55 0.22 40 40
mol color Element
mol selection {all}
mol addrep top
mol representation CPK 0.85 0.0 40 40
mol selection {name Fe}
mol color ColorID 3
mol addrep top
color Element C gray
color Element H white
display resetview
rotate x by -75
scale by 1.35
render TachyonInternal ligand.tga
quit
