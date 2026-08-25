axes location off
display resize 880 880
display projection orthographic
display depthcue off
display rendermode GLSL
color Display Background white
set POCKET "resid 150 222 238 257 260 261 264 286 287 288 290 291"

mol new rec.pdb type pdb waitfor all
mol delrep 0 top
# whole protein, quiet
mol representation NewCartoon 0.28 12 2.5 0
mol color ColorID 6
mol selection {protein}
mol material Transparent
mol addrep top
# the pocket lining, in full
mol representation Licorice 0.16 24 24
mol color Name
mol selection "$POCKET and noh"
mol material Opaque
mol addrep top

mol new pose1.pdb type pdb waitfor all
mol delrep 0 top
mol representation CPK 0.62 0.25 32 32
mol color Element
mol selection {all}
mol material Opaque
mol addrep top
mol representation VDW 0.55 32
mol color ColorID 3
mol selection {name Fe}
mol addrep top
color Element C gray

# frame the pocket, not the whole protein
mol top 1
display resetview
scale by 0.30
rotate y by 25
rotate x by -12
render TachyonInternal pocket.tga
quit
