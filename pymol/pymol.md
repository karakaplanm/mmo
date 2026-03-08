# Using PyMol


fetch 1bna
remove solvent

hide everthing
show sticks
show spheres
set sphere_scale, 0.25

ray 2400, 1800
png dna_modeli.png, dpi=300


set ray_trace_mode, 1
ray 2400, 1800
png dna_modeli.png, dpi=300

color red, resn DA
color blue, resn DC


hide everthing
show cartoon

reinitialize
fetch 1bna
show surface, 1bna
set transparency, 0.5



reinitialize
fetch 261d


reinitialize
fetch 1d12
select ilac, resn DM2
show spheres, ilac
color yellow, ilac


