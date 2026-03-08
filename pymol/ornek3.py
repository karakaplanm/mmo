import pymol
from pymol import cmd

# PyMOL'u arayüzsüz (headless) veya arayüzlü başlat
# 'gui' dersen pencere açılır, 'quiet' dersen sadece arka planda çalışır
#pymol.finish_launching(['pymol', '-qc']) 
pymol.finish_launching(['pymol', '-A3']) 

cmd.fetch("1d12")
cmd.select("ilac", "resn DM2")
cmd.show("spheres", "ilac")

# Görüntüyü kaydet ve çık
cmd.ray(1200, 1200)
cmd.png("dna_cikti.png")
#cmd.quit()
