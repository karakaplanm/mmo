# pymol ornek2.py
# veya
# PyMol > run ornek2.py

from pymol import cmd

def dna_setup():
    cmd.reinitialize()
    cmd.fetch("1d12")
    cmd.show("sticks", "resn DM2")
    cmd.color("red", "resn DM2")

dna_setup()
