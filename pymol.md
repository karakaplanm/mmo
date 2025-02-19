# PyMol

PyMOL is an open-source molecular visualization software widely used in structural biology, computational chemistry, and drug design. It allows users to visualize molecular structures, analyze interactions, and create high-quality molecular graphics.
Key Features of PyMOL

- **3D Molecular Visualization:** Displays proteins, nucleic acids, and small molecules in various representations (stick, cartoon, surface, etc.).
- **Publication-Quality Images:** Generates high-resolution molecular graphics.
- **Molecular Editing:** Modifies structures by adding or removing atoms, bonds, or ligands.
- **Protein-Ligand Interaction Analysis:** Highlights hydrogen bonds, hydrophobic contacts, and binding pockets.
- **Trajectory Visualization:** Supports molecular dynamics trajectory visualization from tools like GROMACS, AMBER, or CHARMM.
- **Scripting and Automation:** Uses Python scripting for batch processing and automation.
- **Surface and Electrostatic Potential Mapping:** Visualizes molecular surfaces and electrostatic potentials.

## Download
https://pymol.org

## Install
   `$ sudo apt install pymol`

Get Licence with Student/Teacher
https://pymol.org/buy.html

## Basic PyMOL Interface and Commands
https://www.rcsb.org/

```
PyMOL> fetch 1crn

> color red, 1crn
> show stick
> hide cartoon
> show surface
```

## Rotating, Zooming, and Translating:

Use the mouse left button to rotate, middle button to zoom, and right button to translate


## Train this commands
```
fetch 1crn
show cartoon
color blue, 1crn

Show disulfide bonds

select disulfides, resn CYS
show sticks, disulfides
color yellow, disulfides

Show Surface 

show surface
set transparency, 0.3
```
