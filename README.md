# Biomolecular Modelling


## Course Content

@ Basic foundations and applications of molecular modeling.
@ Visualization of Biomolecules. Calculating physical properties of biomolecules. Theoretical background of molecular simulations.
@ Computational methods and algorithms of molecular computations. Current trends and software in molecular modeling.
@ Structure prediction of biomolecules and drug design.

## Weekly Contents
+ Basic foundations and applications of molecular modeling
+ Physical foundations of molecular modeling. Statistical mechanics. Quantum Mechanics
+ Introduction to Nucleic Acids Structure. Visualization of Nucleic Acids. Basic Calculations on Nucleic Acids
+ Structure of Proteins and Visualization of a Protein. Classification and basic calculations.
+ Structure prediction: Protein Folding. Homology Modeling
+ Molecular mechanics. Basic assumptions. Force field models.
+ Potential Energy Surfaces: Saddle Points, First-Order Methods, Second-Order Methods
+ Molecular Mechanics Examples. Geometry Optimization. Amino Acids
+ Electrostatics & Solvation in Biomolecules
+ Introduction to molecular editing and visualization software
+ Introduction to Drug Design
+ Introduction to Monte Carlo Methods
+ Ab Initio Methods
+ Applying biomolecular modeling to ongoing research
    

## Week 1

+ Do PyMol Installation <a href=https://pymol.org/> https://pymol.org/</a></li>
+ Get GitHub Repository account like <a href=https://github.com/karakaplanm> https://github.com/karakaplanm</a></li>
+ Visit the github page of this lecture <a href=https://github.com/karakaplanm/mmo>https://github.com/karakaplanm/mmo</a></li>
+ Get Gromacs <a href=https://gromacs.org>https://gromacs.org</a></li>
+ Get AutoDock Vina <a href=https://vina.scripps.edu>https://vina.scripps.edu/</a></li>
+ Install WSL (Windows Subsystem for Linux) to Windows

## Week 2

### PyMol

Download
https://pymol.org


Get Licence with Student/Teacher
https://pymol.org/buy.html

#### Basic PyMOL Interface and Commands
https://www.rcsb.org/

```
PyMOL> fetch 1crn

> color red, 1crn
> show stick
> hide cartoon
> show surface
```

### Rotating, Zooming, and Translating:

Use the mouse left button to rotate, middle button to zoom, and right button to translate


### Train this commands
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
