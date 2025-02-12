# Biomolecular Modelling

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
