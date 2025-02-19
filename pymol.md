# PyMol


## Download
https://pymol.org

## Install
   `$ sudo apt install pymol`

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
