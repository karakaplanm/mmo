# Folding Simulation with GROMACS

## Simulate a Small Protein (e.g., Villin Headpiece)


### Download a PDB file ( e.g., 2F4K.pdb from <a href="https://www.rcsb.org">https://www.rcsb.org</a> ).

```
PyMol> fetch 1AKI
PyMol> remove solvent
PyMol > save 1AKI.pdb, 1AKI
```

### Prepare input files:

```console
# Convert PDB to GROMACS format
$ gmx pdb2gmx -f 1AKI.pdb -o protein.gro -water spc
# Enter 15

# Generate simulation box
$ gmx editconf -f protein.gro -o box.gro -c -d 1.0

# Energy minimization
$ wget http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
$ gmx grompp -f minim.mdp -c box.gro -p topol.top -o em.tpr
$ gmx mdrun -v -deffnm em
```
