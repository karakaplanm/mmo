# Folding Simulation with GROMACS

## Simulate a Small Protein (e.g., Villin Headpiece)


### Download a PDB file ( e.g., 2F4K.pdb from <a href="https://www.rcsb.org">https://www.rcsb.org</a> ).

```
PyMol> fetch 2F4K
PyMol> remove solvent
PyMol> select yoket, resi 70
PyMol> remove yoket 
PyMol > save 2F4K.pdb, 2F4K
```

### Prepare input files:

```console
# Convert PDB to GROMACS format
$ gmx pdb2gmx -f 2F4K.pdb -o protein.gro -water spc

# Generate simulation box
$ gmx editconf -f protein.gro -o box.gro -c -d 1.0

# Energy minimization
$ gmx grompp -f em.mdp -c box.gro -p topol.top -o em.tpr
$ gmx mdrun -v -deffnm em
```
