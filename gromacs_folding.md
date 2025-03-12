# Folding Simulation with GROMACS

## Simulate a Small Protein (e.g., Villin Headpiece)


### Download a PDB file (e.g., 2F4K.pdb from <a href="https://www.rcsb.org">https://www.rcsb.org</a>).

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
