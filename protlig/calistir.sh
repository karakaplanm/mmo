#!/usr/bin/bash

printf "1\n1\n1\n0\n" | gmx pdb2gmx -f 3HTB_clean.pdb -o 3HTB_processed.gro -ter

# open JZ4.pdb  with Avogadro and add Build -> Hydrogens -> Add Hydrogens -> save as JZ4.mol2

# Run cgenff_charmm2gmx_py3_nx2.py
# python3 cgenff_charmm2gmx_py3_nx2.py JZ4.mol2 charmm36_tip3p.prm topol.top -o JZ4.itp
#python3 cgenff_charmm2gmx_py3_nx2.py JZ4 JZ4.mol2 jz4.str charmm36-feb2026_cgenff-5.0.ff

gmx editconf -f JZ4_gmx.pdb -o jz4.gro

# Properly merge protein and ligand into complex.gro
N_PROT=$(sed -n '2p' 3HTB_processed.gro)
N_LIG=$(sed -n '2p' jz4.gro)
N_TOT=$((N_PROT + N_LIG))

head -n 1 3HTB_processed.gro > complex.gro
echo "$N_TOT" >> complex.gro
tail -n +3 3HTB_processed.gro | head -n -1 >> complex.gro
tail -n +3 jz4.gro | head -n -1 >> complex.gro
tail -n 1 3HTB_processed.gro >> complex.gro

# Modify topol.top correctly
if ! grep -q "JZ4_gmx.itp" topol.top; then
    # Add ligand dihedral types right after forcefield
    sed -i '/#include ".*forcefield\.itp"/a \
\n[ dihedraltypes ]\n;      i        j        k        l  func         phi0         kphi  mult\n  CG2R61    CG321    CG321    CG331     9     0.000000     0.167360     3\n' topol.top

    # Insert ligand topology include
    sed -i 's/; Include water topology/; Include ligand topology\n#include "JZ4_gmx.itp"\n\n; Include water topology/g' topol.top
    
    # Append the ligand to the molecules list
    echo "JZ4                 1" >> topol.top
fi

gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Energy Minimization
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em


# Step Six: Equilibration
printf "0 & ! a H*\nq\n" | gmx make_ndx -f jz4.gro -o index_jz4.ndx


echo "3" | gmx genrestr -f jz4.gro -n index_jz4.ndx -o posre_jz4.itp -fc 1000 1000 1000

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

gmx mdrun -deffnm nvt


# Step Seven: Equilibration, Part 2
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr

gmx mdrun -deffnm npt

#Step Eight: Production MD

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_10.tpr

gmx mdrun -deffnm md_0_10



printf "1\n1\n1\n0\n" | gmx trjconv -f md_0_10.xtc -s md_0_10.tpr -o trajectory.pdb -pbc mol -center
pymol trajectory.pdb