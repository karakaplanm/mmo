#!/usr/bin/bash

echo "Step 1: Prepare Lysozyme protein molecule"
echo "15" | gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce

echo "Step 2: Defining the Unit Cell & Adding Solvent"
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top

echo "Step 3: Examine the Topology"
pymol 1AKI_processed.gro
pymol 1AKI_solv.gro

echo "Step 4: Add Ions for neutralize"
#wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
echo "13" | gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

echo "Step 5: Energy Minimization"
#wget http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

echo "10 0" | gmx energy -f em.edr -o potential.xvg
xmgrace potential.xvg

echo "Step 6: Equilibration, Equilibration, Part 1"

#wget http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -nb gpu

echo "16 0" | gmx energy -f nvt.edr -o temperature.xvg
xmgrace temperature.xvg


echo "Step 7: Equilibration, Part 2"
#wget http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -nb gpu

echo "18 0" | gmx energy -f npt.edr -o pressure.xvg
xmgrace pressure.xvg

echo "24 0" | gmx energy -f npt.edr -o density.xvg
xmgrace density.xvg



echo "Step 8: Production MD"

#wget http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1 -nb gpu

echo "Step 9: Analysis"
echo "1 0 " | gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center

echo "4 4" | gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
xmgrace rmsd.xvg
echo "1 1" | gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
xmgrace rmsd_xtal.xvg

echo "1" | gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg
xmgrace gyrate.xvg

echo "Step 10: Visualize vith PyMol"
echo "1" | gmx trjconv -f md_0_1.xtc  -s md_0_1.gro  -o trajectory.pdb
pymol trajectory.pdb














