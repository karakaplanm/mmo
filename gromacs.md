# Gromacs

## Lysozyme In Water

### Step 1: Prepare Lysozyme protein molecule
- Open PyMOL and give fetch command
  ```console
  PyMOL > fetch 1AKI
  ```
- Remove water molecules<br>
  ```console
  PyMOL > remove solvent
  ```
- Save clean 1AKI molecule
  ```console
  PyMOL> save 1AKI_clean.pdb, 1AKI
  ```
- Execute pdb2gmx by issuing the following command:
  ```console
  $ gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
  ```
  type 15 at the command prompt, followed by 'Enter'.
  You have now generated three new files: 1AKI_processed.gro, topol.top, and posre.itp.
  You can check the files with;
  ```console
  $ ls
  ```
  
### Step 2: Examine the Topology

### Step 3: Defining the Unit Cell & Adding Solvent
- Creating a 1.0 nm box
  ```console
  $ gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
  ```
- Fill with solvent water
  ```console
  $ gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
  ```
### Step 4: Add Ions for neutralize
  ```console
  $ wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
  $ gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
  ```
  When prompted, choose group 13 "SOL" for embedding ions
  
### Step 5: Energy Minimization
  ```console
  $ wget http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
  $ gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
  $ gmx mdrun -v -deffnm em
  ```
  Plotting
  ```console
  $ gmx energy -f em.edr -o potential.xvg
  ```
  Type "10 0"
  ```console
  $ xmgrace potential.xvg
  
  ```
### Step 6: Equilibration
  ```console
  $ wget http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
  $ gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
  $ gmx mdrun -deffnm nvt
  ```
  Plot the results
  ```console
  $ gmx energy -f nvt.edr -o temperature.xvg
  $ xmgrace temperature.xvg
  ```
  
  
  
