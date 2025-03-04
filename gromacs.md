# Gromacs

## Installation

### with apt from Ubuntu repos
```console
$ sudo apt update
$ sudo apt install gromacs
```

### Installation from Source

```console
$ tar xfz gromacs-2025.0.tar.gz
$ cd gromacs-2025.0
$ mkdir build
$ sudo apt install cmake build-essential
$ cd build
$ cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
$ make
$ make check
$ sudo make install
$ source /usr/local/gromacs/bin/GMXRC
```

Add following options for typical installation

`-DGMX_MPI=on` to build using MPI support

`-DGMX_GPU=CUDA` to build with NVIDIA CUDA support enabled.

Example cmake command
```console
$ cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DCMAKE_CXX_COMPILER=gcc-12 -DCMAKE_C_COMPILER=gcc-12
$ make -j8
```

## Lysozyme In Water

### Step 1: Prepare Lysozyme protein molecule
Open PyMOL and give fetch command
```console
PyMOL > fetch 1AKI
```

Remove water molecules<br>
```console
PyMOL > remove solvent
```
Save clean 1AKI molecule
```console
PyMOL> save 1AKI_clean.pdb, 1AKI
```
Execute pdb2gmx by issuing the following command:
```console
$ gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
```
type 15 at the command prompt, followed by 'Enter'.
You have now generated three new files: 1AKI_processed.gro, topol.top, and posre.itp.
You can check the files with;
```console
$ ls
```

### Step 2: Defining the Unit Cell & Adding Solvent
Creating a 1.0 nm box
```console
$ gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
```
Fill with solvent water
```console
$ gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```
### Step 3: Examine the Topology
Open `.gro` files
```console
$ pymol 1AKI_processed.gro
$ pymol 1AKI_solv.gro
$ pymol 1AKI_solv_ions.gro
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

### Step 6: Equilibration, Equilibration, Part 1
```console
$ wget http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
$ gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$ gmx mdrun -deffnm nvt
```
Plot the results
```console
$ gmx energy -f nvt.edr -o temperature.xvg
```
Type "16 0" at the prompt
  
```console
$ xmgrace temperature.xvg
```

### Step 7: Equilibration, Part 2
```console
$ wget http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
$ gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
$ gmx mdrun -deffnm npt
```
Plot
```console
$ gmx energy -f npt.edr -o pressure.xvg
Type "18 0"

$ xmgrace pressure.xvg

$ gmx energy -f npt.edr -o density.xvg
Type "24 0"
$ xmgrace density.xvg
```

### Step 8: Production MD

```console
$ wget http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp
$ gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
$ gmx mdrun -deffnm md_0_1
```

### Step 9: Analysis

```console
$ gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
```

```console
$ gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
```
Choose 4 ("Backbone")
```console
$ xmgrace rmsd.xvg
```

```console
$ gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
```

```console
$ xmgrace rmsd_xtal.xvg
```

```console
$ gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg
```
Choose group 1 (Protein) for analysis.

  
  
