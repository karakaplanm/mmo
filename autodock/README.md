# AutoDock


https://ccsb.scripps.edu/mgltools/downloads/


AutoDock is a suite of automated docking tools designed to predict how small molecules, such as substrates or drug candidates, bind to a receptor of known 3D structure. It is widely used in the field of molecular modeling and drug design. The software is particularly valuable for virtual screening, where it helps researchers identify potential drug candidates by simulating how they interact with target proteins.

## Download and Install

https://ccsb.scripps.edu/mgltools/downloads/

```
$ cd mgltools_x86_64Linux2_X.X.X
$ /install.sh
$ export PATH=~/mgltools_x86_64Linux2_1.5.7/bin:$PATH
$ cd mmo/autodock/
$ adt
```


## Ligand Preparation Steps

- File -> Read Molecule -> load jz4.pdb
- Ligand -> Input -> Choose -> jz4.pdb -> Open
- Ligand -> Torsion Tree -> Detect Root
- Ligand -> Torsion Tree -> Choose Torsions
- Ligand -> Output -> Save pdbqt (ligand.pdbqt)
- Close the adt GUI


## Receptor Steps

- Open adt
- File -> Read Molecule -> load 3HTB_clean.pdb
- Edit -> Hydrogens -> Add -> Polar Only 
- Grid -> Macromolecule -> Choose -> protein -> Select Molecule -> Save As ->protein.pdbqt
- Grid -> Set Map Types -> Open Ligand -> ligand.pdbqt -> Open 
- Grid -> Grid Box Where we select the docking area
- Grid Box -> Close Saving Curent 
- Grid -> Output -> Save mygrid.gpf 
- Docking -> Macromolecule -> Set Rigid Filename -> protein.pdbqt
- Docking -> Ligand -> Choose -> ligand.pdbqt -> Accept
- Docking -> Search Parameter -> Genetic Algorithm -> Accept
- Docking -> Docking Parameters -> Accept
- Docking -> Other Options -> AutoDock 4.2 Parameters -> Accept
- Docking -> Ouput -> Lamarkian GA -> Docking results

## Docking steps Command line

```
autogrid4 -p mydock.gpf -l protein.dlg
autodock4 -p mydock.dpf -l ligand.dlg
```

## Analysis Docking results

```
pymol ligand.dlg
```
