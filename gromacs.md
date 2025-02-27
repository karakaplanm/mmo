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

  
