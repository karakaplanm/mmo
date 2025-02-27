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
  PyMOL> save clean_1AKI.pdb, 1AKI
  ```

  
