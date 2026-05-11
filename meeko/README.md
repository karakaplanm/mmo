# Molecular Docking with Meeko and AutoDock Vina

## Installation
```bash
pip install meeko --break-system-packages
pip install gemmi --break-system-packages
pip install pyproj --break-system-packages
pip install prody --break-system-packages
pip install pandas-stubs --break-system-packages
```

## Ligand Preparation
The ligand `jz4.mol2` is prepared using Meeko:
```bash
mk_prepare_ligand.py -i jz4.mol2 -o jz4.pdbqt
```

## Receptor Preparation
The unmodified receptor structure was downloaded directly from RCSB PDB. To prevent the ligand from docking into an already occupied pocket, we remove the native JZ4 ligand and other heteroatoms (like water and buffers) by keeping only the `ATOM` records. Finally, alternate locations are addressed automatically:
```bash
curl -s -o 3htb_rcsb.pdb https://files.rcsb.org/download/3HTB.pdb
# Clean the PDB to remove the native JZ4 ligand and other HETATM records
python -c "
with open('3htb_rcsb.pdb') as f, open('3htb_clean.pdb', 'w') as w:
    for line in f:
        if line.startswith('ATOM') or line.startswith('TER') or line.startswith('END'):
            w.write(line)
"
mk_prepare_receptor.py -i 3htb_clean.pdb -o 3htb -p -a --default_altloc A
```
This generates the `3htb.pdbqt` file required by AutoDock Vina.

## Molecular Docking (AutoDock Vina)
The grid box was centered around the original crystallographic coordinates of the JZ4 ligand (`x=22.7142`, `y=-25.2627`, `z=-3.283`), and an exploration box size of 20x20x20 Ångströms was used:
```bash
vina --receptor 3htb.pdbqt --ligand jz4.pdbqt \
     --center_x 22.7142 --center_y -25.2627 --center_z -3.283 \
     --size_x 20 --size_y 20 --size_z 20 \
     --out jz4_docked.pdbqt
```

The resulting top 9 predicted binding poses and their calculated affinities (top score: -6.948 kcal/mol) are saved in `jz4_docked.pdbqt`.

## Visualization
To visualize the docked complex in PyMOL, you can load both the receptor and the docking output together:
```bash
pymol 3htb.pdbqt jz4_docked.pdbqt
```
*(Alternatively, you can also load `3htb_clean.pdb` instead of `3htb.pdbqt` if you prefer to see the standard PDB format for the receptor.)*
In PyMOL, you can use the left/right arrow keys at the bottom right corner of the screen to cycle through the 9 different docked poses (states) of the ligand.