# Install psi4: !pip install psi4
import psi4
import matplotlib.pyplot as plt
import numpy as np

# Set up H₂ molecule
psi4.set_memory('2 GB')
psi4.set_num_threads(2)
h2 = psi4.geometry("""
0 1
H
H 1 0.74
""")

# Compute orbitals with Hartree-Fock
psi4.set_options({'basis': 'sto-3g', 'scf_type': 'pk'})
e, wfn = psi4.energy('scf', return_wfn=True)

# Plot the first molecular orbital
psi4.cubeprop(wfn)
grid = np.loadtxt('Psi_a_1_1-A.cube', skiprows=6)  # Load cube file
plt.contourf(grid[:, 0], grid[:, 1], grid[:, 3], cmap='RdBu')
plt.colorbar()
plt.title('H₂ Molecular Orbital (σ bonding)')
plt.show()
