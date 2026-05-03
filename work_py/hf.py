# Install PySCF first: !pip install pyscf
from pyscf import gto, scf

# Define the H₂ molecule
mol = gto.M(
    atom = 'H 0 0 0; H 0 0 0.74',  # Bond length = 0.74 Å
    basis = 'sto-3g'                # Minimal basis set
)

# Run Hartree-Fock
mf = scf.RHF(mol)
mf.kernel()

print(f"Hartree-Fock Energy: {mf.e_tot:.4f} Hartree")
