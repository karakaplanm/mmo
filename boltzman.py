import numpy as np
import matplotlib.pyplot as plt

# Constants
kB = 8.617e-5  # eV/K (Boltzmann constant)
T = 300         # Temperature (K)
energies = np.linspace(0, 0.5, 100)  # Energy states (eV)

# Boltzmann probabilities
probabilities = np.exp(-energies / (kB * T)) / np.sum(np.exp(-energies / (kB * T)))

# Plot
plt.plot(energies, probabilities, 'r-', linewidth=2)
plt.xlabel('Energy (eV)')
plt.ylabel('Probability')
plt.title('Boltzmann Distribution at T = 300 K')
plt.grid(True)
plt.show()
