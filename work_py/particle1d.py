import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e

# Box parameters
L = 1e-9  # Box length (1 nm)
n_states = 4  # Number of quantum states to plot
x = np.linspace(0, L, 1000)  # Position grid

def wavefunction(n, x, L):
    return np.sqrt(2/L) * np.sin(n * np.pi * x / L)

def energy(n, L, m=m_e):
    return (n**2 * np.pi**2 * hbar**2) / (2 * m * L**2) / 1.6e-19  # Convert to eV

# Plot wavefunctions and energies
plt.figure(figsize=(10, 6))
for n in range(1, n_states + 1):
    psi = wavefunction(n, x, L)
    E = energy(n, L)
    plt.plot(x, psi + E, label=f'n = {n}, E = {E:.2f} eV')
    
plt.xlabel('Position (m)')
plt.ylabel('Energy (eV) / Ψ(x)')
plt.title('Wavefunctions for a Particle in a 1D Box')
plt.legend()
plt.grid(True)
plt.show()
