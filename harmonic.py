import numpy as np
import matplotlib.pyplot as plt

# Parameters
k = 1.0        # Force constant (eV/Å²)
mass = 12.0     # Mass (u)
dt = 0.1        # Time step (fs)
n_steps = 1000  # Simulation steps

# Initial conditions
x0 = 0.5        # Initial displacement (Å)
v0 = 0.0        # Initial velocity (Å/fs)

# Arrays to store results
positions = np.zeros(n_steps)
velocities = np.zeros(n_steps)

# Velocity Verlet algorithm
x = x0
v = v0
for i in range(n_steps):
    a = -k * x / mass  # Acceleration (F = -kx)
    x += v * dt + 0.5 * a * dt**2
    a_new = -k * x / mass
    v += 0.5 * (a + a_new) * dt
    positions[i] = x
    velocities[i] = v

# Plot trajectory
plt.plot(positions)
plt.xlabel('Time Step')
plt.ylabel('Position (Å)')
plt.title('1D Harmonic Oscillator MD Simulation')
plt.grid(True)
plt.show()
