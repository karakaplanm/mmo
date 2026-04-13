import numpy as np
import matplotlib.pyplot as plt

# Sabitler
r_e = 1.012  # Denge bağ uzunluğu (Angstrom)
D_e = 450.0  # N-H bağ enerjisi (kJ/mol)
a = 2.2      # Potansiyel eğrilik parametresi
barrier_height = 24.2 # Amonyak inversiyon bariyeri (kJ/mol)

# 1. Morse Potansiyeli (Bağ Gerilmesi)
def morse_potential(r):
    return D_e * (1 - np.exp(-a * (r - r_e)))**2

# 2. Inversiyon Potansiyeli (Çift Kuyu Modeli)
# Azotun H3 düzlemine olan mesafesi x
def inversion_potential(x):
    x_0 = 0.38  # Denge mesafesi (Angstrom)
    return barrier_height * ((x/x_0)**2 - 1)**2

# Veri setleri oluşturma
r_range = np.linspace(0.5, 3.0, 500)
x_range = np.linspace(-0.8, 0.8, 500)

energy_morse = morse_potential(r_range)
energy_inv = inversion_potential(x_range)

# Grafikleştirme
plt.figure(figsize=(12, 5))

# Sol Grafik: Bağ Gerilmesi
plt.subplot(1, 2, 1)
plt.plot(r_range, energy_morse, 'b-', linewidth=2)
plt.axvline(r_e, color='r', linestyle='--', label=f'Denge ($r_e$={r_e}Å)')
plt.title('Amonyak N-H Bağ Gerilme Enerjisi')
plt.xlabel('Bağ Uzunluğu (Å)')
plt.ylabel('Potansiyel Enerji (kJ/mol)')
plt.legend()
plt.grid(alpha=0.3)

# Sağ Grafik: Şemsiye İnversiyonu
plt.subplot(1, 2, 2)
plt.plot(x_range, energy_inv, 'g-', linewidth=2)
plt.axhline(barrier_height, color='orange', linestyle='--', label=f'Bariyer ({barrier_height} kJ/mol)')
plt.title('Amonyak Şemsiye İnversiyonu Enerjisi')
plt.xlabel('Azotun Düzleme Mesafesi (Å)')
plt.ylabel('Potansiyel Enerji (kJ/mol)')
plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
#plt.savefig('amonyak_enerji_grafigi.png')

