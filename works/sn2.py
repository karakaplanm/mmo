import numpy as np
import matplotlib.pyplot as plt

# Reaksiyon koordinatı (-1'den 1'e, girenlerden ürünlere)
x = np.linspace(-1, 1, 500)

# Potansiyel Enerji Eğrisi Oluşturma (Aktivasyon enerjisi ve Reaksiyon ısısı modellemesi)
# Girenler ve ürünler arasındaki enerji farkı (Ekzotermik reaksiyon varsayımı)
delta_h = -20  # Reaksiyon entalpisi (kJ/mol)
activation_barrier = 80 # Aktivasyon enerjisi (kJ/mol)

# Enerji profili için matematiksel fonksiyon (Asimetrik bir çan eğrisi)
def sn2_energy_profile(x, barrier, dh):
    # Temel çan eğrisi
    energy = barrier * np.exp(-4 * x**2)
    # Ürün tarafını (x > 0) aşağı çekerek ekzotermik karakter kazandırma
    energy += (dh / 2) * (1 + np.tanh(2 * x))
    return energy

y = sn2_energy_profile(x, activation_barrier, delta_h)

# Grafikleştirme
plt.figure(figsize=(10, 6))
plt.plot(x, y, color='darkred', linewidth=3)

# Önemli noktaları işaretleme
plt.annotate('Girenler\n(Cl⁻ + CH₃Br)', xy=(-0.9, y[25]), xytext=(-1.2, y[25]+10),
             arrowprops=dict(arrowstyle="->", color='black'))

plt.annotate('Geçiş Hali (TS)\n[Cl···CH₃···Br]⁻', xy=(0, activation_barrier-5), xytext=(-0.2, activation_barrier+10),
             weight='bold', color='blue')

plt.annotate('Ürünler\n(CH₃Cl + Br⁻)', xy=(0.9, y[-25]), xytext=(0.7, y[-25]-20),
             arrowprops=dict(arrowstyle="->", color='black'))

# Eksenleri ve Görünümü Düzenleme
plt.title('$S_N2$ Reaksiyonu Enerji Profili: $Cl^- + CH_3Br$', fontsize=14)
plt.xlabel('Reaksiyon Koordinatı (İlerleme)', fontsize=12)
plt.ylabel('Potansiyel Enerji (kJ/mol)', fontsize=12)
plt.xticks([]) # Koordinat birimsizdir
plt.grid(alpha=0.2)
plt.ylim(min(y)-30, activation_barrier+30)

plt.tight_layout()
plt.show()

