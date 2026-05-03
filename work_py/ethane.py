import numpy as np
import matplotlib.pyplot as plt

# Dihedral açı 0'dan 360 dereceye (Tam tur)
theta = np.linspace(0, 360, 500)
theta_rad = np.deg2rad(theta)

# Etan için Enerji Bariyeri (yaklaşık 12.5 kJ/mol veya 2.9 kcal/mol)
V0 = 12.5

# Enerji fonksiyonu: E(theta) = (V0 / 2) * (1 + cos(3 * theta))
# Bu formül 0, 120, 240, 360 derecelerde max (eclipsed) 
# 60, 180, 300 derecelerde min (staggered) değer verir.
energy = (V0 / 2) * (1 + np.cos(3 * theta_rad))

plt.figure(figsize=(10, 6))
plt.plot(theta, energy, color='teal', linewidth=2.5)

# Önemli Noktaların İşaretlenmesi
plt.scatter([0, 120, 240, 360], [V0]*4, color='red', label='Eclipsed (Çakışık)')
plt.scatter([60, 180, 300], [0]*3, color='blue', label='Staggered (Çapraz)')

# Grafik Düzenlemeleri
plt.title('Etanda C-C Bağı Etrafında Dönme ve Konformasyonel Enerji Analizi', fontsize=14)
plt.xlabel('Dihedral Açı (Derece $^\circ$)', fontsize=12)
plt.ylabel('Potansiyel Enerji (kJ/mol)', fontsize=12)
plt.xticks(np.arange(0, 361, 60))
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

# Açıklamalar
plt.annotate('Maksimum Enerji (Eclipsed)', xy=(0, V0), xytext=(25, V0+0.5),
             arrowprops=dict(facecolor='black', shrink=0.05, width=1))
plt.annotate('Minimum Enerji (Staggered)', xy=(60, 0), xytext=(85, 1),
             arrowprops=dict(facecolor='black', shrink=0.05, width=1))

plt.tight_layout()
plt.show()
#plt.savefig('etan_konformasyon_analizi.png')

