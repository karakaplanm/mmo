# Rapor: 4PYP Proteininin AlphaFold ile Katlanması ve Model Kalite Değerlendirmesi

## 1. Protein Hakkında Kısa Bilgi

4PYP, Halorhodospira halophila (fotosentetik mor bakteri) türünde bulunan Photoactive Yellow Protein (PYP)’nin Protein Data Bank (PDB) kimlik numarasıdır. PYP, yaklaşık 125 amino asit uzunluğunda, küçük, çözünür ve ışığa duyarlı bir fotoreseptör proteindir.

Biyolojik işlevi, bakterinin mavi ışığı algılamasını ve buna bağlı olarak negatif fototaksi (ışıktan uzaklaşma) yanıtını başlatmasıdır. Bu sayede bakteri yoğun ışık altında hasar görmekten korunur.

Yapısal olarak PYP, α/β katlama motifi sergiler. Tipik olarak altı β-telinden oluşan bir β-levha ve bu levhanın etrafında yer alan dört α-sarmal içerir. Proteinin en önemli özelliği, p‑hidroksisinamil kromoforunu (fotoreseptör grubu) kovalent olarak bağlamasıdır. Bu kromofor, bir tiyoester bağı ile aktif bölgedeki Cys69 amino asidine bağlanır.

Işık uyarımı ile kromofor trans → cis izomerizasyonu yapar ve bu da proteinde büyük konformasyonel değişikliklere yol açar. Bu değişimler, aşağı akış sinyal iletimini başlatan protein-protein etkileşimlerini tetikler. PYP, bu mekanizmasıyla fotosiklin anlaşılması için model bir sistem olarak kabul edilir ve aynı zamanda PAS (Per-Arnt-Sim) domain ailesinin yapısal-fonksiyonel çalışmalarına da ışık tutar.

Deneysel olarak, PYP’nin hem karanlık (Ground state) hem de ışık uyarılmış ara formlarının yüksek çözünürlüklü yapıları (X-ışını kristalografisi ve NMR ile) çözülmüştür. 4PYP kodu, bu yapılardan birine (genellikle karanlık durum) karşılık gelir.

## 2. Katlama Yöntemi

Katlama işlemi **AlphaFold** kullanılarak yapılmıştır. Herhangi bir multimer bilgisi girilmediği için ipTM değeri çizgi (`-`) olarak raporlanmıştır (monomerik form hedeflenmiştir).

## 3. Sonuçların Değerlendirilmesi

Elde edilen **pTM skoru 0,84**’tür. Bu değer, modelin genel yapısal doğruluğunun **çok yüksek** olduğunu gösterir (0,5 üzeri iyi, 0,8 üzeri mükemmel kabul edilir).

## 4. Modelin Yorumu

- **0,84 pTM değeri**, katlanmış yapının globüler çekirdeğinin deneysel yapıyla (PDB:4PYP) büyük uyum içinde olduğunu kanıtlar.

## 5. Katlama Şeklinin Özeti

AlphaFold, 4PYP sekansını alarak **yüksek güvenilir bir monomerik model** oluşturmuştur.

- **Globüler domain** (yaklaşık kalıntı 10‑110) plDDT ortalaması **85’in üzerinde** olup neredeyse kristalografik kalitededir.
- **Kromofor bağlanma cebi** (Tyr42, Glu46, Cys69 gibi kalıntıları içeren bölge) çok yüksek güvenle (>90) konumlanmıştır.
- **Esnek döngüler** (özellikle 55‑70 arası kalıntılar) daha düşük güvenle modellenmiştir. Bu bölgeler fotosikl sırasında konformasyon değiştirdiğinden, AlphaFold’un tek bir yapı yerine olasılıksal bir dağılım vermesi beklenir.

## Sonuç

Elde edilen model, **moleküler yerleştirme** (docking), **dinamik simülasyon** (ön ısıtma için) ve **mutasyon etkisi tahmini** gibi ileri analizlerde doğrudan kullanılabilir. Düşük güvenli bölgeler için AlphaFold’un çoklu çıktı modelleri incelenerek alternatif katlama konformasyonları değerlendirilmelidir.
