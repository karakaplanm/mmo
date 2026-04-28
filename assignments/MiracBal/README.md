# Human CTLA-4 (3OSK) Protein Model (AlphaFold Prediction)

**Yazar:** Miraç Bal
**Kurum:** İnönü Üniversitesi - Moleküler Biyoloji ve Genetik Bölümü
**Ders:** Biyomoleküler Modelleme

---

## 🔬 Genel Bakış
Bu dosya, **AlphaFold** algoritması kullanılarak tahmin edilen **CTLA-4** (Sitotoksik T-Lenfosit İlişkili Protein 4) proteininin üç boyutlu yapısını içermektedir [cite: 1]. AlphaFold, amino asit dizilerinden protein yapılarını yüksek doğrulukla tahmin eden yapay zeka tabanlı bir sistemdir [cite: 2].

## 🧬 CTLA-4 Nedir?
* **Gen Adı:** CTLA4 [cite: 3]
* **Protein Tipi:** İmmün Kontrol Noktası Reseptörü [cite: 3]
* **Ekspresyon:** Öncelikle aktive edilmiş T hücrelerinde, özellikle düzenleyici T hücrelerinde (Tregs) bulunur [cite: 3].

### Temel Fonksiyonlar
* T hücresi aktivasyonunu inhibe eder [cite: 3].
* Aşırı immün yanıtları önleyerek immün homeostazı korur [cite: 3].

### Ligand Etkileşimi ve CD28 İlişkisi
* **Ligandlar:** CD80 (B7-1) ve CD86 (B7-2) [cite: 3].
* CTLA-4, CD28 ile aynı ligandlara bağlanır ancak CD28 aktivasyonu teşvik ederken, CTLA-4 inhibitör sinyal sağlar [cite: 3].
* CTLA-4, bu ligandlara karşı CD28'den daha yüksek afiniteye sahiptir [cite: 3].

---

## 🛠 Metodoloji
1. Protein seçildi ve Protein Data Bank (PDB) üzerinden takma adı ile aratıldı [cite: 3].
2. Protein **Fasta** formatında indirildi ve Google AlphaFold sistemine yüklendi [cite: 4].
3. Üretilen 5 farklı konformasyon arasından, en yüksek küresel güven sıralamasına sahip olan **model_0** referans yapı olarak seçildi [cite: 5].

---

## 📊 Yapısal Analiz ve Görselleştirme (pLDDT)
AlphaFold modelleri, güven skorlarına göre renklendirilmiştir [cite: 6]:

* 🔵 **Mavi (>90):** Yüksek doğruluk, stabil çekirdek yapılar [cite: 6].
* 🟡 **Siyan/Sarı (90-50):** Düşük güvenilirlik, esnek bağlayıcı (linker) bölgeler [cite: 6].
* 🟠 **Turuncu (<50):** %24 oranında düzensiz (disordered) uç bölgeler [cite: 7].

---

## 📂 Dosya İçeriği
* `fold_2026_03_30_11_26_model_0.cif`: Google AlphaFold'dan alınan 3D model [cite: 8].
* `fold_2026_03_30_11_26_summary_confidences_0.json`: Güven metrikleri [cite: 8].
* `shh_pymol_render.png`: PyMOL görselleştirmesi [cite: 8].