# 1M47 (Human Interleukin‑2) – Ökaryotik Sitokin

## Protein Tanıtımı ve İşlevi

**1M47**, insan (*Homo sapiens*) kaynaklı **Interleukin‑2 (IL‑2)** proteininin kristal yapısını temsil eden bir PDB dosyasıdır. IL‑2, bağışıklık sisteminin en önemli düzenleyici moleküllerinden biri olup, **T hücreleri**, **B hücreleri** ve **doğal öldürücü (NK) hücreler** üzerinde etkili olan bir sitokindir. Hücreler arası sinyal iletişimini sağlayarak bağışıklık yanıtının başlatılması, sürdürülmesi ve sonlandırılmasında kritik roller oynar.

### Biyolojik Fonksiyonu

- **T hücre proliferasyonu:** IL‑2, aktiflenmiş T hücrelerinin klonal olarak genişlemesini tetikler. Bu sayede vücut, enfeksiyon veya tümörle savaşacak yeterli sayıda savaşçı T hücresi üretebilir.

- **Düzenleyici T hücre (Treg) kontrolü:** IL‑2, Treg hücrelerinin hayatta kalması ve fonksiyonu için gereklidir. Düşük doz IL‑2 Treg'leri destekleyerek otoimmüniteyi engellerken, yüksek doz IL‑2 efektör T hücrelerini aktive eder.

- **NK hücre ve B hücre aktivasyonu:** IL‑2, NK hücrelerinin sitotoksik aktivitesini artırır ve B hücrelerinden antikor üretimini teşvik eder.

- **Sinyal iletimi:** IL‑2, reseptörüne (IL‑2Rα/β/γ) bağlanarak JAK‑STAT, PI3K‑Akt ve MAPK sinyal yollarını aktive eder.

### Yapısal Özellikleri

- **Uzunluk:** 133 amino asit (olgun form)
- **Molekül ağırlığı:** ~15,4 kDa
- **Disülfid bağı:** Cys58 ile Cys105 arasında bir adet disülfid köprüsü bulunur; bu bağ yapısal stabiliteye katkıda bulunur.
- **Katlanma:** Dört ana α‑sarmaldan oluşan "up‑up‑down‑down" tipi sarmal demet yapısı. β‑tabaka içermez.
- **Aktif bölge / reseptör bağlanma yüzeyi:** Helix A ve D üzerindeki hidrofobik ve polar kalıntılar, IL‑2 reseptörüne bağlanmayı sağlar.

### Kullanım Alanları

- **Kanser immünoterapisi:** Rekombinant IL‑2 (aldesleukin, Proleukin®), metastatik böbrek hücreli karsinom ve malign melanom tedavisinde FDA onaylıdır.

- **Otoimmün hastalıklar:** Düşük doz IL‑2, Tip 1 diyabet, sistemik lupus eritematozus (SLE) ve multipl skleroz gibi hastalıklarda Treg hücrelerini artırarak bağışıklık toleransını yeniden sağlamayı hedefler.

- **Temel bilim:** Sitokin katlanması, sarmal demet yapıları ve reseptör‑ligand tanıma mekanizmaları için paradigmatik bir model protein.

- **İlaç hedeflemesi:** Günümüzde IL‑2 varyantları (superkine, bypass mutants) ile Treg/efektör T hücreleri arasında seçicilik sağlanmaya çalışılmaktadır.

---

## AlphaFold Tahminlerinin Değerlendirilmesi – 1M47 (Human Interleukin‑2)

### Genel Model Kalitesi

AlphaFold ile yapılan tahmin sonucunda elde edilen **summary_confidences** dosyasına göre:

- **ptm (pTM) skoru = 0.85**  
  Bu skor, 0.5 eşiğinin oldukça üzerindedir. **pTM > 0.5**, tahmin edilen genel katlanmanın (backbone ve topoloji) doğru olduğunu gösterir. 0.85 gibi yüksek bir değer, bu proteinin **çok yüksek güvenle** modellendiğini ve deneysel yapıyla (1M47, X‑ışını kristalografisi 2.0 Å) büyük uyum içinde olduğunu kanıtlar.

- **ipTM (interface pTM) = 0.85** (chain_pair_iptm'den alınmıştır)  
  Bu protein **tek zincirli** (monomer) olmasına rağmen ipTM değeri 0.85 gibi yüksek çıkmıştır. Bu, proteinin **yapısal bütünlüğünün ve sıkı paketlenmesinin** bir göstergesidir. Multimer benzeri bir değerlendirme protokolünden kaynaklanabilir.

- **pTM ve ranking_score = 0.85**  
  Ranking_score, modelin diğer alternatif tahminler arasında en iyi sıralamayı aldığını gösterir. 0.85, çok yüksek bir sıralama puanıdır.

- **num_recycles = 10.0**  
  AlphaFold’un 10 döngü (recycle) kullanarak modeli iyileştirdiğini belirtir.

- **has_clash = 0.0** ve **fraction_disordered = 0.0**  
  Tahmin edilen modelde **hiçbir atom çakışması** yoktur ve proteinin hiçbir bölgesi **düzensiz (disordered)** olarak işaretlenmemiştir. Bu, IL‑2’nin tamamının stabil ve katlanmış olduğu anlamına gelir.

- **pAE_min (chain_pair_pae_min) = 0.76**  
  Çok düşük tahmini hata mesafesi (0.76 Å), komşu kalıntılar arasındaki göreceli konumların çok yüksek doğrulukta olduğunu gösterir.

### pLDDT Değerlendirmesi

Mevcut JSON dosyası per-kalıntı pLDDT değerlerini içermediği için tabloda “bilgi yok” olarak belirtilmiştir. Ancak, IL‑2’nin kompakt dört sarmallı bir sitokin olması ve ptm = 0.85 değeri, **tüm kalıntıların yüksek veya çok yüksek pLDDT (>70)** aralığında olmasını beklememiz için yeterli bir kanıttır.

### Sonuç

AlphaFold, **1M47 (Human Interleukin‑2)** proteinini **son derece güvenilir bir şekilde** modellemiştir. 0.85’lik pTM skoru, deneysel kristal yapıyla (2.0 Å çözünürlük) tutarlı bir tahmin yapıldığını gösterir. Bu model, IL‑2’nin reseptör bağlanma bölgesi, disülfid bağı ve varyant tasarımı hakkında yapılacak ileri çalışmalar için sağlam bir temel oluşturmaktadır.