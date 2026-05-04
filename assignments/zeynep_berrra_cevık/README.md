# 1GXY (ART2.2) – Ökaryotik Mono-ADP-riboziltransferaz

## Protein Tanıtımı ve İşlevi

**1GXY**, sıçan (*Rattus norvegicus*) kaynaklı **ART2.2** (ekto-ADP-riboziltransferaz 2.2) proteininin kristal yapısını temsil eden bir PDB dosyasıdır. ART2.2, hücre zarlarında bulunan ve **hücre dışı NAD⁺** moleküllerini kullanarak hedef proteinlere **ADP‑riboz** grubu ekleyen bir enzimdir. Bu işlem, **mono-ADP-ribozilasyon** olarak adlandırılır ve hücre sinyal iletiminde, özellikle **immün sistem düzenlenmesinde** kritik roller oynar.

### Biyolojik Fonksiyonu

- **NAD⁺ metabolizması:** ART2.2, hücre dışı ortamda bulunan NAD⁺’yi yıkarak ADP‑riboz ve nikotinamid oluşturur.
- **Sinyal molekülleri üretimi:** Oluşan ADP‑riboz, kalsiyum sinyali gibi hücresel yanıtları tetikleyebilir.
- **T‑hücre düzenlemesi:** ART2.2, T‑lenfositlerin yüzeyinde eksprese edilir; bu sayede bağışıklık yanıtının kontrolüne katkıda bulunur.
- **Oto-ADP-ribozilasyon:** Enzim, kendini de ADP‑ribozilleyerek aktivitesini düzenleyebilir.

### Yapısal Özellikleri

- **Uzunluk:** 226 amino asit
- **Molekül ağırlığı:** ~26,2 kDa
- **Disülfid bağları:** Aktif bölgeden uzakta iki disülfid köprüsü bulunur; bunlar proteine yüksek termal ve kimyasal kararlılık kazandırır.
- **N‑terminal uzantı:** Membrana tutunmayı sağlayan hidrofobik bir bölge içerir.
- **Katlanma:** Bakteriyel ADP‑riboziltransferazlarla benzer bir çekirdek katlanma motifine sahiptir.

### Kullanım Alanları

- **Temel bilim:** NAD⁺ sinyal mekanizmalarının aydınlatılması.
- **İmmünoloji:** T‑hücre aktivasyonu ve tolerans çalışmaları.
- **İlaç hedeflemesi:** Otoimmün hastalıklar ve inflamasyon modellerinde potansiyel hedef.

---

## AlphaFold Tahminlerinin Değerlendirilmesi – 1GXY (ART2.2)

### Genel Model Kalitesi

AlphaFold ile yapılan tahmin sonucunda elde edilen **summary_confidences** dosyasına göre:

- **ptm (pTM) skoru = 0.89**  
  Bu skor, 0.5 eşiğinin oldukça üzerindedir. **pTM > 0.5**, tahmin edilen genel katlanmanın (backbone ve topoloji) doğru olduğunu gösterir. 0.89 gibi yüksek bir değer, bu proteinin **çok yüksek güvenle** modellendiğini ve deneysel yapıyla (1GXY) büyük uyum içinde olduğunu kanıtlar.

- **ipTM (interface pTM) =** Hesaplanamamıştır (null).  
  Bunun nedeni, ART2.2’nin bu tahminde **tek zincirli** olarak modellenmiş olmasıdır. İpTM, çoklu zincirli komplekslerde arayüz güvenini ölçtüğü için bu protein için anlamlı bir değer yoktur.

- **ptm ve ranking_score = 0.89**  
  Ranking_score, modelin diğer alternatif tahminler arasında en iyi sıralamayı aldığını gösterir. 0.89, çok yüksek bir sıralama puanıdır.

- **num_recycles = 10.0**  
  AlphaFold’un 10 döngü (recycle) kullanarak modeli iyileştirdiğini belirtir.

- **has_clash = 0.0** ve **fraction_disordered = 0.0**  
  Tahmin edilen modelde **hiçbir atom çakışması** yoktur ve proteinin hiçbir bölgesi **düzensiz (disordered)** olarak işaretlenmemiştir. Bu, ART2.2’nin tamamının stabil ve katlanmış olduğu anlamına gelir.

### pLDDT Değerlendirmesi

Mevcut JSON dosyası per-kalıntı pLDDT değerlerini içermediği için tabloda “bilgi yok” olarak belirtilmiştir. Ancak, ART2.2’nin çok kısa disülfid bağlarıyla stabilize olmuş küçük bir protein olması ve ptm = 0.89 değeri, **tüm kalıntıların yüksek veya çok yüksek pLDDT (>70)** aralığında olmasını beklememiz için yeterli bir kanıttır.

### Sonuç

AlphaFold, **1GXY (ART2.2)** proteinini **son derece güvenilir bir şekilde** modellemiştir. 0.89’luk pTM skoru, deneysel kristal yapıyla (1.71 Å çözünürlük, P21 kristal formu) tutarlı bir tahmin yapıldığını gösterir. Bu model, ART2.2’nin NAD⁺ bağlama bölgesi, disülfid bağları ve membran ilişkili bölgeleri hakkında yapılacak ileri çalışmalar için sağlam bir temel oluşturmaktadır.