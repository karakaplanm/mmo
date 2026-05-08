# Transmembran Protein Modeli (AlphaFold Tahmini)

## File Information
- **File Name:** `unknown_TM_protein_pTM0.85_model.cif`  
- **Format:** CIF (Crystallographic Information File)  
- **Model Type:** In silico prediction (AlphaFold2 monomer)

## Overview
Bu dosya, AlphaFold2 algoritması kullanılarak tahmin edilmiş, bilinmeyen işleve sahip bir transmembran proteinin üç boyutlu yapısını içerir.

AlphaFold, amino asit dizisinden yüksek doğrulukla protein yapısı tahmin eden bir yapay zeka sistemidir. Bu model **tamamen hesaplamalıdır** ve deneysel olarak (X‑ışını kristalografisi, NMR, Cryo‑EM) belirlenmemiştir.

## Biyolojik Bağlam (Çıkarımsal)
Dizinin hidropatisite profili, 12–14 aday transmembran heliksine işaret etmektedir. Bu durum, aşağıdaki ailelerle benzerlik gösterebilir:

- **Major Facilitator Superfamily (MFS)** taşıyıcıları  
- Mitokondriyal taşıyıcı ailesi  
- Antiporter veya simporter özelliğinde bir membran proteini  

Bununla birlikte, **deneysel doğrulama yapılmadığından kesin işlev bilinmemektedir**.

## Yapısal Özellikler (Tahmini)
- Birden fazla transmembran heliks (≥12)  
- Kısa sitoplazmik ve ekstraselüler halkalar (loop)  
- C‑terminalinde 10×His saflaştırma etiketi (esnek bölge)  
- Potansiyel substrat bağlanma bölgeleri (örneğin `GSLQFGY` motifi)

## Model Güvenilirliği (AlphaFold pLDDT)
AlphaFold, her kalıntı için **pLDDT** skoru verir:

| pLDDT aralığı | Güven | Yorum |
|----------------|-------|-------|
| >90 | Çok yüksek | Helikslerin çoğu bu aralıkta (yeşil/mavi) |
| 70–90 | Yüksek | Transmembran bölgeleri güvenilir |
| 50–70 | Düşük | Bazı halkalar veya etiket kısmı |
| <50 | Çok düşük | Muhtemelen düzensiz/esnek |

- **pTM (küresel TM‑skor) = 0.85**  
  (>0.5 → katlanma doğru; 0.85 çok yüksek güven)

> Not: His etiketi ve terminal bölgeler düşük pLDDT içerebilir; bu beklenen bir durumdur.

## Kullanım Alanları (Kısıtlı Deneysel Doğrulama ile)
- 3B yapının görselleştirilmesi (PyMOL, ChimeraX, VMD)  
- Ön ligand yerleştirme çalışmaları (sanal tarama için **dikkatli** olunmalı)  
- Evrimsel korunum analizi yardımıyla fonksiyon tahmini  
- Mutasyon etkilerinin modellenmesi  

## Sınırlamalar
- Deneysel olarak doğrulanmamıştır  
- Protein dinamiklerini yakalamaz (tek bir konformer)  
- Ligand bağlı konformasyon garanti edilmez  
- Esnek döngüler (özellikle büyük loop’lar) yanlış olabilir  
- His etiketi nedeniyle C‑terminal katlanması doğal durumu yansıtmayabilir  

## Referans
- Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596, 583–589.  
- pTM skoru yorumu: AlphaFold EBI FAQ (2026)
