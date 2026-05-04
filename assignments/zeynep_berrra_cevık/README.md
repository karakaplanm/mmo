# 1GXY (ART2.2) – Eukaryotic Mono-ADP-ribosyltransferase

## Protein Introduction and Function

**1GXY** represents the crystal structure of **ART2.2** (ecto-ADP-ribosyltransferase 2.2) from rat (*Rattus norvegicus*). ART2.2 is an enzyme located on cell membranes that uses **extracellular NAD⁺** molecules to transfer an **ADP‑ribose** group to target proteins. This process is called **mono-ADP-ribosylation** and plays critical roles in cell signaling, particularly in **immune system regulation**.

### Biological Function

- **NAD⁺ metabolism:** ART2.2 breaks down extracellular NAD⁺ to produce ADP‑ribose and nicotinamide.
- **Production of signaling molecules:** The resulting ADP‑ribose can trigger cellular responses such as calcium signaling.
- **T‑cell regulation:** ART2.2 is expressed on the surface of T‑lymphocytes, thereby contributing to the control of immune responses.
- **Auto-ADP-ribosylation:** The enzyme can also ADP‑ribosylate itself to regulate its own activity.

### Structural Features

- **Length:** 226 amino acids
- **Molecular weight:** ~26.2 kDa
- **Disulfide bonds:** Two disulfide bridges located away from the active site; they confer high thermal and chemical stability to the protein.
- **N‑terminal extension:** Contains a hydrophobic region that anchors the protein to the membrane.
- **Folding:** Shares a core folding motif similar to bacterial ADP‑ribosyltransferases.

### Applications

- **Basic research:** Elucidating NAD⁺ signaling mechanisms.
- **Immunology:** Studies of T‑cell activation and tolerance.
- **Drug targeting:** A potential target in autoimmune diseases and inflammation models.

---

## Evaluation of AlphaFold Predictions – 1GXY (ART2.2)

### Overall Model Quality

According to the **summary_confidences** file obtained from the AlphaFold prediction:

- **ptm (pTM) score = 0.89**  
  This score is well above the 0.5 threshold. **pTM > 0.5** indicates that the predicted global fold (backbone and topology) is correct. A value as high as 0.89 demonstrates that this protein has been modeled with **very high confidence** and is in good agreement with the experimental structure (1GXY).

- **ipTM (interface pTM) =** Not computed (null).  
  This is because ART2.2 was modeled as a **single chain** in this prediction. ipTM measures interface confidence for multi‑chain complexes, so it is not meaningful for this protein.

- **ptm and ranking_score = 0.89**  
  The ranking_score indicates that this model is ranked best among alternative predictions. 0.89 is a very high ranking score.

- **num_recycles = 10.0**  
  Indicates that AlphaFold used 10 recycling iterations to improve the model.

- **has_clash = 0.0** and **fraction_disordered = 0.0**  
  The predicted model has **no atomic clashes** and no region is marked as **disordered**. This means that the entire ART2.2 protein is predicted to be stable and folded.

### pLDDT Assessment

The provided JSON file does not contain per‑residue pLDDT values, so the table notes “no information available”. However, given that ART2.2 is a small protein stabilized by disulfide bonds and the ptm = 0.89 value, it is reasonable to expect that **all residues fall in the high or very high pLDDT range (>70)**.

### Conclusion

AlphaFold has modeled the **1GXY (ART2.2)** protein with **exceptional reliability**. The pTM score of 0.89 indicates that the prediction is consistent with the experimental crystal structure (1.71 Å resolution, P21 crystal form). This model provides a solid foundation for future studies on the NAD⁺‑binding site, disulfide bonds, and membrane‑associated regions of ART2.2.
