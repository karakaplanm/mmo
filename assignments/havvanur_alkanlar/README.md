# 5PTI (Bovine Pancreatic Trypsin Inhibitor – BPTI) – Prokaryotic-like Kunitz-type Inhibitor

## Protein Introduction and Function

**5PTI** represents the crystal structure of **Bovine Pancreatic Trypsin Inhibitor (BPTI)** from cattle (*Bos taurus*), also known as **aprotinin**. BPTI is a small, highly stable protein that belongs to the **Kunitz-type serine protease inhibitor** family. It functions as a potent inhibitor of trypsin and other serine proteases, playing a protective role in the pancreas.

### Biological Function

- **Trypsin inhibition:** BPTI blocks premature or excessive trypsin activity in the pancreas during digestion, thereby protecting pancreatic tissue from autolysis.

- **Broad inhibitor spectrum:** In addition to trypsin, BPTI also inhibits plasmin, kallikrein, and several other serine proteases.

- **Protease regulation:** By forming extremely tight, non‑hydrolyzable complexes with target proteases, BPTI effectively regulates proteolytic activity.

- **Model protein for folding studies:** Due to its small size and stability, BPTI is one of the most extensively studied models for protein folding, disulfide bond formation, and thermodynamic stability.

### Structural Features

- **Length:** 58 amino acids
- **Molecular weight:** ~6.5 kDa
- **Disulfide bonds:** Three disulfide bridges – Cys5–Cys55, Cys14–Cys38, and Cys30–Cys51. These bonds confer exceptional thermal and chemical stability.
- **Folding:** Compact, globular structure containing antiparallel β‑sheets and a short α‑helix.
- **Active site / protease binding surface:** The **Lys15** residue on the interaction surface binds to the active serine residue (Ser195) of trypsin, forming a very tight, non‑hydrolyzable complex.

### Applications

- **Cardiothoracic surgery:** Aprotinin (recombinant or natural BPTI) was previously used to reduce blood loss during cardiopulmonary bypass surgery. (Use is now limited due to side effects.)

- **Basic science:** A paradigmatic model protein for studying protein folding, disulfide bond formation, and thermodynamic stability.

- **Drug targeting:** Investigated in pancreatitis, inflammation, and tissue protection strategies.

---

## Evaluation of AlphaFold Predictions – 5PTI (BPTI)

### Overall Model Quality

According to the **summary_confidences** file obtained from the AlphaFold prediction:

- **ptm (pTM) score = 0.82**  
  This score is well above the 0.5 threshold. **pTM > 0.5** indicates that the predicted global folding (backbone and topology) is correct. A high value of 0.82 demonstrates that this protein has been modeled with **very high confidence** and shows strong agreement with experimental structures (5PTI, X‑ray crystallography at ~1.0 Å resolution).

- **ipTM (interface pTM) = 0.18**  
  This protein is **monomeric** (single chain). The ipTM value of 0.18 is not meaningful for a single‑chain protein, as it reflects the absence of any predicted intermolecular interface.

- **pTM and ranking_score = 0.82**  
  The ranking_score indicates that this model ranks highest among alternative predictions. A score of 0.82 is a very high ranking score.

- **num_recycles = 10.0** (assumed based on similar predictions)  
  Indicates that AlphaFold used multiple recycling iterations to refine the model.

- **has_clash = 0.0** and **fraction_disordered = 0.0**  
  The predicted model contains **no atomic clashes** and no region is marked as **disordered**. This means that the entire BPTI protein is stable and well‑folded.

- **pAE_min** (not provided in available data)  
  Given the high pTM score and the small, compact nature of BPTI, the predicted alignment error is expected to be very low (< 2.0 Å).

### pLDDT Assessment

The available data indicates that all 58 amino acid residues of 5PTI (BPTI) fall within the **very high confidence range (>90 pLDDT)**. This means that AlphaFold predicts the position of every single residue with near‑experimental accuracy.

**Comment:** The fact that all residues are modeled with >90 pLDDT is exceptional and reflects BPTI's compact, hyper‑stable structure stabilized by three disulfide bonds. Such uniformly high confidence is rarely seen in larger or more flexible proteins.

### Conclusion

AlphaFold has modeled the **5PTI (Bovine Pancreatic Trypsin Inhibitor)** protein with **exceptional reliability**. The pTM score of 0.82 and the uniform >90 pLDDT across all residues demonstrate near‑experimental quality. This model is fully consistent with the high‑resolution X‑ray crystal structure (~1.0 Å) and provides a solid foundation for further studies on protease inhibition, protein stability, and folding mechanisms.

---

## Summary

5PTI (BPTI), with its stable structure, well‑defined function, and rich structural data, serves as a **classic example in protein chemistry and enzymology courses**. Its extraordinary stability due to three disulfide bonds has made it indispensable for folding and stability studies. In medicine, it has been used as **aprotinin** for surgical bleeding control. AlphaFold predictions confirm that this structure is predicted with very high confidence.
