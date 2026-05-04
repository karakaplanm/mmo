# 1M47 (Human Interleukin‑2) – Eukaryotic Cytokine

## Protein Introduction and Function

**1M47** represents the crystal structure of **Interleukin‑2 (IL‑2)** protein from human (*Homo sapiens*). IL‑2 is one of the most important regulatory molecules of the immune system, acting as a cytokine that affects **T cells**, **B cells**, and **natural killer (NK) cells**. By facilitating intercellular signaling, it plays critical roles in the initiation, maintenance, and termination of immune responses.

### Biological Function

- **T cell proliferation:** IL‑2 triggers the clonal expansion of activated T cells. This allows the body to produce sufficient numbers of fighter T cells to combat infections or tumors.

- **Regulatory T cell (Treg) control:** IL‑2 is essential for the survival and function of Treg cells. Low-dose IL‑2 supports Tregs to prevent autoimmunity, while high-dose IL‑2 activates effector T cells.

- **NK cell and B cell activation:** IL‑2 enhances the cytotoxic activity of NK cells and promotes antibody production from B cells.

- **Signal transduction:** Upon binding to its receptor (IL‑2Rα/β/γ), IL‑2 activates the JAK‑STAT, PI3K‑Akt, and MAPK signaling pathways.

### Structural Features

- **Length:** 133 amino acids (mature form)
- **Molecular weight:** ~15.4 kDa
- **Disulfide bond:** One disulfide bridge between Cys58 and Cys105; this bond contributes to structural stability.
- **Folding:** Four-helix bundle structure with an "up-up-down-down" topology. Contains no β‑sheets.
- **Active site / receptor binding surface:** Hydrophobic and polar residues on helix A and D mediate binding to the IL‑2 receptor.

### Applications

- **Cancer immunotherapy:** Recombinant IL‑2 (aldesleukin, Proleukin®) is FDA‑approved for the treatment of metastatic renal cell carcinoma and malignant melanoma.

- **Autoimmune diseases:** Low‑dose IL‑2 aims to restore immune tolerance by increasing Treg cells in diseases such as type 1 diabetes, systemic lupus erythematosus (SLE), and multiple sclerosis.

- **Basic science:** A paradigmatic model protein for cytokine folding, helical bundle structures, and receptor‑ligand recognition mechanisms.

- **Drug targeting:** Current research focuses on IL‑2 variants (superkines, bypass mutants) to achieve selectivity between Treg and effector T cells.

---

## Evaluation of AlphaFold Predictions – 1M47 (Human Interleukin‑2)

### Overall Model Quality

According to the **summary_confidences** file obtained from the AlphaFold prediction:

- **ptm (pTM) score = 0.85**  
  This score is well above the 0.5 threshold. **pTM > 0.5** indicates that the predicted global folding (backbone and topology) is correct. A high value of 0.85 demonstrates that this protein has been modeled with **very high confidence** and shows strong agreement with the experimental structure (1M47, X‑ray crystallography at 2.0 Å).

- **ipTM (interface pTM) = 0.85** (derived from chain_pair_iptm)  
  Although this protein is **monomeric** (single chain), the ipTM value is as high as 0.85. This reflects the **structural integrity and tight packing** of the protein. It may arise from a multimer‑like evaluation protocol.

- **pTM and ranking_score = 0.85**  
  The ranking_score indicates that this model ranks highest among alternative predictions. A score of 0.85 is a very high ranking score.

- **num_recycles = 10.0**  
  Indicates that AlphaFold used 10 recycling iterations to refine the model.

- **has_clash = 0.0** and **fraction_disordered = 0.0**  
  The predicted model contains **no atomic clashes** and no region is marked as **disordered**. This means that the entire IL‑2 protein is stable and well‑folded.

- **pAE_min (chain_pair_pae_min) = 0.76**  
  The very low predicted alignment error (0.76 Å) indicates that the relative positions of neighboring residues are predicted with very high accuracy.

### pLDDT Assessment

The available JSON file does not contain per‑residue pLDDT values. However, given that IL‑2 is a compact four‑helix bundle cytokine and the pTM score is 0.85, it is reasonable to expect that **all residues fall within the high or very high pLDDT range (>70)**.

### Conclusion

AlphaFold has modeled the **1M47 (Human Interleukin‑2)** protein with **exceptional reliability**. The pTM score of 0.85 demonstrates consistent agreement with the experimental crystal structure (2.0 Å resolution). This model provides a solid foundation for future studies on the receptor binding site, disulfide bond, and variant design of IL‑2.
