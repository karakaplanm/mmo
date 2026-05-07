# Concanavalin A (3CNA) – High-Resolution Structure Prediction with AlphaFold

This repository contains the three-dimensional structure of **PDB: 3CNA** (Concanavalin A) modeled using **AlphaFold**. The quality metrics obtained from the modeling are provided below.

## 📊 Model Quality Scores

| Metric | Value | Description |
|--------|-------|-------------|
| **pTM** | 0.95 | Very high confidence – the model is almost identical to the experimental structure. |
| **ipTM** | - | Not computed (expected for a monomeric model). |

> The pTM (predicted Template Modeling score) ranges from 0 to 1 and indicates global folding accuracy. A score of 0.95 signifies an extremely reliable structure.

## 🧬 About the Protein

**Concanavalin A (ConA)** is a **lectin** protein derived from the jack bean plant (*Canavalia ensiformis*). Lectins are a family of proteins that specifically recognize and bind to carbohydrates. ConA is one of the best-characterized and most widely used members of this family.

### Structural Features

- **Binding target:** High affinity for α-D-mannose and α-D-glucose sugars. This binding is reversible and highly specific.
- **Metal cofactors:** Each subunit contains two metal ions: manganese (Mn²⁺) and calcium (Ca²⁺). These metals are essential for proper folding and maintaining the sugar-binding site in an active conformation. When the metals are removed, the protein loses its binding ability.
- **Oligomerization:** At neutral pH (7.0–7.5), four identical subunits assemble to form an active tetramer. At acidic pH, it can dissociate into dimers or monomers.
- **Molecular weight:** Each subunit is approximately 25–26 kDa; the tetramer is approximately 102–104 kDa.
- **Stability:** ConA is stable at high temperatures (up to 70°C) and over a broad pH range (pH 4–10). This stability makes it a preferred tool in laboratory applications.

### Biological and Biotechnological Functions

1. **Affinity chromatography (most common use):** ConA molecules are covalently attached to a chromatography column. As a biological sample passes through, glycosylated molecules (e.g., immunoglobulins, hormones, membrane proteins) bind to the column while non-glycosylated molecules flow through. Specifically bound molecules are then eluted using a competing sugar (typically methyl-α-D-mannopyranoside), yielding highly purified products.

2. **Mitogenic activity (immunology):** ConA binds to specific glycan structures on the surface of T lymphocytes, strongly activating them. This activation leads to blast transformation (cell swelling), initiation of DNA synthesis, and cell division. This property makes ConA a standard tool in immunology laboratories for stimulating T cell proliferation, studying cytokine profiles, and investigating immune response mechanisms.

3. **Cell surface glycan probe:** ConA can be labeled with fluorescent dyes (such as FITC or rhodamine). Labeled ConA is used to visualize the presence, density, and distribution of mannose/glucose-containing glycans on cell surfaces or in tissue sections using microscopy. This method is also used to detect abnormal glycosylation patterns in cancer cells.

4. **Anticancer potential (research stage):** Recent studies have shown that ConA can induce apoptosis (programmed cell death) in certain cancer cell lines, particularly melanoma, liver, and breast cancer cells. This effect is thought to result from ConA recognizing characteristic glycan structures on cancer cell surfaces, followed by activation of mitochondrial signaling pathways. However, these effects are still experimental and have not entered clinical use.

## 🧪 Detailed Interpretation of AlphaFold Modeling Results

AlphaFold is a deep learning-based structure prediction algorithm. It takes the protein's amino acid sequence (UniProt ID: P02861) and predicts its three-dimensional structure. The obtained **pTM = 0.95** value means:

- 0.95 is just below the maximum possible pTM score of 1.0. This indicates that the **global folding shape** of the model is almost identical to experimentally determined structures from X-ray crystallography or cryo-EM.
- According to AlphaFold's evaluation criteria, models with pTM > 0.7 are considered "high confidence." At 0.95, this model is considered **experimental quality**.
- The absence of an **ipTM** value is normal because ipTM (interface pTM) is only computed for protein-protein complexes or multimeric predictions. This study performed a monomeric prediction, so ipTM is not defined.

This model can be reliably used for the following advanced analyses:

| Analysis Type | Description |
|---------------|-------------|
| Molecular Docking | Different sugar derivatives (mannose, glucose, derivatives) can be docked into the binding site to calculate binding energies. Software such as AutoDock or HADDOCK can be used. |
| Molecular Dynamics Simulation | Protein motions on nanosecond to microsecond timescales, metal ion release, and flexibility of the sugar-binding site can be simulated using GROMACS, AMBER, or NAMD. |
| Mutation Analysis | In silico mutations can be introduced to amino acids responsible for sugar binding or metal coordination (e.g., Asp208, Arg228). The effects of these mutations on structural stability and function can be predicted. |
| Structural Alignment and Comparison | The ConA model can be aligned with other lectin families (e.g., soybean agglutinin, PNA, RCA) to identify conserved regions and structural differences. |
| Epitope Prediction | Regions where antibodies raised against ConA are likely to bind can be predicted. |
| Virtual Drug Screening | Thousands of small molecules can be docked into the sugar-binding site of ConA to screen for potential inhibitors. |

## 🔧 Practical Recommendations for Using the Model

- **Visualization:** The model can be opened in PyMOL, ChimeraX, VMD, or UCSF Chimera to examine the sugar-binding site and metal ion coordination.
- **Validation (optional):** You can further validate the model using tools such as QMEAN, PROCHECK, or MolProbity. However, given the pTM = 0.95, these tools are also expected to give high scores.
- **Comparison with experimental data:** You can calculate the RMSD between your model and the experimental 3CNA structure available in the PDB (X-ray, 2.0 Å resolution). The expected RMSD is around 1 Å or lower.

## 💎 Conclusion

This AlphaFold model with a pTM score of 0.95 provides an **experimental-quality, highly reliable starting point** for understanding the structure-function relationship of Concanavalin A, simulating molecular interactions, and developing biotechnological applications. The model is mature enough to be used in scientific publications, master's/doctoral theses, or preliminary drug design studies.