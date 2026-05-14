# Ab Initio Structure Prediction of Klotho Beta (KLB)

**Author:** Kübra NAZLI  
**Institution:** Inonu University – Department of Molecular Biology and Genetics  
**Course:** Molecular Modeling  

---

## Project Overview

This project focuses on the *ab initio* three-dimensional (3D) structure prediction of the Klotho Beta (KLB) protein using AlphaFold 3. KLB is a transmembrane co-receptor protein that plays critical roles in metabolic regulation, bile acid homeostasis, and fibroblast growth factor (FGF) signaling pathways.

Because experimentally resolved structural information for KLB is limited, computational modeling provides an important approach for understanding its structural organization and potential functional regions.

---

## Methodology

The structural prediction was performed using **Google DeepMind AlphaFold 3**.

- Target Protein: **Klotho Beta (KLB)**
- UniProt ID: **Q86Z14**
- Prediction Type: **Monomeric structure**
- Total Generated Models: **5**

Among the generated conformations, **model_0** was selected as the reference structure due to its highest confidence ranking and overall structural consistency.

Structural confidence analysis was performed using:
- pTM score
- pLDDT confidence distribution
- Predicted disorder region analysis
- PAE (Predicted Aligned Error) matrices

---

## Structural Analysis and Visualization

The predicted `.cif` structure was imported into **PyMOL** for visualization and structural inspection.

The model was colored according to per-residue pLDDT confidence scores:

| Confidence Level | Color | Interpretation |
|---|---|---|
| pLDDT > 90 | Blue | Highly reliable structured domains |
| 90 > pLDDT > 50 | Cyan / Yellow | Flexible regions with moderate confidence |
| pLDDT < 50 | Orange | Low-confidence or intrinsically disordered regions |

The extracellular domains displayed high structural confidence, while terminal regions showed increased flexibility and disorder propensity.

---

## Highlighted Functional Region

In the PyMOL render (`klb_pymol_render.png`), residues predicted to participate in the **FGF21 interaction interface** were highlighted as **red spheres** to emphasize the probable ligand-binding region.

This region may contribute to receptor specificity and downstream signaling interactions.

---

## Repository Contents

| File Name | Description |
|---|---|
| `fold_klb_model_model_0.cif` | Reference 3D structural model predicted by AlphaFold |
| `fold_klb_model_summary_confidences_0.json` | Confidence metrics including pLDDT values and PAE matrices |
| `klb_pymol_render.png` | PyMOL visualization with confidence coloring and highlighted interaction regions |

---

## Conclusion

The AlphaFold-based *ab initio* prediction successfully generated a plausible structural model of Klotho Beta (KLB).

The predicted structure revealed:
- Stable extracellular core domains
- Flexible terminal regions
- Potential ligand interaction interfaces

These findings are consistent with the dynamic signaling role of KLB in FGF-mediated pathways.

This project demonstrates how modern computational structural biology tools can provide valuable structural insight for proteins lacking experimentally resolved structures.
