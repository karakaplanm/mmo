# 🧬 Casein Protein Structure Prediction with AlphaFold

## 📌 Project Overview

This project focuses on predicting the three-dimensional (3D) structure of a casein protein using state-of-the-art deep learning approaches. The protein sequence was obtained in FASTA format and analyzed using AlphaFold-based tools. The resulting structural model was exported in `.cif` format and shared in this repository for reproducibility and further analysis.

---

## 🧪 Protein Information

* **Protein Family:** Casein
* **Selected Type:** β-casein (example; update if different)
* **Organism:** *Bos taurus* (cow)
* **Database Source:** UniProt
* **UniProt ID:** P02666 *(update if needed)*

Caseins are a group of milk proteins known for their flexible and partially disordered structures. This characteristic plays a critical role in their biological function and affects computational structure prediction accuracy.

---

## 🔬 Methodology

### 1. Sequence Retrieval

* The amino acid sequence was retrieved from the UniProt database in FASTA format.
* The sequence header and residues were preserved as required for structure prediction.

### 2. Structure Prediction

* The FASTA sequence was submitted to an AlphaFold-based prediction pipeline (ColabFold recommended).
* Default parameters were used unless otherwise specified.
* Multiple sequence alignment (MSA) was generated automatically.

### 3. Output Files

The prediction pipeline generated several output files:

* `.pdb` → Standard protein structure format
* `.json` → Confidence metrics (pLDDT, PAE)
* `.cif` → **Final structure file (used in this repository)**

---

## 📂 Repository Contents

* `casein_structure.cif` → Predicted 3D structure
* `confidence_metrics.json` → Model confidence scores
* `sequence.fasta` → Input amino acid sequence
* `README.md` → Project documentation

---

## 📊 Model Confidence & Interpretation

AlphaFold provides confidence metrics:

* **pLDDT Score (per-residue confidence):**

  * > 90 → High confidence
  * 70–90 → Medium confidence
  * <70 → Low confidence

⚠️ **Important Note:**
Casein proteins are classified as **intrinsically disordered proteins (IDPs)**.
This means:

* They may not adopt a single stable 3D structure
* Lower confidence scores are expected
* Structural flexibility is biologically relevant

---

## 🧠 Scientific Considerations

* The predicted structure should be interpreted cautiously due to the disordered nature of caseins.
* Regions with low pLDDT may correspond to functionally important flexible domains.
* Experimental validation (e.g., NMR, SAXS) is recommended for high-accuracy studies.

---

## 🛠 Tools & Resources Used

* AlphaFold / ColabFold (structure prediction)
* UniProt (sequence database)
* PyMOL or ChimeraX (structure visualization)
* GitHub (data sharing)

---

## 🚀 How to Reproduce

1. Retrieve the FASTA sequence from UniProt
2. Open ColabFold notebook
3. Paste the sequence
4. Run all cells
5. Download output files
6. Upload `.cif` to this repository

---

## 📈 Future Work

* Compare different casein isoforms (αS1, αS2, β, κ)
* Perform disorder prediction analysis
* Conduct molecular dynamics simulations
* Analyze protein-protein interactions

---

## 📚 References

* Jumper et al., 2021 — Highly accurate protein structure prediction with AlphaFold
* UniProt Consortium — Protein sequence database

---

## 👤 Author

Prepared by a Molecular Biology & Genetics student interested in protein structure, evolution, and computational biology.

---
