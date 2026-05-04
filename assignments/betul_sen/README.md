# 🧬 Protein Structure Prediction and Visualization (1CDH)

## 📌 Project Overview

This project aims to predict and visualize the three-dimensional (3D) structure of the protein with PDB ID **1CDH** using computational biology tools. The workflow includes sequence retrieval, structure prediction, and visualization.

---

## 🔬 Methodology

### 1. Sequence Retrieval

* The protein ID **1CDH** was searched in the RCSB Protein Data Bank
* The FASTA sequence of the protein was obtained

### 2. Structure Prediction

* The retrieved FASTA sequence was submitted to the AlphaFold Server
* The predicted structure files were downloaded as a ZIP archive

### 3. File Processing

* The ZIP archive was extracted
* The `.cif` structure files were obtained from the output

### 4. Visualization

* The predicted structure was visualized using PyMOL
* A high-quality PNG image of the protein structure was generated

---

## 📁 Repository Contents

* `structure.png` → 3D visualization of the protein
* `fold_2026_05_04_10_41_model_0.cif` → Selected predicted structure model
* `fold_2026_05_04_10_41_full_data_0.json` → Detailed prediction data
* `fold_2026_05_04_10_41_summary_confidence_0.json` → Confidence scores of the model

---

## 📊 Confidence Evaluation

AlphaFold provides confidence metrics for predicted structures:

* **pLDDT (Predicted Local Distance Difference Test):**

  * 90–100 → Very high confidence
  * 70–90 → Reliable
  * 50–70 → Low confidence
  * <50 → Very low confidence

In this project, **Model 0** was selected for analysis as it represents one of the top-ranked predicted structures.

---

## 📷 Protein Structure

![Protein Structure](structure.png)

---

## 🧪 Tools and Resources

* RCSB Protein Data Bank
* AlphaFold Protein Structure Prediction Server
* PyMOL Molecular Visualization Tool

---

## 🧠 Conclusion

The AlphaFold-based prediction successfully generated a reliable 3D model of the selected protein. Visualization in PyMOL allowed structural interpretation, while confidence scores provided insight into prediction accuracy.

---

## ✍️ Author

**Betül ŞEN**
Molecular Biology and Genetics Student
İnönü University
