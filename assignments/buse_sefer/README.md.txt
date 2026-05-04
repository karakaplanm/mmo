# 🧬 Protein Structure Prediction and Visualization (3-BIK – Human PD-1)

## 📌 Project Overview

This project focuses on predicting and visualizing the three-dimensional (3D) structure of the human protein **PD-1 (Programmed Cell Death Protein 1)** using computational tools. The workflow includes sequence retrieval, structure prediction, and structural visualization.

---

## 🔬 Methodology

### 1. Sequence Retrieval

* The protein ID **3-BIK** was searched in the RCSB Protein Data Bank
* The FASTA sequence of Human PD-1 was obtained

### 2. Structure Prediction

* The FASTA sequence was submitted to the AlphaFold Server
* Predicted structure files were downloaded as a ZIP archive

### 3. File Processing

* The ZIP archive was extracted
* The `.cif` structure files were obtained from the output

### 4. Visualization

* The predicted structure was visualized using PyMOL
* A PNG image of the protein structure was generated

---

## 📁 Repository Contents

* `structure.png` → 3D visualization of Human PD-1
* `fold_..._model_2.cif` → Selected predicted structure (Model 2)
* `fold_..._full_data_2.json` → Detailed prediction data
* `fold_..._summary_confidence_2.json` → Confidence scores of the model

---

## 📊 Confidence Evaluation

AlphaFold provides confidence metrics for predicted structures:

* **pLDDT (Predicted Local Distance Difference Test):**

  * 90–100 → Very high confidence
  * 70–90 → Reliable
  * 50–70 → Low confidence
  * <50 → Very low confidence

In this project, **Model 2** was selected for analysis based on the available output files.

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

The AlphaFold-based prediction successfully generated a 3D structural model of Human PD-1. Visualization in PyMOL enabled structural analysis, while confidence scores provided insight into the reliability of the predicted regions.

---

## ✍️ Author

**Buse SEFER**
Molecular Biology and Genetics Student
İnönü University
