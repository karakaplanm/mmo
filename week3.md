
# Week 3: Introduction to Nucleic Acids Structure, Visualization, and Basic Calculations

---

## **I. Introduction to Nucleic Acids**

### **1.1 What are Nucleic Acids?**
- **Definition**: Biopolymers essential for life; store, transmit, and express genetic information.
- **Types**:
  - **DNA (Deoxyribonucleic Acid)**: Double-stranded, stores genetic code.
  - **RNA (Ribonucleic Acid)**: Single-stranded, involved in protein synthesis (e.g., mRNA, tRNA, rRNA).

### **1.2 Structure of Nucleic Acids**
#### **Nucleotides: Building Blocks**
- **Components**:
  1. **Phosphate Group**: Negatively charged backbone.
  2. **Sugar**:
     - DNA: Deoxyribose (lacks -OH at 2' position).
     - RNA: Ribose (has -OH at 2' position).
  3. **Nitrogenous Base**:
     - **Purines**: Adenine (A), Guanine (G).
     - **Pyrimidines**: Cytosine (C), Thymine (T) [DNA], Uracil (U) [RNA].

#### **DNA Double Helix**
- **Key Features**:
  - Antiparallel strands (5' → 3' and 3' → 5').
  - Complementary base pairing: **A-T** (2 hydrogen bonds), **C-G** (3 hydrogen bonds).
  - Major and minor grooves: Sites for protein interactions.

#### **RNA Structure**
- Single-stranded but forms secondary structures (e.g., hairpins, stem-loops).
- **Key Differences from DNA**:
  - Uracil replaces Thymine.
  - 2' -OH group makes RNA more reactive.

---

## **II. Visualization of Nucleic Acids**

### **2.1 Experimental Techniques**
- **X-ray Crystallography**: Historic method (e.g., Franklin’s Photo 51 for DNA).
- **Cryo-Electron Microscopy (Cryo-EM)**: Visualizes large complexes (e.g., ribosomes).
- **NMR Spectroscopy**: For small, dynamic structures (e.g., RNA aptamers).

### **2.2 Computational Tools**
- **Software**:
  - **PyMOL**: Renders 3D structures; highlights base pairs, grooves.
  - **ChimeraX**: Analyzes hydrogen bonds, electrostatic surfaces.
  - **Jmol**: Web-based viewer for PDB files.
- **Databases**:
  - **Protein Data Bank (PDB)**: Repository for 3D structures (e.g., [PDB ID 1BNA](https://www.rcsb.org/)).

#### **Example: Visualizing DNA in PyMOL**
1. Load PDB file.
2. Color strands differently.
3. Highlight base pairs (`select bp`).
4. Measure distances (e.g., H-bond lengths).

---

## **III. Basic Calculations on Nucleic Acids**

### **3.1 Concentration and Purity**
- **UV Absorbance**:
  - **OD₆₂₀**: Measures nucleic acid concentration.
  - **Beer-Lambert Law**:  
    \[
    \text{Concentration (µg/µL)} = \text{OD}_{260} \times \text{Dilution Factor} \times 50 \, (\text{for DNA}) \, \text{or} \, 40 \, (\text{for RNA})
    \]
  - **Purity Check**: OD₂₆₀/OD₂₈₀ ratio (~1.8 for pure DNA; ~2.0 for RNA).

### **3.2 Melting Temperature (Tₘ)**
- **Definition**: Temperature at which 50% of DNA strands dissociate.
- **Formula for Short DNA**:
  \[
  T_m = 2^\circ C \times (A+T) + 4^\circ C \times (G+C)
  \]
- **GC Content**: Higher GC content → higher Tₘ (more H-bonds).

### **3.3 Molecular Weight (MW)**
- **Formula**:
  \[
  \text{MW} = (n_A \times 313.2) + (n_T \times 304.2) + (n_C \times 289.2) + (n_G \times 329.2) + 79.0 \, (\text{for terminal phosphates})
  \]

### **3.4 Molar Conversions**
- **Example**: Convert µg of DNA to pmol:
  \[
  \text{pmol} = \frac{\text{µg} \times 10^6}{\text{MW (g/mol)}}
  \]

### **3.5 Primer Design**
- **Key Parameters**:
  - Length: 18–30 bases.
  - GC Content: 40–60%.
  - Avoid self-complementarity (secondary structures).

---

## **IV. Summary**
- Nucleic acids are polymers of nucleotides with distinct DNA/RNA structures.
- Visualization tools (PyMOL, Cryo-EM) reveal 3D architecture.
- Basic calculations (concentration, Tₘ, MW) are critical for lab work.