# Week3: Structure of Proteins and Visualization of a Protein. Classification and Basic Calculations

---

## 1. Introduction to Proteins
- **Definition**: Proteins are large, complex biomolecules composed of amino acids, essential for the structure, function, and regulation of the body's tissues and organs.
- **Functions**:
  - Enzymes (catalyze biochemical reactions).<br>
    Well-Known Enzymes and Their RCSB PDB Codes
    1. **DNA Polymerase** – [1TAQ](https://www.rcsb.org/structure/1TAQ) (Taq DNA Polymerase)
    2. **RNA Polymerase** – [4JKR](https://www.rcsb.org/structure/4JKR) (Bacterial RNA Polymerase)
    3. **Amylase** – [1SMD](https://www.rcsb.org/structure/1SMD) (Human Salivary Amylase)
    4. **Lipase** – [1LPB](https://www.rcsb.org/structure/1LPB) (Pancreatic Lipase)
    5. **Protease (Trypsin)** – [2PTN](https://www.rcsb.org/structure/2PTN) (Bovine Trypsin)
    6. **Lactase** – [3OBA](https://www.rcsb.org/structure/3OBA) (Beta-Galactosidase, Lactase-related)
    7. **Catalase** – [1DGF](https://www.rcsb.org/structure/1DGF) (Bovine Catalase)
    8. **Hexokinase** – [1HKC](https://www.rcsb.org/structure/1HKC) (Rat Hexokinase)
    9. **Trypsin** – [3TGI](https://www.rcsb.org/structure/3TGI) (Human Trypsin)
    10. **Carbonic Anhydrase** – [1CA2](https://www.rcsb.org/structure/1CA2) (Human Carbonic Anhydrase II)
    11. **ATP Synthase** – [1QO1](https://www.rcsb.org/structure/1QO1) (F1 ATP Synthase)

  - Structural components (e.g., collagen, keratin).<br>
    Structural Components of Proteins and Their PDB Codes
    1. **Alpha-Helix (Myoglobin as an example)** – [1MBN](https://www.rcsb.org/structure/1MBN) (Sperm Whale Myoglobin)
    2. **Beta-Sheet (Concanavalin A as an example)** – [3CNA](https://www.rcsb.org/structure/3CNA) (Jack Bean Concanavalin A)
    3. **Beta-Turn (Bovine Pancreatic Trypsin Inhibitor as an example)** – [5PTI](https://www.rcsb.org/structure/5PTI) (BPTI)
    4. **Coiled-Coil (GCN4 Leucine Zipper as an example)** – [2ZTA](https://www.rcsb.org/structure/2ZTA) (GCN4 Leucine Zipper)
    5. **Collagen Triple Helix** – [1BKV](https://www.rcsb.org/structure/1BKV) (Collagen-like Peptide)
    6. **Zinc Finger (Zif268 as an example)** – [1ZAA](https://www.rcsb.org/structure/1ZAA) (Zif268 Zinc Finger Protein)
    7. **Immunoglobulin Fold (Antibody Fab Fragment)** – [1IGT](https://www.rcsb.org/structure/1IGT) (Human Immunoglobulin G1 Fab Fragment)
    8. **TIM Barrel (Triosephosphate Isomerase as an example)** – [1TIM](https://www.rcsb.org/structure/1TIM) (Triosephosphate Isomerase)
    9. **Rossmann Fold (Lactate Dehydrogenase as an example)** – [6LDH](https://www.rcsb.org/structure/6LDH) (Lactate Dehydrogenase)
    10. **Ferritin (Iron Storage Protein with Alpha-Helices)** – [1FHA](https://www.rcsb.org/structure/1FHA) (Human Ferritin)
  - Transport and storage (e.g., hemoglobin, ferritin).
  - Signaling (e.g., hormones like insulin).
  - Immune response (e.g., antibodies).

---

## 2. Structure of Proteins
Proteins have four levels of structural organization:

### A. Primary Structure
- Linear sequence of amino acids linked by peptide bonds.
- Determined by the genetic code.
- Example: Ala-Gly-Val-Leu-Ser.

### B. Secondary Structure
- Local folding of the polypeptide chain into regular structures.
- Common types:
  - **Alpha-helix**: Spiral structure stabilized by hydrogen bonds.
  - **Beta-sheet**: Pleated sheets stabilized by hydrogen bonds between adjacent strands.
  - **Turns and loops**: Irregular structures connecting helices and sheets.

### C. Tertiary Structure
- Three-dimensional folding of the entire polypeptide chain.
- Stabilized by:
  - Hydrophobic interactions.
  - Hydrogen bonds.
  - Ionic bonds.
  - Disulfide bridges (covalent bonds between cysteine residues).
- Example: Myoglobin.

### D. Quaternary Structure
- Arrangement of multiple polypeptide chains (subunits) into a functional protein complex.
- Example: Hemoglobin (4 subunits).

---

## 3. Visualization of Proteins
- **Importance**: Visualization helps understand protein structure, function, and interactions.
- **Tools and Techniques**:
  - **X-ray Crystallography**: High-resolution 3D structures.
  - **Nuclear Magnetic Resonance (NMR)**: Solution-state structures.
  - **Cryo-Electron Microscopy (Cryo-EM)**: Large protein complexes.
  - **Computational Tools**:
    - **PyMOL**: Molecular visualization software.
    - **Chimera**: Interactive visualization and analysis.
    - **RasMol**: Simple visualization tool.

---

## 4. Classification of Proteins
Proteins can be classified based on:

### A. Shape
- **Fibrous Proteins**: Long, elongated structures (e.g., collagen, keratin).
- **Globular Proteins**: Compact, spherical structures (e.g., enzymes, hemoglobin).

### B. Composition
- **Simple Proteins**: Composed only of amino acids (e.g., albumin).
- **Conjugated Proteins**: Contain non-protein groups (prosthetic groups):
  - Glycoproteins (carbohydrates).
  - Lipoproteins (lipids).
  - Metalloproteins (metal ions).
  - Hemoproteins (heme group).

### C. Function
- Enzymes, structural proteins, transport proteins, etc.

---

## 5. Basic Calculations in Protein Science

### A. Molecular Weight Calculation
- Formula: Sum of the molecular weights of all amino acids in the protein.
- Example: Ala (89) + Gly (75) + Val (117) = 281 g/mol.

### B. Amino Acid Composition
- Calculate the percentage of each amino acid in the protein.
- Formula: `(Number of a specific amino acid / Total number of amino acids) × 100`.

### C. Isoelectric Point (pI)
- The pH at which a protein has no net charge.
- Calculated using the pKa values of ionizable groups (e.g., amino and carboxyl groups).

### D. Extinction Coefficient
- Measures how much light a protein absorbs at a specific wavelength (usually 280 nm).
- Calculated based on the number of tryptophan, tyrosine, and cysteine residues.

---

## 6. Applications of Protein Structure and Visualization
- Drug design (e.g., targeting active sites of enzymes).
- Understanding disease mechanisms (e.g., misfolded proteins in Alzheimer’s).
- Protein engineering (e.g., designing enzymes for industrial use).

---

## 7. Summary
- Proteins are essential biomolecules with complex structures.
- Visualization tools help us understand their 3D organization.
- Classification and calculations provide insights into their properties and functions.

---

## 8. Questions for Discussion
1. How does the primary structure determine the tertiary structure of a protein?
2. What are the advantages and limitations of different protein visualization techniques?
3. How can protein classification aid in understanding their biological roles?

---

## 9. Further Reading
- "Introduction to Protein Structure" by Carl Branden and John Tooze.
- "Protein Structure and Function" by Gregory Petsko and Dagmar Ringe.
- Online resources: [Protein Data Bank (PDB)](https://www.rcsb.org/), [UniProt](https://www.uniprot.org/).

---
