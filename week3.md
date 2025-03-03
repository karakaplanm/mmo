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
  - Transport and storage (e.g., hemoglobin, ferritin)<br>
    Transport Proteins and Their PDB Codes<br>
    1. **Hemoglobin (Oxygen Transporter in Blood)** – [1GZX](https://www.rcsb.org/structure/1GZX) (Human Deoxyhemoglobin)
    2. **Myoglobin (Oxygen Storage in Muscles)** – [1MBN](https://www.rcsb.org/structure/1MBN) (Sperm Whale Myoglobin)
    3. **Serum Albumin (Carrier of Fatty Acids, Hormones, and Drugs)** – [1E7I](https://www.rcsb.org/structure/1E7I) (Human Serum Albumin)
    4. **Transferrin (Iron Transport in Blood)** – [1H76](https://www.rcsb.org/structure/1H76) (Human Serum Transferrin)
    5. **Ferritin (Iron Storage in Cells)** – [1FHA](https://www.rcsb.org/structure/1FHA) (Human Ferritin)
    6. **Aquaporin (Water Transport Across Membranes)** – [1J4N](https://www.rcsb.org/structure/1J4N) (Aquaporin-1)
    7. **Glucose Transporter (GLUT1, Facilitates Glucose Uptake)** – [4PYP](https://www.rcsb.org/structure/4PYP) (Human GLUT1)
    8. **ABC Transporter (ATP-Binding Cassette Transporter, Moves Molecules Across Membranes)** – [3G5U](https://www.rcsb.org/structure/3G5U) (Bacterial ABC Transporter)
    9. **Lactoferrin (Iron-Binding and Immune Response Protein)** – [1LFG](https://www.rcsb.org/structure/1LFG) (Human Lactoferrin)
    10. **Vitamin B12 Transporter (Binds and Delivers Vitamin B12 in the Body)** – [2BBC](https://www.rcsb.org/structure/2BBC) (Vitamin B12-Binding Protein)
  
  Storage Proteins<br>
  1. **Ferritin (Stores and Releases Iron in Cells)** – [1FHA](https://www.rcsb.org/structure/1FHA) (Human Ferritin)
  2. **Casein (Milk Protein, Stores Nutrients for Growth)** – [1GXY](https://www.rcsb.org/structure/1GXY) (Casein Micelle Model)
  3. **Ovalbumin (Egg White Protein, Nutrient Storage)** – [1OVA](https://www.rcsb.org/structure/1OVA) (Ovalbumin)
  4. **Leghemoglobin (Stores Oxygen in Leguminous Plants' Root Nodules)** – [1BIN](https://www.rcsb.org/structure/1BIN) (Soybean Leghemoglobin)
  5. **Zeins (Corn Seed Storage Protein)** – [2GBY](https://www.rcsb.org/structure/2GBY) (Gamma-Zein)
  6. **Glycogenin (Protein Core of Glycogen, a Glucose Storage Molecule)** – [1LL0](https://www.rcsb.org/structure/1LL0) (Glycogenin-1)
  7. **Phytochelatin Synthase (Detoxification and Storage of Heavy Metals in Plants)** – [5TP9](https://www.rcsb.org/structure/5TP9) (Phytochelatin Synthase)
     
  - Signaling (e.g., hormones like insulin).<br>
    Signaling Proteins and Their PDB Codes
    - G Protein-Coupled Receptors (GPCRs)
      1. **Beta-2 Adrenergic Receptor (GPCR for Epinephrine Signaling)** – [2RH1](https://www.rcsb.org/structure/2RH1)
      2. **Rhodopsin (Light-Sensitive GPCR in Vision)** – [1F88](https://www.rcsb.org/structure/1F88)
      3. **Muscarinic Acetylcholine Receptor (GPCR for Neurotransmission)** – [4MQS](https://www.rcsb.org/structure/4MQS)
      4. **Adenosine A2A Receptor (GPCR for Cardiovascular Function)** – [2YDV](https://www.rcsb.org/structure/2YDV)
    - Kinases (Protein Phosphorylation Signaling)
      - **Protein Kinase A (PKA, cAMP-Dependent Signaling)** – [1ATP](https://www.rcsb.org/structure/1ATP)
      - **Tyrosine Kinase (Src Kinase, Signal Transduction in Growth Control)** – [2SRC](https://www.rcsb.org/structure/2SRC)
      - **MAP Kinase (Mitogen-Activated Protein Kinase, Cell Growth and Survival)** – [2OJG](https://www.rcsb.org/structure/2OJG)
      - **EGF Receptor Tyrosine Kinase (Cell Growth and Proliferation)** – [1M17](https://www.rcsb.org/structure/1M17)
    - Small GTPases (Molecular Switches in Signaling)
      - **Ras Protein (Cell Growth and Proliferation Signaling)** – [5P21](https://www.rcsb.org/structure/5P21)
      - **RhoA GTPase (Cytoskeletal Regulation)** – [1A2B](https://www.rcsb.org/structure/1A2B)
      - **Rab5 (Endocytosis and Vesicle Transport)** – [1R2Q](https://www.rcsb.org/structure/1R2Q)
    - Signal Transducers and Adapters
      - **STAT3 (Signal Transducer and Activator of Transcription, JAK-STAT Pathway)** – [1BG1](https://www.rcsb.org/structure/1BG1)
      - **14-3-3 Protein (Regulatory Adapter Protein in Cell Signaling)** – [1QJA](https://www.rcsb.org/structure/1QJA)
    - Calcium and Second Messenger Signaling
      - **Calmodulin (Calcium Sensor in Cell Signaling)** – [1CLL](https://www.rcsb.org/structure/1CLL)
      - **Protein Kinase C (PKC, Calcium and DAG-Dependent Kinase)** – [2A0M](https://www.rcsb.org/structure/2A0M)
      - **Adenylyl Cyclase (cAMP Producer in GPCR Signaling)** – [1CJU](https://www.rcsb.org/structure/1CJU)
    - Nuclear Receptors (Hormone Signaling)
      - **Glucocorticoid Receptor (Steroid Hormone Receptor)** – [1M2Z](https://www.rcsb.org/structure/1M2Z)
      - **Estrogen Receptor (Nuclear Receptor for Estrogen Signaling)** – [1G50](https://www.rcsb.org/structure/1G50)
      - **Retinoic Acid Receptor (Vitamin A Derivative Signaling)** – [1XAP](https://www.rcsb.org/structure/1XAP)


  - Immune response (e.g., antibodies).<br>
    Immune Response Proteins and Their PDB Codes
    - Antibodies and Immunoglobulins
      - **Immunoglobulin G (IgG) Antibody** – [1IGT](https://www.rcsb.org/structure/1IGT) (Human IgG1 Fab Fragment)
      - **Immunoglobulin A (IgA) Antibody** – [4N02](https://www.rcsb.org/structure/4N02) (Human Secretory IgA)
      - **Immunoglobulin E (IgE) Antibody (Allergy-Related)** – [4J4P](https://www.rcsb.org/structure/4J4P) (Human IgE-Fc)
      - **T-Cell Receptor (TCR) Complex** – [1TCR](https://www.rcsb.org/structure/1TCR) (TCR Bound to Peptide-MHC)
      - **Major Histocompatibility Complex I (MHC I, Antigen Presentation)** – [1HLA](https://www.rcsb.org/structure/1HLA) (HLA-A2 MHC Class I)
      - **Major Histocompatibility Complex II (MHC II, Antigen Presentation)** – [1DLH](https://www.rcsb.org/structure/1DLH) (HLA-DR1 MHC Class II)
    - Complement System Proteins
      - **Complement C3 (Central Component of Complement System)** – [2A73](https://www.rcsb.org/structure/2A73) (Human Complement C3)
      - **Complement C5 (Involved in MAC Complex Formation)** – [3CU7](https://www.rcsb.org/structure/3CU7) (Human Complement C5)
      - **Membrane Attack Complex (MAC) Component C9** – [6DLR](https://www.rcsb.org/structure/6DLR) (C9 MAC Pore)
    - Cytokines and Interleukins
      - **Interleukin-2 (IL-2, T-Cell Growth Factor)** – [1M47](https://www.rcsb.org/structure/1M47) (Human IL-2)
      - **Interleukin-6 (IL-6, Inflammatory Response Regulator)** – [1ALU](https://www.rcsb.org/structure/1ALU) (Human IL-6)
      - **Tumor Necrosis Factor Alpha (TNF-α, Pro-Inflammatory Cytokine)** – [2AZ5](https://www.rcsb.org/structure/2AZ5) (Human TNF-α)
      - **Interferon-Gamma (IFN-γ, Antiviral and Immune Response Mediator)** – [1HIG](https://www.rcsb.org/structure/1HIG) (Human IFN-γ)
    - Immune Cell Surface Receptors
      - **CD4 (T-Helper Cell Receptor, HIV Binding Target)** – [1CDH](https://www.rcsb.org/structure/1CDH) (Human CD4)
      - **CD8 (Cytotoxic T-Cell Co-Receptor)** – [1BQH](https://www.rcsb.org/structure/1BQH) (Human CD8α)
      - **PD-1 (Programmed Cell Death Protein 1, Immune Checkpoint)** – [3BIK](https://www.rcsb.org/structure/3BIK) (Human PD-1)
      - **CTLA-4 (Cytotoxic T-Lymphocyte-Associated Protein 4, Immune Regulation)** – [3OSK](https://www.rcsb.org/structure/3OSK) (Human CTLA-4)
    - Innate Immunity and Pattern Recognition Receptors
      - **Toll-Like Receptor 4 (TLR4, Pathogen Recognition Receptor)** – [3FXI](https://www.rcsb.org/structure/3FXI) (Human TLR4-MD-2-LPS Complex)
      - **NOD2 (Intracellular Pattern Recognition Receptor)** – [5IRN](https://www.rcsb.org/structure/5IRN) (Human NOD2)
      - **C-Reactive Protein (CRP, Acute Phase Inflammation Marker)** – [1B09](https://www.rcsb.org/structure/1B09) (Human CRP)


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
