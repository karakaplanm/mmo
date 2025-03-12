# Lecture 2: Protein Structure Prediction — Protein Folding and Homology Modeling

---

## **I. Introduction to Protein Structure Prediction**
- **Why predict structure?**  
  - Determines function, enables drug design, and elucidates disease mechanisms.
- **Key challenges**:  
  - Computational complexity, conformational flexibility, and limited experimental data.

---

## **II. Protein Folding**

### **2.1 Principles of Protein Folding**
- **Anfinsen's Dogma**:  
  - *"Sequence dictates structure."* Native 3D structure is thermodynamically favored.
- **Levinthal's Paradox**:  
  - Proteins fold in seconds, despite ~10³⁰⁰ possible conformations.  
  - Implies folding follows a directed pathway or "funnel-like" energy landscape.

#### **Stages of Folding**  
1. **Hydrophobic collapse**: Rapid burial of hydrophobic residues.  
2. **Molten globule**: Partially folded intermediate.  
3. **Native state**: Final functional conformation.

### **2.2 Computational Approaches**
- **Molecular Dynamics (MD) Simulations**:  
  - Numerically solves Newton’s equations to track atomic movements.  
  - Limited by timescale (microseconds to milliseconds).  
- **Monte Carlo Methods**:  
  - Random conformational sampling guided by energy functions.  
- **Ab Initio Prediction**:  
  - Predicts structure from sequence alone (no templates).  
  - Tools: **Rosetta**, **I-TASSER**.  
  - Energy function example:  
    \[
    E_{\text{total}} = E_{\text{van der Waals}} + E_{\text{electrostatic}} + E_{\text{solvation}} + \dots
    \]

---

## **III. Homology (Comparative) Modeling**

### **3.1 Core Principle**  
- **"Evolution preserves structure more than sequence."**  
  - If two proteins share >25% sequence identity, they likely share a similar structure.

### **3.2 Step-by-Step Workflow**  
1. **Template Identification**:  
   - Use **BLAST** or **HMMER** to find homologous PDB structures.  
2. **Sequence Alignment**:  
   - Align target and template sequences (**Clustal Omega**, **MAFFT**).  
   - Critical step: Errors here propagate to the model.  
3. **Model Building**:  
   - Copy conserved regions (e.g., α-helices, β-sheets) from the template.  
   - Tools: **MODELLER**, **SWISS-MODEL**.  
4. **Loop Modeling**:  
   - Predict flexible regions (loops) with **Rosetta**, **FALC-Loop**.  
5. **Side-Chain Placement**:  
   - Rotamer libraries optimize side-chain conformations.  
6. **Refinement**:  
   - Energy minimization to resolve steric clashes.  
7. **Validation**:  
   - **Ramachandran plot**: Checks backbone torsion angles.  
   - **MolProbity**: Assesses steric clashes and H-bond geometry.  
   - **RMSD**: Compares model to template (lower = better).

### **3.3 Case Study: Modeling a Kinase**  
- **Target**: Human kinase (unknown structure).  
- **Template**: PDB ID 1ATP (ATP-bound kinase, 35% identity).  
- **Result**: Model identifies ATP-binding pocket and catalytic residues.  

### **3.4 Limitations**  
- Fails for proteins with no homologs in PDB.  
- Low sequence identity (<20%) → unreliable models.  

---

## **IV. Advances and Tools**  
- **AlphaFold2**: Deep learning predicts structures with atomic accuracy (no templates needed).  
- **SWISS-MODEL Repository**: Automated homology modeling pipeline.  
- **CASP (Critical Assessment of Structure Prediction)**: Biennial competition driving innovation.  

---

## **V. Summary**  
- **Protein folding** relies on energy landscapes and computational brute force.  
- **Homology modeling** leverages evolutionary relationships to build 3D models.  
- **Validation** ensures models are biophysically plausible.  
