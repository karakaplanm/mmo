# **Week 13: Introduction to Ab Initio Methods**

## **1. Overview of Ab Initio Methods**
**Ab initio** (Latin for "from first principles") methods solve quantum mechanical equations without empirical parameters. Key features:

- Based on **fundamental physical laws** (Schrödinger equation)
- Require only:
  - Atomic numbers
  - Positions of nuclei
- Provide **high accuracy** for:
  - Molecular structures
  - Electronic properties
  - Reaction mechanisms

---

## **2. Fundamental Theory**

### **A. Schrödinger Equation**
The cornerstone of ab initio methods:

$$
\hat{H}\Psi(\mathbf{r},\mathbf{R}) = E\Psi(\mathbf{r},\mathbf{R})
$$

where:
- $\hat{H}$ = Hamiltonian operator
- $\Psi$ = Wavefunction
- $\mathbf{r}$ = Electron coordinates
- $\mathbf{R}$ = Nuclear coordinates

### **B. Born-Oppenheimer Approximation**
Separates electronic and nuclear motion:
1. Fix nuclear positions
2. Solve electronic Schrödinger equation
3. Repeat for different nuclear configurations

---

## **3. Hierarchy of Ab Initio Methods**

| Method | Description | Accuracy | Cost |
|--------|-------------|----------|------|
| **Hartree-Fock (HF)** | Mean-field approximation | Moderate | Low |
| **Post-HF Methods** | Include electron correlation | High | Very High |
| **DFT** | Electron density functional | Good | Medium |
| **Coupled Cluster** | Gold standard for small systems | Very High | Extreme |

---

## **4. Key Methods Explained**

### **A. Hartree-Fock Theory**
- Solves for molecular orbitals via SCF procedure
- Limitations:
  - No electron correlation
  - Overestimates energies

**SCF Cycle**:
1. Guess initial orbitals
2. Build Fock matrix
3. Solve Roothaan equations
4. Check convergence
5. Repeat until converged

### **B. Post-HF Methods**
Include electron correlation:
- **MP2**: Møller-Plesset perturbation (2nd order)
- **CCSD**: Coupled Cluster with Singles/Doubles
- **CI**: Configuration Interaction

### **C. Density Functional Theory (DFT)**
- Uses electron density instead of wavefunction
- Popular functionals:
  - B3LYP (hybrid)
  - PBE (GGA)
  - ωB97X-D (dispersion-corrected)

---

## **5. Computational Implementation**

### **Typical Workflow**
1. **Geometry Optimization**
   - Find minimum energy structure
2. **Frequency Calculation**
   - Verify minima (no imaginary frequencies)
3. **Property Calculation**
   - Energies, spectra, orbitals

### **Common Software**
- **Gaussian**
- **ORCA**
- **Q-Chem**
- **PySCF** (Python-based)

---

## **6. Applications**

### **A. Molecular Properties**
- Dipole moments
- Vibrational frequencies
- NMR chemical shifts

### **B. Chemical Reactions**
- Reaction pathways
- Transition states
- Activation energies

### **C. Materials Science**
- Band structures
- Defect properties
- Catalytic surfaces

---

## **7. Advantages and Limitations**

### **Advantages**
✔ Parameter-free  
✔ Systematic improvability  
✔ Accurate for small systems  
✔ Can predict new phenomena  

### **Limitations**
✖ High computational cost  
✖ Scaling with system size (O(N⁴) to O(N⁷))  
✖ Basis set convergence issues  

---

## **8. Current Developments**

- **Linear-scaling methods**
- **Embedding techniques** (QM/MM)
- **Machine learning potentials**
- **Quantum computing applications**

---

## **9. Practical Example: Water Molecule**

```python
# PySCF example for H2O calculation
from pyscf import gto, scf

mol = gto.M(
    atom = 'O 0 0 0; H 0 1 0; H 0 0 1',
    basis = 'cc-pVDZ'
)

mf = scf.RHF(mol)
mf.kernel()
print(f"Total energy: {mf.e_tot:.6f} Hartree")
