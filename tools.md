# Open-Source Molecular Mechanics Tools for Ubuntu

Molecular mechanics (MM) tools allow simulation and analysis of molecular systems using classical force fields. Below is a list of popular **open-source** MM tools that run on **Ubuntu/Linux**.

---

## 1. GROMACS

- **Language**: C/C++
- **Features**:
  - Molecular mechanics
  - Energy minimization
  - Molecular dynamics (MD)
  - Solvation
- **Force Fields**: AMBER, CHARMM, OPLS-AA, GROMOS
- **Best For**: High-performance MD of proteins, membranes, and biomolecular complexes
- **Installation**:

  ```bash
  sudo apt install gromacs
  ```

## 2. OpenMM

- Language: Python/C++
- Features:
    - GPU-accelerated MM and MD
    - Custom force fields via XML
    - Real-time and scripted simulations
- Force Fields: AMBER, CHARMM, GAFF
- Best For: Python-based flexible simulations and machine learning integration
- Installation (via conda):
  ```
  conda install -c conda-forge openmm
  ```

## Avogadro (with Open Babel)

- Type: GUI application for molecular modeling
- Features:
    - Interactive molecule builder
    - MM geometry optimization
- Force Fields: MMFF94, UFF, GAFF (via Open Babel)
- Best For: Education, quick MM modeling of small molecules
- Installation:
  ```bash
  sudo apt install avogadro openbabel
  ```
