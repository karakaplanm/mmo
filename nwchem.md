# NWChem

NWChem is an open-source computational chemistry software package designed to run on high-performance parallel computing systems. Developed at the Pacific Northwest National Laboratory (PNNL) under the U.S. Department of Energy, it is widely used for quantum chemistry, molecular dynamics, and other simulations in chemical and material sciences.
Features of NWChem:

+ Quantum Chemistry Methods:
  - Hartree-Fock (HF)
  - Density Functional Theory (DFT)
  - Møller–Plesset Perturbation Theory (MP2)
  - Coupled Cluster (CC) methods
  - Configuration Interaction (CI)
  - Multireference methods like CASSCF

+ Molecular Dynamics (MD):
  - Classical molecular dynamics simulations using force fields
  - QM/MM (Quantum Mechanics/Molecular Mechanics) hybrid simulations

+ Molecular and Periodic System Calculations:
 - Supports molecular and periodic boundary conditions (PBC) for solid-state simulations

+ Scalability and Parallel Performance:
  - Designed to scale efficiently on parallel architectures
  - Uses MPI for distributed computing

+ Advanced Computational Capabilities:
  - Excited-state calculations (TDDFT, EOM-CC)
  - Relativistic effects (ZORA, DKH)
  - Solvation models (COSMO, PCM)

+ Extensive Basis Set Library:
  - Includes a comprehensive set of Gaussian basis sets
  - Users can define custom basis sets

## Applications:
- Computational studies in catalysis, materials science, and spectroscopy
- Biomolecular simulations (e.g., protein-ligand interactions)
- Drug discovery research
- Electronic structure calculations of complex systems
    

## Installation 

`$ sudo apt install nwchem`

## Optimization Example

Create a file named water.nw and paste the following code.
```
start water_opt

echo

geometry units angstroms
  O  0.000000  0.000000  0.000000
  H  0.757000  0.586000  0.000000
  H -0.757000  0.586000  0.000000
end

basis
  * library 6-31G**
end

task scf optimize
```

```
$ nwchem water.nw
```

or write output to a file

```
$ nwchem water.nw > water.out
```

## DFT Calculation

```
start water
memory 500 mb
geometry
  O  0.000  0.000  0.000
  H  0.757  0.586  0.000
  H -0.757  0.586  0.000
end
basis
  O library 6-31G
  H library 6-31G
end
dft
  xc b3lyp
end
task dft energy
```
