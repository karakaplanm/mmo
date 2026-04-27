# CXCR4 Protein Model (AlphaFold Prediction)

## File Information
- **File Name:** UgurCemYildiz_cxcr4_model.cif  
- **Format:** CIF (Crystallographic Information File)  
- **Model Type:** In silico prediction (AlphaFold)

## Overview
This file contains the predicted three-dimensional structure of the CXCR4 (C-X-C Chemokine Receptor Type 4) protein generated using the AlphaFold algorithm.

AlphaFold is an artificial intelligence–based system that predicts protein structures from amino acid sequences with high accuracy. This model is computational and not experimentally determined (e.g., X-ray crystallography, NMR, Cryo-EM).

## Biological Context
CXCR4 is a G protein–coupled receptor (GPCR) involved in several key biological processes:
- Cell migration  
- Immune system regulation  
- Atherosclerosis development  
- Cancer metastasis  

Primary ligand: **CXCL12 (SDF-1)**

## Structural Features (Expected)
- Seven transmembrane helices (typical GPCR architecture)  
- Extracellular ligand-binding regions  
- Intracellular domains for G-protein interaction  

## Model Confidence
AlphaFold provides confidence scores using **pLDDT**:
- >90 → Very high confidence  
- 70–90 → Reliable  
- 50–70 → Low confidence  
- <50 → Likely disordered / unreliable regions  

> Note: Loop regions and extracellular domains often have lower confidence.

## Use Cases
This model can be used for:
- Structural visualization (PyMOL, ChimeraX, VMD)  
- Preliminary ligand docking (with caution)  
- Mutation impact analysis  
- Domain organization studies  

## Limitations
- Not experimentally validated  
- Does not capture protein dynamics  
- Ligand-bound conformation is not guaranteed  
- Flexible regions (especially in GPCRs) may be inaccurate  

## Reference
- Jumper et al., 2021, Nature — AlphaFold
