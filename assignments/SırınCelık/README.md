# 3OSK - CTLA-4 Homodimer Structure (AlphaFold Prediction)

This repository contains the structural model of human CTLA-4 (Cytotoxic T-Lymphocyte Associated Protein 4) folded using AlphaFold. The model is based on the crystallographic structure deposited under PDB ID 3OSK, which represents the apo (ligand-free) homodimer form of the protein. The folding run resulted in a pTM score of 0.82. This value indicates a high-confidence prediction and confirms that the dimeric organization of the protein has been correctly modeled by AlphaFold.

## Biological Background

CTLA-4 is an inhibitory receptor that plays a central role in the regulation of the immune system and is expressed on the surface of T cells. Its primary function is to downregulate T cell activation by binding to the costimulatory ligands B7-1 (CD80) and B7-2 (CD86), which are present on the surface of antigen-presenting cells. Through this mechanism, the body prevents excessive immune responses and avoids autoimmune reactions. Although CTLA-4 is structurally similar to CD28, the two receptors have opposing functions: CD28 transmits an activating signal, whereas CTLA-4 delivers an inhibitory one. Both receptors bind to the same ligands, but CTLA-4 does so with significantly higher affinity than CD28, making it the dominant negative regulator under physiological conditions.

## Structural Features

CTLA-4 belongs to the immunoglobulin superfamily and exhibits an immunoglobulin-like beta-sandwich fold. In solution and on the cell surface, the protein exists as a non-covalently stabilized homodimer. Structural studies on the 3OSK crystal form have demonstrated that CTLA-4 does not undergo large conformational changes upon ligand binding. Instead, it operates through a rigid-body recognition mechanism, meaning the receptor possesses a pre-organized binding surface and does not require structural rearrangement to engage its ligands. The pTM value of 0.82 obtained here confirms that the AlphaFold model successfully captures this dimeric architecture.

## Biotechnological and Clinical Relevance

The significance of CTLA-4 in biotechnology and clinical medicine is substantial. In cancer immunotherapy, monoclonal antibodies targeting CTLA-4 (such as ipilimumab) are employed as immune checkpoint inhibitors to reactivate the immune system against tumor cells. Conversely, for the treatment of autoimmune diseases, the inhibitory function of CTLA-4 is leveraged therapeutically. CTLA-4-Fc fusion proteins (such as abatacept) are used to suppress immune responses in conditions like rheumatoid arthritis. Accurate structural models of CTLA-4 are therefore critically important for drug discovery and development efforts.

## Repository Contents

This dataset includes the following files:

- **3osk_alphafold.pdb** - The raw AlphaFold-predicted structural model
- **3osk_relaxed.pdb** - Energy-minimized structure
- **metrics.json** - Computed pTM and ipTM scores
- **README.md** - This documentation file

## Usage Notes

This model can serve as a starting point for advanced molecular dynamics simulations, molecular docking studies, or structure-based drug design analyses. It should be noted that CTLA-4 is not a monomeric protein but a functional homodimer, and all downstream structural analyses must take this dimeric organization into account. The AlphaFold-generated model also provides a valuable resource for comparative analyses with the experimentally determined 3OSK crystal structure.

## References

- Yu, C., et al. Rigid-body ligand recognition drives CTLA-4 function. (2010)
- PDB Entry: [3OSK](https://www.rcsb.org/structure/3OSK)