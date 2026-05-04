# Ab Initio Structure Prediction of Hemoglobin α2β2 Tetramer

**Author:** Helin Yeliz Duman

**Institution:** Inonu University - Department of Molecular Biology and Genetics

**Course:** Molecular Modeling

## Project Overview
In this project, the 3D structure of human Hemoglobin, a critical metalloprotein responsible for oxygen transport in vertebrates, was modeled ab initio. The functional hemoglobin complex is a heterotetramer consisting of two alpha (α) and two beta (β) subunits (Chains A, C and B, D, respectively). Hemoglobin is a cornerstone of structural biology, illustrating key concepts such as allosteric regulation, cooperativity in ligand binding, and the transition between Tense (T) and Relaxed (R) conformational states.

## Methodology
* The structural modeling was performed using Google DeepMind AlphaFold 3.
* The target sequences (PDB ID: 1GZX) were submitted as a multimeric assembly to predict the native tetrameric state.
* Among the 5 generated conformations, model_0 was selected as the reference structure. Although model_3 showed a slightly higher global pLDDT, model_0 was chosen for its optimal ranking in the overall thermodynamic scoring pipeline.
* Structural validation metrics, including pLDDT and Predicted Aligned Error (PAE), were retrieved from the summary_confidences_0.json file.

## Structural Analysis and Visualization
The predicted .cif model was visualized in PyMOL and analyzed based on per-residue confidence scores (pLDDT):
* **Blue (pLDDT > 90):** Indicates high-accuracy structural prediction. The majority of the α-helical globin domains, which constitute the core of the protein, fall into this category.
* **Cyan/Yellow (90 > pLDDT > 50):** Represents lower confidence regions, primarily localized in the flexible loop segments connecting the helical bundles.
* **Highlighted Active Center:** The heme-binding pockets, crucial for oxygen coordination, were specifically examined. The proximal and distal Histidine residues (His58, His87 for α; His63, His92 for β) were highlighted as red sticks/spheres to showcase the architecture of the oxygen-binding site.

## RMSD Validation
To evaluate the accuracy of the ab initio prediction, the generated model_0 was superimposed onto the experimental crystal structure of human hemoglobin (PDB ID: 1GZX). The alignment yielded an RMSD of 0.583 Å (over 1868 atoms).This exceptionally low RMSD value confirms that AlphaFold 3 successfully captured the complex quaternary arrangement and side-chain orientations of the hemoglobin tetramer.

## Repository Contents
* **fold_2026_05_03_19_25_model_0.cif:** The reference 3D structural model.
* **fold_2026_05_03_19_25_summary_confidences_0.json:** Detailed confidence and error matrices.
* **hemoglobin_pymol_render.png:** High-resolution visualization showcasing pLDDT coloring and highlighted active sites.
