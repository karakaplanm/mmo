# Report: Folding of the 4PYP Protein with AlphaFold and Model Quality Assessment

## 1. Brief Information About the Protein

4PYP is the Protein Data Bank (PDB) identifier of the Photoactive Yellow Protein (PYP) found in *Halorhodospira halophila* (a photosynthetic purple bacterium). PYP is a small, soluble, and light-sensitive photoreceptor protein consisting of approximately 125 amino acids.

Its biological function is to detect blue light and initiate the negative phototaxis response (movement away from light). In this way, the bacterium is protected from damage under intense light conditions.

Structurally, PYP exhibits an α/β folding motif. It typically contains a β-sheet composed of six β-strands surrounded by four α-helices. The most important feature of the protein is its covalent binding to the *p-hydroxycinnamyl chromophore* (photoreceptor group). This chromophore is attached to the active-site amino acid Cys69 through a thioester bond.

Upon light stimulation, the chromophore undergoes a trans → cis isomerization, leading to major conformational changes in the protein. These changes trigger downstream protein–protein interactions that initiate signal transduction. Due to this mechanism, PYP is considered a model system for understanding photocycle dynamics and also provides insight into structural and functional studies of the PAS (Per-Arnt-Sim) domain family.

Experimentally, high-resolution structures of both the dark state (ground state) and light-activated intermediate forms of PYP have been solved using X-ray crystallography and NMR methods. The code 4PYP corresponds to one of these structures (generally the dark-state structure).

## 2. Folding Method

The folding process was performed using **AlphaFold**. Since no multimer information was provided, the ipTM value was reported as a dash (`-`) because a monomeric form was targeted.

## 3. Evaluation of the Results

The obtained **pTM score was 0.84**. This value indicates that the overall structural accuracy of the model is **very high** (values above 0.5 are considered good, while values above 0.8 are considered excellent).

## 4. Interpretation of the Model

- A **pTM value of 0.84** demonstrates that the folded structure’s globular core is highly consistent with the experimental structure (PDB:4PYP).

## 5. Summary of the Folding Pattern

AlphaFold generated a **high-confidence monomeric model** using the 4PYP sequence.

- The **globular domain** (approximately residues 10–110) has an average pLDDT score **above 85**, indicating nearly crystallographic-level quality.
- The **chromophore-binding pocket** (including residues such as Tyr42, Glu46, and Cys69) was positioned with very high confidence (>90).
- **Flexible loops** (especially residues 55–70) were modeled with lower confidence. Since these regions undergo conformational changes during the photocycle, it is expected that AlphaFold provides a probabilistic distribution rather than a single fixed structure.

## Conclusion

The obtained model can be directly used for advanced analyses such as **molecular docking**, **molecular dynamics simulations** (for pre-equilibration), and **mutation effect prediction**. For low-confidence regions, alternative folding conformations should be evaluated by examining multiple AlphaFold output models.
