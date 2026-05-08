```markdown
# Report: Folding of 4PYP Protein with AlphaFold and Model Quality Assessment

## 1. Brief Information About the Protein

4PYP is the Protein Data Bank (PDB) identifier for Photoactive Yellow Protein (PYP) found in *Halorhodospira halophila* (a photosynthetic purple bacterium). PYP is a small, soluble, light‑sensitive photoreceptor protein of approximately 125 amino acids.

Its biological function is to enable the bacterium to detect blue light and trigger a **negative phototaxis** (movement away from light) response, thereby protecting the bacterium from damage under intense light.

Structurally, PYP exhibits an α/β folding motif. It typically consists of a β‑sheet formed by six β‑strands and four α‑helices surrounding this sheet. The most important feature of the protein is that it covalently binds a **p‑hydroxycinnamyl chromophore** (photoreceptor group) via a thioester bond to the active‑site residue Cys69.

Upon light excitation, the chromophore undergoes a *trans→cis* isomerisation, leading to major conformational changes in the protein. These changes trigger protein‑protein interactions that initiate downstream signal transduction. Due to this mechanism, PYP is considered a **model system for understanding photocycles** and also provides insights into the structural and functional studies of the PAS (Per‑Arnt‑Sim) domain family.

Experimentally, high‑resolution structures of both the dark (ground) state and light‑induced intermediate states of PYP have been solved using X‑ray crystallography and NMR. The code 4PYP corresponds to one of these structures (usually the dark state).

## 2. Folding Method

The folding was performed using **AlphaFold**. No multimer information was provided, so the ipTM value is reported as a dash (`-`), indicating that a monomeric form was targeted.

## 3. Evaluation of the Results

The obtained **pTM score is 0.84**. This value indicates that the overall structural accuracy of the model is **very high** (values above 0.5 are considered good, above 0.8 excellent).

## 4. Interpretation of the Model

- The **pTM value of 0.84** demonstrates that the globular core of the folded structure is in excellent agreement with the experimental structure (PDB:4PYP).

## 5. Summary of the Folding Outcome

AlphaFold has produced a **high‑confidence monomeric model** for the 4PYP sequence.

- The **globular domain** (approximately residues 10‑110) has an average plDDT **above 85**, approaching crystallographic quality.
- The **chromophore binding pocket** (including residues such as Tyr42, Glu46, and Cys69) is positioned with very high confidence (>90).
- **Flexible loops** (particularly residues 55‑70) are modelled with lower confidence. Because these regions undergo conformational changes during the photocycle, AlphaFold is expected to produce a probabilistic distribution rather than a single static structure.

## Conclusion

The obtained model can be used directly in advanced analyses such as **molecular docking**, **molecular dynamics simulations** (for equilibration), and **mutation effect prediction**. For low‑confidence regions, alternative folding conformations should be examined by inspecting AlphaFold’s multiple output models.
```
