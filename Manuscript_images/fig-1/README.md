# Figure 1 (Atomic Structure + Band Structure)

This folder contains the files used to generate Figure 1 of the manuscript, which combines the atomic structure of monolayer Sb$_2$Te$_3$ and its corresponding band structure.

Purpose

The figure is designed to provide a clear and consistent visualization of:

	1.	The atomic structure of Sb$_2$Te$_3$ is displayed as a 5×5×1 supercell, with reference axes provided to clarify the orientation of the applied uniaxial strain.
	2.	The electronic band structure, computed with DFT and plotted (with Julia/PGFPlotsX.jl) for visual clarity.
	3.	A final composition that joins both panels into a single figure, suitable for publication.

Method

* **Atomic structure**  
  - Generated using the Atomic Simulation Environment (ASE) in Python.  
  - A 5×5×1 supercell was created from the unit cell.  
  - The image was exported directly from ASE visualization tools.  

* **Band structure**  
  - Numerical data were processed and plotted in Julia using the PGFPlotsX package, which produces LaTeX-quality vector graphics.  
  - No smoothing or fitting was applied; the plots are a direct rendering of the calculated eigenvalues.  

* **Final composition**  
  - The structure image and the band-structure plots were combined in a single panel using LaTeX with TikZ.  
  - Reference axes, labels, and annotations indicating the direction of the applied strain were added within TikZ.  
  - This ensured consistent fonts, vector quality, and a unified figure layout across all subpanels.  

LaTeX workflow

	1. Sb2Te3-sc.tex 
	Generates a PDF containing the structure image rendered with ASE.
	2. fig_center.tex
	Takes the PDF output from Sb2Te3-sc.tex and overlays axis references and labels.
	3. fig-1.tex
	Combines the processed structure figure (fig_center.tex) with the band structure plots generated in Julia/PGFPlotsX (with and without SOC) to produce the final Figure 1 included in the manuscript.



> :warning: **Important note**  
> No generative AI tools were used in this figure. All visualizations come from ASE, Julia/PGFPlotsX, and LaTeX/TikZ.