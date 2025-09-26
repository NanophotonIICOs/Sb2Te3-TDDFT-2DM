<img align="left" src="https://github.com/NanophotonIICOs/.github/blob/main/profile/images/Escudo.png">
<img align="right" src="https://github.com/NanophotonIICOs/.github/blob/main/profile/images/uaslp-1.png">
<br><br><br><br><br><br>

# Sb2Te3-TDDFT-2DM  

Research data and scripts repository for the paper:  
**“Linear optical response of monolayer Sb₂Te₃ under uniaxial strain assessed by Time-Dependent Density Functional Theory” (2D Materials, IOP Publishing 2025).**  

This repository provides **reproducibility assets**: input files (Quantum ESPRESSO & YAMBO), representative outputs, and the scripts used to process and plot the figures included in the article.  

---

## Repository contents  

- **Input files**  
  Quantum ESPRESSO (`.in`) and YAMBO (`yambo.in`) input files corresponding to the simulations reported in the paper.  

- **Representative outputs**  
  Selected processed data (e.g., band structures, dielectric function, absorption spectra) required to reproduce the figures.  
  Large raw files (`.save/`, `tmp/`, YAMBO databases) are **not included** due to size, but can be provided upon request.  

- **Processing and plotting scripts**  
  - **Python (ASE)** → for structure manipulation and atomic visualizations.  
  - **Julia (PGFPlotsX)** → for high-quality plotting of band structures and optical spectra.  
    - PGFPlotsX: [PGFPlotsX.jl](https://kristofferc.github.io/PGFPlotsX.jl/stable/)  
  - **LaTeX/TikZ** → for figure composition, axis labeling, and final layout.  

- **Institutional tools**  
  Pre-processing of Quantum ESPRESSO outputs was done using scripts from the institutional repository [AtomistIICO](https://github.com/NanophotonIICOs/AtomistIICO).  

---


## Data policy  

- Only **input files, selected outputs, and scripts** are included here.  
- Full raw datasets from QE/YAMBO are excluded for size reasons but can be made available upon request.  

---

## Statement on AI usage  

> ⚠️ **No generative AI tools were used in the production of this repository or its figures.**  
> All results and visualizations are based on first-principles calculations (DFT, TDDFT, GW/BSE) and processed with standard scientific workflows (Python/ASE, Julia/PGFPlotsX, LaTeX/TikZ).  

---


<h2>Acknowledgments</h2>

<p>We gratefully acknowledge the support provided by:</p>

<div style="display: flex; align-items: center; gap: 40px; flex-wrap: wrap;">
  <div style="text-align: center;">
    <img src="https://secihti.mx/wp-content/uploads/2024/12/logotipo_SCyT_color_803x97px_v02.svg" alt="SECTIHI Logo" style="height: 80px;">
    <p><strong>Secretaría de Ciencia, Humanidades, Tecnología e Innovación</strong></p>
  </div>
  
  <div style="text-align: center;">
    <img src="http://registro.lnsa.buap.mx/imagenes/LNS.png" alt="LNS Logo" style="height: 80px;">
    <p><strong>Laboratorio Nacional de Supercómputo del Sureste de México</strong></p>
  </div>
</div>

