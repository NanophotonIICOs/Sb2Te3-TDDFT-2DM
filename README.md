<table width="100%" border="0" cellspacing="0" cellpadding="0">
  <tr>
    <td align="left" width="20%">
      <img src="https://github.com/NanophotonIICOs/.github/blob/main/profile/images/Escudo.png" width="100%">
    </td>
    <td align="center" width="60%">
      <img src="figure-workflow/Sb2Te3-sc.gif" width="50%">
    </td>
    <td align="right" width="20%">
      <img src="https://github.com/NanophotonIICOs/.github/blob/main/profile/images/uaslp-1.png" width="100%">
    </td>
  </tr>
</table>

<br><br>
# Sb2Te3-TDDFT-2DM  

Research data and scripts repository for the paper:  
**“Linear optical response of monolayer Sb2Te3 under uniaxial strain assessed by Time-Dependent Density Functional Theory” (2D Materials, IOP Publishing 2025).**  

This repository provides **reproducibility resources**: input files (Quantum ESPRESSO & YAMBO), representative outputs, and the scripts used to process and plot the figures included in the article.  

---

## Repository contents  

- **Input files**  
  Quantum ESPRESSO (`.in`) and YAMBO (`yambo.in`) input files corresponding to the simulations reported in the paper.  

- **Representative outputs**  
  Selected processed data (e.g., band structures, dielectric function, absorption spectra) required to reproduce the figures.  
  Large raw files (`.save/`, `tmp/`, YAMBO databases) are **not included** due to size, but can be provided upon request.  

- **Processing and plotting scripts**  
  - **Python (ASE)** &#8594; for structure manipulation and atomic visualizations.  
  - **Julia (PGFPlotsX)** &#8594; for high-quality plotting of band structures and optical spectra.  
    - PGFPlotsX: [PGFPlotsX.jl](https://kristofferc.github.io/PGFPlotsX.jl/stable/)  
  - **LaTeX/TikZ** &#8594; for figure composition, axis labeling, and final layout.  

- **Institutional tools**  
  Post-processing of Quantum ESPRESSO outputs was done using scripts from the institutional repository [AtomistIICO](https://github.com/NanophotonIICOs/AtomistIICO).  

---


## Data policy  

- Only the files required to reproduce the results presented in the manuscript are included here:  
**input files, selected representative outputs, and the corresponding processing/plotting scripts.**  
- Large intermediate data (e.g. Quantum ESPRESSO `.save/` directories, temporary files, or YAMBO databases) are excluded due to size, but can be made available upon reasonable request.

---

## Statement on AI usage  

> ⚠️ **No generative AI tools were used in the production of this repository or its figures.**  
> All results and visualizations are based on first-principles calculations (DFT, TDDFT, GW/BSE) and processed with standard scientific workflows (Python/ASE, Julia/PGFPlotsX, LaTeX/TikZ).  

---

## Contact  

For questions, clarifications, or requests regarding the data and scripts:  

- **Dr. Oscar Ruiz Cigarrillo**  
  Instituto de Investigación en Comunicación Óptica (IICO)  
  Universidad Autónoma de San Luis Potosí (UASLP), México  
  📧 Email: [oscar.ruiz@uaslp.mx](mailto:oscar.ruiz@uaslp.mx)  
- **Dr. Raul Eduardo Balderas Navarro**  
  Instituto de Investigación en Comunicación Óptica (IICO)  
  Universidad Autónoma de San Luis Potosí (UASLP), México  
  📧 Email: [raul.balderas@uaslp.mx;](mailto:raul.balderas@uaslp.mx;)  

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

