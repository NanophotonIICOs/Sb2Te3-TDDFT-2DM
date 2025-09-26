# YAMBO Calculations (TDDFT, G0W0, BSE)

This folder contains the inputs and representative outputs of the **YAMBO** calculations performed for the manuscript.  
The simulations include **Time-Dependent Density Functional Theory (TDDFT)**, **G0W0 quasiparticle corrections**, and **Bethe–Salpeter Equation (BSE)** optical spectra.  

---

## Cluster accounts

Due to restrictions on the maximum number of nodes per job in the cluster, the calculations were distributed across **two different user accounts**:

- **acct-1** &#8594; TDDFT calculations, including the $\alpha$-convergence tests, as well as the G0W0 and BSE runs  
- **acct-2** &#8594; complementary TDDFT calculations  

This splitting was purely a **technical requirement of the HPC allocation policy**.  
All calculations used the **same computational parameters** (cutoffs, number of bands, k-point meshes, etc.) to ensure consistency of the results.  

---

## Contents

- `acct-1/` &#8594; TDDFT (including $\alpha$-convergence tests), G0W0, and BSE runs performed under account 1  
- `acct-2/` &#8594; TDDFT runs performed under account 2  
- Representative outputs (`ndb.QP`, `o-eh.eps`, dielectric function data) and plotting scripts  

---

## TDDFT $\alpha$-convergence tests

Within `acct-1/`, a dedicated set of calculations was performed to verify the convergence of the **$\alpha$ parameter** in TDDFT.  
This set includes:  

- YAMBO input files with varying $\alpha$ values  
- Output data used to assess the stability of the optical response with respect to $\alpha$  
- Python scripts for parsing and analyzing the YAMBO outputs, as well as generating convergence plots  

These tests guided the choice of the representative $\alpha$ value employed in the optical spectra reported in the manuscript.  

---

## Reproducibility

The use of two accounts does **not affect the reproducibility or validity of the results**.  
All runs were performed with the same code version (YAMBO 5.3.0), compiled with the same libraries and cluster modules.  
The split into two accounts is only reflected in the directory organization of this repository.  

---