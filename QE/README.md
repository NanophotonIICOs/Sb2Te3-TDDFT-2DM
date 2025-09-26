# Quantum ESPRESSO calculations

This folder contains the scripts and input files used to run **Quantum ESPRESSO** calculations on the LNS cluster.  
The main goal of these calculations is to generate the `SAVE/` directory required by **YAMBO**, and to provide the electronic structure data (with and without SOC) that is later used for GW, BSE, and TDDFT simulations.

---

## Subfolders

- `bands-pbe/`  
  Contains the SLURM job script and input generators for **SCF + NSCF runs** using the GGA-PBE functional.  
  Band structures are calculated **with and without spin–orbit coupling (SOC)**.  
  At the end of the NSCF run, a `prefix.save/` directory (e.g. `Sb2Te3.save/`) is created.  

- `pbe-yambo-inputs/`  
  Contains the scripts used to generate **YAMBO-ready inputs** from Quantum ESPRESSO outputs.  
  These scripts prepare the required files for:  
  - **G0W0 calculations**  
  - **Bethe–Salpeter Equation (BSE)**  
  - **Time-Dependent Density Functional Theory (TDDFT)**  

---

## Purpose

* The `pbe-yambo-inputs/` scripts convert these results into the proper `SAVE/` directory structure that YAMBO requires.  
* This workflow ensures that YAMBO can directly initialize and run G0W0, BSE, and TDDFT calculations based on consistent DFT inputs.  

---

## Notes on reproducibility

- The scripts are tailored to the **Huapáctic (LNS) cluster**, using SLURM job submission and environment modules (GCC, OpenMPI, FFTW, etc.).  
- Users on other HPC systems should adapt the `module load` section and SLURM headers as needed.  
- ***No post-processing or artificial smoothing is applied: the scripts simply run QE and prepare the data for YAMBO.*** 

---