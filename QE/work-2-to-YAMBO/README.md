This folder contains the SLURM job script used on the LNS cluster to run Quantum ESPRESSO calculations required to generate the SAVE directory.

Purpose

The script runs both the SCF and NSCF stages of Quantum ESPRESSO using the provided input files.
At the end of the NSCF run, a prefix.save/ directory (e.g. Sb2Te3.save/) is created.

This SAVE/ directory is essential for YAMBO, as it stores the eigenvalues, wavefunctions, and electronic structure information needed to initialize GW and BSE calculations.