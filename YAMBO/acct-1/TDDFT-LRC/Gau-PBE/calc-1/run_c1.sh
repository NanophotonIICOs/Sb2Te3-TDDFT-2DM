#!/bin/bash
#SBATCH -J tddft             # Nombre del trabajo
#SBATCH -o yo.o-%j            # Archivo de salida (%j se expande al jobID)
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=24 # number of tasks per node (maximum 24)
#SBATCH --exclude=cn0605-3
#SBATCH -p normal # SLURM queue (partition)
#SBATCH -t 120:00:00 # run time (hh:mm:ss)

module purge
#------------ GCC
module load GCC/12.2.0  
module load OpenMPI/4.1.4-GCC-12.2.0 
module load FFTW/3.3.10-GCC-12.2.0
module load ScaLAPACK/2.2.0-gompi-2022b-fb
module load OpenBLAS/0.3.21-GCC-12.2.0
# export LD_LIBRARY_PATH=/
#------------ INTEL
# module load intel-compilers/2021.4.0
# module load impi/2021.4.0-intel-compilers-2021.4.0
# module load imkl-FFTW/2021.4.0-iimpi-2021b
# module load HDF5/1.12.1-iimpi-2021b
# module load libxc/5.1.6-intel-compilers-2021.4.0 
#------------
export OMP_NUM_THREADS=1
source ~/scripts/files.sh

#---------------------------------JOb names and Dirs
cd $SLURM_SUBMIT_DIR
#---------------------------------Print Parallel Params
RJOB=$SLURM_JOB_NAME
NTASKS=$SLURM_NTASKS
#NTASKS=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))
NNODES=$SLURM_NNODES
NTASKPN=$SLURM_TASKS_PER_NODE
CALC="yambo"
print_parallel_params "$RJOB" "$NTASKS" "$NNODES" "$NTASKPN" "$CALC"

for alpha in -0.4 -0.3 -0.2 -0.1; do
cat > ox_${alpha}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
FFTGvecs= 50             Ry    # [FFT] Plane-waves
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= ${alpha}              # [TDDFT] LRC alpha factor
NGsBlkXd= 3               Ry    # [Xd] Response block size
% QpntsRXd
   1 | 1 |                         # [Xd] Transferred momenta
%
% BndsRnXd
   1 | 60 |                         # [Xd] Polarization function bands
%
% EnRngeXd
  0.00000 | 5.00000 |         eV    # [Xd] Energy range
%
% DmRngeXd
 0.050000 | 0.050000 |         eV    # [Xd] Damping range
%
ETStpsXd= 1500                    # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="ox_${alpha}"
mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J ${prefix}_run -C ${prefix}_out
echo "Alpha=$alpha terminado: $(date)"
done
echo "End job:$(date)"
