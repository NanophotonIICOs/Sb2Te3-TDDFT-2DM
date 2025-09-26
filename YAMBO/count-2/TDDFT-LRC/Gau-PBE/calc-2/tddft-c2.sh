#!/bin/bash
#SBATCH -J tddft-c1             # Job name             
#SBATCH -o yambo-tddft.o-%j           
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12 # number of tasks per node (maximum 24)
#SBATCH --cpus-per-task=2       #  threads  
#SBATCH -p normal # SLURM queue (partition)
#SBATCH --exclude=cn0605-3
#SBATCH -t 120:00:00 # run time (hh:jobmm:ss)

module purge
# module load GCC/12.2.0  
# module load OpenMPI/4.1.4-GCC-12.2.0 
# module load FFTW/3.3.10-GCC-12.2.0
# module load ScaLAPACK/2.2.0-gompi-2022b-fb
# module load OpenBLAS/0.3.21-GCC-12.2.0
# export LD_LIBRARY_PATH=/zfs-home/202401017n/software/libxc-6.2.2/lib:$LD_LIBRARY_PATH

module load intel-compilers/2021.4.0
module load impi/2021.4.0-intel-compilers-2021.4.0
module load imkl-FFTW/2021.4.0-iimpi-2021b
module load HDF5/1.12.1-iimpi-2021b
module load libxc/5.1.6-intel-compilers-2021.4.0 

#--------------------------------- OpenMPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
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
echo "OMP_NUM_THREADS              : ${OMP_NUM_THREADS}"
echo "YAMBOROOT                    : ${YAMBOINTELROOT}"
echo " "

for ngsBlkXd in 3 4 5 6 7; do
prefix="ox_${ngsBlkXd}_Ry"
echo "Input file                   : ${prefix}.in"
echo "Job directory                : ${prefix}_run"
echo "Output directory             : ${prefix}_out"
cat > ${prefix}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4              # [TDDFT] LRC alpha factor
NGsBlkXd= ${ngsBlkXd}            Ry    # [Xd] Response block size
% QpntsRXd
   1 | 1 |                         # [Xd] Transferred momenta
%
% BndsRnXd
   1 | 100 |                         # [Xd] Polarization function bands
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
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J ${prefix}_run -C ${prefix}_out
echo "ngsBlkXd=${ngsBlkXd} Ry terminado: $(date)"
done
echo "End job:$(date)"
