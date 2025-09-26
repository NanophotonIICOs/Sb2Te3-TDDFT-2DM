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
#--------------------------------- Job names and Dirs
cd ${SLURM_SUBMIT_DIR}
#--------------------------------- Print Job
source ~/scripts/pms2slurm.sh
RJOB=${SLURM_JOB_NAME}
CALC="yambo"
print_job "$RJOB" "$CALC" "${SLURM_SUBMIT_DIR}"

#--------------------------------- Memory Calculation Function
calculate_memory_per_core() {
    local nodes=$1
    local ntasks_per_node=$2
    local cpus_per_task=$3
    local memory_per_node_gb=${4:-128}  # Default 128GB if not specified
    
    # Calculate total tasks, cores and memory
    local total_tasks=$((nodes * ntasks_per_node))
    local total_cores=$((total_tasks * cpus_per_task))
    local total_memory_gb=$((nodes * memory_per_node_gb))
    
    # Memory calculations
    local memory_per_core_gb=$((total_memory_gb / total_cores))
    local memory_per_core_mb=$((memory_per_core_gb * 1024))
    local memory_per_task_gb=$((total_memory_gb / total_tasks))
    echo "== MEMORY ALLOCATION SUMMARY =="
    echo "Nodes                        : $nodes"
    echo "Tasks per node               : $ntasks_per_node"
    echo "CPUs per task                : $cpus_per_task"
    echo "Memory per node              : ${memory_per_node_gb} GB"
    echo "Total tasks                  : $total_tasks"
    echo "Total cores                  : $total_cores"
    echo "Total memory                 : ${total_memory_gb} GB"
    echo "Memory per task              : ${memory_per_task_gb} GB"
    #echo "Memory per core              : ${memory_per_core_gb} GB (${memory_per_core_mb} MB)"
    echo "==============================="
    
    # Export for use in calculations if needed
    export MEMORY_PER_CORE_GB=$memory_per_core_gb
    export MEMORY_PER_CORE_MB=$memory_per_core_mb
    export MEMORY_PER_TASK_GB=$memory_per_task_gb
    export TOTAL_TASKS=$total_tasks
    export TOTAL_CORES=$total_cores
    export TOTAL_MEMORY_GB=$total_memory_gb
}

#--------------------------------- Print nodes configuration
calculate_memory_per_core $SLURM_NNODES $SLURM_NTASKS_PER_NODE $SLURM_CPUS_PER_TASK $memory_per_node


for alpha in -0.7 -0.6 -0.5 -0.4 -0.3 ; do
prefix="tddft_eps_gaup_${alpha}"
echo "Input file                   : ${prefix}.in"
echo "Job directory                : ${prefix}_run"
echo "Output directory             : ${prefix}_out"
cat > ${prefix}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= ${alpha}              # [TDDFT] LRC alpha factor
NGsBlkXd= 3               Ry    # [Xd] Response block size
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
echo "Alpha                        : ${alpha} terminado"
done
echo "End job:$(date)"
