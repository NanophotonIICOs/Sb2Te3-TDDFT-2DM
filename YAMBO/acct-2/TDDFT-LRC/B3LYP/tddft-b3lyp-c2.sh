#!/bin/bash
#SBATCH -J tddft-b3lyp-c3             # Job name             
#SBATCH -o yambo-tddft.o-%j           
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12 # number of tasks per node (maximum 24)
#SBATCH --cpus-per-task=2       #  threads  
#SBATCH -p normal # SLURM queue (partition)
#SBATCH --exclude=cn0605-3
#SBATCH -t 120:00:00 # run time (hh:jobmm:ss)

module purge
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
    local first_node=$(scontrol show hostnames $SLURM_NODELIST | head -1)
    local mem_kb=$(srun -N1 -n1 --nodelist=$first_node grep MemTotal /proc/meminfo 2>/dev/null | awk '{print $2}')
    local memory_per_node_gb=$((mem_kb / 1024 / 1024))
    # local memory_per_node_gb=${SLURM_MEM_PER_NODE} 
    
    # Calculate total tasks, cores and memory
    local total_tasks=$((nodes * ntasks_per_node))
    local total_cores=$((total_tasks * cpus_per_task))
    local total_memory_gb=$((nodes * memory_per_node_gb))
    
    # Memory calculations
    local memory_per_core_gb=$((total_memory_gb / total_cores))
    local memory_per_core_mb=$((memory_per_core_gb * 1024))
    local memory_per_task_gb=$((total_memory_gb / total_tasks))
    
    echo "=== RESOURCE ALLOCATION SUMMARY ==="
    echo "Cores-Threads        : ${total_cores}(CPU)-${cpus_per_task}(threads)"
    echo "MPI Cores            : ${total_tasks}"
    echo "Threads per core     : ${cpus_per_task}"
    echo "Threads total        : ${total_cores}"
    echo "Nodes Computing      : ${nodes}"
    echo "Memory per task      : ${memory_per_task_gb} GB"
    echo "Total memory         : ${total_memory_gb} GB"
    echo "====================================="
    echo " "
    
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

calc0="eps-b3lyp-0.0"
calc1="eps-b3lyp-1.0"
calc2="eps-b3lyp-2.0"
calc3="eps-b3lyp-3.0"
calc4="eps-b3lyp-4.0"
calc5="eps-b3lyp-5.0"

check_and_create_dir "./$calc0"
check_and_create_dir "./$calc1"
check_and_create_dir "./$calc2"
check_and_create_dir "./$calc3"
check_and_create_dir "./$calc4"
check_and_create_dir "./$calc5"

cd ./$calc0
mpirun -np ${SLURM_NTASKS} p2y -I ../../../../../../QE/b3lyp/calc-2/Sb2Te3-uatsx-c1-0.0.save
mpirun -np ${SLURM_NTASKS} yambo 

cat > ox-${calc0}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4              # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF

prefix="ox-${calc0}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc0} x-dir"

cat > oy-${calc0}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="oy-${calc0}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc0} y-dir"

#=========================================================== uatsx-1.0=====================================================
cd ..
cd ./$calc1
mpirun -np ${SLURM_NTASKS} p2y -I ../../../../../../QE/b3lyp/calc-2/Sb2Te3-uatsx-c1-1.0.save
mpirun -np ${SLURM_NTASKS} yambo 

cat > ox-${calc1}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="ox-${calc1}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc1} x-dir"

# # #---------------------------------------------ydir
cat > oy-${calc1}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="oy-${calc1}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc1} y-dir"

#=========================================================== uatsx-2.0=====================================================
cd ..
cd ./$calc2
mpirun -np ${SLURM_NTASKS} p2y -I ../../../../../../QE/b3lyp/calc-2/Sb2Te3-uatsx-c1-2.0.save
mpirun -np ${SLURM_NTASKS} yambo 

cat > ox-${calc2}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="ox-${calc2}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc2} x-dir"

# #---------------------------------------------ydir
cat > oy-${calc2}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="oy-${calc2}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc2} y-dir"



#=========================================================== uatsx-3.0=====================================================
cd ..
cd ./$calc3
mpirun -np ${SLURM_NTASKS} p2y -I ../../../../../../QE/b3lyp/calc-2/Sb2Te3-uatsx-c1-3.0.save
mpirun -np ${SLURM_NTASKS} yambo 

cat > ox-${calc3}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="ox-${calc3}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc3} x-dir"

#---------------------------------------------ydir
cat > oy-${calc3}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="oy-${calc3}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc3} y-dir"


#=========================================================== uatsx-4.0=====================================================
cd ..
cd ./$calc4
mpirun -np ${SLURM_NTASKS} p2y -I ../../../../../../QE/b3lyp/calc-2/Sb2Te3-uatsx-c1-4.0.save
mpirun -np ${SLURM_NTASKS} yambo 

cat > ox-${calc4}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="ox-${calc4}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc4} x-dir"

# #---------------------------------------------ydir
cat > oy-${calc4}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}     # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="oy-${calc4}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc4} y-dir"



#=========================================================== uatsx-5.0=====================================================
cd ..
cd ./$calc5
mpirun -np ${SLURM_NTASKS} p2y -I ../../../../../../QE/b3lyp/calc-2/Sb2Te3-uatsx-c1-5.0.save
mpirun -np ${SLURM_NTASKS} yambo 

cat > ox-${calc5}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="ox-${calc5}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc5} x-dir"

# # #---------------------------------------------ydir
cat > oy-${calc5}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
X_Threads=${OMP_NUM_THREADS}                   # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.4             # [TDDFT] LRC alpha factor
NGsBlkXd= 5              Ry    # [Xd] Response block size
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
  0.050000 | 0.050000 |          eV    # [Xd] Damping range
%
ETStpsXd= 1501                    # [Xd] Total Energy steps
% LongDrXd
 0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
%
EOF
prefix="oy-${calc5}"
mpirun -np ${SLURM_NTASKS} yambo -F ${prefix}.in -J "${prefix}_run" -C "${prefix}_out"
echo "Done ${calc5} y-dir"

echo "End job:$(date)"