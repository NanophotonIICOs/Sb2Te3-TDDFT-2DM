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
export LD_LIBRARY_PATH=/zfs-home/202401017n/software/libxc-6.2.2/lib:$LD_LIBRARY_PATH
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

calc0="o1-0.0"
calc1="o1-1.0"
calc2="o1-2.0"
calc3="o1-3.0"
calc4="o1-4.0"
calc5="o1-5.0"

check_and_create_dir "./$calc0"
check_and_create_dir "./$calc1"
check_and_create_dir "./$calc2"
check_and_create_dir "./$calc3"
check_and_create_dir "./$calc4"
check_and_create_dir "./$calc5"

# cd ./$calc0
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/p2y -I ../../../../../QE-work-2/b3lyp/calc-6/tmp/Sb2Te3-uatsx-c1-0.0.save
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo 

# cat > ox-${calc0}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry        # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1              # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF

# prefix="ox-${calc0}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc0} x dir"

# cat > oy-${calc0}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry    # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF
# prefix="oy-${calc0}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc0} y dir"

#=========================================================== uatsx-1.0=====================================================
# cd ..
# cd ./$calc1
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/p2y -I ../../../../../QE-work-2/b3lyp/calc-6/tmp/Sb2Te3-uatsx-c1-1.0.save
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo 

# cat > ox-${calc1}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry    # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF
# prefix="ox-${calc1}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc1} x dir"

# # #---------------------------------------------ydir
# cat > oy-${calc1}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry    # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF
# prefix="oy-${calc1}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc1} y dir"

#=========================================================== uatsx-2.0=====================================================
# cd ..
# cd ./$calc2
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/p2y -I ../../../../../QE-work-2/b3lyp/calc-6/tmp/Sb2Te3-uatsx-c1-2.0.save
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo 

# cat > ox-${calc2}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry    # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF
# prefix="ox-${calc2}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc2} x dir"

# # #---------------------------------------------ydir
# cat > oy-${calc2}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry    # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF
# prefix="oy-${calc2}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc2} y dir"



#=========================================================== uatsx-3.0=====================================================
#cd ..
cd ./$calc3
mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/p2y -I ../../../../../QE-work-2/b3lyp/calc-6/tmp/Sb2Te3-uatsx-c1-3.0.save
mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo 

cat > ox-${calc3}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
FFTGvecs=  50          Ry    # [FFT] Plane-waves
X_Threads=12                     # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
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
mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
echo "Done ${calc3} x dir"

#---------------------------------------------ydir
cat > oy-${calc3}.in << EOF
optics                           # [R] Linear Response optical properties
chi                              # [R][CHI] Dyson equation for Chi.
tddft                            # [R][K] Use TDDFT kernel
FFTGvecs=  50          Ry    # [FFT] Plane-waves
X_Threads=12                     # [OPENMP/X] Number of threads for response functions
Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
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
mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
echo "Done ${calc3} y dir"


#=========================================================== uatsx-4.0=====================================================
# cd ..
# cd ./$calc4
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/p2y -I ../../../../../QE-work-2/b3lyp/calc-6/tmp/Sb2Te3-uatsx-c1-4.0.save
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo 

# cat > ox-${calc4}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry    # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF
# prefix="ox-${calc4}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc4} x dir"

# # #---------------------------------------------ydir
# cat > oy-${calc4}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry    # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF
# prefix="oy-${calc4}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc4} y dir"



#=========================================================== uatsx-5.0=====================================================
# cd ..
# cd ./$calc5
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/p2y -I ../../../../../QE-work-2/b3lyp/calc-6/tmp/Sb2Te3-uatsx-c1-5.0.save
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo 

# cat > ox-${calc5}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry    # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  1.000000 | 0.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF
# prefix="ox-${calc5}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc5} x dir"

# # #---------------------------------------------ydir
# cat > oy-${calc5}.in << EOF
# optics                           # [R] Linear Response optical properties
# chi                              # [R][CHI] Dyson equation for Chi.
# tddft                            # [R][K] Use TDDFT kernel
# FFTGvecs=  50          Ry    # [FFT] Plane-waves
# X_Threads=12                     # [OPENMP/X] Number of threads for response functions
# Chimod= "LRC"                    # [X] IP/Hartree/ALDA/LRC/PF/BSfxc
# LRC_alpha= -0.1             # [TDDFT] LRC alpha factor
# NGsBlkXd= 5              Ry    # [Xd] Response block size
# % QpntsRXd
#    1 | 1 |                         # [Xd] Transferred momenta
# %
# % BndsRnXd
#    1 | 60 |                         # [Xd] Polarization function bands
# %
# % EnRngeXd
#   0.00000 | 5.00000 |         eV    # [Xd] Energy range
# %
# % DmRngeXd
#   0.050000 | 0.050000 |          eV    # [Xd] Damping range
# %
# ETStpsXd= 1501                    # [Xd] Total Energy steps
# % LongDrXd
#  0.000000 | 1.000000 | 0.000000 |        # [Xd] [cc] Electric Field
# %
# EOF
# prefix="oy-${calc5}"
# mpirun -np ${SLURM_NTASKS} $YAMBOGCCROOT/yambo -F ${prefix}.in -J "${prefix}-run" -C "${prefix}-out"
# echo "Done ${calc5} y dir"

echo "End job:$(date)"
