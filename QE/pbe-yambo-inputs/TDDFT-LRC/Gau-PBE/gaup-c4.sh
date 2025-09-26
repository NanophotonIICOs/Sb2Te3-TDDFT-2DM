#!/bin/bash
#SBATCH -J QE-GAUP # job name
#SBATCH -o qe.o-%j # output and error file name (%j expands to jobID)
#SBATCH --nodes=7 # total of nodes (N,nodes) or  total number of MPI tasks requested (n)
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
#---------------------------------Functions
# Variables
RJOB=$SLURM_JOB_NAME
NTASKS=$SLURM_NTASKS
NNODES=$SLURM_NNODES
NTASKPN=$SLURM_TASKS_PER_NODE
CALC="pw.x"
TMP_DIR="tmp"
SYS0="Sb2Te3-uatsx-c4-0.0"
SYS1="Sb2Te3-uatsx-c4-1.0"
SYS2="Sb2Te3-uatsx-c4-2.0"
SYS3="Sb2Te3-uatsx-c4-3.0"
SYS4="Sb2Te3-uatsx-c4-4.0"
SYS5="Sb2Te3-uatsx-c4-5.0"

#---------------------------------Print Parallel Params
print_parallel_params "$RJOB" "$NTASKS" "$NNODES" "$NTASKPN" "$CALC"
check_and_create_dir "./$SCF"
check_and_create_dir "./$NSCF"
check_and_create_dir "./$BANDS"
check_and_create_dir "./$TMP_DIR"
#==================================== Uniaxial-tensile strain at 0% ======================================================

#-----------------------------------SCF
scf_in_0="./$SCF/$SYS0.$SCF.$PW.in"
scf_out_0="./$SCF/$SYS0.$SCF.$PW.out"
print_calc_params  "$SYS0" "$SCF" "$scf_in_0" "$scf_out_0"

cat > $scf_in_0 << EOF
&CONTROL
calculation   = '$SCF'
prefix        = '$SYS0'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/

&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
4.310483557   0.000000000   0.000000000
-2.155241778   3.734391482   0.000000000
0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.4788874857        1.8207748514       13.0864344015
Te               2.4791053292        1.8195801340        7.3723391896
Sb               0.3229804050        0.5757787051        9.0763801732
Te               0.3239477704        0.5746866929       14.7905057904
Te               0.3240262797        3.0640956665       11.0813769653

K_POINTS (automatic)
13 13 1 0 0 0
EOF
check_input_file "$scf_in_0"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x < $scf_in_0 > $scf_out_0
echo "$SCF calculation is end "
echo " "

#-----------------------------------NSCF
nscf_in_0="./$NSCF/$SYS0.$NSCF.$PW.in"
nscf_out_0="./$NSCF/$SYS0.$NSCF.$PW.out"
print_calc_params  "$SYS0" "$NSCF" "$nscf_in_0" "$nscf_out_0"

cat > $nscf_in_0 << EOF
&CONTROL
calculation   = '$NSCF'
prefix        = '$SYS0'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/

&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
4.310483557   0.000000000   0.000000000
-2.155241778   3.734391482   0.000000000
0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.4788874857        1.8207748514       13.0864344015
Te               2.4791053292        1.8195801340        7.3723391896
Sb               0.3229804050        0.5757787051        9.0763801732
Te               0.3239477704        0.5746866929       14.7905057904
Te               0.3240262797        3.0640956665       11.0813769653

K_POINTS (automatic)
13 13 1 0 0 0
EOF
check_input_file "$nscf_in_0"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x < $nscf_in_0 > $nscf_out_0
echo "$NSCF calculation is end "
echo " "


#==================================== Uniaxial-tensile strain at 1% ======================================================

#-----------------------------------SCF
scf_in_1="./$SCF/$SYS1.$SCF.$PW.in"
scf_out_1="./$SCF/$SYS1.$SCF.$PW.out"
print_calc_params  "$SYS1" "$SCF" "$scf_in_1" "$scf_out_1"

cat > $scf_in_1 << EOF
&CONTROL
calculation   = '$SCF'
prefix        = '$SYS1'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/

&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
4.353588393   0.000000000   0.000000000
-2.155241778   3.732952155   0.000000000
0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.5200981363        1.8318200729       13.0785054255
Te               2.5044351829        1.8017863878        7.3893288453
Sb               0.3235158827        0.5652608372        9.0842898074
Te               0.3405583053        0.5927577591       14.7734419753
Te               0.3449800532        3.0632909930       11.0814704664

K_POINTS (automatic)
13 13 1 0 0 0
EOF
check_input_file "$scf_in_1"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $scf_in_1 > $scf_out_1
echo "$SCF calculation is end "
echo " "

# #-------------------------------NSCF

nscf_in_1="./$NSCF/$SYS1.$NSCF.$PW.in"
nscf_out_1="./$NSCF/$SYS1.$NSCF.$PW.out"
print_calc_params  "$SYS1" "$NSCF" "$nscf_in_1" "$nscf_out_1"

cat > $nscf_in_1 << EOF
&CONTROL
calculation   = '$NSCF'
prefix        = '$SYS1'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/
&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
4.353588393   0.000000000   0.000000000
-2.155241778   3.732952155   0.000000000
0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.5200981363        1.8318200729       13.0785054255
Te               2.5044351829        1.8017863878        7.3893288453
Sb               0.3235158827        0.5652608372        9.0842898074
Te               0.3405583053        0.5927577591       14.7734419753
Te               0.3449800532        3.0632909930       11.0814704664

K_POINTS (automatic)
13 13 1 0 0 0

EOF
check_input_file "$nscf_in_1"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $nscf_in_1 > $nscf_out_1
echo "$NSCF calculation is end "
echo " "


#==================================== Uniaxial-tensile strain at 2% ======================================================

#-----------------------------------SCF
scf_in_2="./$SCF/$SYS2.$SCF.$PW.in"
scf_out_2="./$SCF/$SYS2.$SCF.$PW.out"
print_calc_params  "$SYS2" "$SCF" "$scf_in_2" "$scf_out_2"

cat > $scf_in_2 << EOF
&CONTROL
calculation   = '$SCF'
prefix        = '$SYS2'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/

&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
   4.396693228   0.000000000   0.000000000
   -2.155241778   3.732952155   0.000000000
   0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.5631547487        1.8420419422       13.0690161419
Te               2.5290548815        1.7840313771        7.4086200306
Sb               0.3236276847        0.5528005097        9.0938387521
Te               0.3569472771        0.6115407227       14.7544304125
Te               0.3654432586        3.0645014983       11.0811311828

K_POINTS (automatic)
13 13 1 0 0 0
EOF
check_input_file "$scf_in_2"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $scf_in_2 > $scf_out_2
echo "$SCF calculation is end "
echo " "

# #-------------------------------NSCF

nscf_in_2="./$NSCF/$SYS2.$NSCF.$PW.in"
nscf_out_2="./$NSCF/$SYS2.$NSCF.$PW.out"
print_calc_params  "$SYS2" "$NSCF" "$nscf_in_2" "$nscf_out_2"

cat > $nscf_in_2 << EOF
&CONTROL
calculation   = '$NSCF'
prefix        = '$SYS2'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/
&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
   4.396693228   0.000000000   0.000000000
   -2.155241778   3.732952155   0.000000000
   0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.5631547487        1.8420419422       13.0690161419
Te               2.5290548815        1.7840313771        7.4086200306
Sb               0.3236276847        0.5528005097        9.0938387521
Te               0.3569472771        0.6115407227       14.7544304125
Te               0.3654432586        3.0645014983       11.0811311828

K_POINTS (automatic)
13 13 1 0 0 0

EOF
check_input_file "$nscf_in_2"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $nscf_in_2 > $nscf_out_2
echo "$NSCF calculation is end "
echo " "

#==================================== Uniaxial-tensile strain at 3% ======================================================

#-----------------------------------SCF
scf_in_3="./$SCF/$SYS3.$SCF.$PW.in"
scf_out_3="./$SCF/$SYS3.$SCF.$PW.out"
print_calc_params  "$SYS3" "$SCF" "$scf_in_3" "$scf_out_3"

cat > $scf_in_3 << EOF
&CONTROL
calculation   = '$SCF'
prefix        = '$SYS3'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/

&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
   4.439798064   0.000000000   0.000000000
   -2.155241778   3.732952155   0.000000000
   0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.6051902369        1.8527206076       13.0607506976
Te               2.5548625461        1.7669293459        7.4260643044
Sb               0.3231821485        0.5417762286        9.1021363550
Te               0.3732530860        0.6290395471       14.7367093177
Te               0.3863801236        3.0644503208       11.0813758453

K_POINTS (automatic)
13 13 1 0 0 0
EOF
check_input_file "$scf_in_3"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $scf_in_3 > $scf_out_3
echo "$SCF calculation is end "
echo " "

# #-------------------------------NSCF

nscf_in_3="./$NSCF/$SYS3.$NSCF.$PW.in"
nscf_out_3="./$NSCF/$SYS3.$NSCF.$PW.out"
print_calc_params  "$SYS3" "$NSCF" "$nscf_in_3" "$nscf_out_3"

cat > $nscf_in_3 << EOF
&CONTROL
calculation   = '$NSCF'
prefix        = '$SYS3'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/
&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
   4.439798064   0.000000000   0.000000000
   -2.155241778   3.732952155   0.000000000
   0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.6051902369        1.8527206076       13.0607506976
Te               2.5548625461        1.7669293459        7.4260643044
Sb               0.3231821485        0.5417762286        9.1021363550
Te               0.3732530860        0.6290395471       14.7367093177
Te               0.3863801236        3.0644503208       11.0813758453

K_POINTS (automatic)
13 13 1 0 0 0

EOF
check_input_file "$nscf_in_3"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $nscf_in_3 > $nscf_out_3
echo "$NSCF calculation is end "
echo " "


#==================================== Uniaxial-tensile strain at 4% ======================================================

#-----------------------------------SCF
scf_in_4="./$SCF/$SYS4.$SCF.$PW.in"
scf_out_4="./$SCF/$SYS4.$SCF.$PW.out"
print_calc_params  "$SYS4" "$SCF" "$scf_in_4" "$scf_out_4"

cat > $scf_in_4 << EOF
&CONTROL
calculation   = '$SCF'
prefix        = '$SYS4'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/

&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
   4.482902899   0.000000000   0.000000000
   -2.155241778   3.732952155   0.000000000
   0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.6469954338        1.8637789214       13.0522023544
Te               2.5800441047        1.7503250968        7.4440806968
Sb               0.3232546401        0.5311881700        9.1105763920
Te               0.3898357779        0.6453247907       14.7187870447
Te               0.4073784748        3.0642990711       11.0813900320

K_POINTS (automatic)
13 13 1 0 0 0
EOF
check_input_file "$scf_in_4"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $scf_in_4 > $scf_out_4
echo "$SCF calculation is end "
echo " "

# #-------------------------------NSCF

nscf_in_4="./$NSCF/$SYS4.$NSCF.$PW.in"
nscf_out_4="./$NSCF/$SYS4.$NSCF.$PW.out"
print_calc_params  "$SYS4" "$NSCF" "$nscf_in_4" "$nscf_out_4"

cat > $nscf_in_4 << EOF
&CONTROL
calculation   = '$NSCF'
prefix        = '$SYS4'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/
&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
   4.482902899   0.000000000   0.000000000
   -2.155241778   3.732952155   0.000000000
   0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.6469954338        1.8637789214       13.0522023544
Te               2.5800441047        1.7503250968        7.4440806968
Sb               0.3232546401        0.5311881700        9.1105763920
Te               0.3898357779        0.6453247907       14.7187870447
Te               0.4073784748        3.0642990711       11.0813900320

K_POINTS (automatic)
13 13 1 0 0 0

EOF
check_input_file "$nscf_in_4"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $nscf_in_4 > $nscf_out_4
echo "$NSCF calculation is end "
echo " "

#==================================== Uniaxial-tensile strain at 5% ======================================================

#-----------------------------------SCF
scf_in_5="./$SCF/$SYS5.$SCF.$PW.in"
scf_out_5="./$SCF/$SYS5.$SCF.$PW.out"
print_calc_params  "$SYS5" "$SCF" "$scf_in_5" "$scf_out_5"

cat > $scf_in_5 << EOF
&CONTROL
calculation   = '$SCF'
prefix        = '$SYS5'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/

&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
   4.526007735   0.000000000   0.000000000
   -2.155241778   3.732952155   0.000000000
   0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.6891559420        1.8740387011       13.0438468704
Te               2.6059906632        1.7342622830        7.4621461653
Sb               0.3229747280        0.5209365354        9.1189378471
Te               0.4057659489        0.6614010248       14.7006895716
Te               0.4282614395        3.0642775056       11.0814160656

K_POINTS (automatic)
13 13 1 0 0 0
EOF
check_input_file "$scf_in_5"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $scf_in_5 > $scf_out_5
echo "$SCF calculation is end "
echo " "

# #-------------------------------NSCF

nscf_in_5="./$NSCF/$SYS5.$NSCF.$PW.in"
nscf_out_5="./$NSCF/$SYS5.$NSCF.$PW.out"
print_calc_params  "$SYS5" "$NSCF" "$nscf_in_5" "$nscf_out_5"

cat > $nscf_in_5 << EOF
&CONTROL
calculation   = '$NSCF'
prefix        = '$SYS5'
pseudo_dir    = '$PSEUDO_DIR'
verbosity     = 'high' 
outdir        = '$TMP_DIR'
/
&SYSTEM 
ibrav            = 0                  
nat              = 5          
ntyp             = 2          
ecutwfc          = 70
input_dft        = 'gaup'
nbnd             = 100
nqx1             = 2
nqx2             = 2
nqx3             = 1 
vdw_corr         = 'DFT-D'
exxdiv_treatment='none'
x_gamma_extrapolation = .false.
force_symmorphic =.true.
/

&ELECTRONS
conv_thr         = 1.0D-7
diago_full_acc = .TRUE.
/

ATOMIC_SPECIES
Sb  121.760  Sb_ONCV_PBE-1.2.upf
Te  127.600  Te_ONCV_PBE-1.2.upf 

CELL_PARAMETERS (angstrom)
   4.526007735   0.000000000   0.000000000
   -2.155241778   3.732952155   0.000000000
   0.000000000   0.000000000   22.162984280

ATOMIC_POSITIONS (angstrom)
Sb               2.6891559420        1.8740387011       13.0438468704
Te               2.6059906632        1.7342622830        7.4621461653
Sb               0.3229747280        0.5209365354        9.1189378471
Te               0.4057659489        0.6614010248       14.7006895716
Te               0.4282614395        3.0642775056       11.0814160656

K_POINTS (automatic)
13 13 1 0 0 0
EOF
check_input_file "$nscf_in_5"
echo " "
mpirun -np $SLURM_NTASKS $QEGCCROOT/pw.x  < $nscf_in_5 > $nscf_out_5
echo "$NSCF calculation is end "
echo " "

echo "End job:$(date)"