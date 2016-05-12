#!/bin/bash
#SBATCH -J cellsim_34_vcl
#SBATCH -A nesi00119          # Project Account
#SBATCH --time=0:15:00       # Walltime HH:MM:SS
#SBATCH --mem-per-cpu=16G     # Memory
#SBATCH --ntasks=1            # number of tasks
#SBATCH --cpus-per-task=1     # number of threads
#SBATCH -C avx                 # sb=Sandybridge wm=Westmere avx=Sandybridge or Ivybridge

# output some information
echo $HOSTNAME

# load module(s)
module load intel/2015a
module load Python/3.5.0-intel-2015a

# number of MKL/OMP threads
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

# copy mesh and parameters
cp -f ../../test-ref/cell01m_HARMONIC_100p.msh cs.msh
cp -f ../../test-ref/generic3d_03-cs.dat cs.dat
#cp -f ../../test-ref/generic3d_04-cs.dat

# run the job
srun ../generic3d_03
#srun ../generic3d_04

# compare the result
python cs_reduce.py
diff cR.bin generic3d_03-cR.bin
#diff cR.bin generic3d_04-cR.bin
