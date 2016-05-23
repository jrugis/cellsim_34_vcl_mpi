#!/bin/bash
#SBATCH -J cellsim_34_vcl_cuda
#SBATCH -A nesi00119          # Project Account
#SBATCH --time=0:15:00       # Walltime HH:MM:SS
#SBATCH --mem-per-cpu=32G     # Memory
#SBATCH --ntasks=1            # number of tasks
#SBATCH --cpus-per-task=1     # number of threads
#SBATCH --gres=gpu:1
#SBATCH -C sb&kepler

# output some information
echo $HOSTNAME

# load module(s)
module load intel/2015a
module load Python/3.5.0-intel-2015a
module load CUDA/7.5.18

# number of MKL/OMP threads
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

# copy mesh and parameters
cp -f ../../test-ref/cell01m_HARMONIC_100p.msh cs.msh
cp -f ../../test-ref/generic3d_03-cs.dat cs.dat
#cp -f ../../test-ref/generic3d_04-cs.dat cs.dat

# run the job
srun ../generic3d_03_cuda
#srun ../generic3d_04_cuda

# compare the result
python cs_reduce.py
diff cR.bin generic3d_03_cuda-cR.bin
#diff cR.bin generic3d_04_cuda-cR.bin
