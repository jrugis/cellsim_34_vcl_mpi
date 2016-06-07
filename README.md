# cellsim_34_vcl

##cmake

###pan (or anywhere you want to use the Intel compilers)
1. load required modules
  1. module load intel/2015a
  2. module load Python/3.5.0-intel-2015a
  3. module load CUDA/7.5.18
  4. module load CMake/3.4.1-GCC-4.9.2
2. git clone https://github.com/jrugis/cellsim_34_vcl.git
3. build CPU versions (CUDA version must be built separately due to not supporting Intel compiler currently)
  1. mkdir -p cellsim_34_vcl/build
  2. cd cellsim_34_vcl/build
  3. CC=icc CXX=icpc cmake .. -DBUILD_OPENMP=ON -DBUILD_MKL=ON
  4. make -j6
4. build CUDA versions
  1. mkdir -p ../build_cuda
  2. cd ../build_cuda
  3. CC=gcc CXX=g++ cmake .. -DBUILD_SERIAL=OFF -DBUILD_CUDA=ON
  4. make
5. run a simulation
  1. cd ../pan/test
  2. edit run_test.sl (e.g. to point to one of the binaries you just built)
  3. sbatch run_test.sl
6. check the results (this is also performed at the end of run_test.sl)
  1. python cs_reduce.py
  2. python cs_results_r.py
  3. diff cR.bin generic3d_03-cR.bin

###linux (GCC compilers)
1. install Boost, Eigen and ViennaCL (optionally MKL and CUDA)
2. git clone https://github.com/jrugis/cellsim_34_vcl.git
3. build required version (set BUILD_* flags to OFF if not required)
  1. mkdir -p cellsim_34_vcl/build
  2. cd cellsim_34_vcl/build
  3. CC=gcc CXX=g++ cmake .. -DBUILD_SERIAL=ON -DBUILD_OPENMP=ON -DBUILD_MKL=ON -DBUILD_CUDA=ON -DTHREE_VARIABLES=ON -DFOUR_VARIABLES=ON
  4. you may require extra options if the dependencies are not in default locations
    1. add -DEIGEN3_INCLUDE_DIR=/path/to/eigen3 if Eigen cannot be located
    2. add -DBOOST_ROOT=/path/to/boost if Boost cannot be located
    3. add -DVIENNACL_INCLUDE_DIR=/path/to/viennacl if ViennaCL cannot be located
    4. make sure the MKLROOT environment variable is set properly to point to the MKL location
  5. make
4. run a simulation (e.g. generic3d_03)
  1. cd ../linux/test
  2. cp ../../test-ref/cell01m_HARMONIC_100p.msh cs.msh
  3. cp ../../test-ref/generic3d_03-cs.dat cs.dat
  4. ../../build/generic3d_03
5. check the results
  1. python cs_reduce.py
  2. python cs_results_r.py
  3. diff cR.bin ../../test-ref/generic3d_03-cR.bin

##linux
1. install eigen and ViennaCL
2. git clone https://github.com/jrugis/cellsim_34_vcl.git
3. build the simulation code
  1. cd cellsim_34_vcl/linux/build
  2. edit Makefile
  3. make
4. run a simulation (e.g. generic3d_03)
  1. cd ../test
  2. cp ../../test-ref/cell01m_HARMONIC_100p.msh cs.msh
  3. cp ../../test-ref/generic3d_03-cs.dat cs.dat
  4. ../generic3d_03
5. check the results
  1. python cs_reduce_min-max.py
  2. python cs_results_r.py
  3. diff cR.bin ../../test-ref/generic3d_03-cR.bin

##pan
1. load required modules
  1. module load intel/2015a
  2. module load Python/3.5.0-intel-2015a
  3. module load CUDA/7.5.18 (for the CUDA version only)
2. git clone https://github.com/jrugis/cellsim_34_vcl.git
3. build the simulation code
  1. cd cellsim_34_vcl/pan/build
  2. edit Makefile
  3. make
4. run a simulation
  1. cd ../test
  2. edit run_test.sl
  3. sbatch run_test.sl
5. check the results
  1. python cs_reduce.py
  2. python cs_results_r.py
  3. diff cR.bin generic3d_03-cR.bin
