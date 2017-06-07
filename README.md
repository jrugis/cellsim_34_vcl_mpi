# cellsim_34_vcl
[![DOI](https://zenodo.org/badge/23023/jrugis/cellsim_34_vcl.svg)](https://zenodo.org/badge/latestdoi/23023/jrugis/cellsim_34_vcl)
[![Build Status](https://travis-ci.org/jrugis/cellsim_34_vcl.svg?branch=master)](https://travis-ci.org/jrugis/cellsim_34_vcl)

## pan (cmake)
1. load required modules
  1. `module load CMake/3.6.1`
  2. `module load Python/3.5.1-intel-2015a`
2. `git clone https://github.com/jrugis/cellsim_34_vcl.git`
3. build the simulation code
  1. `mkdir cellsim_34_vcl/build && cd cellsim_34_vcl/build`
  2. `CXX=icpc cmake .. -DBUILD_OPENMP=ON -DBUILD_MKL=ON -DBUILD_CUDA=ON -DCMAKE_BUILD_TYPE=RELEASE`
  3. `make` (the executables will be located in the "src" directory underneath the build directory)
4. run a test simulation
  * `ctest --output-on-failure` (will run all test cases in serial and print a summary)
  * `ctest -N` (will print the names of all available test cases)
  * `ctest -V -R generic3d_03_pan_vcl-serial-intel$` (will run the specified test case with verbose output)
  * `ctest -j4` (will run all test cases in parallel, using 4 threads)

## linux
1. install boost, eigen and ViennaCL
2. git clone https://github.com/jrugis/cellsim_34_vcl.git
3. build the simulation code
  1. cd cellsim_34_vcl/linux/build
  2. edit Makefile
  3. make
4. run a simulation (e.g. generic3d_03)
  1. cd ../test
  2. cp ../../test/cell01m_HARMONIC_100p.msh cs.msh
  3. cp ../../test/generic3d_03-cs.dat cs.dat
  4. ../generic3d_03
5. check the results
  1. python ../../test/cs_reduce_min-max.py (only if not already reduced)
  2. python ../../test/cs_results_r.py
  3. python ../../test/cs_compare_peaks.py cR.bin ../../test/generic3d-cR.bin
