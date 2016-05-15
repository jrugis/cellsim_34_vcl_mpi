# cellsim_34_vcl

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
  1. python cs_reduce.py
  2. python cs_results_r.py
  3. diff cR.bin ../../test-ref/generic3d_03-cR.bin

##pan
1. load required modules
  1. module load intel/2015a
  2. module load Python/3.5.0-intel-2015a
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
