
# Required module on pan:
#     module load intel/2015a
#     module load Python/3.5.0-intel-2015a
#     module load CUDA/7.5.18

# NOTE: uncomment lines for the 3 or 4 variable model (CPU or CUDA version) respectively
#
EXE = generic3d_03
#EXE = generic3d_03_cuda
#EXE = generic3d_03_debug
CPP_DEFS =
#
#EXE = generic3d_04
#EXE = generic3d_04_cuda
#EXE = generic3d_04_debug
#CPP_DEFS = -DFOUR_VARIABLES

# NOTE: uncomment lines for the CPU or CUDA version respectively
#
SUFFIX = cpp
CXX = icpc
CXXFLAGS = -O3 -xSSE4.2 -axAVX -fp-model precise
LDFLAGS = 
#
#SUFFIX = cu
#CXX = nvcc
#CXXFLAGS = -arch sm_35 -O3 -Xcompiler -march=native,-mtune=native -DVIENNACL_WITH_CUDA
#LDFLAGS = 

# NOTE: uncomment lines for the OpenMP version (must be used with CPU version above only)
#CXXFLAGS += -openmp -DVIENNACL_WITH_OPENMP
#LDFLAGS += -openmp

# library include paths
EIGEN_INC = -I/share/easybuild/RHEL6.3/sandybridge/software/Eigen/3.2.4-ictce-5.4.0/include
BOOST_INC = -I/share/easybuild/RHEL6.3/sandybridge/software/Boost/1.55.0-goolf-1.5.14-Python-2.7.8/include
VCL_INC = -I/projects/nesi00119/code/CS_dev/viennacl/ViennaCL-1.7.1
# uncomment to disable debugging for release version
DEFS = -DNDEBUG

# nothing should need to change below here

# add defs/includes
CXXFLAGS += $(CPP_DEFS) $(DEFS) $(EIGEN_INC) $(BOOST_INC) $(VCL_INC)

SRC = ../../src
OBJS = cell_3d.o cGeneric3dModel.o cCellMesh.o cVCLSolver.o
EXEDIR = ..

all: $(EXEDIR)/$(EXE)

cell_3d.o: $(SRC)/cell_3d.$(SUFFIX)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

cGeneric3dModel.o: $(SRC)/cGeneric3dModel.$(SUFFIX) $(SRC)/cGeneric3dModel.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

cCellMesh.o: $(SRC)/cCellMesh.$(SUFFIX) $(SRC)/cCellMesh.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

cVCLSolver.o: $(SRC)/cVCLSolver.$(SUFFIX) $(SRC)/cVCLSolver.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(EXEDIR)/$(EXE): $(OBJS)
	$(CXX) -o $@ $+ $(LDFLAGS)

.PHONY: clean distclean

clean:
	rm -f *.o

distclean: clean
	rm -f $(EXEDIR)/generic3d_*
