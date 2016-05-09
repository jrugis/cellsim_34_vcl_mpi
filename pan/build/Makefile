
# Required module on pan:
#     module load intel/2015a
#     module load Python/3.5.0-intel-2015a

# NOTE: uncomment lines for the 3 or 4 variable model respectively
#
EXE = generic3d_03
#EXE = generic3d_03_debug
CPP_DEFS =
#CPP_DEFS = -DVIENNACL_WITH_OPENMP -fopenmp
#
#EXE = generic3d_04
#EXE = generic3d_04_debug
#CPP_DEFS = -DFOUR_VARIABLES

# library include paths
EIGEN_INC = -I/share/easybuild/RHEL6.3/sandybridge/software/Eigen/3.2.4-ictce-5.4.0/include
BOOST_INC = -I/share/easybuild/RHEL6.3/sandybridge/software/Boost/1.55.0-goolf-1.5.14-Python-2.7.8/include
VCL_INC = -I/projects/nesi00119/code/CS_dev/viennacl/ViennaCL-1.7.1
# uncomment to disable debugging for release version
DEFS = -DNDEBUG

CPP = icpc
CPP_FLAGS = -std=c++11 -O3 -xSSE4.2 -axAVX -fp-model precise $(CPP_DEFS) $(DEFS) $(EIGEN_INC) $(BOOST_INC) $(VCL_INC)
LINK_FLAGS = 

# nothing should need to change below here

SRC = ../../src
OBJS = cell_3d.o cGeneric3dModel.o cCellMesh.o cVCLSolver.o
EXEDIR = ..

all: $(EXEDIR)/$(EXE)

cell_3d.o: $(SRC)/cell_3d.cpp
	$(CPP) $(CPP_FLAGS) -o $@ -c $<

cGeneric3dModel.o: $(SRC)/cGeneric3dModel.cpp $(SRC)/cGeneric3dModel.h
	$(CPP) $(CPP_FLAGS) -o $@ -c $<

cCellMesh.o: $(SRC)/cCellMesh.cpp $(SRC)/cCellMesh.h
	$(CPP) $(CPP_FLAGS) -o $@ -c $<

cVCLSolver.o: $(SRC)/cVCLSolver.cpp $(SRC)/cVCLSolver.h
	$(CPP) $(CPP_FLAGS) -o $@ -c $<

$(EXEDIR)/$(EXE): $(OBJS)
	$(CPP) -o $@ $+ $(LINK_FLAGS)

.PHONY: clean distclean

clean:
	rm -f *.o

distclean: clean
	rm -f $(EXEDIR)/generic3d_*
