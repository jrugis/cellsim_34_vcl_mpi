
# NOTE: uncomment lines for the 3 or 4 variable model respectively
#
EXE = generic3d_03
#EXE = generic3d_03_debug
CPP_DEFS =
#
#EXE = generic3d_04
#EXE = generic3d_04_debug
#CPP_DEFS = -D FOUR_VARIABLES

#NOTE: uncomment for debug or optimised respectively  
#
#CPP_FLAGS =  $(CPP_DEFS) $(CPP_INCLUDES) -m64 -Og -Wall -Wextra
CPP_FLAGS =  $(CPP_DEFS) $(CPP_INCLUDES) -m64 -O3 

CPP = g++
CPP_INCLUDES = -I/usr/include/eigen3

LINK = $(CPP)
LINK_FLAGS = -m64

OBJS = cell_3d.o cGeneric3dModel.o cCellMesh.o cVCLSolver.o
SRC = ../../src

all: $(EXE)

cell_3d.o: $(SRC)/cell_3d.cpp
	$(CPP) $(CPP_FLAGS) -o $@ -c $<

cGeneric3dModel.o: $(SRC)/cGeneric3dModel.cpp $(SRC)/cGeneric3dModel.h
	$(CPP) $(CPP_FLAGS) -o $@ -c $<

cCellMesh.o: $(SRC)/cCellMesh.cpp $(SRC)/cCellMesh.h
	$(CPP) $(CPP_FLAGS) -o $@ -c $<

cVCLSolver.o: $(SRC)/cVCLSolver.cpp $(SRC)/cVCLSolver.h
	$(CPP) $(CPP_FLAGS) -o $@ -c $<

$(EXE): $(OBJS)
	$(LINK) $(LINK_FLAGS) -o ../$@ $+

clean:
	rm *.o


