/*
 ============================================================================
 Name        : cell_3d.cpp
 Author      : J.Rugis
 Version     :
 Copyright   : (c) 2017 J.Rugis
 Description : Cell Simulation with ViennaCL and MPI
 ============================================================================
 */

#define MESH_FILE "cs.msh"
#include <iostream>
#include <sys/time.h>
#include <unistd.h>
#include <mpi.h>

#include "cCellMesh.h"
#include "cGeneric3dModel.h"

#define MPI_NODES 7 // one per cell, DO NOT CHANGE
#define MPI_ERROR_ABORT -101
#define MPI_NODES_ABORT -102
#define TEMP_SIZE 40

// shutdown mpi on error
void mpi_abort(int err)
{
  std::cout << "MPI shutdown: " << err << std::endl;
  MPI_Abort(MPI_COMM_WORLD, err);
}

// mpi error checking macro
#define MPI_CHECK(call) \
  if((call) != MPI_SUCCESS) { \
    std::cout << "MPI error calling: " << call << std::endl; \
    mpi_abort(MPI_ERROR_ABORT); }

// the main program function for each mpi node (cell)
int main(int argc,char **args){
  int commSize, commRank;
  cCellMesh *mesh_00;
  cGeneric3dModel *model;
  struct timeval start, end;
  int duration;
  std::string host_name;

  gettimeofday(&start, NULL); // get the time

  // get the hostname
  char temp[TEMP_SIZE];
  gethostname(temp, TEMP_SIZE);
  host_name = temp;

  // initialize mpi
  MPI_CHECK(MPI_Init(&argc, &args));
  MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &commSize));
  MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &commRank));
  if(commSize != MPI_NODES) mpi_abort(MPI_NODES_ABORT); // check mpiexec node count
  std::cout << "<MAIN " << commRank << "/" << commSize << "> hostname: " << host_name << std::endl;

  // read in the mesh file and display some mesh information
  std::cout << "<MAIN " << commRank << "/" << commSize << "> creating mesh object..." << std::endl;
  //mesh_00 = new cCellMesh(MESH_FILE, commRank);
  //mesh_00->print_info();

  // setup and run a model
  //  - reads in a model parameter file
  std::cout << "<MAIN " << commRank << "/" << commSize << "> creating model object..." << std::endl;
  //model = new cGeneric3dModel(mesh_00, commRank);
  //model->run();

  // save the results
  //model->save_results();
  gettimeofday(&end, NULL);
  duration = end.tv_sec - start.tv_sec;
  std::cout << "<MAIN " << commRank << "/" << commSize << "> execution time: " << duration << " sec" << std::endl;

  //delete model;
  //delete mesh_00;

  // shutdown mpi
  MPI_CHECK(MPI_Finalize());
  return 0;
}
