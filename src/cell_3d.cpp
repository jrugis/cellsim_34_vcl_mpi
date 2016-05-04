/*
 ============================================================================
 Name        : cell_3d.cpp
 Author      : J.Rugis
 Version     :
 Copyright   : (c) 2015 J.Rugis
 Description : Cell Simulation with PETSc and MPI
 ============================================================================
 */

#define MESH_FILE "cs.msh"
#include <iostream>
#include <sys/time.h>
#include <unistd.h>

#include "cCellMesh.h"
#include "cGeneric3dModel.h"

#define TEMP_SIZE 40

int main(int argc,char **args){
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
	std::cout << "<MAIN> hostname:" << host_name<< std::endl;

	// read in the mesh file and display some mesh information
	std::cout << "<MAIN> creating mesh object..." << std::endl;
	mesh_00 = new cCellMesh(MESH_FILE);
	mesh_00->print_info();

	// setup and run a model
	//  - reads in a model parameter file
	std::cout << "<MAIN> creating model object..." << std::endl;
	model = new cGeneric3dModel(mesh_00);
	model->run();

	// save the results
	model->save_results();
	gettimeofday(&end, NULL);
	duration = end.tv_sec - start.tv_sec;
	std::cout << "<MAIN> execution time: " << duration << " sec" << std::endl;

	delete model;
	delete mesh_00;

	return 0;
}
