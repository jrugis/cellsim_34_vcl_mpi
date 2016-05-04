/*
 * cCellMesh.cpp
 *
 *  Created on: 14/04/2015
 *      Author: jrugis
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cmath>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "cCellMesh.h"

cCellMesh::cCellMesh(std::string file_name){
	// initialise member variables
	nodes_count = 0;
	total_elements_count = 0;
	surface_elements_count = volume_elements_count = 0;

	get_mesh(file_name); /* rank 0 only??? */
	calc_dist();
}

cCellMesh::~cCellMesh(){
}

void cCellMesh::calc_dist(){
	std::cout << "<MESH> calculating node distance to surface..." << std::endl;
	Eigen::Matrix<tCoord,1,3> v1, v2;
	for(tElement n = 0; n < nodes_count; n++){
		if(surface_node(n)){
			node_data(n, dist_surface) = 0.0;
			continue;
		}
		v1 = coordinates.block<1,3>(n, 0);
		node_data(n, dist_surface) = std::numeric_limits<tCalcs>::max();
		for(tElement s = 1; s < nodes_count; s++){
			if(!surface_node(s)) continue;
			v2 = coordinates.block<1,3>(s, 0);
			tCalcs d = (v1 - v2).squaredNorm();
			if(d < node_data(n, dist_surface)) node_data(n, dist_surface) = d;
		}
		node_data(n, dist_surface) = std::sqrt(tCalcs(node_data(n, dist_surface)));
	}
}

void cCellMesh::get_mesh(std::string file_name){
    // local variables
	std::ifstream cell_file(file_name.c_str()); // open the mesh file
	std::string line;                           // file line buffer
    std::vector <std::string> tokens;           // tokenized line

    // get the mesh nodes
	std::cout << "<MESH> getting the mesh nodes..." << std::endl;
	while(getline(cell_file, line)){
		if(line != "$Nodes") continue;
		getline(cell_file, line);
		nodes_count = std::atof(line.c_str());
		coordinates.resize(nodes_count, Eigen::NoChange);
		for(tElement n = 0; n < nodes_count; n++){
			getline(cell_file, line);
			boost::split(tokens, line, boost::is_any_of(", "), boost::token_compress_on);
			for(int m = 0; m < 3; m++) coordinates(n, m) = atof(tokens[m + 1].c_str());
		}
		break;
	}
	surface_node.resize(nodes_count, Eigen::NoChange);
	surface_node.Constant(nodes_count, false);

	// get the mesh elements
	std::cout << "<MESH> getting the mesh elements..." << std::endl;
	while(getline(cell_file, line)){
		if(line != "$Elements") continue;
		getline(cell_file, line);
		total_elements_count = std::atof(line.c_str());
		surface_elements.resize(total_elements_count, Eigen::NoChange);   // overkill for now
		volume_elements.resize(total_elements_count, Eigen::NoChange);    //
		for(tElement n = 0; n < total_elements_count; n++){
			getline(cell_file, line);
			boost::split(tokens, line, boost::is_any_of(", "), boost::token_compress_on);
			int element_type = atoi(tokens[1].c_str());
			if(element_type == 2){  // surface triangles
				for(int m = 0; m < 3; m++){                                         // index from zero
					tElement i = atol(tokens[m + 5].c_str()) - 1;
					surface_elements(surface_elements_count, m) = i;
					surface_node(i) = true;
				}
				surface_elements_count++;
			}
			if(element_type == 4){  // volume tetrahedrons
				for(int m = 0; m < 4; m++)                                          // index from zero
					volume_elements(volume_elements_count, m) = atol(tokens[m + 5].c_str()) - 1;
				volume_elements_count++;
			}
		}
		surface_elements.conservativeResize(surface_elements_count, Eigen::NoChange);  // correct the size
		volume_elements.conservativeResize(volume_elements_count, Eigen::NoChange);    //
		break;
	}
	// get the node data
	std::cout << "<MESH> getting the mesh node data..." << std::endl;
	while(getline(cell_file, line)){
		if(line != "\"distance to nearest lumen\"") continue;
		else for(int i = 0; i < 6; i++) getline(cell_file, line); // skip six lines
		node_data.resize(nodes_count, Eigen::NoChange);
		for(tElement n = 0; n < nodes_count; n++){
			getline(cell_file, line);
			boost::split(tokens, line, boost::is_any_of(", "), boost::token_compress_on);
			node_data(n, dist_lumen) = atof(tokens[1].c_str());
		}
		break;
	}
	cell_file.close();
}

void cCellMesh::print_info(){
	std::cout << "<MESH> nodes_count: " << nodes_count << std::endl;
	std::cout << "<MESH> total_elements_count: " << total_elements_count << std::endl;
	std::cout << "<MESH> surface_elements_count: " << surface_elements_count << std::endl;
	std::cout << "<MESH> volume_elements_count: " << volume_elements_count << std::endl;
}
