/*
 * cCellMesh.h
 *
 *  Created on: 14/04/2015
 *      Author: jrugis
 */

#ifndef CCELLMESH_H_
#define CCELLMESH_H_

#include <string>
#include <Eigen/Dense>

#include "cVCLSolver.h"

typedef double tCoord;
typedef long tElement;

enum mesh_node_values{dist_lumen, dist_surface, MESHNCOUNT};

class cCellMesh {
public:
	cCellMesh(std::string filename);
	virtual ~cCellMesh();
	void print_info();

	tElement nodes_count;
	tElement total_elements_count, surface_elements_count, volume_elements_count;
	Eigen::Array<tCoord, Eigen::Dynamic, 3, Eigen::RowMajorBit> coordinates;
	Eigen::Array<tElement, Eigen::Dynamic, 3, Eigen::RowMajorBit> surface_elements;
	Eigen::Array<tElement, Eigen::Dynamic, 4, Eigen::RowMajorBit> volume_elements;
	Eigen::Array<tCalcs, Eigen::Dynamic, MESHNCOUNT, Eigen::RowMajorBit> node_data;
	Eigen::Array<bool, Eigen::Dynamic, 1> surface_node;

private:
	void get_mesh(std::string file_name);
	void calc_dist();
};

#endif /* CCELLMESH_H_ */
