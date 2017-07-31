/*
 * cCellMesh.h
 *
 * Author: jrugis
 */

#ifndef CCELLMESH_H_
#define CCELLMESH_H_

#include <fstream>
#include <string>
#include <Eigen/Dense>

#ifdef MKL_SOLVER
#include "cMKLSolver.h"
#else
#include "cVCLSolver.h"
#endif

typedef double tCoord;
typedef long tElement;

enum mesh_node_values{dist_lumen, dist_surface, MESHNCOUNT};

class cCellMesh {
public:
  cCellMesh(int commRank, std::fstream *fstr);
  virtual ~cCellMesh();
  void print_info();

  int rank;
  tElement nodes_count;
  tElement total_elements_count, surface_elements_count, volume_elements_count;
  Eigen::Array<tCoord, Eigen::Dynamic, 3, Eigen::RowMajorBit> coordinates;
  Eigen::Array<tElement, Eigen::Dynamic, 3, Eigen::RowMajorBit> surface_elements;
  Eigen::Array<tElement, Eigen::Dynamic, 4, Eigen::RowMajorBit> volume_elements;
  Eigen::Array<tCalcs, Eigen::Dynamic, MESHNCOUNT, Eigen::RowMajorBit> node_data;
  Eigen::Array<bool, Eigen::Dynamic, 1> surface_node;

private:
  std::fstream *fs;
  std::string filename;
  void get_mesh(std::string file_name);
  void calc_dist();
};

#endif /* CCELLMESH_H_ */
