/*
 * cVCLSolver.h
 *
 * Author: jrug001
 */

#ifndef CVCLSOLVER_H_
#define CVCLSOLVER_H_

typedef double tCalcs;

#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/ilu.hpp>

typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, Eigen::Dynamic> MatrixXXC;
typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, 1> MatrixX1C;
typedef viennacl::compressed_matrix<tCalcs> vcl_sparse_t;
typedef viennacl::linalg::ilut_precond<vcl_sparse_t> vcl_precond_t;
typedef Eigen::SparseMatrix<tCalcs> SparseMatrixTCalcs;

class cVCLSolver {
public:
  cVCLSolver(SparseMatrixTCalcs &sparseA, std::fstream *fstr);
  virtual ~cVCLSolver();
  void step(MatrixX1C &solvec, MatrixX1C &rhsvec);

private:
  std::fstream *fs;
  vcl_sparse_t vcl_sparseA;    // sparse matrix for passing to ViennaCL
  int size;                    // number of columns
  vcl_precond_t  *vcl_precond; // preconditioner for passing to gmres
};

#endif /* CVCLSOLVER_H_ */
