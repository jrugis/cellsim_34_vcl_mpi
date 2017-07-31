/*
 * cVCLSolver.cpp
 *
 * Author: jrug001
 */

#ifndef MKL_SOLVER

#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Sparse>

// enable Eigen wrappers within ViennaCL
#define VIENNACL_HAVE_EIGEN
#include "viennacl/linalg/gmres.hpp"

#include "cVCLSolver.h"
#include "cCellMesh.h"

cVCLSolver::cVCLSolver(SparseMatrixTCalcs &sparseA, std::fstream *fstr){
  fs = fstr;
  *fs << "<SOLVER> initialising the VCL solver..." << std::endl;
    
  // copy A matrix into vcl matrix
  viennacl::copy(sparseA, vcl_sparseA);
  #ifdef VIENNACL_WITH_CUDA
  viennacl::context cuda_context(viennacl::CUDA_MEMORY);
  viennacl::switch_memory_context(vcl_sparseA, cuda_context);
  #endif
 
  // number of columns
  size = sparseA.cols();
    
  // preconditioner configuration
  viennacl::linalg::ilut_tag precond_conf(40, 1e-5);
    
  // precondition the matrix
  vcl_precond = new vcl_precond_t(vcl_sparseA, precond_conf);
}

cVCLSolver::~cVCLSolver(){
  delete vcl_precond;
}

void cVCLSolver::step(MatrixX1C &solvec, MatrixX1C &rhsvec){
  // cast to VectorXd for vcl copy
  Eigen::VectorXd evec_rhs = static_cast<Eigen::VectorXd>(rhsvec);
  Eigen::VectorXd evec_sol = static_cast<Eigen::VectorXd>(solvec);

  // copy data to viennacl vector
  #ifdef VIENNACL_WITH_CUDA
  viennacl::context cuda_context(viennacl::CUDA_MEMORY);
  viennacl::vector<tCalcs> vcl_sol(size, cuda_context);
  viennacl::vector<tCalcs> vcl_rhs(size, cuda_context);
  #else
  viennacl::vector<tCalcs> vcl_sol;
  viennacl::vector<tCalcs> vcl_rhs(size);
  #endif
  viennacl::copy(rhsvec, vcl_rhs);
    
  // solve
  viennacl::linalg::gmres_tag my_gmres_tag(1e-8, 500, 15);
  vcl_sol = viennacl::linalg::solve(vcl_sparseA, vcl_rhs, my_gmres_tag, *vcl_precond);
    
  // copy/cast data back to Eigen vector
  viennacl::copy(vcl_sol, evec_sol);
  solvec = static_cast<MatrixX1C>(evec_sol);

  *fs << my_gmres_tag.iters() << " " << my_gmres_tag.error() << std::endl;
}

#endif
