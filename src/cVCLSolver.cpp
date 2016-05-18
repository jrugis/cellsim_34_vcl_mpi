/*
 * cVCLSolver.cpp
 *
 *  Created on: Mar 7, 2016
 *      Author: jrug001
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Sparse>

// enable Eigen wrappers within ViennaCL
#define VIENNACL_HAVE_EIGEN
#include "viennacl/linalg/gmres.hpp"

#include "cVCLSolver.h"
#include "cCellMesh.h"

cVCLSolver::cVCLSolver(MatrixXXC &Amat){
	std::cout << "<SOLVER> initialising the solver..." << std::endl;
    
    // convert A matrix to sparse format
    Eigen::SparseMatrix<tCalcs, Eigen::RowMajor> sparseA;
    sparseA = Amat.sparseView();
    
    // copy A matrix into vcl matrix
    viennacl::copy(sparseA, vcl_sparseA);
    
    // number of columns
    size = Amat.cols();
    
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
    viennacl::vector<tCalcs> vcl_sol;
    viennacl::vector<tCalcs> vcl_rhs(size);
    viennacl::copy(rhsvec, vcl_rhs);
    
    // solve
    viennacl::linalg::gmres_tag my_gmres_tag(1e-8, 500, 15);
    vcl_sol = viennacl::linalg::solve(vcl_sparseA, vcl_rhs, my_gmres_tag, *vcl_precond);
    
    // copy/cast data back to Eigen vector
    viennacl::copy(vcl_sol, evec_sol);
    solvec = static_cast<MatrixX1C>(evec_sol);
    
    std::cout << my_gmres_tag.iters() << " " << my_gmres_tag.error() << std::endl;
}
