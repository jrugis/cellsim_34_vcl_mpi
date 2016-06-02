/*
 * cMKLSolver.h
 *
 *  Created on: May 10, 2016
 *      Author: Chris Scott
 */

#ifndef CMKLSOLVER_H_
#define CMKLSOLVER_H_

typedef double tCalcs;

#include <string>
#include <Eigen/Dense>
#include <mkl.h>


typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, Eigen::Dynamic> MatrixXXC;
typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, 1> MatrixX1C;
typedef Eigen::SparseMatrix<tCalcs> SparseMatrixTCalcs;


class cMKLSolver {
public:
	cMKLSolver(MatrixXXC &Amat);
	virtual ~cMKLSolver();
	void step(MatrixX1C &solvec, MatrixX1C &rhsvec);

private:
    void fatal_error(std::string msg);
    
    // storage for sparse matrix
    double *Acsr;
    MKL_INT *Aj, *Ai;
    
    // storage for preconditioner
    double *bilut;
    MKL_INT *ibilut, *jbilut;
    
    // size of the system
    MKL_INT size;
    
    // gmres parameters
    MKL_INT gmres_restarts;
    double gmres_relative_tol, gmres_absolute_tol;
    
    // gmres workspace
    double *gmres_tmp;
    double *gmres_trvec;
};

#endif /* CMKLSOLVER_H_ */
