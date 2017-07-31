/*
 * cGeneric3dModel.h
 *
 * Author: jrugis
 */

#ifndef CGENERIC3DMODEL_H_
#define CGENERIC3DMODEL_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include "cCellMesh.h"
#ifdef MKL_SOLVER
#include "cMKLSolver.h"
#else
#include "cVCLSolver.h"
#endif

#ifndef FOUR_VARIABLES
#define VARIABLES 3 // the number of model variables
#else
#define VARIABLES 4 // the number of model variables
#endif

#define REF_MASS_SIZE 4   // reference mass dimension

// the generic 3D model parameters
#ifndef FOUR_VARIABLES
enum model_parameters{delt, tend, reduce, \
	PLCsrt, PLCfin, \
	IPRdn, IPRdf, IPRmin, PLCds, PLCdl, \
	c0, ip0, h0, ct, gama, \
	kIPR, Kc, Kp, Ki, tau, Dc, \
	kRyR, KRyR, \
	VS, KS, kleak, \
	VPLC, Vdeg, K3K, Dp, \
	PCOUNT};
#else
enum model_parameters{delt, tend, reduce, \
	PLCsrt, PLCfin, \
	IPRdn, IPRdf, IPRmin, PLCds, PLCdl, \
	c0, ip0, h0, ce0, gama, \
	kIPR, Kc, Kp, Ki, tau, Dc, Dce, \
	kRyR, KRyR, \
	VS, KS, kleak, \
	VPLC, Vdeg, K3K, Dp, \
	PCOUNT};
#endif

enum model_node_values{IPR_n, PLC_n, MODELNCOUNT};            // node spatial factors
enum model_element_values{VOL_e, IPR_e, PLC_e, MODELECOUNT};  // element volume and spatial factors

// some convenience typedefs
typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, Eigen::Dynamic> MatrixXXC;
typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, 1> MatrixX1C;
typedef Eigen::Array<tCalcs, Eigen::Dynamic, 1> ArrayX1C;
typedef Eigen::Array<tCalcs, 1, VARIABLES> Array1VC;
typedef Eigen::Array<tCalcs, REF_MASS_SIZE, REF_MASS_SIZE> ArrayRefMass;
typedef Eigen::Triplet<tCalcs> Triplet;

class cGeneric3dModel {
public:
  cGeneric3dModel(cCellMesh *mesh, std::fstream *fstr);
  virtual ~cGeneric3dModel();
  void run();
  void save_results();

  MatrixXXC u; // solution matrix
  SparseMatrixTCalcs sparseA, sparseMass; // A and mass matrices

private:
  std::fstream *fs;
  void get_parameters();
  void init_u();
  MatrixX1C make_load(long i);
  ArrayRefMass make_ref_mass();
#ifndef FOUR_VARIABLES
  Array1VC getbodyreactions(tCalcs c, tCalcs ip, tCalcs h, tCalcs ipr_f, tCalcs plc_f);
#else
  Array1VC getbodyreactions(tCalcs c, tCalcs ip, tCalcs h, tCalcs ce, tCalcs ipr_f, tCalcs plc_f);
#endif
  tCalcs getboundaryflux(tCalcs c);
  void make_matrices();
  void load_node_data(std::string file_name, int dindex);
  void save_matrix(std::string file_name, MatrixXXC mat);
  void save_matrix_reduce(std::string file_name, MatrixXXC mat);
  void fatal_error(std::string msg);

  cCellMesh *mesh;
#ifdef MKL_SOLVER
  cMKLSolver *solver;
#else
  cVCLSolver *solver;
#endif
  tCalcs p[PCOUNT]; // the model parameters array
  long numt, plc_st, plc_ft; // number of time steps, PLC start and finish time steps
  Eigen::Array<tCalcs, Eigen::Dynamic, MODELNCOUNT> node_data;
  Eigen::Array<tCalcs, Eigen::Dynamic, MODELECOUNT> element_data;
};

#endif /* CGENERIC3DMODEL_H_ */

