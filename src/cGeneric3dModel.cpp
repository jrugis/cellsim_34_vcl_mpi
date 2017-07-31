/*
 * cGeneric3dModel.cpp
 *
 * Author: jrugis
 */

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#ifdef MKL_SOLVER
#include "cMKLSolver.h"
#else
#include "cVCLSolver.h"
#endif
#include "cGeneric3dModel.h"

cGeneric3dModel::cGeneric3dModel(cCellMesh *m, std::fstream *fstr) {
  mesh = m;
  fs = fstr;

  get_parameters();
  make_matrices();  // create the constant matrices
  init_u();
  *fs << "<MODEL> creating solver object..." << std::endl;
  #ifdef MKL_SOLVER
  solver = new cMKLSolver(sparseA, fs);
  #else
  solver = new cVCLSolver(sparseA, fs);
  #endif
}

cGeneric3dModel::~cGeneric3dModel() {
  delete solver;
}

void cGeneric3dModel::get_parameters(){
  std::ifstream model_file("cs"+std::to_string(mesh->rank+1)+".dat"); // open the model parameters file
  std::string line;                    // file line buffer
  std::vector <std::string> tokens;    // tokenized line

  // is the file is open?
  if (not model_file.is_open()) {
    fatal_error("the model parameters file could not be opened");
  }

  *fs << "<MODEL> reading model parameters..." << std::endl;
  int n = 0;   // read in the model parameters
  while(getline(model_file, line)){
    if(line.data()[0] == '%') continue;
    boost::split(tokens, line, boost::is_any_of(", "), boost::token_compress_on);
    if((n + tokens.size()) > PCOUNT) fatal_error("too many parameters in the model parameters file");
    for(unsigned int m = 0; m < tokens.size(); m++) p[n++] = atof(tokens[m].c_str());
  }
  model_file.close();
  if(n != PCOUNT) fatal_error("too few parameters in the model parameters file");
  numt = p[tend] / p[delt];
  plc_st = p[PLCsrt] / p[delt];
  plc_ft = p[PLCfin] / p[delt];
}

void cGeneric3dModel::init_u(){
  *fs << "<MODEL> initialising the solution matrix..." << std::endl;
  tElement np = mesh->nodes_count;
  u.resize(VARIABLES * np, numt); // NOTE: the variables ordering is c, ip, h, ce
  u.block(0, 0, np, 1) = MatrixXXC().Constant(np, 1, p[c0]);
  u.block(np, 0, np, 1) = MatrixXXC().Constant(np, 1, p[ip0]);
  u.block(2 * np, 0, np, 1) = MatrixXXC().Constant(np, 1, p[h0]);
  #ifdef FOUR_VARIABLES
  u.block(3 * np, 0, np, 1) = MatrixXXC().Constant(np, 1, p[ce0]);
  #endif
}

#ifndef FOUR_VARIABLES
Array1VC cGeneric3dModel::getbodyreactions(tCalcs c, tCalcs ip, tCalcs h, tCalcs ipr_f, tCalcs plc_f){
#else
Array1VC cGeneric3dModel::getbodyreactions(tCalcs c, tCalcs ip, tCalcs h, tCalcs ce, tCalcs ipr_f, tCalcs plc_f){
#endif
  tCalcs phi_c = pow(c, 3) / (pow(p[Kc], 3) + pow(c, 3));
  tCalcs phi_p = pow(ip, 4) / (pow(p[Kp], 4) + pow(ip, 4));
  tCalcs PO = phi_c * phi_p * h;
  tCalcs hinf = pow(p[Ki], 2) / (pow(p[Ki], 2) + pow(c, 2));
  tCalcs P_RyR = pow(c, 2) / (pow(p[KRyR], 2) + pow(c, 2));

  tCalcs JIPR = ipr_f * p[kIPR] * PO;
  tCalcs JS = p[VS] * pow(c, 2) / (pow(p[KS], 2) + pow(c, 2));
  #ifndef FOUR_VARIABLES
  tCalcs ce = (p[ct] - c) / p[gama];
  #endif
  tCalcs erc = (JIPR + p[kleak] + p[kRyR] * P_RyR) * (ce - c) - JS; // calcium from ER
  Array1VC reactions;
  reactions(0) = erc;
  reactions(1) = (plc_f * p[VPLC]) - (p[Vdeg] * ip * pow(c, 2) / (pow(p[K3K], 2) + pow(c, 2)));
  reactions(2) = (hinf - h) / p[tau];
  #ifdef FOUR_VARIABLES
  reactions(3) = -1.0 * erc / p[gama];
  #endif
  return reactions;
}

tCalcs cGeneric3dModel::getboundaryflux(tCalcs c){
  //tCalcs Jinflux = 0.005 + (p[mu] * 0.01);
  //tCalcs Jpm =  0.00005 * c;
  //return(Jinflux - Jpm);
  return(0.0); // NO BOUNDRY FLUX !!!
}

MatrixX1C cGeneric3dModel::make_load(long i){
  tElement np = mesh->nodes_count;
  ArrayX1C c, ip, h;
  ArrayX1C load_c, load_ip;
  #ifdef FOUR_VARIABLES
  ArrayX1C ce;
  ArrayX1C load_ce;
  #endif
  MatrixX1C load;

  c = u.block(0, i, np, 1);
  ip = u.block(np, i, np, 1);
  h = u.block(2 * np, i, np, 1);
  #ifdef FOUR_VARIABLES
  ce = u.block(3 * np, i, np, 1);
  #endif

  load_c = load_c.Zero(np, 1);
  load_ip = load_ip.Zero(np, 1);
  load.resize(VARIABLES * np, Eigen::NoChange);
  #ifdef FOUR_VARIABLES
  load_ce = load_ce.Zero(np, 1);
  #endif

  // reaction terms
  for(tElement n = 0; n < (mesh->volume_elements_count); n++){ // for each volume element...
    Eigen::Array<tElement,1,4> vi;      // tetrahedron vertex indices
    vi = mesh->volume_elements.block<1,4>(n, 0);

    Eigen::Array<tCoord,4,3> vert;      // tetrahedron vertex coordinates
    for(int i = 0; i < 4; i++)
      vert.block<1,3>(i, 0) = mesh->coordinates.block<1,3>(tElement(vi(i)), 0);

    tCalcs cav  = 0.25 * (c(vi(0))  + c(vi(1))  + c(vi(2))  + c(vi(3)));
    tCalcs ipav = 0.25 * (ip(vi(0)) + ip(vi(1)) + ip(vi(2)) + ip(vi(3)));
    tCalcs hav  = 0.25 * (h(vi(0)) + h(vi(1)) + h(vi(2)) + h(vi(3)));
    #ifdef FOUR_VARIABLES
    tCalcs ceav = 0.25 * (ce(vi(0)) + ce(vi(1)) + ce(vi(2)) + ce(vi(3)));

    Array1VC reactions = getbodyreactions(cav, ipav, hav, ceav,
        tCalcs(element_data(n, IPR_e)), tCalcs(element_data(n, PLC_e)));
    #else
    Array1VC reactions = getbodyreactions(cav, ipav, hav,
        tCalcs(element_data(n, IPR_e)), tCalcs(element_data(n, PLC_e)));
    #endif

    for(int i = 0; i < 4; i++){ // for each tetrahedron vertex
      load_c(vi(i)) += element_data(n, VOL_e) * 0.25 * reactions(0); // reaction terms, scaled by one-quarter volume
      load_ip(vi(i)) += element_data(n, VOL_e) * 0.25 * reactions(1);
      #ifdef FOUR_VARIABLES
      load_ce(vi(i)) += element_data(n, VOL_e) * 0.25 * reactions(3);
      #endif
    }
  }

  // add calcium boundary flux
  for(tElement n = 0; n < (mesh->surface_elements_count); n++){ // for each surface element...
    Eigen::Array<tElement,1,3> vi;      // triangle vertex indices
    vi = mesh->surface_elements.block<1,3>(n, 0);

    Eigen::Array<tCoord,3,3> vert;      // triangle vertex coordinates
    for(int i = 0; i < 3; i++)
      vert.block<1,3>(i, 0) = mesh->coordinates.block<1,3>(tElement(vi(i)), 0); // why is typecast needed???

    Eigen::Matrix<tCoord,1,3> side1 = vert.block<1,3>(0, 0) - vert.block<1,3>(1, 0);
    Eigen::Matrix<tCoord,1,3> side2 = vert.block<1,3>(0, 0) - vert.block<1,3>(2, 0);
    tCalcs triangle_area = 0.5 * (side1.cross(side2)).norm();

    tCalcs cav  = (c(vi(0)) + c(vi(1)) + c(vi(2))) / 3.0;
    tCalcs bndy_flux = getboundaryflux(cav);
    for(int i = 0; i < 3; i++) // for each triangle vertex
      load_c(vi(i)) += triangle_area / 3.0 * bndy_flux; // integrate over triangle surface
  }

  // the diffusing variables
  load.block(0, 0, np, 1) = load_c;
  load.block(np, 0, np, 1) = load_ip;
  #ifdef FOUR_VARIABLES
  load.block(3*np, 0, np, 1) = load_ce;
  #endif

  // the non-diffusing variables
  for(tElement n = 0; n < np; n++){ // for each node...
    Array1VC reactions;
    #ifndef FOUR_VARIABLES
    reactions = getbodyreactions(tCalcs(c(n)), tCalcs(ip(n)), tCalcs(h(n)),
      tCalcs(node_data(n, IPR_n)), tCalcs(node_data(n, PLC_n)));
    #else
    reactions = getbodyreactions(tCalcs(c(n)), tCalcs(ip(n)), tCalcs(h(n)), tCalcs(ce(n)),
      tCalcs(node_data(n, IPR_n)), tCalcs(node_data(n, PLC_n)));
    #endif
    load((2 * np) + n) = reactions(2); // h
  }
  //save_matrix("load.txt", load);
  return load;
}

void cGeneric3dModel::make_matrices(){
  tElement np = mesh->nodes_count;

  node_data.resize(np, Eigen::NoChange);
  std::string f1 = "ipr_REF"+std::to_string(mesh->rank+1)+".bin"; // optional pre-calculated spatial reference data
  std::string f2 = "plc_REF"+std::to_string(mesh->rank+1)+".bin"; // optional pre-calculated spatial reference data

  // EITHER read in the optional node spatial factor array
  if((access(f1.c_str(), F_OK) != -1) and (access(f2.c_str(), F_OK) != -1)){ 
    *fs << "<MODEL> reading in optional spatial factors..." << std::endl;
    load_node_data(f1, IPR_n);
    load_node_data(f2, PLC_n);
  }
  // OR calculate the node spatial factor array
  else{
    *fs << "<MODEL> calculating the spatial factors..." << std::endl;
    for(tElement n = 0; n < (np); n++){ // for each node...
      tCalcs d = mesh->node_data(n, dist_lumen);
      if(d <= p[IPRdn]) node_data(n, IPR_n) =  1.0; // linear gradient between d1 and d2
      else if(d >= p[IPRdf]) node_data(n, IPR_n) =  p[IPRmin];
      else node_data(n, IPR_n) = 1.0 - (1.0 - p[IPRmin]) * ((d - p[IPRdn]) / (p[IPRdf] - p[IPRdn]));
      node_data(n, PLC_n) = ((mesh->node_data(n, dist_surface) < p[PLCds])
        and (mesh->node_data(n, dist_lumen) >= p[PLCdl])) ? 1.0 : 0.0;
    }
    save_matrix("ipr"+std::to_string(mesh->rank+1)+".bin", MatrixXXC(node_data.block(0, IPR_n, np, 1)));
    save_matrix("plc"+std::to_string(mesh->rank+1)+".bin", MatrixXXC(node_data.block(0, PLC_n, np, 1)));
  }

  // make the reference mass matrix
  ArrayRefMass ref_mass;
  ref_mass = make_ref_mass();

  // make the mass and stiffness matrices
  *fs << "<MODEL> calculating the constant matrices..." << std::endl;
  MatrixXXC stiffc, stiffp, small_mass;
  #ifdef FOUR_VARIABLES
  MatrixXXC stiffce;
  #endif
  stiffc = stiffc.Zero(np, np);
  stiffp = stiffp.Zero(np, np);
  #ifdef FOUR_VARIABLES
  stiffce = stiffce.Zero(np, np);
  #endif
  small_mass = small_mass.Zero(np, np);
  element_data.resize(mesh->volume_elements_count, Eigen::NoChange);

  for(tElement n = 0; n < (mesh->volume_elements_count); n++){ // for each volume element...
    Eigen::Matrix<tElement,1,4> vi;      // tetrahedron vertex indices
    vi = mesh->volume_elements.block<1,4>(n, 0);

    element_data(n, IPR_e) = 0.25 *     // element ipr spatial factor
      (node_data(vi[0], IPR_n) + node_data(vi[1], IPR_n) + node_data(vi[2], IPR_n) + node_data(vi[3], IPR_n));
    element_data(n, PLC_e) = 0.25 *     // element plc spatial factor
      (node_data(vi[0], PLC_n) + node_data(vi[1], PLC_n) + node_data(vi[2], PLC_n) + node_data(vi[3], PLC_n));

    Eigen::Matrix<tCoord,4,3> vert; // tetrahedron vertex coordinates
    for(int i = 0; i < 4; i++)
      vert.block<1,3>(i, 0) = mesh->coordinates.block<1,3>(tElement(vi(i)), 0); // why is typecast needed???

    Eigen::Matrix<tCoord,3,3> J;    // tetrahedron edge vectors
    for(int i = 0; i < 3; i++)
      J.block<1,3>(i, 0) = vert.block<1,3>(i + 1, 0) - vert.block<1,3>(0, 0);
    tCalcs V, Vx6;                  // tetrahedron volume, (6x) volume
    Vx6 = J.determinant();
    V = Vx6 / 6.0;
    element_data(n, VOL_e) = V;   // save the tetrahedron volume

    tCalcs Ic = V * p[Dc]; // diffusion coefficients
    tCalcs Ip = V * p[Dp];
    #ifdef FOUR_VARIABLES
    tCalcs Ice = V * p[Dce];
    #endif

    Eigen::Matrix<tCoord,4,4> M, C, G;
    M.col(0) << 1, 1, 1, 1;
    M.block<4,3>(0, 1) = vert;
    C = M.inverse();
    G = C.block<3,4>(1, 0).transpose() * C.block<3,4>(1, 0); // gradients of the basis functions

    // construct the mass and stiffness matrix components
    for(int i = 0; i < 4; i++){
      stiffc(vi(i), vi(i)) += G(i, i) * Ic;
      stiffp(vi(i), vi(i)) += G(i, i) * Ip;
      #ifdef FOUR_VARIABLES
      stiffce(vi(i), vi(i)) += G(i, i) * Ice;
      #endif
      small_mass(vi(i), vi(i)) += ref_mass(i, i) * Vx6;
    }
    for(int i = 0; i < 3; i++){
      stiffc(vi(0), vi(i + 1)) += G(0, i + 1) * Ic;
      stiffp(vi(0), vi(i + 1)) += G(0, i + 1) * Ip;
      #ifdef FOUR_VARIABLES
      stiffce(vi(0), vi(i + 1)) += G(0, i + 1) * Ice;
      #endif
      small_mass(vi(0), vi(i + 1)) += ref_mass(0, i + 1) * Vx6;
      stiffc(vi(i + 1), vi(0)) = stiffc(vi(0), vi(i + 1));
      stiffp(vi(i + 1), vi(0)) = stiffp(vi(0), vi(i + 1));
      #ifdef FOUR_VARIABLES
      stiffce(vi(i + 1), vi(0)) = stiffce(vi(0), vi(i + 1));
      #endif
      small_mass(vi(i + 1), vi(0)) = small_mass(vi(0), vi(i + 1));
    }
    for(int i = 0; i < 2; i++){
      stiffc(vi(1), vi(i + 2)) += G(1, i + 2) * Ic;
      stiffp(vi(1), vi(i + 2)) += G(1, i + 2) * Ip;
      #ifdef FOUR_VARIABLES
      stiffce(vi(1), vi(i + 2)) += G(1, i + 2) * Ice;
      #endif
      small_mass(vi(1), vi(i + 2)) += ref_mass(1, i + 2) * Vx6;
      stiffc(vi(i + 2), vi(1)) = stiffc(vi(1), vi(i + 2));
      stiffp(vi(i + 2), vi(1)) = stiffp(vi(1), vi(i + 2));
      #ifdef FOUR_VARIABLES
      stiffce(vi(i + 2), vi(1)) = stiffce(vi(1), vi(i + 2));
      #endif
      small_mass(vi(i + 2), vi(1)) = small_mass(vi(1), vi(i + 2));
    }
    stiffc(vi(2), vi(3)) += G(2, 3) * Ic;
    stiffp(vi(2), vi(3)) += G(2, 3) * Ip;
    #ifdef FOUR_VARIABLES
    stiffce(vi(2), vi(3)) += G(2, 3) * Ice;
    #endif
    small_mass(vi(2), vi(3)) += ref_mass(2, 3) * Vx6;
    stiffc(vi(3), vi(2)) = stiffc(vi(2), vi(3));
    stiffp(vi(3), vi(2)) = stiffp(vi(2), vi(3));
    #ifdef FOUR_VARIABLES
    stiffce(vi(3), vi(2)) = stiffce(vi(2), vi(3));
    #endif
    small_mass(vi(3), vi(2)) = small_mass(vi(2), vi(3));
  }

  // construct the mass matrix from a list of triplets (non zero elements)
  std::vector<Triplet> triplet_list;
    
  // construct list of triplets
  int np2 = np * 2;
  #ifdef FOUR_VARIABLES
  int np3 = np * 3;
  #endif
  for (int j = 0; j < np; j++) {
    for (int i = 0; i < np; i++) {
      double v_ij = small_mass(i, j);
      // add non zeros in first, second and fourth (if required) blocks
      if (v_ij != 0) {
        // mass.block(0, 0, np, np) = small_mass;
        triplet_list.push_back(Triplet(i, j, v_ij));
                
        // mass.block(np, np, np, np) = small_mass;
        triplet_list.push_back(Triplet(np + i, np + j, v_ij));
                
        #ifdef FOUR_VARIABLES
        // mass.block(3*np, 3*np, np, np) = small_mass;
        triplet_list.push_back(Triplet(np3 + i, np3 + j, v_ij));
        #endif
      }
            
      if (i == j) {
        // set identity in third block
        triplet_list.push_back(Triplet(np2 + i, np2 + j, 1.0));
      }
    }
  }
    
  // create sparse mass matrix from triplets
  sparseMass.resize(VARIABLES * np, VARIABLES * np);
  sparseMass.setFromTriplets(triplet_list.begin(), triplet_list.end());
    
  // construct stiff matrix from list of triplets
  triplet_list.clear();
    
  // construct list of triplets
  for (int j = 0; j < np; j++) {
    for (int i = 0; i < np; i++) {
      // first block
      // stiff.block(0, 0, np, np) = stiffc;
      double v_stiffc = stiffc(i, j);
      if (v_stiffc != 0) {
        triplet_list.push_back(Triplet(i, j, v_stiffc));
      }
            
      // second block
      // stiff.block(np, np, np, np) = stiffp;
      double v_stiffp = stiffp(i, j);
      if (v_stiffp != 0) {
        triplet_list.push_back(Triplet(np + i, np + j, v_stiffp));
      }
            
      #ifdef FOUR_VARIABLES
      // fourth block
      // stiff.block(3*np, 3*np, np, np) = stiffce;
      double v_stiffce = stiffce(i, j);
      if (v_stiffce != 0) {
        triplet_list.push_back(Triplet(np3 + i, np3 + j, v_stiffce));
      }
      #endif
    }
  }

  // create sparse stiff matrix from triplets
  SparseMatrixTCalcs sparseStiff(VARIABLES * np, VARIABLES * np);
  sparseStiff.setFromTriplets(triplet_list.begin(), triplet_list.end());
    
  // make the A matrix
  sparseA.resize(VARIABLES * np, VARIABLES * np);
  sparseA = sparseMass + (p[delt] * sparseStiff);
}

ArrayRefMass cGeneric3dModel::make_ref_mass(){
  ArrayRefMass ref_mass;
  tCalcs v = (1.0 / 6.0) * 0.25 * 0.25;
  for(int i = 0; i < REF_MASS_SIZE; i++){
    for(int j = 0; j < REF_MASS_SIZE; j++){
      ref_mass(i, j) = v;
    }
  }
  return ref_mass;
}

void cGeneric3dModel::run(){
  *fs << "<MODEL> running the model..." << std::endl;
  MatrixX1C solvec, rhsvec; // the solution and right-hand side vectors
  solvec.resize(VARIABLES * mesh->nodes_count, Eigen::NoChange);
  rhsvec.resize(VARIABLES * mesh->nodes_count, Eigen::NoChange);
  solvec = u.col(0);
  tCalcs plc = p[VPLC];
  for(long i = 1; i < numt; i++){
    if(mesh->rank == 0){
      std::cout << std::fixed << std::setprecision(3) << i * p[delt] << std::endl;
      std::cout.flush();
    }
    *fs << std::fixed << std::setprecision(3) << i * p[delt] << " ";
    p[VPLC] = ((i < plc_st) or (i > plc_ft)) ? 0.0 : plc;
    rhsvec = (sparseMass * solvec) + (p[delt] * make_load(i - 1));
    //*********************************************************
    //solvec = Amat.llt().solve(rhsvec);
    //solvec = Amat.ldlt().solve(rhsvec);
    //solvec = Amat.partialPivLu().solve(rhsvec);
    //solvec = Amat.fullPivHouseholderQr().solve(rhsvec);
    solver->step(solvec, rhsvec);
    //*********************************************************
    u.col(i) = solvec;
  }
  //save_matrix("mass_ei.bin", mass);
  //save_matrix("rhsvec_ei.bin", rhsvec);
  //save_matrix("Amat_ei.bin", Amat);
}

void cGeneric3dModel::load_node_data(std::string file_name, int dindex){
  std::ifstream file(file_name.c_str(), std::ios::binary); // open the file
  tElement rows, cols;
  float f;
  int esize = sizeof(f);
  file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
  file.read(reinterpret_cast<char*>(&cols), sizeof(cols));
  if((rows != mesh->nodes_count) or (cols != 1)){
    *fs << "<MODEL> ERROR: load_node_data size mismatch" << std::endl;
    exit(0);
  }
  for(int i = 0; i < rows; i++){
    file.read(reinterpret_cast<char*>(&f), esize); // column order
    node_data(i, dindex) = tCalcs(f);
  }
  file.close();
}

void cGeneric3dModel::save_matrix_reduce(std::string file_name, MatrixXXC mat){
  // peak value for each row
  MatrixX1C max_per_row = mat.rowwise().maxCoeff();
    
  // which rows have the highest and lowest peaks
  int max_index, min_index;
  max_per_row.maxCoeff(&max_index);
  max_per_row.minCoeff(&min_index);
  *fs << "<MODEL> Reducing results in " << file_name << ": " << max_index << " and " << min_index << std::endl;
    
  // write data
  std::ofstream file(file_name.c_str(), std::ios::binary); // create the file
  tElement rows = 2;
  tElement cols = mat.cols();
  float f;
  int esize = sizeof(f);
  file.write(reinterpret_cast<char*>(&rows), sizeof(rows));
  file.write(reinterpret_cast<char*>(&cols), sizeof(cols));
  for(int j = 0; j < cols; j++){
    f = mat(max_index, j);  // convert to float for smaller file size
    file.write(reinterpret_cast<char*>(&f), esize); // column order
    f = mat(min_index, j);  // convert to float for smaller file size
    file.write(reinterpret_cast<char*>(&f), esize); // column order
  }
  file.close();
}

void cGeneric3dModel::save_matrix(std::string file_name, MatrixXXC mat){
  std::ofstream file(file_name.c_str(), std::ios::binary); // create the file
  tElement rows = mat.rows();
  tElement cols = mat.cols();
  float f;
  int esize = sizeof(f);
  file.write(reinterpret_cast<char*>(&rows), sizeof(rows));
  file.write(reinterpret_cast<char*>(&cols), sizeof(cols));
  for(int j = 0; j < cols; j++){
    for(int i = 0; i < rows; i++){
      f = mat(i, j);  // convert to float for smaller file size
      file.write(reinterpret_cast<char*>(&f), esize); // column order
    }
  }
  file.close();
}

void cGeneric3dModel::save_results(){
  tElement np = mesh->nodes_count;
    
  // reduce first or output full results
  if (p[reduce]) {
    save_matrix_reduce("cR"+std::to_string(mesh->rank)+".bin", MatrixXXC(u.block(0, 0, np, numt)));     // calcium
    save_matrix_reduce("ip3R"+std::to_string(mesh->rank)+".bin", MatrixXXC(u.block(np, 0, np, numt)));  // ip3
  }
  else {
    save_matrix("c"+std::to_string(mesh->rank)+".bin", MatrixXXC(u.block(0, 0, np, numt)));     // calcium
    save_matrix("ip3"+std::to_string(mesh->rank)+".bin", MatrixXXC(u.block(np, 0, np, numt)));  // ip3
  }

  //save_matrix("d"+std::to_string(mesh->rank)+".bin", MatrixXXC(u.block(2 * np, 0, np, numt)));  // d
  //save_matrix("ce"+std::to_string(mesh->rank)+".bin", MatrixXXC(u.block(3 * np, 0, np, numt)));  // ce
}

void cGeneric3dModel::fatal_error(std::string msg){
  *fs << "<MODEL> ERROR: " << msg << std::endl;
  exit(0);
}
