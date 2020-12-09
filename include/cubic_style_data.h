#ifndef CUBIC_SYTLE_DATA_H
#define CUBIC_SYTLE_DATA_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/min_quad_with_fixed.h>

struct cubic_style_data
{
    // Declare the members that can be precomputed
    // Write some brief documentation for each component you have defined
    
    // Initilize some parameters according to section 3.1
    double lambda = 0.0; // control the cubeness
    double rho = 1e-4; // update value for ADMM
    double epsabs = 1e-5; // absolute criterion for ADMM
    double epsrel = 1e-3; // relative criterion for ADMM
    double mu = 10; // mu value for penalty rho updata
    double tauincr = 2; // tao value for penalty rho updata
    double taudecr = 2; // tao value for penalty rho updata
    Eigen::MatrixXd zs, us; // update value for ADMM
    Eigen::VectorXd rhos; // rho values
    
    // For ASAP
    Eigen::MatrixXd N, a; // n denotes the unit area-weighted normal vecotr
                          // a is the barycentric area
    Eigen::SparseMatrix<double> L, K; // L is familiar cotangent discrete Laplacian matrix
                                      // K contains cotangents multiplied against differences across edges in the rest mesh
    igl::min_quad_with_fixed_data<double> data; // data

    std::vector<Eigen::MatrixXi> he; // Store the half-dege information
	std::vector<Eigen::MatrixXd> D; // Dij is edge vector of vertices i,j
	std::vector<Eigen::VectorXd> W; // Wij is the cotangent weight
    
    Eigen::VectorXi b;
    Eigen::MatrixXd bc;
};

#endif
