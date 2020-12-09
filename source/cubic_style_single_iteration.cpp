#include "cubic_style_single_iteration.h"
#include <igl/slice.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/columnize.h>
#include <math.h>
#include <iostream>

void cubic_style_single_iteration(
    cubic_style_data & data,
    Eigen::MatrixXd & U) 
{
    // [1] ADMM paper: https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf 
    // [2] ADMM lasso matlob code example: https://web.stanford.edu/~boyd/papers/admm/lasso/lasso.html     
    //      The code follows these two reference a lot!
    // [3] Cubic Stylization paper: https://www.dgp.toronto.edu/projects/cubic-stylization/ 
    // Local step for R
    Eigen::MatrixXd Rs(3,U.rows()*3);

    int MAX_ITER = 1000; 
    for (int i=0; i<U.rows(); i++) {
        // Initial values
        Eigen::Matrix3d Ri;
        Eigen::VectorXd z = data.zs.col(i);
        Eigen::VectorXd u = data.us.col(i);
        double rho = data.rhos(i);

        // compute Mi
        // Mi = [Di ni][Wi 0, 0 rho][Dtildei' (z-u)'] [3]
        Eigen::MatrixXd Di = data.D[i];
        Eigen::VectorXd ni = data.N.row(i).transpose();
        Eigen::VectorXd Wi = data.W[i];
        
        // Get the distance of new positives by half edges
        Eigen::MatrixXd hei, hej;
        igl::slice(U, data.he[i].col(0), 1, hei);
        igl::slice(U, data.he[i].col(1), 1, hej);
        Eigen::MatrixXd Dtildei(3, data.he[i].rows());
        Dtildei = (hej - hei).transpose();

        // Hint: Mi = D * W * Dtildei' + ni * rho * (z-u)
        // Atb in [2]
        Eigen::Matrix3d firstterm = Di * Wi.asDiagonal() * Dtildei.transpose();

        // ADMM
        for (int k=0; k<MAX_ITER; k++) {
            // update Ri = vi ui' (svd of Mi) [3]
            Eigen::Matrix3d Mi = firstterm + (ni * rho * (z-u).transpose());
            // SVD of Mi
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(Mi, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Matrix3d v_svd = svd.matrixV();
            Eigen::Matrix3d u_svd = svd.matrixU();
            Ri = v_svd * u_svd.transpose();
            // Ensure that the det(Ri) > 0 
            if (Ri.determinant() < 0){
                u_svd.col(2) = -u_svd.col(2);
                Ri = v_svd * u_svd.transpose();
            }

            // update z lasso problem
            // See sec 4.4.3 [1] & z-update with relaxatio [2]
            Eigen::VectorXd zold = z;
            Eigen::VectorXd x = Ri * ni + u;
            double kappa = data.lambda * data.a(i)/rho;
            // shrinkage
            auto shrinkage = [](Eigen::VectorXd & x, double & k){
                Eigen::VectorXd z = (x.array()-k).array().max(0.0) - (-x.array()-k).array().max(0.0);
                return z;
            };
            z = shrinkage(x, kappa);
            
            // update u = u + Ri*ni - z [3]
            u.noalias() += Ri * ni - z; // From fit_rotations_l1.cpp in CubicStylization code, noalias() more efficient

            // update rho, u
            // See sec 3.4.1 [1]
            double r = (Ri*ni-z).norm();
            double s = (-rho*(z - zold)).norm();

            if (r > data.mu * s) {
                rho = data.taoincr * rho;
                u = u/data.taoincr;
            } else if (s > data.mu * r) {
                rho = rho / data.taodecr;
                u = u * data.taodecr;
            }

            // Early stop
            // See sec 3.4.1 [1] & termination checks [2]
            double n = z.size();
            double eps_pri = sqrt(1.0*n) * data.epsabs + data.epsrel * std::max((Ri*ni).norm(), z.norm());
            double eps_dual = sqrt(1.0*n) * data.epsabs + data.epsrel * ((rho * u).norm());

            if ((r < eps_pri) && (s < eps_dual)) {
                // stop here
                data.zs.col(i) = z;
                data.us.col(i) = u;
                data.rhos(i) = rho;
                Rs.block(0, 3*i, 3, 3) = Ri;
                break;
            }
        }
    }

    // Global step
    // std::cout << 'V ' << U.rows() << std::endl;
    // => n
    // std::cout << 'K ' << data.K.rows() << ' ' << data.K.cols() << std::endl;
    // => 3n * 3(3n)
    // std::cout << 'r ' << Rs.rows() << ' ' << Rs.cols() << std::endl;  
    // => 3 * 3n
    Eigen::VectorXd Rcol;
    Eigen::VectorXd Beq;
    igl::columnize(Rs, U.rows(), 2, Rcol);
    Eigen::VectorXd B = data.K * Rcol; // 3(3n) * 1
    Eigen::MatrixXd BT(B.size()/3, 3);
    BT.col(0) = B.block(0, 0, U.rows(), 1);
    BT.col(1) = B.block(1*U.rows(), 0, U.rows(), 1);
    BT.col(2) = B.block(2*U.rows(), 0, U.rows(), 1);
    igl::min_quad_with_fixed_solve(data.data,BT,data.bc,Beq,U);
}
    
