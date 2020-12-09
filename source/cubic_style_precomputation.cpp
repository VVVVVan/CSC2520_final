#include "cubic_style_precomputation.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/arap_linear_block.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/massmatrix.h>
#include <igl/arap_rhs.h>
#include <igl/vertex_triangle_adjacency.h>
#include <vector>
#include <igl/slice.h>

void cubic_style_precomputation(
	const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    cubic_style_data & data)
{
    // Initial z, u
    data.zs.resize(3, V.rows());
    data.zs.setZero();
    data.us.resize(3, V.rows());
    data.us.setZero();
    data.rhos.resize(V.rows());
    data.rhos.setConstant(data.rho);


    // Pre-compute variables in ARAP term
    // tr(V'LV) + tr(V'KR) according to asap deformation in class
    // compute L
    igl::cotmatrix(V,F,data.L);
    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute(data.L, data.b, Aeq, false, data.data);
    
    // compute K <= N(i) 'spokes and rims' edges of the ith vertex
    // https://libigl.github.io/libigl-python-bindings/igl_docs/#arap_rhs 
    igl::arap_rhs(V, F, V.cols(), igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS, data.K);
    
    // compute Di, Wi <= "spokes and rims"
    // Get the faces of a vertex
    std::vector<std::vector<int>> VF, VI;
    igl::vertex_triangle_adjacency(V.rows(),F,VF,VI);

    data.he.resize(V.rows());
    data.D.resize(V.rows());
    data.W.resize(V.rows());
    for (int i=0; i<V.rows(); i++) {
        // Get the adjacent faces
        std::vector<int> adjVF = VF[i];
        // Half edge is stored to compute Dtilta in iteration
        data.he[i].resize(adjVF.size()*3, 2);
        data.D[i].resize(3, adjVF.size()*3); // According to M in paper, Di houuld be 3*3adjVF
        data.W[i].resize(adjVF.size()*3);
        
        for (int j=0; j<adjVF.size(); j++) {
            // Get the index of each face
            int v0 = F(adjVF[j],0);
            int v1 = F(adjVF[j],1);
            int v2 = F(adjVF[j],2);

            // compute half edge for i
            Eigen::MatrixXi tmp1(3,2);
            tmp1 << v0, v1,
                    v1, v2,
                    v2, v0;
            data.he[i].block(3*j, 0, 3, 2) = tmp1;

            // // compute Di also works
            // Eigen::MatrixXd tmp2(3, 3);
            // tmp2.col(0) = V.row(v1) - V.row(v0);
            // tmp2.col(1) = V.row(v2) - V.row(v1);
            // tmp2.col(2) = V.row(v0) - V.row(v2);
            // data.D[i].block(0, 3*j, 3, 3) = tmp2;

            // compute w for i
            Eigen::VectorXd tmp3(3);
            tmp3 << data.L.coeff(v0, v1),
                    data.L.coeff(v1, v2),
                    data.L.coeff(v2, v0);
            data.W[i].block(3*j, 0, 3, 1) = tmp3; 
        }
        // compute distance between each pair of half edges
        // According to cube_style_precomputation.cpp in CubicStylization code
        Eigen::MatrixXd hei, hej;
        igl::slice(V, data.he[i].col(0), 1, hei);
        igl::slice(V, data.he[i].col(1), 1, hej);
        data.D[i] = (hej - hei).transpose();
    }
    
    // Pre-compute variables in cubeness term
    // compute N, unit area-weighted normal vector
    igl::per_vertex_normals(V,F,data.N);

    // compute a, which is the barycentric area of vertex *i* 
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V,F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    data.a = M.diagonal();
}