#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
    Eigen::SparseMatrix<double> L, M, Aeq;
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(M);
    // Q = (M.inverse()*L).transpose() * L
    Eigen::SparseMatrix<double> Q = solver.solve(L).transpose() * L;
    //This line of code is adapted from the libigl tutorial
    igl::min_quad_with_fixed_precompute(Q,b,Aeq,false,data);
}

