#include "biharmonic_precompute.h"

#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <iostream>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
    data.n = V.rows();
    const int n = V.rows();

    // create the bi-laplacian matrix, Q
    Eigen::SparseMatrix<double> L( n, n );
    igl::cotmatrix( V, F, L );

    Eigen::SparseMatrix<double> M( n, n );
    igl::massmatrix( V, F, igl::MassMatrixType::MASSMATRIX_TYPE_VORONOI, M );

    // Ax = b -> Mx = I
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double > > solver;
    solver.compute( M );
    if( solver.info() != Eigen::Success )
        std::cout<<"SimplicialLDLT failed to prefactor M, might crash..."<<std::endl;

    Eigen::SparseMatrix<double> I( n, n );
    I.setIdentity();
    auto Minverse = solver.solve( I );
    if( solver.info() != Eigen::Success )
        std::cout<<"SimplicialLDLT failed to SOLVE M, might crash..."<<std::endl;
    
    
    Eigen::SparseMatrix<double> Qpartial( n, n );
    Qpartial = L.transpose() * Minverse;
    Eigen::SparseMatrix<double> Q( n, n );
    Q = Qpartial * L;

    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute( Q, b, Aeq, false, data );
}
