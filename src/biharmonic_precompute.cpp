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
    // REPLACE WITH YOUR CODE
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
    else
        std::cout<<"SimplicialLDLT prefactor succeeded!"<<std::endl;

    Eigen::SparseMatrix<double> I( n, n );
    I.setIdentity();
    auto Minverse = solver.solve( I );
    if( solver.info() != Eigen::Success )
        std::cout<<"SimplicialLDLT failed to SOLVE M, might crash..."<<std::endl;
    else
        std::cout<<"SimplicialLDLT solve succeeded!"<<std::endl;
    

    std::cout<<"M is: "<<M.rows()<<"x"<<M.cols()<<std::endl;
    std::cout<<"Minverse is: "<<Minverse.rows()<<"x"<<Minverse.cols()<<std::endl;
    std::cout<<"L is: "<<L.rows()<<"x"<<L.cols()<<std::endl;
    
    
    Eigen::SparseMatrix<double> Qpartial( n, n );
    Qpartial = L.transpose() * Minverse;
    Eigen::SparseMatrix<double> Q( n, n );
    Q = Qpartial * L;
    std::cout<<"Q is: "<<Q.rows()<<"x"<<Q.cols()<<std::endl;
    std::cout<<"b is: "<<b<<std::endl;

    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute( Q, b, Aeq, false, data );

    
    // std::cout<<"Minverse is:"<<Minverse<<std::endl;
    //std::cout<<"Product is: " << (L.transpose() * Minverse )<< std::endl;
    // std::cout<<"L.transpose() is:"<<L.transpose()<<std::endl;
    // std::cout<<"L is:"<<L<<std::endl;
    
    // //igl::min_quad_with_fixed( )
    // Eigen::VectorXd B = Eigen::VectorXd::Zero( V.rows(), 1 ); // linear term
    // Eigen::VectorXd Beq; // empty constraints
    // Eigen::SparseMatrix<double> Aeq;
    // igl::min_quad_with_fixed_precompute( (-L).eval(), b, Aeq, true, data );

    // Eigen::VectorXd U1 = U.col(0), U2 = U.col(1);
    // igl::min_quad_with_fixed_solve( mqwf, B, boundaryUV.col(0), Beq, U1 );
    // igl::min_quad_with_fixed_solve( mqwf, B, boundaryUV.col(1), Beq, U2 );

    // U.col(0) = U1;
    // U.col(1) = U2;

}
