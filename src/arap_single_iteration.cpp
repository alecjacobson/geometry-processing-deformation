#include "arap_single_iteration.h"

#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

#include <ostream>

void arap_single_iteration(
    const igl::min_quad_with_fixed_data<double> & data,
    const Eigen::SparseMatrix<double> & K,
    const Eigen::MatrixXd & bc,
    Eigen::MatrixXd & U)
{
    const int n = data.n;

    // local - solve for the rotations...
    Eigen::MatrixXd R( 3*n, 3 );
    R = Eigen::MatrixXd::Zero( 3*n, 3 );
    Eigen::MatrixXd C( 3*n, 3 );
    C = (U.transpose() * K).transpose();
    //std::cout<<"C is: "<<C.rows()<<"x"<<C.cols()<<std::endl;
    
    // solve min tr( Ct R ), walk each vertexes (k) rotation Rk (3,3)
    // and weighted covariance Ck (3,3)
    for (int k = 0; k < n; ++k)
    {
        //std::cout<<"k is:"<<k<<std::endl;
        Eigen::Matrix3d Rk = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d Ck = C.block( 3*k, 0, 3, 3 );
        //std::cout<<"Ck is: "<<Ck.rows()<<"x"<<Ck.cols()<<std::endl<<Ck<<std::endl;

        
        igl::polar_svd3x3( Ck, Rk );
        //std::cout<<"Rk after svd is: "<<Rk.rows()<<"x"<<Rk.cols()<<std::endl<<Rk<<std::endl;
        //std::cout<<"attempting to write to "<<3*k<<","<<0<<std::endl;
        //std::cout<<"R size is: "<<R.rows()<<"x"<<R.cols()<<std::endl;
        //std::cout<<"That part of R is:"<<std::endl<<R.block( 3*k, 0, 3, 3 )<<std::endl;
        //std::cout<<"proceeding..."<<std::endl;
        // write into R - I love blocks!
        //Rk = Eigen::Matrix3d::Identity();
        R.block( 3*k, 0, 3, 3 ) = Rk;
        //std::cout<<"Wrote to R..."<<std::endl;
    }

    //std::cout<<"Done with SVD"<<std::endl;

    //std::cout<<"K is: "<<K.rows()<<"x"<<K.cols()<<std::endl;
    //std::cout<<"R is: "<<R.rows()<<"x"<<R.cols()<<std::endl;
    //std::cout<<"U is: "<<U.rows()<<"x"<<U.cols()<<std::endl;    

    // global - use the rotations and solve the quadratic
    Eigen::MatrixXd B = K * R;
    //std::cout<<"B is: "<<B.rows()<<"x"<<B.cols()<<std::endl;
    
    Eigen::VectorXd Beq; // empty equality constraints
    //std::cout<<"Beq is: "<<Beq.rows()<<"x"<<Beq.cols()<<std::endl;
    
    igl::min_quad_with_fixed_solve( data, B, bc, Beq, U );}
