#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  /*
    min_x ∫ ‖∆d‖² dA ≈ tr(DᵀLᵀM⁻ᵀMM⁻¹LD) = tr(DᵀLᵀM⁻¹LD) = tr(DᵀQD)
    
    For which Q has already been prefactorized for our problem, and so we just need to perform the solve with our given handle locations

    The bc specifies the (x,y,z) location of the point which is fixed, and due to the minimization of the trace, solve each X,Y,Z indep:
   */

  D = Eigen::MatrixXd::Zero(data.n, 3);
  
  // size of the system is data.n
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(data.n, 1);
  Eigen::MatrixXd Beq = Eigen::MatrixXd::Zero(0,0);

  // solve for the (X, Y, Z)'s separately
  for(int32_t i = 0; i < 3; i++)
  {
    Eigen::VectorXd dcol;
    if(!igl::min_quad_with_fixed_solve(data, B, bc.col(i), Beq, dcol))
      printf("ERROR mqwf compute\n");
    D.col(i) = dcol;
  }
}

