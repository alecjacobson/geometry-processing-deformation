#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // REPLACE WITH YOUR CODE
  D = Eigen::MatrixXd::Zero(data.n,3);



  // Sparse type not working with igl.
  //Eigen::SparseMatrix<double> B, Beq;
  //B.resize(data.n, 3);
  //Beq.resize(data.n, 3);

  // Sparse isn't working. So use a usual matrix.
  Eigen::MatrixXd B, Beq;

  // B needs to be initialized for igl to work without assertion fail.
  B = Eigen::MatrixXd::Zero(data.n, 3);

  
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);
}



