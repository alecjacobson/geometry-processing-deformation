#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/polar_dec.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>

#include "fit_rotation.h"

using namespace Eigen;

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  std::ostream & stream,
  Eigen::MatrixXd & U)
{
  // two level optimization: local and global
  // construct C first
  MatrixXd C = K.transpose() * U;
  // construct R: a stack of rotation matrices
  MatrixXd R(3 * data.n, 3); // a stack of rotation 3x3 matrix for each vertex
  Matrix3d Ck, Rk, Rk2, T;
  for (int k = 0; k < data.n; k++) {
      Ck = C.block(3*k, 0, 3, 3);
      // stream << Ck.determinant() << std::endl;

      // // svd
      // igl::polar_svd3x3(Ck, Rk2);

      // igl::polar_dec(Ck, Rk, T);
      // // Rk = Rk.transpose().eval();

      // double diff = abs((Ck*Rk).eval().trace() - (Ck*Rk2).eval().trace());

      // if (diff > 3e-2) {
      //     std::cout << (Ck*Rk).eval().trace() << std::endl;
      //     std::cout << (Ck*Rk2.transpose()).eval().trace() << std::endl;
      // }
      
      // fit rotation
      fit_rotation(Ck, 2.2, false, false, Rk, Rk);
      Rk = Rk.transpose().eval();

      R.block(3 * k, 0, 3, 3) = Rk;
      
  }
  // solve for new U
  igl::min_quad_with_fixed_solve(data, K * R, bc, MatrixXd(), U);

}
