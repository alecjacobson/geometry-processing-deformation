#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/polar_dec.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>
#include <chrono>

#include "fit_rotation.h"

#include "minitrace.h"

using namespace Eigen;


const auto rotationMatrixToAngles = [](Eigen::Matrix3d & R) {
    double sy = sqrt(R(0,0)*R(0,0)+R(1,0)*R(1,0));
    double theta_x, theta_y, theta_z;
    if (sy > 1e-6) {
        theta_x = atan2(R(2, 1), R(2, 2));
        theta_y = atan2(-R(2, 0), sy);
        theta_z = atan2(R(1,0), R(0,0));
    }
    else {
        theta_x = atan2(-R(1, 2), R(1, 1));
        theta_y = atan2(-R(2, 0), sy);
        theta_z = 0;
    }
    return RowVector3d(theta_x, theta_y, theta_z);
};


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
  int iters = 0;
  Matrix3d Q = MatrixXd::Identity(3, 3);
  using nano = std::chrono::nanoseconds;

  for (int k = 0; k < data.n; k++) {
      Ck = C.block(3*k, 0, 3, 3);

      // stream << Ck.determinant() << std::endl;


      // // MTR_BEGIN("C++", "svd");
      // // svd
      // auto start = std::chrono::high_resolution_clock::now();
      // igl::polar_svd3x3(Ck, Rk);
      // auto finish = std::chrono::high_resolution_clock::now();
      // std::cout << "polar_svd3x3() took "
      //     << std::chrono::duration_cast<nano>(finish - start).count()
      //     << " nanoseconds\n";
      // // MTR_END("C++", "svd");



      auto start = std::chrono::high_resolution_clock::now();
      // fit rotation
      // iters = fit_rotation(Ck, 2.2, false, true, Rk, Q);
      igl::polar_dec(Ck, Rk, T);
      auto finish = std::chrono::high_resolution_clock::now();
      std::cout << "polar_dec() took "
          << std::chrono::duration_cast<nano>(finish - start).count()
          << " nanoseconds\n";



      // double diff = abs((Ck*Rk).eval().trace() - (Ck*Rk2).eval().trace());

      // if (diff > 3e-2) {
      //     std::cout << (Ck*Rk).eval().trace() << std::endl;
      //     std::cout << (Ck*Rk2.transpose()).eval().trace() << std::endl;
      // }

  
      // auto start = std::chrono::high_resolution_clock::now();
      // // fit rotation
      // // iters = fit_rotation(Ck, 2.2, false, true, Rk, Q);
      // iters = fit_rotation(Ck, 2.2, false, false, Rk, Rk);
      // auto finish = std::chrono::high_resolution_clock::now();
      // std::cout << "fit_rotation() took "
      //     << std::chrono::duration_cast<nano>(finish - start).count()
      //     << " nanoseconds\n";
      // Rk = Rk.transpose().eval();


      // stream << rotationMatrixToAngles(Rk)*180/M_PI << std::endl;
      // stream << iters << std::endl;

      R.block(3 * k, 0, 3, 3) = Rk;

  }
  // solve for new U
  igl::min_quad_with_fixed_solve(data, K * R, bc, MatrixXd(), U);

}
