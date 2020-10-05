#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/polar_dec.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>
#include <chrono>


using namespace Eigen;


Vector3d arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  std::ostream & stream,
  Eigen::MatrixXd & R_last,
  Eigen::MatrixXd & U,
  float* Rf,
  float* Mf,
  int num_of_group)
{

  // two level optimization: local and global
  // construct C first
  MatrixXd C = K.transpose().eval() * U;
  MatrixXd R(3 * data.n, 3); // a stack of rotation 3x3 matrix for each vertex
  Matrix3d Ck, Rk, T, R_lastk;
  Matrix3d Q = MatrixXd::Identity(3, 3);



//   auto start1 = std::chrono::high_resolution_clock::now();
//   //   new version
//   for (int i = 0; i < num_of_group-1; i++) {
//       for (int j = 0; j < 9; j++) {
//           Mf[i*72+j*8+0] = (float)C(i*24+0*3+j/3, j%3);
//           Mf[i*72+j*8+1] = (float)C(i*24+1*3+j/3, j%3);
//           Mf[i*72+j*8+2] = (float)C(i*24+2*3+j/3, j%3);
//           Mf[i*72+j*8+3] = (float)C(i*24+3*3+j/3, j%3);
//           Mf[i*72+j*8+4] = (float)C(i*24+4*3+j/3, j%3);
//           Mf[i*72+j*8+5] = (float)C(i*24+5*3+j/3, j%3);
//           Mf[i*72+j*8+6] = (float)C(i*24+6*3+j/3, j%3);
//           Mf[i*72+j*8+7] = (float)C(i*24+7*3+j/3, j%3);
//       }
//   }

// //   auto start = std::chrono::high_resolution_clock::now();
//   fit_rotation_small_avx(Mf, Rf, num_of_group-1);
// //   auto finish = std::chrono::high_resolution_clock::now();
// //   std::cout << "caylay avx took "
// //       << std::chrono::duration_cast<nano>(finish - start).count()
// //       << std::endl;

//   // new version
//   for (int i = 0; i < num_of_group-1; i++) {
//       for (int j = 0; j < 9; j++) {
//           R(i*24+0*3+j/3, j%3) = (double)Rf[i*72+j*8+0];
//           R(i*24+1*3+j/3, j%3) = (double)Rf[i*72+j*8+1];
//           R(i*24+2*3+j/3, j%3) = (double)Rf[i*72+j*8+2];
//           R(i*24+3*3+j/3, j%3) = (double)Rf[i*72+j*8+3];
//           R(i*24+4*3+j/3, j%3) = (double)Rf[i*72+j*8+4];
//           R(i*24+5*3+j/3, j%3) = (double)Rf[i*72+j*8+5];
//           R(i*24+6*3+j/3, j%3) = (double)Rf[i*72+j*8+6];
//           R(i*24+7*3+j/3, j%3) = (double)Rf[i*72+j*8+7];
//       }
//   }
//   auto finish1 = std::chrono::high_resolution_clock::now();
//     // std::cout << "fit rotation overall: "
//     // << std::chrono::duration_cast<micro>(finish1 - start1).count()
//     // << " microseconds\n";
//   auto local_time = std::chrono::duration_cast<micro>(finish1 - start1).count();



  auto start1 = std::chrono::high_resolution_clock::now();
  // old version
  for (int k = 0; k < data.n; k++) {

      Ck = C.block(3*k, 0, 3, 3);

      // fit rotation
      R_lastk = R_last.block(3*k, 0, 3, 3);

      Ck = R_lastk.transpose().eval() * Ck;

      igl::polar_svd3x3(Ck, Rk);

      R_lastk = R_lastk * Rk;
      R.block(3 * k, 0, 3, 3) = R_lastk.eval();
      R_last.block(3 * k, 0, 3, 3) = R_lastk.eval();

  }


  // solve for new U
  igl::min_quad_with_fixed_solve(data, K * R, bc, MatrixXd(), U);

  return Vector3d(0., 0., 0.);

}




