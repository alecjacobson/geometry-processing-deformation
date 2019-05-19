#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/polar_dec.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>
#include <chrono>

#include "utils.h"
#include "fit_rotation.h"
#include "fit_rotation_sse.h"
#include "fit_rotation_avx.h"

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
  using micro = std::chrono::microseconds;
  MatrixXd C = K.transpose().eval() * U;
  MatrixXd R(3 * data.n, 3); // a stack of rotation 3x3 matrix for each vertex
  Matrix3d Ck, Rk, Rk2, T, R_lastk;
  int iters = 0;
  Matrix3d Q = MatrixXd::Identity(3, 3);



  auto start1 = std::chrono::high_resolution_clock::now();
  //   new version
  for (int i = 0; i < num_of_group-1; i++) {
      for (int j = 0; j < 9; j++) {
          Mf[i*72+j*8+0] = (float)C(i*24+0*3+j/3, j%3);
          Mf[i*72+j*8+1] = (float)C(i*24+1*3+j/3, j%3);
          Mf[i*72+j*8+2] = (float)C(i*24+2*3+j/3, j%3);
          Mf[i*72+j*8+3] = (float)C(i*24+3*3+j/3, j%3);
          Mf[i*72+j*8+4] = (float)C(i*24+4*3+j/3, j%3);
          Mf[i*72+j*8+5] = (float)C(i*24+5*3+j/3, j%3);
          Mf[i*72+j*8+6] = (float)C(i*24+6*3+j/3, j%3);
          Mf[i*72+j*8+7] = (float)C(i*24+7*3+j/3, j%3);
      }
  }


//   auto start = std::chrono::high_resolution_clock::now();
  fit_rotation_small_avx(Mf, Rf, num_of_group-1);
//   auto finish = std::chrono::high_resolution_clock::now();
//   std::cout << "caylay avx took "
//       << std::chrono::duration_cast<nano>(finish - start).count()
//       << std::endl;


  // new version
  for (int i = 0; i < num_of_group-1; i++) {
      for (int j = 0; j < 9; j++) {
          R(i*24+0*3+j/3, j%3) = (double)Rf[i*72+j*8+0];
          R(i*24+1*3+j/3, j%3) = (double)Rf[i*72+j*8+1];
          R(i*24+2*3+j/3, j%3) = (double)Rf[i*72+j*8+2];
          R(i*24+3*3+j/3, j%3) = (double)Rf[i*72+j*8+3];
          R(i*24+4*3+j/3, j%3) = (double)Rf[i*72+j*8+4];
          R(i*24+5*3+j/3, j%3) = (double)Rf[i*72+j*8+5];
          R(i*24+6*3+j/3, j%3) = (double)Rf[i*72+j*8+6];
          R(i*24+7*3+j/3, j%3) = (double)Rf[i*72+j*8+7];
      }
  }
  auto finish1 = std::chrono::high_resolution_clock::now();
    // std::cout << "fit rotation overall: "
    // << std::chrono::duration_cast<micro>(finish1 - start1).count()
    // << " microseconds\n";


  auto local_time = std::chrono::duration_cast<micro>(finish1 - start1).count();






//   auto start1 = std::chrono::high_resolution_clock::now();
//   // old version
//   for (int k = 0; k < data.n; k++) {

//       Ck = C.block(3*k, 0, 3, 3);

//       // fit rotation
//       R_lastk = R_last.block(3*k, 0, 3, 3);

//       Ck = R_lastk.transpose().eval() * Ck;

//     //   iters = fit_rotation_small_no_sse(Ck, Rk);
//       igl::polar_svd3x3(Ck, Rk);
//     //   stream << rotationMatrixToAxisAngle(Rk) << std::endl;

//       R_lastk = R_lastk * Rk;

//       R.block(3 * k, 0, 3, 3) = R_lastk.eval();
//       R_last.block(3 * k, 0, 3, 3) = R_lastk.eval();

//   }
//   auto finish1 = std::chrono::high_resolution_clock::now();
//     std::cout << "sifakis simd overall: "
//     << std::chrono::duration_cast<nano>(finish1 - start1).count()
//     << " nanoseconds\n";



  // global step
  auto start2 = std::chrono::high_resolution_clock::now();
  // solve for new U
  igl::min_quad_with_fixed_solve(data, K * R, bc, MatrixXd(), U);
  auto finish2 = std::chrono::high_resolution_clock::now();
    // std::cout << "min quad: "
    // << std::chrono::duration_cast<micro>(finish2 - start2).count()
    // << " microseconds\n";

  auto global_time = std::chrono::duration_cast<micro>(finish2 - start2).count();




  return Vector3d(local_time, global_time, 0.);

}




// iters = fit_rotation(Ck, 2.2, false, true, Rk, Q);
// iters = fit_rotation(Ck, 2.2, false, true, Rk, R_lastK_T);


// if (k % 100 == 0) {
//   stream << Ck << std::endl;
// }


// // svd
// auto start = std::chrono::high_resolution_clock::now();
// igl::polar_svd3x3(Ck, Rk);
// auto finish = std::chrono::high_resolution_clock::now();
// std::cout << "polar_svd3x3() took "
//     << std::chrono::duration_cast<nano>(finish - start).count()
//     << " nanoseconds\n";
// R.block(3 * k, 0, 3, 3) = Rk.eval();



//  // new version
//   for (int i = 0; i < num_of_group; i++) {
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


//   auto start = std::chrono::high_resolution_clock::now();
//   fit_rotation_small_no_avx(Mf, Rf, num_of_group);
//   auto finish = std::chrono::high_resolution_clock::now();
//   std::cout << "fit_rotation() took "
//       << std::chrono::duration_cast<nano>(finish - start).count()/data.n
//       << " nanoseconds\n";

//   // new version
//   for (int i = 0; i < num_of_group; i++) {
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

//   stream << R.block(0, 0, 3, 3) << std::endl;