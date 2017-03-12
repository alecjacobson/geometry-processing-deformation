#ifndef CLOSEST_ROTATION_H
#define CLOSEST_ROTATION_H
#include <Eigen/Core>
// Given a 3×3 matrix `M`, find the closest rotation matrix `R`.
//
// Inputs:
//   M  3x3 matrix
// Outputs:
//   R  3x3 rotation matrix
//
void closest_rotation(
	const Eigen::Matrix3d & M,
	Eigen::Matrix3d & R);
#endif