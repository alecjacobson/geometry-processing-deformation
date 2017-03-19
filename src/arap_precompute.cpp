#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/min_quad_with_fixed.h>

typedef Eigen::Triplet<double> tri;

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  /*
    Adding a non-linear rotation into the optimization allows for better behaviour for large displacements:

    min_{x,R} ∫‖∇x - R∇x̑‖² dA = 1/6*∑ ∑cᵢⱼ‖(vᵢ - vⱼ) - R_k(v̑ᵢ - v̑ⱼ)‖²

    where the sums are over all of the vertices (k=1 to n), and all of the faces attached to the vertex (i,jϵF(k))
    The summations can be massaged into a form:

    1/6*∑ ∑cᵢⱼ(vᵢ - vⱼ)ᵀ(vᵢ - vⱼ) + 1/6*∑ ∑cᵢⱼ(vᵢ - vⱼ)ᵀR_k(v̑ᵢ - v̑ⱼ)

    which can be rewritten as:

    tr(VᵀLV) + tr(VᵀKR)

    where L is the cotangent laplacian, Rϵℝ^{3nxn} which consists of all of the induvidual rotation matrices R_k stacked, and 
    Kϵℝ^{nx3n} containing the cotangent weights multiplied with the differences in the rest mesh.

    K is formed by taking all of the entries which would lie within the cotangent matrix (igl::cotmatrix_entries) and their corresponding
    edges, and multiplying them against eachother.

    The K and L matrices can be precomputed, and L prefactorized for efficiency.
   */

  // Precompute and factorize the tr(VᵀLV) solve, where the fixed vertices are given by b:
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  // in order to solve with MQWF, need to modify the arguments slightly as it aims to minimize:
  // 0.5*zᵀAz + zᵀB
  L = 2*L;

  Eigen::SparseMatrix<double> Aeq(0,0);
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

  // Precompute the K matrix:
  std::vector<tri> kVal;
  // we loop through all faces on the mesh and add entries for cᵢⱼ*(v̑ᵢ - v̑ⱼ).
  // cotmatrix_entries returns a list of cotmatrix values where it corresponds to the vertex opposite the edge [1,2], [2,0], [0,1]
  Eigen::MatrixXd cotmat_entries;
  igl::cotmatrix_entries(V, F, cotmat_entries);

  for(int32_t face = 0; face < F.rows(); face++)
  {
    Eigen::RowVector3d c = cotmat_entries.row(face);
    // loop through the vertices on the face:
    for(int32_t vertex = 0; vertex < 3; vertex++)
    {
      // need to line up with the entry in cᵢⱼ, so need order [1,2], [2,0], [0,1]
      int32_t i = F(face, (vertex + 1) % 3);
      int32_t j = F(face, (vertex + 2) % 3);

      Eigen::RowVector3d edgeij = V.row(i) - V.row(j);

      // loop through the indices of the vertices of interest to determine the position to add into the matrix
      for(int32_t k = 0; k < 3; k++)
      {
	int32_t vert_ind = F(face, k);
	// add this value in to the matrix for the xyz's
	for(int32_t dim = 0; dim < 3; dim++)
	{
	  // the stride is 3, due to the rotations being 3x3. The sign changes because of the half edge
	  kVal.push_back(tri(i, 3*vert_ind + dim,  c(vertex)*edgeij(dim)));
	  kVal.push_back(tri(j, 3*vert_ind + dim, -c(vertex)*edgeij(dim)));
	}
      }      
    }
  }
  K.resize(V.rows(), 3*V.rows());
  K.setFromTriplets(kVal.begin(), kVal.end());
}
