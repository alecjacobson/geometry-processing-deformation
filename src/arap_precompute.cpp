#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

typedef Eigen::Triplet<double> T;

void arap_precompute(
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		const Eigen::VectorXi & b,
		igl::min_quad_with_fixed_data<double> & data,
		Eigen::SparseMatrix<double> & K) {
  
	int vertexCount = V.rows();
	int faceCount = F.rows();

	// First we deal with "data". We prefactor transpose(V) * L * V.
	Eigen::SparseMatrix<double> L, Aeq;
	igl::cotmatrix(V, F, L);
	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

	// Now we build K as in the readme.

	Eigen::MatrixXd cotangents;
	igl::cotmatrix_entries(V, F, cotangents);

	std::vector<T> triplets;
	for (int faceIndex = 0; faceIndex < faceCount; faceIndex++) {

		for (int faceVertex = 0; faceVertex < 3; faceVertex++) {

			// Where b, c are in {0, 1, 2} and indexA, indexB, indexC are 
			// the actual vertex indices.
			int b = (faceVertex + 1) % 3;
			int c = (faceVertex + 2) % 3;
			int indexA = F(faceIndex, faceVertex);
			int indexB = F(faceIndex, b);
			int indexC = F(faceIndex, c);

			// This is the cotangent for the angle opposite to the edge.
			double cotangent = cotangents(faceIndex, faceVertex);

			// Wasn't sure whether I should half these because the README
			// says "half-edges". Doing so made the Mesh skinnier, so I
			// did not.
			Eigen::RowVector3d difference = V.row(indexB) - 
					V.row(indexC);
			Eigen::RowVector3d entryVector = difference * cotangent;
			for (int k = 0; k < 3; k++) {
				for (int axis = 0; axis < 3; axis++) {

					int kVertexIndex = F(faceIndex, k);
					triplets.push_back(T(indexB, 3 * kVertexIndex + axis,
							entryVector(axis)));
					triplets.push_back(T(indexC, 3 * kVertexIndex + axis,
							-entryVector(axis)));
				}
			}
		}

		K.resize(vertexCount, 3 * vertexCount);
		K.setFromTriplets(triplets.begin(), triplets.end());

		// The issue in the github suggests doing this. Also, not doing it 
		// causes the mesh to become bigger...
		K /= 3.0;
	}




}
