#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/arap.h>
#include <igl/arap_linear_block.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/cotmatrix_entries.h>

using namespace Eigen;
using namespace igl;


/*

***************FIRST ATTEMPT******************
//following implementation details in Alec's paper - https://igl.ethz.ch/projects/fast/fast-automatic-skinning-transformations-siggraph-2012-jacobson-et-al.pdf

struct DirectEdge 
{
	int v1, v2;
	double cotW;

	double getCoef(int v) 
	{
		if (v1 == v) 
			return 1;
		else if (v2 == v)
			return -1;
		else
			return 0;
	}
};


void build_adjacency(std::vector<int>& facering, const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, MatrixXd& C, MatrixXd& result)
{	
	
	int nVert = V.rows();
	const int nEdges = 3 * facering.size();
	
	std::vector<DirectEdge> edges;
	//       for triangles, columns correspond to edges [1,2],[2,0],[0,1]
	
	// Find all halfedges for ridges and spokes of the ring around each vertex
	for (int m = 0; m < facering.size(); ++m)
	{
		int fj = facering[m];		

		int i0 = F(fj, 0), i1 = F(fj, 1), i2 = F(fj, 2);
		
		double cot12 = C(fj, 0), cot20 = C(fj, 1), cot01 = C(fj, 2);
				
		edges.push_back({ i1, i2, cot12 });
		edges.push_back({ i2, i0, cot20 });
		edges.push_back({ i0, i1, cot01 });
	}
	
	std::vector< Eigen::Triplet<double> > tripletList;
	tripletList.reserve(nEdges);

	for (int m = 0; m < facering.size(); ++m)
	{
		int fj = facering[m];

		int i0 = F(fj, 0),  i1 = F(fj, 1), i2 = F(fj, 2);

		for (int k = 0; k < nEdges; ++k)
		{
	     	tripletList.push_back(Triplet<double>(i0, k, edges[k].getCoef(i0)));
			tripletList.push_back(Triplet<double>(i1, k, edges[k].getCoef(i1)));
			tripletList.push_back(Triplet<double>(i2, k, edges[k].getCoef(i2)));
		}	
	}
	
	VectorXd Ctemp(nEdges);
	for (int k = 0; k < nEdges; ++k)
	{
		Ctemp(k) = edges[k].cotW;
	}

	SparseMatrix<double> A;
	A.resize(nVert, nEdges);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	
	result = V.transpose()*A*Ctemp.asDiagonal()*A.transpose();	
}

void arap_precompute(
const Eigen::MatrixXd & V,
const Eigen::MatrixXi & F,
const Eigen::VectorXi & b,
igl::min_quad_with_fixed_data<double> & data,
Eigen::SparseMatrix<double> & K)
{
	int nVert = V.rows();

	SparseMatrix<double> L;
	cotmatrix(V, F, L);

	SparseMatrix<double> Aeq;
	bool result = min_quad_with_fixed_precompute(L, b, Aeq, false, data);
	assert(result);

	//std::vector<std::vector<int>> M;
	//adjacency_list(F, M, true);
	std::vector<std::vector<int>> VF;
	std::vector<std::vector<int>> VFi;
	vertex_triangle_adjacency(V.rows(), F, VF, VFi);

	//     C  #F by 3 list of 1/2*cotangents corresponding angles
	//       for triangles, columns correspond to edges [1,2],[2,0],[0,1]

	MatrixXd C;
	cotmatrix_entries(V, F, C);

	MatrixXd Kfinal;
	Kfinal.resize(3 * nVert, nVert);

	for (int fi = 0; fi < VF.size(); ++fi)
	{
	MatrixXd result;
	build_adjacency(VF[fi], V, F, C, result);
	assert(result.rows() == 3 && result.cols() == nVert);
	Kfinal.block(3 * fi, 0, 3, nVert) = result;
	}

	K = Kfinal.transpose().sparseView();
}
*/


void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
	int nVert = V.rows();

	SparseMatrix<double> L;
	cotmatrix(V, F, L);

	SparseMatrix<double> Aeq;
	L = 2*L; //because of how minimization is formulated
	bool result = min_quad_with_fixed_precompute(L, b, Aeq, false, data);
	assert(result);
	
	MatrixXd C;
	cotmatrix_entries(V, F, C);
	
	std::vector< Eigen::Triplet<double>> tripletList;
	
	int nFaces = F.rows();
	
	for (int i = 0; i < nFaces; ++i) 
	{		
		for (int k = 0; k < 3; ++k) // since every edge will apear 3 times for every vertex of the face
		{
			int r = F(i, k);

			for (int v = 0; v < 3; ++v)
			{				
				int vi = F(i, (v + 1) % 3); //Since cotmatrix_entries gives edges [1,2],[2,0],[0,1]	
				int vj = F(i, (v + 2) % 3);
				double c_ij = C(i, v);
				VectorXd e_ij = c_ij * (V.row(vi) - V.row(vj));

				for (int d = 0; d < 3; ++d)
				{
					int w = 3 * r + d;
					tripletList.push_back(Triplet<double>(vi, w, e_ij(d)));
					tripletList.push_back(Triplet<double>(vj, w, -e_ij(d)));
				}
			}
		}
	}
	
	K.resize(nVert, 3 * nVert);
	K.setFromTriplets(tripletList.begin(), tripletList.end());
}
