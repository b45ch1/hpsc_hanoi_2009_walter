#ifndef DENSE_QR_DECOMPOSITION_BY_GIVENS_ROTATIONS_HPP
#define DENSE_QR_DECOMPOSITION_BY_GIVENS_ROTATIONS_HPP

#include <iostream>
#include <cstdlib>
#include <vector>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>

/** \brief Computes the Q^T R decomposition
Warning: this algorithm returns Q^T and not Q.
*/
template <class Tdouble>
int qr(const mtl::dense2D<Tdouble> &in_A, mtl::dense2D<Tdouble> &QT, mtl::dense2D<Tdouble> &R){
	/* check inputs */
	const int N = in_A.num_rows();
	const int M = in_A.num_cols();
	
	if(N != M){
		return -1;
	}
	
	/* prepare Q and R */
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			QT[n][m] = n==m;
		}
	}
	
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			R[n][m] = in_A[n][m];
		}
	}

	/* MAIN ALGORITHM */
	for(int n = 0; n!=N; ++n){
		for(int m = n+1; m!=N; ++m){
			/* defining coefficients of the Givens rotation */
			const Tdouble a = R[n][n];
			const Tdouble b = R[m][n];
			const Tdouble r = sqrt(a*a + b*b);
			const Tdouble c = a/r;
			const Tdouble s = b/r;
			
			for(int k = 0; k!=N; ++k){
				/* update R */
				const Tdouble Rnk = R[n][k];
				R[n][k] = c*Rnk + s*R[m][k];
				R[m][k] =-s*Rnk + c*R[m][k];
				/* update Q */
				const Tdouble QTnk = QT[n][k];
				QT[n][k] = c*QTnk + s*QT[m][k];
				QT[m][k] =-s*QTnk + c*QT[m][k];
			}
		}
	}
	return 0;
} 
#endif
