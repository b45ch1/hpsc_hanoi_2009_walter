#ifndef DENSE_INVERSE_HPP
#define DENSE_INVERSE_HPP

#include <iostream>
#include <cstdlib>
#include <vector>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>

#include "dense_qr_decomposition_by_givens_rotations.hpp"

/** \brief Compute the inverse of a matrix with givens rotations
*/
template <class Tdouble> int inv(const mtl::dense2D<Tdouble> &A, mtl::dense2D<Tdouble> &QT){
	const int N = A.num_rows();
	mtl::dense2D<Tdouble> R(N,N);
	qr(A,QT,R);


	/*
	Solve now the extended linear system
			(Q R | I) = ( R | QT )
			i.e.
			/R_11 R_12 R_13 ... R_1M | 1 0 0 0 ... 0 \
			| 0   R_22 R_23 ... R_2M | 0 1 0 0 ... 0 |
			| 0         ... ... .... |               |
			\                   R_NM | 0 0 0 0 ... 1 /
		
	*/
	
	for(int n = N-1; n>=0; --n){
		const Tdouble Rnn = R[n][n];
		for(int m = 0; m!=N; ++m){
			R[n][m] = R[n][m] / Rnn;
			QT[n][m] = QT[n][m] / Rnn;
		}
		for(int m = n+1; m<N; ++m){
			const Tdouble Rnm = R[n][m];
			R[n][m] = 0;
			for(int k = 0; k != N; ++k){
				QT[n][k] -= QT[m][k]*Rnm;
			}
		}
	}
	return 0;
}

#endif

