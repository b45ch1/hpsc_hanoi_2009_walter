#ifndef DENSE_INVERSE_HPP
#define DENSE_INVERSE_HPP

#include <iostream>
#include <cstdlib>
#include <vector>
#include "dense_qr_decomposition_by_givens_rotations.hpp"

/** \brief Compute the inverse of a matrix with Givens rotations
*/
template <class Tdouble>
void inv(Tdouble *a, Tdouble *qt, Tdouble *r, int na){
	int n,m,k;
	Tdouble rnn, rnm;
	qr(a,qt,r,na);
	
	for(n = na-1; n>=0; --n){
		rnn = *(r +myindex(n,n,na)) ;
		for(m = 0; m!=na; ++m){
			*(r +myindex(n,m,na))   = *(r  + myindex(n,m,na)) / rnn;
			*(qt + myindex(n,m,na)) = *(qt + myindex(n,m,na)) / rnn;
		}
		for(m = n+1; m<na; ++m){
			rnm = *(r + myindex(n,m,na));
			*(r +  myindex(n,m,na)) = 0;
			for(k = 0; k != na; ++k){
				*(qt + myindex(n,k,na))  = *(qt + myindex(n,k,na))  - *(qt + myindex(m,k,na))*rnm;
			}
		}
	}
}

#endif

