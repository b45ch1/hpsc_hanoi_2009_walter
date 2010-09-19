#include "math.h"
#include "dense_inverse.h"

/** \brief Computes the Q^T R decomposition
Warning: this algorithm returns Q^T and not Q.
A na x na array
Q na x na array
r na x na array
Computes:
				A = Q R
				========
*/

inline int myindex(int n, int m, int na){
	return n*na+m;
}

void qr(double *a, double *qt, double *r, int na){
    double at,bt,rt,ct,st, rnk, qtnk;
	int n,m, k, tmp; 

	/* prepare Q and r */
	for(n = 0; n!=na; ++n){
		for(m = 0; m!=na; ++m){
			if(n==m) tmp = 1;
			else     tmp = 0;
			qt[myindex(n,m,na)] = tmp;
		}
	}
	
	for(n = 0; n!=na; ++n){
		for(m = 0; m!=na; ++m){
			r[myindex(n,m,na)] = a[myindex(n,m,na)];
		}
	}

	/* MAIN ALGORITHM */
	for(n = 0; n!=na; ++n){
		for(m = n+1; m!=na; ++m){
			/* defining coefficients of the Givens rotation */
			at = r[myindex(n,n,na)];
			bt = r[myindex(m,n,na)];
			rt = sqrt(at*at + bt*bt);
			ct = at/rt;
			st = bt/rt;
			
			for(k = 0; k!=na; ++k){
				/* update r */
				rnk = r[myindex(n,k,na)];
				r[myindex(n,k,na)] = ct*rnk   + st* (r[myindex(m,k,na)]);
				r[myindex(m,k,na)] =-st*rnk   + ct* (r[myindex(m,k,na)]);
				/* update Q */
				qtnk = qt[myindex(n,k,na)];
				qt[myindex(n,k,na)] = ct*qtnk + st*(qt[myindex(m,k,na)]);
				qt[myindex(m,k,na)] =-st*qtnk + ct*(qt[myindex(m,k,na)]);
			}
		}
	}
}

void inv(double *a, double *qt, double *r, int na){
	int n,m,k;
	double rnn, rnm;
	qr(a,qt,r,na);
	
	for(n = na-1; n>=0; --n){
		rnn = r[myindex(n,n,na)] ;
		for(m = 0; m!=na; ++m){
			r[myindex(n,m,na)]   = r[myindex(n,m,na)] / rnn;
			qt[myindex(n,m,na)] = qt[myindex(n,m,na)] / rnn;
		}
		for(m = n+1; m<na; ++m){
			rnm = r[myindex(n,m,na)];
			r[ myindex(n,m,na)] = 0;
			for(k = 0; k != na; ++k){
				qt[myindex(n,k,na)]  = qt[myindex(n,k,na)]  - qt[myindex(m,k,na)]*rnm;
			}
		}
	}
}
