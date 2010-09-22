#include <stdio.h>
#include <stdlib.h>

#include "adolc/adolc.h"
#include "dense_inverse.hpp"
#include "dot.hpp"

extern "C"{
#include "clapack.h"
#include "cblas.h"
#include "dense_inverse.h"
#include "inv_b.h"
#include "utpm.h"

}

#include <sys/time.h>
struct timeval tv;
int mtime(void)
{
gettimeofday(&tv,NULL);
return (int)(tv.tv_sec*1000 + (tv.tv_usec / 1000));
}

using namespace std;
/** computes y = trace(A) */
template <class Tdouble>
int trace(Tdouble *A, Tdouble *Ainv, Tdouble *R, Tdouble *y, int N){
    *y = 0;
    inv(A, Ainv, R, N);
    for(int n = 0; n < N; ++n){
        (*y) += Ainv[id(n,n,N)];
    }
}

extern "C"
int function_inv_c(int N, double *runtimes){
    /*
	TIMING MY C-IMPLEMENTATION INVERSE COMPUTATION 
	*/
	int info;
	double start_time, end_time;
	
	double *A = new double[N*N];
	double *Ainv = new double[N*N];
	double *R = new double[N*N];
	double *Id = new double[N*N];
	
	/* fill A with random numbers */
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			A[id(n,m,N)] = rand()/100000000000.;
		}
	}
	
	start_time = mtime();
	inv(A, Ainv, R, N);
	end_time = mtime();
	runtimes[0] = end_time-start_time;
	
	/* check that dot(A,Ainv) = I */
	cblas_dgemm(CblasRowMajor, CblasNoTrans,
				 CblasNoTrans, N, N,
				N, 1., A,
				 N, Ainv, N,
				 0., Id, N);
	{
		double sum = 0;
		for (int n = 0; n != N; ++n){
			for (int m = 0; m != N; ++m){
				sum += Id[id(n,m,N)];
			}
		}
		info = (bool) (abs(sum/N - 1) < N*1e-8 );
	}
	
	delete[] A;
	delete[] Ainv;
	delete[] R;
	delete[] Id;
	
	return 0;
}



extern "C"
int function_inv_adolc(int N, double *runtimes){
    /*
	TIMING TAPING OF THE C++-IMPLEMENTATION INVERSE COMPUTATION 
	*/
	int info;
	double start_time, end_time;
	double *A = new double[N*N];
	double *Ainv = new double[N*N];
	double *R = new double[N*N];
	double *Id = new double[N*N];
	adouble *aA = new adouble[N*N];
	adouble *aR = new adouble[N*N];
	adouble *aAinv = new adouble[N*N];
	
	/* fill A with random numbers */
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			A[id(n,m,N)] = rand()/100000000000.;
		}
	}	
	
	
	start_time = mtime();
	trace_on(0);
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			aA[id(n,m,N)]<<= A[id(n,m,N)];
		}
	}
	inv(aA, aAinv, aR, N);
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			aAinv[id(n,m,N)]>>= Ainv[id(n,m,N)];
		}
	}
	trace_off();
	end_time = mtime();
	runtimes[0] = end_time-start_time;
	
	start_time = mtime();
	function(0,N*N,N*N,A, Ainv);
	end_time = mtime();
	runtimes[1] = end_time-start_time;
	
	/* check that dot(A,Ainv) = I */
	cblas_dgemm(CblasRowMajor, CblasNoTrans,
				 CblasNoTrans, N, N,
				N, 1., A,
				 N, Ainv, N,
				 0., Id, N);
	{
		double sum = 0;
		for (int n = 0; n != N; ++n){
			for (int m = 0; m != N; ++m){
				sum += Id[id(n,m,N)];
			}
		}
		info = (bool) (abs(sum/N - 1) < N*1e-8 );
	}
	
	delete[] aA;
	delete[] aR;
	delete[] aAinv;
	delete[] A;
	delete[] Ainv;
	delete[] R;
	delete[] Id;
	
	
	return 0;
}



extern "C"
int gradient_trace_adolc(int N, double *runtimes){
    /*
	compute gradient of trace(inv(A)) with ADOL-C
	*/
	double start_time, end_time;
	double *A = new double[N*N];
	double *Ainv = new double[N*N];
	double *R = new double[N*N];
	double *Id = new double[N*N];
	adouble *aA = new adouble[N*N];
	adouble *aR = new adouble[N*N];
	adouble *aAinv = new adouble[N*N];
	
	/* fill A with random numbers */
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			A[id(n,m,N)] = rand()/100000000000.;
		}
	}

	start_time = mtime();
	double y;
	trace_on(0);
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			aA[id(n,m,N)]<<= A[id(n,m,N)];
		}
	}
	adouble ay;
	trace(aA, aAinv, aR, &ay, N);
	ay >>= y;
	trace_off();
	end_time = mtime();
	runtimes[0] = end_time-start_time;

	double *g;
	g = myalloc(N*N);
	start_time = mtime();
	gradient(0,N*N,A, g);
	end_time = mtime();
	runtimes[1] = end_time-start_time;
	
	
	delete[] aA;
	delete[] aR;
	delete[] aAinv;
	delete[] A;
	delete[] Ainv;
	delete[] R;
	delete[] Id;
	
	
	return -1;
}


 
