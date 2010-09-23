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
#include "utps.h"


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
	this is equivalent to Univariate Taylor Polynomial Arithmetic of Scalars
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

extern "C"
int gradient_trace_c(int N, double *runtimes){
    /*
    compute gradient of trace(inv(A)) with a C implementation
    this is equivalent to Univariate Taylor Polynomial Arithmetic of Matrices
      
    
    need to compute
    d tr(A^-1) = tr( - (A^-1)^2 dA) and therefore
    \nabla tr(A^-1) = - ((A^-1)^2)^T
    */
    
    double start_time, end_time;
    double *A = new double[N*N];
    double *Ainv = new double[N*N];
    double *R = new double[N*N];
    double *Abar = new double[N*N];
    
    /* fill A with random numbers */
    for(int n = 0; n!=N; ++n){
        for(int m = 0; m!=N; ++m){
            A[id(n,m,N)] = rand()/100000000000.;
        //     Ainv[id(n,m,N)] =  A[id(n,m,N)];
        }
    }
    start_time = mtime();
    
    /* compute A^{-1} */
    inv(A, Ainv, R, N);
     /* compute Abar = - (A^{-1} A^{-1})^T */
    cblas_dgemm(CblasRowMajor, CblasTrans,
                 CblasTrans, N, N,
                N, -1., Ainv,
                 N, Ainv, N,
                 0., Abar, N);
    end_time = mtime();
    runtimes[0] = end_time-start_time;
    
    delete[] A;
    delete[] Ainv;
    delete[] R;
    delete[] Abar;
    
    return -1;
}


extern "C"
int utp_dot_adolc(int P, int D, int N, double *runtimes){
    /* UTP computation of the matrix matrix product using ADOL-C in the forward mode */
    double start_time, end_time;
    adouble *aA = new adouble[(N*N)];
    adouble *aB = new adouble[(N*N)];
    adouble *aC = new adouble[(N*N)];
    double *x = new double[2*N*N];
    double *y = new double[N*N];
    double ***V = myalloc3(2*N*N,P, D-1);
    double ***W = myalloc3(N*N,P, D-1);

    start_time = mtime();
    trace_on(2);
    for(int n = 0; n!=N; ++n){
        for(int m = 0; m!=N; ++m){
            aA[id(n,m,N)]<<= rand()/100000000000.;
            aB[id(n,m,N)]<<= rand()/100000000000.;
        }
    }
    dot(aA, aB, aC, N);
    for(int n = 0; n!=N; ++n){
        for(int m = 0; m!=N; ++m){
            aC[id(n,m,N)]>>= *y;
        }
    }
    trace_off();
    end_time = mtime();
    runtimes[0] = end_time-start_time;
    
    start_time = mtime();
    hov_forward(2,N*N,2*N*N,D-1, P, x, V, y, W);
    end_time = mtime();
    runtimes[1] = end_time-start_time;
    
    delete[] x;
    delete[] y;
    delete[] aA;
    delete[] aB;
    delete[] aC;
    
    return -1;
}

extern "C"
int utp_dot_taylorpoly(int P, int D, int N, double *runtimes){
    /* UTP computation of the matrix matrix product using TAYLORPOLY in the forward mode */
    
    double start_time, end_time;
    double *utpm_A = new double[(1+(D-1)*P)*(N*N)];
    double *utpm_B = new double[(1+(D-1)*P)*(N*N)];
    double *utpm_C = new double[(1+(D-1)*P)*(N*N)];
    int *utpm_strides = new int[3];
    utpm_strides[0] = N*N*sizeof(double);
    utpm_strides[1] = sizeof(double);
    utpm_strides[2] = N*sizeof(double);
    int *utpm_shape = new int[2];
    utpm_shape[0] = N; utpm_shape[1] = N;
    
    /* fill utpm_A and utpm_B with random numbers */
    for(int i = 0; i!=N*N + (D-1)*P*(N*N); ++i){
        utpm_A[i] = rand()/100000000000.;
        utpm_B[i] = rand()/100000000000.;
    }
    
    // print_utpm(P, D, 2, utpm_shape, &utpm_strides[1], utpm_A);
    start_time = mtime();
    utpm_dot(P, D, N, N, N, 1., utpm_A, utpm_strides, utpm_B, utpm_strides, 0., utpm_C, utpm_strides);
    end_time = mtime();
    runtimes[0] = end_time-start_time;
    
    return -1;
}

extern "C"
int utps_dot_taylorpoly(int P, int D, int N, double *runtimes){
    int stride = (1+P*(D-1));
    double start_time, end_time;
    double *x = new double[(1+(D-1)*P)*(N*N)];
    double *y = new double[(1+(D-1)*P)*(N*N)];
    double *z = new double[(1+(D-1)*P)*(N*N)];

    start_time = mtime();
    for(int n = 0; n != N; ++n){ // iterate over all columns
       for(int m = 0; m != N; ++m){ // iterate over all rows
            for(int k = 0; k != N; ++k){
                utps_amul(P, D, x+fid(m,k,N)*stride, y+fid(k,n,N)*stride, z+fid(m,n,N)*stride);
            }
        }
    }
    end_time = mtime();
    runtimes[0] = end_time-start_time;
    return -1;
    
}
    



 
