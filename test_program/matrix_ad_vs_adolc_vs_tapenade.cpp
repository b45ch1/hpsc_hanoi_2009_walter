#include <iostream>
#include <fstream>

// for memory consumption
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include "adolc/adolc.h"
#include "dense_inverse.hpp"

extern "C"{
#include "clapack.h"
#include "cblas.h"
#include "dense_inverse.h"
#include "inv_b.h"

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

int main(int argc, char* argv[]){
	int N = 100;
	if( argc >= 2){
		int tmp = atoi(argv[1]);
		if(tmp <= 300){
			N = tmp;
		}
	}


	cout<<endl<<endl;
	cout<<"================================================"<<endl;
	cout<<"                  RUNNING TESTS                 "<<endl;
	cout<<"================================================"<<endl;
	cout<<endl<<endl;

	cout<<"Setting N = "<<N<<endl;
	cout<<"PID of this process  = "<< getpid() <<endl;
	char mycmd[255];
	snprintf (mycmd, (size_t)255, "echo [N=%d]>> mem_consumption.txt", N);
	system(mycmd);

	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);


	/* write output to runtimes.txt */
	ofstream runtimes_file;
	runtimes_file.open("runtimes.txt", fstream::in | fstream::out | fstream::app);
	runtimes_file<<N<<"\t";
	
	double *A = new double[N*N];
	double *Ainv = new double[N*N];
	double *Ainv2 = new double[N*N];
	double *QT = new double[N*N];
	double *Q = new double[N*N];
	double *B = new double[N*N];
	double *R = new double[N*N];
	double *Id = new double[N*N];
	double *Abar = new double[N*N];
	double *WORK = new double[N*N];
	

	int start_time,end_time;

	/* fill A with random numbers */
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			A[id(n,m,N)] = rand()/100000000000.;
			Ainv[id(n,m,N)] = A[id(n,m,N)];
		}
	}
	/*
	RUNNING TESTS:
	==============
	
		The tests compare:
			* Lapack VS my implementation to compute VS taped version of my implementation to compute A^-1
				* It turns out that my implementation is a less than a factor 2 slower than LAPACK.
				* The taped function evaluation is much slower though.
			* Compare the perfomance to compute the gradient:
				1) Matrix AD
				2) ADOL-C
				3) TAPENADE
			* Check the obtained results for consistency, i.e. check that there are no bugs in the code!
	*/
	
	cout<<endl<<endl;
	cout<<"================================================"<<endl;
	cout<<"           TESTING FUNCTION EVALUATION          "<<endl;
	cout<<"================================================"<<endl;
	cout<<endl<<endl;
	
	
	/* TIMING LAPACK IMPLEMENTATION INVERSE COMPUTATION */
	/* =============================================== */
	/* compute the inverse by combination of dgetrf and dgetri */
	int *ipiv = new int[N];	int info; int LWORK = N*N;
	start_time = mtime();
	clapack_dgetrf(CblasRowMajor, N, N, Ainv, N, ipiv);
	// cout<<"info="<<info<<endl;
	clapack_dgetri(CblasRowMajor, N, Ainv, N, ipiv);
	// cout<<"info="<<info<<endl;
	end_time = mtime();
	printf("normal LAPACK function evaluation of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";


	/* check that dot(A,Ainv) = I */
	double one = 1;
	double zero = 0;
	char notrans = 'n';
	cblas_dgemm(CblasRowMajor, CblasNoTrans,
				 CblasNoTrans, N, N, N, 1., A,
				 N, Ainv, N,
				0., Id, N);
	{
		double sum = 0;
		for (int n = 0; n != N; ++n){
			for (int m = 0; m != N; ++m){
				sum += Id[id(n,m,N)];
			}
		}
		cout<<sum<<endl;
		cout<< "Computation correct? "<< (bool) (abs(sum/N - 1) < N*1e-8 ) <<endl;
	}
	

	/* TIMING MY C-IMPLEMENTATION INVERSE COMPUTATION */
	/* ============================================== */

	start_time = mtime();
	inv(A, Ainv, R, N);
	end_time = mtime();
	printf("normal selfmade function evaluation in C++ of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";
	
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
		cout<< "Computation correct? "<< (bool) (abs(sum/N - 1) < N*1e-8 ) <<endl;
	}
	
	
	/* TIMING MY F77-IMPLEMENTATION INVERSE COMPUTATION */
	/* ================================================ */
	start_time = mtime();
	{
	int Ntmp = N;
	inv(A,Ainv,R, Ntmp);
	}
	end_time = mtime();
	printf("normal selfmade function evaluation in F77 of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";

	
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
		cout<< "Computation correct? "<< (bool) (abs(sum/N - 1) < N*1e-8 ) <<endl;
	}
	
	
	/* TIMING MY C++-IMPLEMENTATION INVERSE COMPUTATION */
	/* ================================================ */
	int Ntmp = N;
	start_time = mtime();
	inv(A,Ainv,R, N);
	end_time = mtime();
	printf("normal selfmade function evaluation in C of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";

	
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
		cout<< "Computation correct? "<< (bool) (abs(sum/N - 1) < N*1e-8 ) <<endl;
	}

	/* taping inv */
	cout<<"start tracing my implementation of inv with ADOL-C"<<endl;
	adouble *aA = new adouble[N*N];
	adouble *aR = new adouble[N*N];
	adouble *aAinv = new adouble[N*N];

	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	
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
	printf("normal selfmade function evaluation in C++ of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";	
	cout<<"finished tracing"<<endl;
	
	
	/* ================================ */
	/* TIMING TAPED INVERSE COMPUTATION */
	/* ================================ */
	start_time = mtime();
	function(0,N*N,N*N,A, Ainv);
	end_time = mtime();
	printf("taped function evaluation of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";

	
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
		cout<< "Computation correct? "<< (bool) (abs(sum/N - 1) < N*1e-8 ) <<endl;
	}
	
// 	/* TIMING MY C++-IMPLEMENTATION f COMPUTATION */
// 	/* ======================================== */
// 	start_time = mtime();
// 	f(A,Ainv);
// 	end_time = mtime();
// 	printf("normal function evaluation of f needs %d ms.\n",(end_time-start_time));
// 	runtimes_file<<end_time-start_time<<"\t";


// 	/* TIMING TAPED f COMPUTATION */
// 	/* ========================== */
// 	start_time = mtime();
// 	function(1,1,N*N,A.data, &y);
// 	end_time = mtime();
// 	printf("taped function evaluation of f needs %d ms.\n",(end_time-start_time));
// 	runtimes_file<<end_time-start_time<<"\t";


	cout<<endl<<endl;
	cout<<"================================================"<<endl;
	cout<<"           TESTING GRADIENT EVALUATION          "<<endl;
	cout<<"================================================"<<endl;
	cout<<endl<<endl;

	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	
	/* TIMING TAPED JACOBIAN COMPUTATION OF INV */
	if(N<=30){
        double **J;
        J=myalloc2(N*N,N*N);
        start_time = mtime();
        jacobian(0,N*N,N*N,A, J);
        end_time = mtime();
        printf("taped jacobian of inv needs %d ms.\n",(end_time-start_time));
    
    }
    
	/* taping inv */
	cout<<"start tracing my implementation of trace(inv(.)) with ADOL-C"<<endl;
	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	double y;
	trace_on(1);
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			aA[id(n,m,N)]<<= A[id(n,m,N)];
		}
	}
	adouble ay;
	trace(aA, aAinv, aR, &ay, N);
	ay >>= y;
	trace_off();
	cout<<"finished tracing"<<endl;
	
	/* ====================================== */
	/* TIMING TAPED GRADIENT COMPUTATION OF f */
	/* ====================================== */
	double *g;
	g = myalloc(N*N);
	start_time = mtime();
	gradient(1,N*N,A, g);
	end_time = mtime();
	printf("ADOL-C gradient of f needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";
	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	
	
	/* ===================================================== */
	/* COMPUTING THE GRADIENT OF f WITH MATRIX AD            */
	/* ===================================================== */
	/*              need to compute
	//              d tr(A^-1) = tr( - (A^-1)^2 dA) and therefore
	//              \nabla tr(A^-1) = - ((A^-1)^2)^T
	*/
	start_time = mtime();
	/* Ainv2 == A here */
	inv(A, Ainv, R, N);

	 /* compute Abar = - (A^{-1} A^{-1})^T */
	cblas_dgemm(CblasRowMajor, CblasTrans,
				 CblasTrans, N, N,
				N, -1., Ainv,
				 N, Ainv, N,
				 0., Abar, N);
	end_time = mtime();
	printf("UTPM gradient evaluation of f needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";
	
	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	
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
				sum += abs( (n==m) - Id[id(n,m,N)]);
			}
		}
		cout<< "Computation correct? "<< (bool) (sum< N*1e-8 ) <<endl;
	}

	/* CHECK THAT TAPED GRADIENT AND LAPACK GRADIENT EVALUATION ARE THE SAME */
	double difference = 0.;
	for(int n = 0; n != N; ++n){
		for(int m = 0; m != N; ++m){
			difference += abs(g[id(n,m,N)] - Abar[id(n,m,N)]);
			// cout<< g[id(n,m,N)]<<" "<<Abar[id(n,m,N)]<<endl;
		}
	} 
    cout<< "Computation correct? "<< (bool) (difference< N*1e-6 ) <<endl;

	
	/* ===================================================== */
	/* COMPUTING THE GRADIENT OF f WITH MATRIX AD via LAPACK */
	/* ===================================================== */
	/*              need to compute
	//              d tr(A^-1) = tr( - (A^-1)^2 dA) and therefore
	//              \nabla tr(A^-1) = - ((A^-1)^2)^T
	*/
	for(int n = 0; n != N; ++n){
		for(int m = 0; m != N; ++m){
		    Ainv[id(n,m,N)] = A[id(n,m,N)];
		}
	}
	
	start_time = mtime();
	/* compute the inverse by combination of dgetrf and dgetri */
	clapack_dgetrf(CblasRowMajor, N, N, Ainv, N, ipiv);
	clapack_dgetri(CblasRowMajor, N, Ainv, N, ipiv);
	
	 /* compute Abar = - (A^{-1} A^{-1})^T */
	cblas_dgemm(CblasRowMajor, CblasTrans,
				 CblasTrans, N, N,
				N, -1., Ainv,
				 N, Ainv, N,
				 0., Abar, N);
	end_time = mtime();
	printf("UTPM lapack gradient evaluation of f needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";

	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	

// 	/* CHECK THAT TAPED GRADIENT AND LAPACK GRADIENT EVALUATION ARE THE SAME */
// 	double difference = 0.;
// 	for(int n = 0; n != N; ++n){
// 		for(int m = 0; m != N; ++m){
// 			difference += abs(g[n*N + m] - Abar[n][m]);
// 		}
// 	} 
// 	assert(difference/(N*N)<=1e-6);
	
	
	/* ============================================= */
	/* COMPUTING THE GRADIENT OF trace WITH TAPENADE */
	/* ============================================= */
	
	/* we have to compute
		fbar df = fbar d tr C
				= tr ( fbar Id d C )
	*/
	

	/* since fbar = 1, fbar Id = Id */
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			Id[id(n,m,N)] = (n == m);
		}
	}
	
	/* now compute Abar */
	double *Ainvbar = new double[N*N];
	double *Rbar = new double[N*N];
	
	start_time = mtime();
	inv_b(A, Abar, Ainv, Ainvbar, R, Rbar, N);
	end_time = mtime();
	printf("TAPENADE gradient evaluation of f using the c code needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<endl;
	
	runtimes_file.close();
	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	
	// // fortran version
	// start_time = mtime();
	// int Ntmp2 = N;
	// inv_b_(A, Abar, Ainv, Id, R, Ntmp2);
	// end_time = mtime();
	// printf("TAPENADE gradient evaluation of f using the fortrang code needs %d ms.\n",(end_time-start_time));
	// runtimes_file<<end_time-start_time<<endl;	

	return 0;
}

