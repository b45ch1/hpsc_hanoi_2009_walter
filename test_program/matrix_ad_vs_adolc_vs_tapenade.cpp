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
}

extern "C"{
#include "cblas.h"
}

#include <sys/time.h>
struct timeval tv;
int mtime(void)
{
gettimeofday(&tv,NULL);
return (int)(tv.tv_sec*1000 + (tv.tv_usec / 1000));
}




extern "C" void inv_(double *A, double *QT, double *R, int *NA);
extern "C" void inv_b__(double *A, double *Ab, double *QT, double *QTb, double *R, int *NA);

extern "C" void dgetrf_( long *M, long *N, double * A, long * LDA, long * IPIV,
	long *INFO );
/**<
*SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
	  INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
	  INTEGER            IPIV( * )
	  DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*/

extern "C" void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK,
	int* LWORK, int* INFO);
///< \brief computes the inverse of a column major matrix
/**<
for example: dgetri_(&n,&BlasA(0,0),&n,&ipiv[0],&work[0],&n,&info);
 *  \param N       (input) INTEGER
 *          The order of the matrix A.  N >= 0.
 *
 *  \param A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
 *          On entry, the factors L and U from the factorization
 *          A = P*L*U as computed by DGETRF.
 *          On exit, if INFO = 0, the inverse of the original matrix A.
 *
 *   \param LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  \param IPIV    (input) INTEGER array, dimension (N)
 *          The pivot indices from DGETRF; for 1<=i<=N, row i of the
 *          matrix was interchanged with row IPIV(i).
 *
 *  \param WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
 *          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
 *
 * \param LWORK   (input) INTEGER
 *          The dimension of the array WORK.  LWORK >= max(1,N).
 *          For optimal performance LWORK >= N*NB, where NB is
 *          the optimal blocksize returned by ILAENV.
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the WORK array, returns
 *          this value as the first entry of the WORK array, and no error
 *          message related to LWORK is issued by XERBLA.
 *
 *  \param INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
 *                singular and its inverse could not be computed.
**/

extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
	double *ALPHA, double *A, int *LDA, double *B, int *LDB, double *BETA,
	double *C, int *LDC);
/**<
	  SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
	  DOUBLE PRECISION ALPHA,BETA
	  INTEGER K,LDA,LDB,LDC,M,N
	  CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
	  DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/

// void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
				// const enum CBLAS_TRANSPOSE TransB, const long M, const long N,
				// const long K, const double alpha, const double *A,
				// const long lda, const double *B, const long ldb,
				// const double beta, double *C, const long ldc);

				// enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
						// AtlasConj=114};



using namespace std;
using namespace mtl;

/// \brief compute the trace of the inverse of a matrix
template<class Tdouble>
Tdouble f(const mtl::dense2D<Tdouble> &A, mtl::dense2D<Tdouble> &Ainv){
	Tdouble retval(0);
	inv(A,Ainv);
	for(int r = 0; r != A.num_rows(); ++r){
		retval += Ainv[r][r];
	}
	return retval;
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

	/*
	SETUP:
	======
	
		We use the mtl library because of the convenient defined ostream operators, so we can do
		cout << A <<endl;
		to output A.
		
	*/

	
	mtl::dense2D<double> A(N,N);
	mtl::dense2D<double> Ainv(N,N);
	mtl::dense2D<double> Ainv2(N,N);
	mtl::dense2D<double> QT(N,N);
	mtl::dense2D<double> Q(N,N);
	mtl::dense2D<double> B(N,N);
	mtl::dense2D<double> R(N,N);
	mtl::dense2D<double> Id(N,N);
	mtl::dense2D<double> dA(N,N);


	int start_time,end_time;

	/* fill A with random numbers */
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			A[n][m] = rand()/100000000000.;
			Ainv2[n][m] = A[n][m];
		}
	}


	/* taping inv */
	mtl::dense2D<adouble> aA(N,N);
	mtl::dense2D<adouble> aAinv(N,N);

	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);

	trace_on(0);
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			aA[n][m]<<= A[n][m];
		}
	}
	inv(aA,aAinv);
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			aAinv[n][m]>>= Ainv[n][m];
		}
	}
	trace_off();

	/* taping  f */
	start_time = mtime();

	adouble ay;
	double y;
	trace_on(1);
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			aA[n][m]<<= A[n][m];
		}
	}
	ay = f(aA,aAinv);
	ay >>= y;
	trace_off();
	end_time = mtime();
	printf("taping of f needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";
	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	
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
	
	
	/* TIMING ATLAS IMPLEMENTATION INVERSE COMPUTATION */
	/* =============================================== */
	/* compute the inverse by combination of dgetrf and dgetri */
	Ainv = A;
	start_time = mtime();
	{
	long ipiv[N];	long info;
	clapack_dgetrf(CblasRowMajor, N, N,
				   Ainv.data, N, ipiv);
	clapack_dgetri(CblasRowMajor, N, Ainv.data, N, ipiv);
	}
	end_time = mtime();
	printf("normal ATLAS function evaluation of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";


	/* check that dot(A,Ainv) = I */
	cblas_dgemm(CblasRowMajor, CblasNoTrans,
				 CblasNoTrans, N, N,
				N, 1., A.data,
				 N, Ainv.data, N,
				 0., Id.data, N);
	{
		double sum = 0;
		for (int n = 0; n != N; ++n){
			for (int m = 0; m != N; ++m){
				sum += Id[n][m];
			}
		}
		cout<< "Computation correct? "<< (bool) (abs(sum/N - 1) < N*1e-8 ) <<endl;
	}
	

	/* TIMING MY C-IMPLEMENTATION INVERSE COMPUTATION */
	/* ============================================== */
	start_time = mtime();
	inv(A,Ainv);
	end_time = mtime();
	printf("normal function evaluation in C++ of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";

	
	/* check that dot(A,Ainv) = I */
	cblas_dgemm(CblasRowMajor, CblasNoTrans,
				 CblasNoTrans, N, N,
				N, 1., A.data,
				 N, Ainv.data, N,
				 0., Id.data, N);
	{
		double sum = 0;
		for (int n = 0; n != N; ++n){
			for (int m = 0; m != N; ++m){
				sum += Id[n][m];
			}
		}
		cout<< "Computation correct? "<< (bool) (abs(sum/N - 1) < N*1e-8 ) <<endl;
	}

	/* TIMING MY F77-IMPLEMENTATION INVERSE COMPUTATION */
	/* ================================================ */
	int Ntmp = N;
	start_time = mtime();
	inv_(A.data,Ainv.data,R.data, &Ntmp);
	end_time = mtime();
	printf("normal function evaluation in F77 of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";

	
	/* check that dot(A,Ainv) = I */
	cblas_dgemm(CblasRowMajor, CblasNoTrans,
				 CblasNoTrans, N, N,
				N, 1., A.data,
				 N, Ainv.data, N,
				 0., Id.data, N);
	{
		double sum = 0;
		for (int n = 0; n != N; ++n){
			for (int m = 0; m != N; ++m){
				sum += Id[n][m];
			}
		}
		cout<< "Computation correct? "<< (bool) (abs(sum/N - 1) < N*1e-8 ) <<endl;
	}
	
	/* TIMING TAPED INVERSE COMPUTATION */
	/* ================================ */
	start_time = mtime();
	function(0,N*N,N*N,A.data, Ainv.data);
	end_time = mtime();
	printf("taped function evaluation of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";

	
	/* check that dot(A,Ainv) = I */
	cblas_dgemm(CblasRowMajor, CblasNoTrans,
				 CblasNoTrans, N, N,
				N, 1., A.data,
				 N, Ainv.data, N,
				 0., Id.data, N);
	{
		double sum = 0;
		for (int n = 0; n != N; ++n){
			for (int m = 0; m != N; ++m){
				sum += Id[n][m];
			}
		}
		cout<< "Computation correct? "<< (bool) (abs(sum/N - 1) < N*1e-8 ) <<endl;
	}
	
	/* TIMING MY C++-IMPLEMENTATION f COMPUTATION */
	/* ======================================== */
	start_time = mtime();
	f(A,Ainv);
	end_time = mtime();
	printf("normal function evaluation of f needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";


	/* TIMING TAPED f COMPUTATION */
	/* ========================== */
	start_time = mtime();
	function(1,1,N*N,A.data, &y);
	end_time = mtime();
	printf("taped function evaluation of f needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";


// 	/* TIMING TAPED JACOBIAN COMPUTATION OF INV */
// 	double **J;
// 	J=myalloc2(N*N,N*N);
// 	start_time = mtime();
// 	jacobian(0,N*N,N*N,A.data, J);
// 	end_time = mtime();
//     printf("taped jacobian of inv needs %d ms.\n",(end_time-start_time));


	////////////////////////////////////////////////////
	////////// TESTING DERIVATIVE EVALUATION ///////////
	////////////////////////////////////////////////////
	cout<<endl<<endl;
	cout<<"================================================"<<endl;
	cout<<"           TESTING GRADIENT EVALUATION          "<<endl;
	cout<<"================================================"<<endl;
	cout<<endl<<endl;

	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	

	/* ====================================== */
	/* TIMING TAPED GRADIENT COMPUTATION OF f */
	/* ====================================== */
	double *g;
	g = myalloc(N*N);
	start_time = mtime();
	gradient(1,N*N,A.data, g);
	end_time = mtime();
	printf("ADOL-C gradient of f needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";
	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);

	
	/* ===================================================== */
	/* COMPUTING THE GRADIENT OF f WITH MATRIX AD via LAPACK */
	/* ===================================================== */
	/*              need to compute
	//              d tr(A^-1) = tr( - (A^-1)^2 dA) and therefore
	//              \nabla tr(A^-1) = - ((A^-1)^2)^T
	*/
	start_time = mtime();
	/* compute the inverse by combination of dgetrf and dgetri */
	/* Ainv2 == A here */
	long ipiv[N];	long info;
	clapack_dgetrf(CblasRowMajor, N, N,
				   Ainv2.data, N, ipiv);
	clapack_dgetri(CblasRowMajor, N, Ainv2.data, N, ipiv);

	 /* compute Abar = - (A^{-1} A^{-1})^T */
	cblas_dgemm(CblasRowMajor, CblasTrans,
				 CblasTrans, N, N,
				N, -1., Ainv2.data,
				 N, Ainv2.data, N,
				 0., dA.data, N);
	end_time = mtime();
	printf("lapack gradient evaluation of f needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";

	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);
	

	/* CHECK THAT TAPED GRADIENT AND LAPACK GRADIENT EVALUATION ARE THE SAME */
	double difference = 0.;
	for(int n = 0; n != N; ++n){
		for(int m = 0; m != N; ++m){
			difference += abs(g[n*N + m] - dA[n][m]);
		}
	} 
	assert(difference/(N*N)<=1e-6);
	
	
	/* ========================================= */
	/* COMPUTING THE GRADIENT OF f WITH TAPENADE */
	/* ========================================= */
	
	/* we have to compute
		fbar df = fbar d tr C
				= tr ( fbar Id d C )
	*/
	
	start_time = mtime();

	/* since fbar = 1, fbar Id = Id */
	for(int n = 0; n!=N; ++n){
		for(int m = 0; m!=N; ++m){
			Id[n][m] = (n == m);
		}
	}
	
	/* now compute Abar */
	int Ntmp2 = N;
	inv_b__(A.data, dA.data, Ainv.data, Id.data, R.data, &Ntmp2);
	
	
	
	end_time = mtime();
	printf("TAPENADE gradient evaluation of f needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<endl;

	runtimes_file.close();
	snprintf (mycmd, (size_t)255, "cat /proc/%d/status | grep VmPeak >> mem_consumption.txt", getpid());
	system(mycmd);







	
	return 0;
}

