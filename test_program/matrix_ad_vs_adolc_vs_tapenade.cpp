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
extern "C" void inv_b_(double *A, double *Ab, double *QT, double *QTb, double *R, int *NA);
extern "C" void inv_b(double *a, double *ab, double *qt, double *qtb, double *r, double *rb, int na) ;


extern "C" void dgetrf_( int *M, int *N, double * A, int * LDA, int * IPIV,
	int *INFO );
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


extern "C" void dgels_(char *TRANS,  int *M, int *N, int *NRHS,
	double *A, int *LDA, double *B, int *LDB, double *WORK, int *LWORK, int *INFO);
/**<
      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
     $                  INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGELS solves overdetermined or underdetermined real linear systems
*  involving an M-by-N matrix A, or its transpose, using a QR or LQ
*  factorization of A.  It is assumed that A has full rank.
*
*  The following options are provided:
*
*  1. If TRANS = 'N' and m >= n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || B - A*X ||.
*
*  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
*     an underdetermined system A * X = B.
*
*  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
*     an undetermined system A**T * X = B.
*
*  4. If TRANS = 'T' and m < n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || B - A**T * X ||.
*
*  Several right hand side vectors b and solution vectors x can be
*  handled in a single call; they are stored as the columns of the
*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
*  matrix X.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          = 'N': the linear system involves A;
*          = 'T': the linear system involves A**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of
*          columns of the matrices B and X. NRHS >=0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit,
*            if M >= N, A is overwritten by details of its QR
*                       factorization as returned by DGEQRF;
*            if M <  N, A is overwritten by details of its LQ
*                       factorization as returned by DGELQF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the matrix B of right hand side vectors, stored
*          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
*          if TRANS = 'T'.
*          On exit, if INFO = 0, B is overwritten by the solution
*          vectors, stored columnwise:
*          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
*          squares solution vectors; the residual sum of squares for the
*          solution in each column is given by the sum of squares of
*          elements N+1 to M in that column;
*          if TRANS = 'N' and m < n, rows 1 to N of B contain the
*          minimum norm solution vectors;
*          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
*          minimum norm solution vectors;
*          if TRANS = 'T' and m < n, rows 1 to M of B contain the
*          least squares solution vectors; the residual sum of squares
*          for the solution in each column is given by the sum of
*          squares of elements M+1 to N in that column.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= MAX(1,M,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          LWORK >= max( 1, MN + max( MN, NRHS ) ).
*          For optimal performance,
*          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
*          where MN = min(M,N) and NB is the optimum block size.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO =  i, the i-th diagonal element of the
*                triangular factor of A is zero, so that A does not have
*                full rank; the least squares solution could not be
*                computed.
*
*  =====================================================================
*/


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
	dgetrf_(&N, &N, Ainv, &N, ipiv, &info);
	// cout<<"info="<<info<<endl;
	dgetri_(&N, Ainv, &N, ipiv, WORK, &LWORK, &info);
	// cout<<"info="<<info<<endl;
	end_time = mtime();
	printf("normal LAPACK function evaluation of inv needs %d ms.\n",(end_time-start_time));
	runtimes_file<<end_time-start_time<<"\t";


	/* check that dot(A,Ainv) = I */
	double one = 1;
	double zero = 0;
	char notrans = 'n';
	dgemm_(&notrans, &notrans, &N, &N, &N, &one, A,
				 &N, Ainv, &N,
				&zero, Id, &N);
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
	inv_(A,Ainv,R, &Ntmp);
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
	
	// cout<<(Ainv[0]*Ainv[0])<<" "<<Abar[0]<<endl;

	/* CHECK THAT TAPED GRADIENT AND LAPACK GRADIENT EVALUATION ARE THE SAME */
	double difference = 0.;
	for(int n = 0; n != N; ++n){
		for(int m = 0; m != N; ++m){
			difference += abs(g[id(n,m,N)] - Abar[id(n,m,N)]);
			// cout<< g[id(n,m,N)]<<" "<<Abar[id(n,m,N)]<<endl;
		}
	} 
    cout<< "Computation correct? "<< (bool) (difference< N*1e-6 ) <<endl;
		// cout<<difference<<endl;

	
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
	// int ipiv[N];	int info;
	
	dgetrf_(&N, &N, Ainv, &N, ipiv, &info);
	// cout<<"info="<<info<<endl;
	dgetri_(&N, Ainv, &N, ipiv, WORK, &LWORK, &info);
	// clapack_dgetrf(CblasRowMajor, N, N,
	// 			   Ainv2, N, ipiv);
	// clapack_dgetri(CblasRowMajor, N, Ainv, N, ipiv);

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
	// inv_b_(A, Abar, Ainv, Id, R, &Ntmp2);
	// end_time = mtime();
	// printf("TAPENADE gradient evaluation of f using the fortrang code needs %d ms.\n",(end_time-start_time));
	// runtimes_file<<end_time-start_time<<endl;	

	return 0;
}

