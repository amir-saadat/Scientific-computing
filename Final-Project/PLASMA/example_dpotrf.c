/**
 *
 * @file example_dpotrf.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example of Cholesky factorization
 *
 * @version 2.4.5
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated d Tue Nov 22 14:35:53 2011
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//#include </home/asaadat/project/plasma-installer_2.4.5/install/include/plasma.h>
//#include </home/asaadat/project/plasma-installer_2.4.5/install/include/cblas.h>
//#include </home/asaadat/project/plasma-installer_2.4.5/install/include/lapacke.h>
//#include </home/asaadat/project/plasma-installer_2.4.5/install/include/core_blas.h>
#include <plasma.h>
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>

#include <sys/time.h>
#include <time.h>

double get_cur_time();

int check_factorization(int, double*, double*, int, int);

int IONE=1;
int ISEED[4] = {0,0,0,1};   /* initial seed for dlarnv() */

int main (int argc, char *argv[])
{

    int cores, N, LDA;   
 
    if(argc > 1){
       cores = atoi(argv[1]);
       N = atoi(argv[2]);
       LDA   = N ;
	}
    else{
       cores = 1;
       N     = 5000 ;
       LDA   = 5000 ;
        }
    int info_factorization;

    double *A1   = (double *)malloc(LDA*N*sizeof(double));
    double *A2   = (double *)malloc(LDA*N*sizeof(double));

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)){
        printf("Out of Memory \n ");
        return EXIT_SUCCESS;
    }

    /* Plasma Initialize */
    PLASMA_Init(cores);
    printf("-- PLASMA is initialized to run on %d cores. \n",cores);

    /* Initialize A1 and A2 for Symmetric Positive Matrix */
    PLASMA_dplgsy( (double)N, N, A1, LDA, 51 );
    PLASMA_dlacpy( PlasmaUpperLower, N, N, A1, LDA, A2, LDA );

    /* Plasma routines */

    double time1=0.0;
    double time2=0.0;
    double flops=0.0;

    time1 = get_cur_time();
    PLASMA_dpotrf(PlasmaUpper, N, A2, LDA);
    time2 = get_cur_time();

    flops = (N/1e3*N/1e3*N/1e3/3 ) / (double)(time2-time1);
    //flops = (double) (N*N*N/3.0 ) / (double)(time2-time1)  ;
    printf("DPOTRF  N:%d time(s): %f GFLOPS:%f \n", N, time2-time1, flops);


    /* Check the factorization */
    info_factorization = check_factorization( N, A1, A2, LDA, PlasmaUpper);

    if ( info_factorization != 0 )
       printf("-- Error in DPOTRF example ! \n");
    else
       printf("-- Run of DPOTRF example successful ! \n");

    free(A1); free(A2);

    PLASMA_Finalize();

    return EXIT_SUCCESS;
}


/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */

double get_cur_time() 
{
  double time;
  struct timeval   tv;
  struct timezone  tz;
  
  gettimeofday(&tv, &tz);
  time = tv.tv_sec + tv.tv_usec / 1000000.0;
  return time;
} 


int check_factorization(int N, double *A1, double *A2, int LDA, int uplo)
{
    double Anorm, Rnorm;
    double alpha;
    int info_factorization;
    int i,j;
    double eps;

    eps = LAPACKE_dlamch_work('e');

    double *Residual = (double *)malloc(N*N*sizeof(double));
    double *L1       = (double *)malloc(N*N*sizeof(double));
    double *L2       = (double *)malloc(N*N*sizeof(double));
    double *work              = (double *)malloc(N*sizeof(double));

    memset((void*)L1, 0, N*N*sizeof(double));
    memset((void*)L2, 0, N*N*sizeof(double));

    alpha= 1.0;

    LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,' ', N, N, A1, LDA, Residual, N);

    /* Dealing with L'L or U'U  */
    if (uplo == PlasmaUpper){
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L1, N);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L2, N);
        cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    }
    else{
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L1, N);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L2, N);
        cblas_dtrmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    }

    /* Compute the Residual || A -L'L|| */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
           Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];

    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Residual, N, work);
    Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A1, LDA, work);

    printf("============\n");
    printf("Checking the Cholesky Factorization \n");
    printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));

    if ( isnan(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 10.0) ){
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else{
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(Residual); free(L1); free(L2); free(work);

    return info_factorization;
}
