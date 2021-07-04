#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <memory.h>
#include <malloc.h>

//for mordor8
//#include "/mnt/scratch/sw/papi-3.6.2/include/papi.h"
#include "papi.h"
#include "/share/apps/mkl/include/mkl_cblas.h"
//#include "mkl_cblas.h"
#include <omp.h>

#define dimension 10000
#define Nthreads 8
#define chunk 1


/*OMP function prototypes*/
void Task_col(double *a, int *ready, int n, int col);
void Col_mod(double *a, int *ready, int n, int col, int prcol);
void Col_div(double *a, int n, int col);
void Task_col_comp(double *a, int *ready, int n, int col);
void Col_mod_comp(double *a, int *ready, int n, int col, int prcol);
void Col_div_comp(double *a, int n, int col);

static void test_fail(char *file, int line, char *call, int retval)
{

	printf("%s\tFAILED\nLine # %d\n", file, line);
	if ( retval == PAPI_ESYS ) {
		char buf[128];
		memset( buf, '\0', sizeof(buf) );
		sprintf(buf, "System error in %s:", call );
		perror(buf);
		}
	else if ( retval > 0 ) {
			printf("Error calculating: %s\n", call );
			}
	else {
		char errstring[PAPI_MAX_STR_LEN];
		PAPI_perror(retval, errstring, PAPI_MAX_STR_LEN );
		printf("Error in %s: %s\n", call, errstring );
	}
	printf("\n");
	exit(1);
}
void show_matrix(double **A, int n,const char *argu) {
     
	 int i,j;
	 for(i = 0; i < n; i++){
         for(j = 0; j < n; j++){
             if (argu=="lower"){
                if(j <= i)
                    printf("%2.5f ", A[i][j]);
                else
                    printf("%2.5f ",0.0);
                    }
             else{
                 if (argu=="upper"){
                    if(j >= i)
                        printf("%2.5f ", A[i][j]);
                    else
                        printf("%2.5f ",0.0);
                }
              }

         }
        printf("\n");
     }

     printf("\n");
}

void show_1dmatrix(const char * desc,double *a, int n){
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < n; i++ ) {
		for( j = 0; j < n; j++ ){
			if (j<=i)
				 printf( " %2.5f", a[i*n+j] );
			else
				 printf( " %2.5f", 0.0);
				}
		    printf( "\n" );		
	}
}

void show_comp_matrix(const char * desc,double *a, int n){

	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < n; i++ ) {
		for( j = 0; j < n; j++ ){
			if (j<=i)
				 printf( " %2.5f", a[i+j*(2*n-j-1)/2] );
			else
				 printf( " %2.5f", 0.0);
				}
		    printf( "\n" );		
	}
}

void chol_col_ori(double **a, int n){

	  int i,j,l,k;
	  for(j = 0; j < n; j++){
		  for (k = 0; k < j; k++)
			  for (i = j; i < n; i++)
				  a[i][j] = a[i][j] - a[i][k]*a[j][k];

	  if (a[j][j]<0){
			printf("Error in Column-Oriented Cholseky Decomposition:\n");
		    printf("Matrix is not positive definite OR dimension of matrix is not set correctly\n");
			exit(EXIT_FAILURE);
				}
	  // Square root of the diagonal.
	  a[j][j] = sqrt(a[j][j]);

	  // Scale the row below.
	  for(k = (j+1); k < n; k++)
		  a[k][j] = a[k][j]/a[j][j];

	  }
	  printf("\nThe lower triangular matrix L in A=L.LT :\n\n");
	  const char *string="lower";
	  //show_matrix(a,n,string);
}


int main (int argc, char *argv[])
{
    
	float real_time, proc_time, mflops;
    	long long flpins;
	int retval;
	
	int i, j, k, l;
    	int n, info;
	double t;
	double *m1, *m, **a,*a_orig,*a_serial,*a_comp_ser,*a_comp;
	int *ready_ser, *ready;		
	// Get the dimensions from the command-line if available.
	if(argc > 1)
		n = atoi(argv[1]);
	else
	    n =4;
	//allocate for ready flag vector
	ready = (int *)calloc(n,sizeof(int));
	ready_ser = (int *)calloc(n,sizeof(int));
	// Allocate the memory for the matrix.
	a = (double **)malloc(n*sizeof(double*));
	a_orig = (double *)malloc(n*n*sizeof(double));
	a_serial = (double *)malloc(n*n*sizeof(double));
	a_comp_ser = (double *)malloc((n*n+n)/2*sizeof(double));
	a_comp = (double *)malloc((n*n+n)/2*sizeof(double));
	m = (double *)malloc(n*n*sizeof(double));
	m1 = (double *)malloc(n*n*sizeof(double));
	
	for(i = 0; i < n; i++){
	//	a_orig[i] = (double *)malloc(n*sizeof(double));
		a[i] = (double *)malloc(n*sizeof(double));
	}
							
										
	/*double m[]={5.0, 1.2, 0.3, -0.6,
				1.2, 6.0, -0.4, 0.9,
				0.3, -0.4, 8.0, 1.7,
			   -0.6, 0.9, 1.7, 10.0};*/

	/* double m[]={1.0, 1.0, 1.0, 1.0, 1.0,
		         1.0, 2.0, 3.0, 4.0, 5.0,
				 1.0, 3.0, 6.0, 10.0,15.0,
				 1.0, 4.0, 10.0, 20.0,35.0,
				 1.0, 5.0, 15.0, 35.0, 70.0};*/


	/*********Constructing Tridiagonal Matrix**********/
	 /*for(i = 0; i < n; i++)
        	for(j = 0; j < n; j++){
                	if(i == j)
                        	m[i*n+j] = 4.0;
                        else if(abs(i - j) == 1.0)
                        	m[i*n+j] = -1.0;
                        else
                        	m[i*n+j] = 0.0;
                          };*/

	/* Generating random positive definite matrix*/
	for (i=0; i<n*n; i++ ){
	  m1[i] = (double)rand() / (double)RAND_MAX;}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, m1, n, m1, n, 0.0, m, n);

    	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			a[i][j] = m[i*n+j];
			a_serial[i*n+j]=m[i*n+j];
			a_orig[i*n+j] = m[i*n+j];
			a_comp[i+(j*(2*n-j-1))/2] = m[i*n+j];
			a_comp_ser[i+(j*(2*n-j-1))/2] = m[i*n+j];
		
		}
	}

    	
	/*********SERIAL CONCEPTUALLY SIMILAR METHOD FOR COMPARISON*********/

	
	/* Setup PAPI library and begin collecting data from the counters */
	if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
		test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
	  
	/* Operation to be timed */
	//chol_col_ori(a,n);  

	/* Collect the data into the variables passed in */
	if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
			    test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
	
	printf("For Basic Serial Without Call to BLAS(General Comparison):\n");	  
	printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
				    real_time, proc_time, flpins, mflops);
	printf("%s\tPASSED\n", __FILE__);
	PAPI_shutdown();
	


        /*************SERIAL UNPACKED VERSION WITH CALLING BLAS FUNCTIONS AND USING FUNCTION*********/
                /*************MAIN SERIAL METHOD FOR COMPARISON WITH OMP CHOLESKY**********/

	if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
		test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
        
	for (j=0; j<n; j++){

                Task_col(a_serial,ready_ser,n,j);
        }
	/* Collect the data into the variables passed in */
	if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
			    test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
	
	printf("\nUnpacked Serial Cholesky With Using Call to BLAS:\n");	  
	printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
				    real_time, proc_time, flpins, mflops);
	printf("%s\tPASSED\n", __FILE__);
    
	double flop = (n/1e3*n/1e3*n/1e3/3 ) / (double)(real_time);
        printf("my Unpacked Serial  GFLOPS %lf \n",flop);
	
	PAPI_shutdown();
	
	//printf("\n");
	//printf("Do you want to see the Matrices L in Serial Cholseky and OpenMP cholseky: (y or n)");
	//char release=fgetc(stdin);
	char release = 'n';
	if (release=='n'){
		//return 0;
	}
	else{
	show_1dmatrix("Unpacked Serail (Matrix is Lower Triangular)",a_serial,n);
	
	}

        /*************SERIAL PACKED VERSION WITH CALLING BLAS FUNCTIONS AND USING FUNCTION*********/
                /*************MAIN SERIAL METHOD FOR COMPARISON WITH OMP CHOLESKY**********/
	
	for (i=0;i<n;i++)
		ready_ser[i]=0;

	/* Setup PAPI library and begin collecting data from the counters */
	if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
		test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
        for (j=0; j<n; j++){

                Task_col_comp(a_comp_ser,ready_ser,n,j);
        }
	/* Collect the data into the variables passed in */
	if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
			    test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
	
	printf("\nPacked Serial Cholesky With Using Call to BLAS:\n");	  
	printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
				    real_time, proc_time, flpins, mflops);
	printf("%s\tPASSED\n", __FILE__);
	
	flop = (n/1e3*n/1e3*n/1e3/3 ) / (double)(real_time);
        printf("my Packed Serial  GFLOPS %lf \n",flop);
	
	PAPI_shutdown();
	
	if (release=='n'){
		//return 0;
	}
	else{
	
	show_comp_matrix("Packed Serial (Matrix is Lower Triangular)",a_comp_ser,n);
	
	}

	/**********OpenMP section************/

 	printf("\nStarting OMP part:\n");
	
	int MaxNt;
	MaxNt = omp_get_max_threads();
	if (Nthreads>MaxNt){
		printf("Warning");
		printf("Number of requested threads are greater than max number of available threads\n");
		printf("Requested %d, Maximum %d\n",Nthreads,MaxNt);
	}
	omp_set_num_threads(Nthreads);
		
	//timing variable
	double start, end;
 	int tid;

        start = omp_get_wtime();
 
        //Starting Parallel Region
	#pragma omp parallel shared(a_orig,ready) private(tid,i,j)
        {
        tid = omp_get_thread_num();
 	#pragma omp for schedule (dynamic)

        for (j=0; j<n; j++){
                //printf("Thread %d did column %d\n",tid,j);

                Task_col(a_orig,ready,n,j);
        }
        }/*** End of parallel region ***/
        
	end = omp_get_wtime();

	printf("\nTime Passed for Unpacked Parallel:%lf\n",(end-start));
	if (release=='n'){
		//return 0;
	}
	else{
	
	show_1dmatrix("Unpacked Parallel (Matrix is Lower Triangular)",a_orig,n);
	}

	flop = (n/1e3*n/1e3*n/1e3/3 ) / (double)(end-start);
        printf("my Packed Serial  GFLOPS %lf \n",flop);


	/****Packed parallel******/
	for(i=0;i<n;i++)
		ready[i]=0;

        start = omp_get_wtime();
 
        //Starting Parallel Region
	#pragma omp parallel shared(a_comp,ready) private(tid,i,j)
        {
        tid = omp_get_thread_num();
 	#pragma omp for schedule (dynamic)

        for (j=0; j<n; j++){
                //printf("Thread %d did column %d\n",tid,j);

                Task_col_comp(a_comp,ready,n,j);
        }
        }/*** End of parallel region ***/
        
	end = omp_get_wtime();

	printf("\nTime Passed for Packed Parallel:%lf\n",(end-start));
	if (release=='n'){
		//return 0;
	}
	else{
	
	show_comp_matrix("Packed Parallel (Matrix is Lower Triangular)",a_comp,n);
	}


	flop = (n/1e3*n/1e3*n/1e3/3 ) / (double)(end-start);
        printf("my Packed Parallel  GFLOPS %lf \n",flop);
	
	/*collecting garbage*/
	
	for (i=0;i<n;i++)
		free(a[i]);
	free(a);free(m);free(a_serial);free(a_comp);free(a_orig);free(a_comp_ser);
	free(ready_ser);free(ready);
	return 0;
}

void Task_col(double *a, int *ready, int n, int col){
        int k;
        //int tid = omp_get_thread_num();
        //printf("tid: %d my column is %d\n",tid,col);

        if (col==0){
                Col_div(a, n, 0);
                ready[0] = 1;
                }
        else{
                for (k = 0; k < col; k++){

                        /*Statement for keeping thread busy in case column k is not still ready*/
                        while (ready[k] == 0){}
                        Col_mod(a, ready, n, col, k);
                        }
                Col_div(a, n, col);
		ready[col] = 1;
        }
}

void Col_mod(double *a, int *ready, int n, int col, int prcol){
        //int i, tid=omp_get_thread_num();
        //printf("hello from %d modifying col: %d with ready[%d]=%d\n",tid,prcol,prcol,ready[prcol]);

        /*performing DAXPY from mkl_cblas*/
        cblas_daxpy(n-col, -a[col*n+prcol], &a[col*n+prcol], n, &a[col*n+col],n);
}

void Col_div(double *a, int n, int col){
        
        a[col*n+col] = sqrt(a[col*n+col]);
        cblas_dscal(n-col-1, 1.0/a[col*n+col], &a[(col+1)*n+col], n);    
}


void Task_col_comp(double *a, int *ready, int n, int col){
        int k;
        if (col==0){
                Col_div_comp(a, n, 0);
                ready[0] = 1;
                }
        else{
                for (k = 0; k < col; k++){

                        /*Statement for keeping thread busy in case column k is not still ready*/
                        while (ready[k] == 0){}
                        Col_mod_comp(a, ready, n, col, k);
                        }
                Col_div_comp(a, n, col);
                ready[col] = 1;
        }
}

void Col_mod_comp(double *a, int *ready, int n, int col, int prcol){

        /*performing DAXPY from mkl_cblas*/
        cblas_daxpy(n-col, -a[col+prcol*(2*n-prcol-1)/2], &a[col+prcol*(2*n-prcol-1)/2], 1, &a[col+col*(2*n-col-1)/2],1);
}

void Col_div_comp(double *a, int n, int col){

        a[col+col*(2*n-col-1)/2] = sqrt(a[col+col*(2*n-col-1)/2]);
        cblas_dscal(n-col-1, 1.0/a[col+col*(2*n-col-1)/2], &a[col+1+col*(2*n-col-1)/2], 1);
}



/*void chol_Serial(double *a_orig,int *ready,int n){
 	
	int i, j, k;
	for (j=0; j<n; j++){
		
		if (j==0){
			a_orig[0] = sqrt(a_orig[0]);
			cblas_dscal(n-1, 1.0/a_orig[0], &a_orig[n], n);
			}
		else{
			for (k = 0; k < j; k++){
				cblas_daxpy(n-j, -a_orig[j*n+k], &a_orig[j*n+k], n, &a_orig[j*n+j],n);
				}	
		
			a_orig[j*n+j] = sqrt(a_orig[j*n+j]);
			cblas_dscal(n-j-1, 1.0/a_orig[j*n+j], &a_orig[(j+1)*n+j], n);
		}
	}
}
*/

/*for packed format*/
void chol_Serial(double *a_orig,int *ready,int n){
 	
	int i, j, k;
	
	for (j=0; j<n; j++){
		
		if (j==0){
			a_orig[0] = sqrt(a_orig[0]);
			cblas_dscal(n-1, 1.0/a_orig[0], &a_orig[1], 1);
			}
		else{
			for (k = 0; k < j; k++){
				cblas_daxpy(n-j, -a_orig[j+k*(2*n-k-1)/2], &a_orig[j+k*(2*n-k-1)/2], 1, &a_orig[j+j*(2*n-j-1)/2],1);
				}	
		
			a_orig[j+j*(2*n-j-1)/2] = sqrt(a_orig[j+j*(2*n-j-1)/2]);
			cblas_dscal(n-j-1, 1.0/a_orig[j+j*(2*n-j-1)/2], &a_orig[j+1+j*(2*n-j-1)/2], 1);
		}
	}
}

