#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <memory.h>
#include <malloc.h>
#include "papi.h"
#include "cblas.h"

 
static void test_fail(char *file, int line, char *call, int retval);

int main(int argc, char** argv) {
	
/*  extern void dummy(void *);*/
  float real_time, proc_time, mflops;
  long long flpins;
  int retval;
  int i,j,k,l;
  int Num_Tot_n = 26;/*Number of different n used here*/
  FILE * InFile;
  FILE * Result2;/*Execuation rate vs Matrix dimension*/
  FILE * ATLAS;/*Final result of Matrix C*/
  
  int dim,dim1; /*Matrix Dimension*/
  double e,f;/*Used for reading data from myfile.txt*/
  /*int n;*/ /*Matrices Dimension*/
  InFile = fopen ( "myfile.txt" , "rb" );
  ATLAS= fopen ( "atlas.txt" , "w+" );
  Result2= fopen ( "result2.txt" , "w+" );
  
  for (l= 0; l<Num_Tot_n; l++)
  {
  fscanf (InFile, "%i", &dim);
  int n = dim;
  printf (" Matrix Dimension: %i\n",n);
  
  double A[n*n], B[n*n], C[n*n], C_[n][n];/*Matrices*/

  
  /*reading the file*/


  if (InFile==NULL) {fputs ("File error",stderr);}
    else
    {
    	/*reading A*/
    	for (i = 0; i < n*n; i++) {
    		fscanf (InFile, "%lf", &e);
    		A[i] = e;
    		printf (" %lf\n",A[i]);
    	}
    	/*reading B*/
    	for (i = 0; i < n*n; i++) {
    		fscanf (InFile, "%lf", &e);
    		B[i] = e;
    		printf (" %lf\n",B[i]);
    	}
    	/*reading C*/
    	for (i = 0; i < n*n; i++) {
    		fscanf (InFile, "%lf", &e);
    		C[i] = e;
    		printf (" %lf\n",C[i]);
    	}
    }
  
  
  /* Setup PAPI library and begin collecting data from the counters */
  if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
    test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
  
  /* Matrix Matrix Multiplication */
  
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 1.0, C, n);

  /* Collect the data into the variables passed in */
  if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
    test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
  
  printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
  real_time, proc_time, flpins, mflops);
  printf("%s\tPASSED\n", __FILE__);
  fprintf(Result2,"%i \t %f\n",n, mflops);
  PAPI_shutdown();

  /* Saving the results */
  
  fprintf (ATLAS,"%i\n", n);
  
  int count = 0;
  for (i = 0; i < n ; i++) {
	  for (j = 0; j < n ; j++){
		  printf("i,j: %i %i", i, j);
		  printf("; C[i,j]: %lf \n", C[count]);
		  C_[i][j] = C[count];
		  fprintf (ATLAS, " %5.20f", C_[i][j]);
		  count = count +1;
	  }
	  fprintf (ATLAS, "\n");
	  }

  
  
  }/* Num_Tot_m loop ends here*/
  
  return 0;
}

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
