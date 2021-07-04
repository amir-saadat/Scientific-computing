#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <memory.h>
#include <malloc.h>
#include "papi.h"

 
static void test_fail(char *file, int line, char *call, int retval);

int main(int argc, char** argv) {
	
/*  extern void dummy(void *);*/
  float real_time, proc_time, mflops;
  long long flpins;
  int retval;
  int i,j,k,l;
  int Num_Tot_n = 26;/*Number of different n used here*/
  FILE * InFile;
  FILE * Result1;/*Final result of Matrix C*/
  FILE * Result2;/*Execuation rate vs Matrix dimension*/
  FILE * Residual;
  FILE * ATLAS;
  
  int dim,dim1; /*Matrix Dimension*/
  double e,f;/*Used for reading data from myfile.txt*/
  /*int n;*/ /*Matrices Dimension*/
  InFile = fopen ( "myfile.txt" , "rb" );
  Result1= fopen ( "result1.txt" , "w+" );
  Result2= fopen ( "result2.txt" , "w+" );
  Residual= fopen ( "residual.txt" , "w+" );
  ATLAS= fopen ( "atlas.txt" , "rb" );
  
  for (l= 0; l<Num_Tot_n; l++)
  {
  fscanf (InFile, "%i", &dim);
  int n = dim;
  printf (" Matrix Dimension: %i\n",n);
  
  double A[n][n], B[n][n], C[n][n], C_[n][n];/*Matrices*/

  
  /*reading the file*/


  if (InFile==NULL) {fputs ("File error",stderr);}
    else
    {
    	/*reading A*/
    	for (i = 0; i < n; i++) {
    		for (k = 0; k < n; k++){
    		fscanf (InFile, "%lf", &e);
    		A[i][k] = e;
    		printf (" %lf\n",A[i][k]);
    	}
    	}
    	/*reading B*/
    	for (i = 0; i < n; i++) {
    		for (k = 0; k < n; k++){
    		fscanf (InFile, "%lf", &e);
    		B[i][k] = e;
    		printf (" %lf\n",B[i][k]);
    	}
    	}
    	/*reading C*/
    	for (i = 0; i < n; i++) {
    		for (k = 0; k < n; k++){
    		fscanf (InFile, "%lf", &e);
    		C[i][k] = e;
    		printf (" %lf\n",C[i][k]);
    	}
    	}
    }
  
  int Nop = 2*n*n*n;/* number of operations*/
  
  /* Setup PAPI library and begin collecting data from the counters */
  if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
    test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
  
  /* Matrix Matrix Multiplication */
  
  for(i=0;i<n;i++)
	for(j=0;j<n;j++)
     for(k=0;k<n;k++)
	C[i][j]=C[i][j]+A[i][k]*B[k][j];

  /* Collect the data into the variables passed in */
  if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
    test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
  
  printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
  real_time, proc_time, flpins, mflops);
  printf("%s\tPASSED\n", __FILE__);
  fprintf(Result2,"%i \t %f\n",n, mflops);
  PAPI_shutdown();

  /* Saving the results */
  
  for (i = 0; i < n ; i++) {
	  for (j = 0; j < n ; j++){
		  printf("i,j: %i %i", i+1, j+1);
		  printf("; C[i,j]: %lf \n", C[i][j]);
		  fprintf (Result1, " %lf", C[i][j]);
	  }
	  fprintf (Result1, "\n");
	  }
  
  /*reporting the residual*/
    double NormC, NormCmC_;
    double r;/*residual*/
    double eps = 1.0e-15;  
    double Catlas[n][n];
    
    fscanf (ATLAS, "%i", &dim1);
    if (dim1==n){
    	for (i = 0; i < n; i++) {
    		for (j = 0; j < n; j++){
    		fscanf (ATLAS, "%lf", &f);
    		Catlas[i][j] = f;
   	   printf("Catlas[i,j]: %lf \n", Catlas[i][j]);
    		}
    	}
  	
    	NormC = 0.0; NormCmC_ = 0.0;
    	
     for (i = 0; i < n; i++) {
   	  for (j = 0; j < n; j++){
   	  	NormC += C[i][j]*C[i][j]; 	  	
   	    NormCmC_ += (C[i][j]-Catlas[i][j])*(C[i][j]-Catlas[i][j]);
   	    printf("DC[i,j]: %lf\n", (C[i][j]-Catlas[i][j]));
   	  }
   	  }

     NormC = sqrt(NormC);
     NormCmC_ = sqrt(NormCmC_);
      }
    r = NormCmC_/NormC/eps;
    
    fprintf(Residual,"%i ", n);
    fprintf(Residual,"%5.20f \n", r);
  
  
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
