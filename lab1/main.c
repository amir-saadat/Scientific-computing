#include <stdio.h>
#include <stdlib.h>
#include "c_timer.h"
#include  <math.h>

int main(int argc, char** argv) {
  int i,j,k,l;
  int Num_Tot_n = 10;/*Number of different n used here*/
  FILE * InFile;
  FILE * Result1;/*Final result of Matrix C*/
  FILE * Result2;/*Execuation rate vs Matrix dimension*/
  FILE * Residual;
  FILE * MATLAB;
  int dim,dim1; /*Matrix Dimension*/
  double e,f;/*Used for reading data from myfile.txt*/
  
  InFile = fopen ( "myfile.txt" , "rb" );
  Result1= fopen ( "result1.txt" , "w+" );
  Result2= fopen ( "result2.txt" , "w+" );
  Residual = fopen ( "residual.txt" , "w+" );
  MATLAB = fopen ( "matlab.txt" , "rb" );
  
  for (l= 0; l<Num_Tot_n; l++)
  {
  fscanf (InFile, "%i", &dim);
  int n = dim;
  printf (" Matrix Dimension: %i\n",n);
  
  double **A, **B, **C;/*Matrices*/
  
  A = malloc(sizeof(double *)*n);
  if (A == NULL) {
	  printf("out of memory\n");
	  exit;
  }
  for (i=0;i<n;i++){
	  A[i] = malloc(n*sizeof(double));
	  if (A[i] == NULL) {
		  printf("out of memory\n");
		  exit;
	  }
  }
  
  B = malloc(sizeof(double *)*n);
  if (B == NULL) {
	  printf("out of memory\n");
	  exit;
  }
  for (i=0;i<n;i++){
	  B[i] = malloc(n*sizeof(double));
	  if (B[i] == NULL) {
		  printf("out of memory\n");
		  exit;
	  }
  }
  
  C = malloc(sizeof(double *)*n);
  if (C == NULL) {
	  printf("out of memory\n");
	  exit;
  }
  for (i=0;i<n;i++){
	  C[i] = malloc(n*sizeof(double));
	  if (C[i] == NULL) {
		  printf("out of memory\n");
		  exit;
	  }
  }
  
  
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
  
  
  double btime, etime;
//  double Norm;
  int Nop = 2*n*n*n;/* number of operations*/
  
  /* Matrix Matrix Multiplication */
  

  btime = get_cur_time();
  
  for (i = 0; i < n ; i++) {
	  for (j = 0; j < n ; j++){
		  for (k = 0; k < n ; k++){
			  C[i][j] = C[i][j] + A[i][k]*B[k][j];
		  }
	  }
	  }
  
  etime = get_cur_time();
  
  for (i = 0; i < n ; i++) {
	  for (j = 0; j < n ; j++){
		  printf("i,j: %i %i", i, j);
		  printf("; C[i,j]: %lf \n", C[i][j]);
		  fprintf (Result1, " %lf", C[i][j]);
	  }
	  fprintf (Result1, "\n");
	  }
  
  /*reporting the residual*/
   double NormC, NormCmC_;
   double r;/*residual*/
   double eps = 1.0e-15;  
   double **Cmatlab;
   
   Cmatlab = malloc(sizeof(double *)*n);
   if (Cmatlab == NULL) {
 	  printf("out of memory\n");
 	  exit;
   }
   for (i=0;i<n;i++){
	   Cmatlab[i] = malloc(n*sizeof(double));
 	  if (Cmatlab[i] == NULL) {
 		  printf("out of memory\n");
 		  exit;
 	  }
   }
   
   fscanf (MATLAB, "%i", &dim1);
   if (dim1==n){
   	for (i = 0; i < n; i++) {
   		for (j = 0; j < n; j++){
   		fscanf (MATLAB, "%lf", &f);
   		Cmatlab[i][j] = f;
  	   printf("Cmatlab[i,j]: %lf \n", Cmatlab[i][j]);
   		}
   	}
 	
   	NormC = 0.0; NormCmC_ = 0.0;
   	
    for (i = 0; i < n; i++) {
  	  for (j = 0; j < n; j++){
  	  	NormC += C[i][j]*C[i][j];
  	    NormCmC_ += (C[i][j]-Cmatlab[i][j])*(C[i][j]-Cmatlab[i][j]);
  	  }
  	  }

    NormC = sqrt(NormC);
    NormCmC_ = sqrt(NormCmC_);
     }
   r = NormCmC_/NormC/eps;
   
   fprintf(Residual,"%i ", n);
   fprintf(Residual,"%5.20f \n", r);
  
  printf("Elapsed time: %f seconds\n", etime-btime);
//  printf("2nd Norm: %f \n", Norm);
  fprintf (Result2, "%i %f\n",n, Nop/(etime-btime));
  
  for(i=0;i<n;i++){
	  free(A[i]);free(B[i]);free(C[i]);free(Cmatlab[i]);
  }
  free(A);free(B);free(C);free(Cmatlab);
  
  }/* Num_Tot_m loop ends here*/
  
  return 0;
}
