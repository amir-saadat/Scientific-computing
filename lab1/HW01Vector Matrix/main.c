#include <stdio.h>
#include <stdlib.h>
#include "c_timer.h"
#include <math.h>

int main(int argc, char** argv) {
  int i,j,k;
  int Num_Tot_n = 50;
  FILE * InFile;
  FILE * Result1;
  FILE * Result2;
  FILE * Residual;
  FILE * MATLAB;
  double e, f;
  int dim1,dim;
  
  InFile = fopen ( "myfile.txt" , "rb" );
  Result1= fopen ( "result1.txt" , "w+" );
  Result2= fopen ( "result2.txt" , "w+" );
  Residual = fopen ( "residual.txt" , "w+" );
  MATLAB = fopen ( "matlab.txt" , "rb" );
  
  for (k= 0; k<Num_Tot_n; k++)
  {
  fscanf (InFile, "%i", &dim);
  int n = dim;
  printf (" Matrix Dimension: %i\n",n);
  
  double y[n], x[n];/*vectors*/
  double A[n][n]; /*Matrix*/
  
  /*reading the file*/


  if (InFile==NULL) {fputs ("File error",stderr);}
    else
    {
    	/*reading y*/
    	for (i = 0; i < n; i++) {
    		fscanf (InFile, "%lf", &e);
    		y[i] = e;
    		printf (" %lf\n",y[i]);
    	}
    	/*reading x*/
    	for (i = 0; i < n; i++) {
    		fscanf (InFile, "%lf", &e);
    		x[i] = e;
    		printf (" %lf\n",x[i]);
    	}
    	/*reading A*/
    	for (i = 0; i < n; i++) {
    		for (j = 0; j < n; j++){
    		fscanf (InFile, "%lf", &e);
    		A[i][j] = e;
    		printf (" %lf\n",A[i][j]);
    	}
    	}
    }
  
  
  double btime, etime;
//  double Norm;
  int Nop = 2*n*n;/* number of operations*/
  
  /* Matrix Vector Multiplication */
  
  fprintf (Result1, "n = %i\n",n);
  
  btime = get_cur_time();
  
  for (i = 0; i < n ; i++) {
	  for (j = 0; j < n ; j++){
	  	y[i] = y[i] + x[j]*A[i][j];
	  }
	  }
  
  etime = get_cur_time();
  
  for (i = 0; i < n ; i++) {
	  fprintf (Result1, "%lf\n",y[i]);
	  printf("i: %i",i);
	  printf("; y[i]: %lf \n", y[i]);
	  }
  
  /*reporting the residual*/
   double Normy, Normymy_;
   double r;/*residual*/
   double eps = 1.0e-15;  
   double ymatlab[n];
   
   fscanf (MATLAB, "%i", &dim1);
   if (dim1==n){
   	for (i = 0; i < n; i++) {
   		fscanf (MATLAB, "%lf", &f);
   		ymatlab[i] = f;
  	   printf("ymatlab[i]: %lf \n", ymatlab[i]);
   	}
 	
   	Normy = 0.0; Normymy_ = 0.0;
   	
    for (i = 0; i < n; i++) 
  	  {
  	  	Normy += y[i]*y[i];
  	    Normymy_ += (y[i]-ymatlab[i])*(y[i]-ymatlab[i]);
  	  }

    Normy = sqrt(Normy);
    Normymy_ = sqrt(Normymy_);

     }
   r = Normymy_/Normy/eps;
   
   fprintf(Residual,"%i ", n);
   fprintf(Residual,"%5.20f \n", r);
   
   
  
  printf("Elapsed time: %f seconds\n", etime-btime);
//  printf("2nd Norm: %f \n", Norm);
  fprintf (Result2, "%i %f\n",n, Nop/(etime-btime));
  
  }/* Num_Tot_m loop ends here*/
  
  return 0;
}
