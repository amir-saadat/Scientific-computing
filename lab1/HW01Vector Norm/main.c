#include <stdio.h>
#include <stdlib.h>
#include "c_timer.h"
#include <math.h>

int main(int argc, char** argv) {
  int i,j;
  int Num_Tot_n = 991;
  FILE * InFile;
  FILE * OutFile;
  FILE * Result1;
  FILE * Result2;
  FILE * Residual;
  FILE * NORM;
  int dim,dim1;
  double f,e;
  
  InFile = fopen ( "myfile.txt" , "rb" );
  Result1= fopen ( "result1.txt" , "w+" );/*calculated Norm*/
  Result2= fopen ( "result2.txt" , "w+" );/*Rate vs Dimension*/
  Residual= fopen ( "residual.txt" , "w+" );
  NORM= fopen ( "norm.txt" , "rb" );
  
  for (j= 0; j<Num_Tot_n; j++)
  {
  OutFile = fopen ( "tempfile.txt" , "wb+" );
  fscanf (InFile, "%i", &dim);
  int n = dim;
  printf (" %i\n",n);
  
  double x[n];/*vector*/
  
  /*reading the file*/


  if (InFile==NULL) {fputs ("File error",stderr);}
    else
    {
    	for (i = 0; i < n; i++) 
    	{
    		fscanf (InFile, "%lf", &f);
    		x[i] = f;
    		printf (" %5.20f\n",x[i]);
    		fprintf (OutFile, "%5.20f\n",x[i]);
    	}
    }
  
  
  long size;
  char * buffer;
  size_t result;

  // obtain file size:
  fseek (OutFile , 0 , SEEK_END);
  size = ftell (OutFile);
  printf ("Size of tempfile.txt: %ld bytes.\n",size);
  rewind (OutFile);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*size);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,size,OutFile);
  if (result != size) {fputs ("Reading error",stderr); exit (3);}
     
  /* the whole file is now loaded in the memory buffer. */
  
  // terminate
  fclose (OutFile);
  free (buffer);
  

  double btime, etime;
  double Norm = 0.0;
  int Nop = 2*n;/* number of operations*/
  
  /* Vector 2nd Norm*/
  

  btime = get_cur_time();
  
  for (i = 0; i < n; i++) 
	  {
	  	Norm += x[i]*x[i];
	  }
  Norm = sqrt(Norm);
  
  etime = get_cur_time();
  
  printf("Elapsed time: %f seconds\n", etime-btime);
  printf("2nd Norm: %5.20f \n", Norm);
  fprintf (Result2, "%i %5.20f\n",n, Nop/(etime-btime));
  fprintf(Result1,"n: %i \n", n);
  fprintf(Result1,"2nd Norm: %5.20f \n", Norm);
  
  /*reporting the residual*/
  double no;
  double r;/*residual*/
  double eps = 1.0e-15;
    
  fscanf (NORM, "%i", &dim1);
  if (dim1==n){
	fscanf (NORM, "%lf", &e);
	no = e;
	printf ("MATLAB Norm: %5.20f \n", no);
    }
  r = fabs(Norm-no)/Norm/eps;
  
  fprintf(Residual,"%i ", n);
  fprintf(Residual,"%5.20f \n", r);
  
  
  }/* Num_Tot_m loop ends here*/
  
  return 0;
}
