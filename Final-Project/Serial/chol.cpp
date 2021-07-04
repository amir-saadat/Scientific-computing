// * Row-wise Cholesky Factorization in C++ with Timing Routines
// 
// * Note: This code will overwrite the input matrix! Create a copy if needed.
//
// * Implemented in ANSI C++. 
//   - To compile with g++             : g++ -o chol_row chol_row.cpp -lm
//   - Intriguing performance results  : g++ -O3 -o chol_row chol_row.cpp -lm
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <cblas.h>
#include "papi.h"

using namespace std;

// For gettimeofday().
#include <sys/time.h>

extern "C"{
void dpotrf_(const char* jobu, int* n,  double*a, int* lda,  int* info );
}


// Timing routine.
double tictoc();


//For demonstrating the matrix
void show_matrix(double **A, int n,const char *argu) {
	 for(int i = 0; i < n; i++){
		 for(int j = 0; j < n; j++){
			 if (argu=="lower"){
				if(j <= i)
					//cout << A[i][j] << "\t\t";
					printf("%2.5f ", A[i][j]);
				else
					//cout << "0.000000\t\t";
					printf("%2.5f ",0.0);
					}
			 else{
				 if (argu=="upper"){
					if(j >= i)
						//cout << A[i][j] << "\t\t";
						printf("%2.5f ", A[i][j]);
					else
						//cout << "0.000000\t\t";
						printf("%2.5f ",0.0);
				}
			  }
				  
		 }
		cout << endl;
     }

	 cout<< endl;
	 /*int i,j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%2.5f ", A[i * n + j]);
			printf("\n");
			}*/
}
 
/*printing 1d matrix with column major format*/
void print_matrix(const char* desc, int m, int n, double* a, int lda ) {
	 int i, j;
	 printf( "\n %s\n", desc );
	 for( i = 0; i < m; i++ ) {
         for( j = 0; j < n; j++ ) printf( " %6.4f", a[i+j*lda] );
			printf( "\n" );
					}
}

void chol_row(double **a, int n){

  for(int i = 0; i < n; i++)
    {
      if (a[i][i]<0){
		  cout<<"Error in Row-Wise Cholseky Decomposition:"<<endl;
		  cout<<"Matrix is not positive definite OR dimension of matrix is not set correctly"<<endl;
		  exit(EXIT_FAILURE);
	  }

      // Square root of the diagonal.
	  a[i][i] = sqrt(a[i][i]);

      // Scale the row below. 
      for(int j = (i+1); j < n; j++)
	      a[i][j] = a[i][j]/a[i][i];
	
      // Subtract the matrix formed by row'*row from the
      // minor matrix.
      for(int l = (i+1); l < n; l++)
	     for(int k = l; k < n; k++)
	        a[l][k] = a[l][k] - a[i][l]*a[i][k];
    }
  //cout<<endl<<"The Upper triangular matrix LT in A=L.LT :"<<endl<<endl;
  const char *string = "upper";
  //show_matrix(a,n,string);
}

void chol_col(double **a, int n){

  for(int i = 0; i < n; i++)
    {
      if (a[i][i]<0){
		  cout<<"Error in Column-Wise Cholseky Decomposition:"<<endl;
		  cout<<"Matrix is not positive definite OR dimension of matrix is not set correctly"<<endl;
		  exit(EXIT_FAILURE);
	  }

      // Square root of the diagonal.
	  a[i][i] = sqrt(a[i][i]);

      // Scale the row below. 
      for(int j = (i+1); j < n; j++)
	      a[j][i] = a[j][i]/a[i][i];
	
      // Subtract the matrix formed by row'*row from the
      // minor matrix.
      for(int l = (i+1); l < n; l++)
	     for(int k = l; k < n; k++)
	        a[k][l] = a[k][l] - a[k][i]*a[l][i];
    }
  //cout<<endl<<"The lower triangular matrix L in A=L.LT :"<<endl<<endl;
  const char *string="lower";
  //show_matrix(a,n,string);
}

int main(int argc, char *argv[]){
	
	double flop;
	int i, j, k, l;
	int n, info;
	double t;
	double **a_row , **a_col , *A, *m1, *m;
  
	// Get the dimensions from the command-line if available.
	if(argc > 1)
		n = atoi(argv[1]);
	else
    n = 4;

	// Allocate the memory for the matrix.
	a_row = new double *[n];
	for(i = 0; i < n; i++)
		a_row[i] = new double [n];
    
	a_col = new double *[n];
	for(i = 0; i < n; i++)
		a_col[i] = new double [n];
	
	A = new double [n*n];
	m1 = new double [n*n];
	m = new double [n*n];

	/* double m[]={1.0, 1.0, 1.0, 1.0, 1.0,
			  1.0, 2.0, 3.0, 4.0, 5.0,
	          1.0, 3.0, 6.0, 10.0,15.0,
			  1.0, 4.0, 10.0, 20.0,35.0,
			  1.0, 5.0, 15.0, 35.0, 70.0};*/

	/* double m[]={5.0, 1.2, 0.3, -0.6,
			  1.2, 6.0, -0.4, 0.9,
	          0.3, -0.4, 8.0, 1.7,
			  -0.6, 0.9, 1.7, 10.0};*/

	for (i=0; i<n*n; i++ ){
		m1[i] = (double)rand() / (double)RAND_MAX;}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, m1, n, m1, n, 0.0, m, n);

	cblas_dcopy( n*n , m , 1 , A , 1);

	for (i = 0; i < n; i++) 
		for (k = 0; k < n; k++){
			a_row[i][k] = m[i*n+k];
			a_col[i][k] = m[i*n+k];
		}



	// Initialize the matrix to a tridiagonal system.
	/*for(i = 0; i < n; i++)
		for(j = 0; j < n; j++){
		 if(i == j)
	        a[i][j] = 4;
		 else if(abs(i - j) == 1)
			a[i][j] = -1;
		 else
			a[i][j] = 0;
      }*/

  // Start the timer.
  t = tictoc(); 
  chol_col(a_col,n);
  // Stop the timer.
  t = tictoc() - t;
  flop = (n/1e3*n/1e3*n/1e3/3 ) / (double)(t);
  printf("\nCol_wise Cholesky GFLOPS %lf \n",flop);


  t = tictoc();  
  chol_row(a_row,n);
  t = tictoc() - t;
  flop = (n/1e3*n/1e3*n/1e3/3 ) / (double)(t);
  printf("\nRow_wise Cholesky GFLOPS %lf \n",flop);


  
  t = tictoc();  
  dpotrf_("L", &n , A , &n , &info);
  t = tictoc() - t;
  flop = (n/1e3*n/1e3*n/1e3/3 ) / (double)(t);
  printf("\nLAPACK dpotrf()GFLOPS %lf \n",flop);
  
  
  if( info > 0 ) {
	    printf( "LAPACK factorization could not be completed.\n" );
		exit( 1 );
	}

  for (j=0;j<n;j++){
	  for (k=0;k<n;k++){
		   if (j<k)
			   A[j+k*n] = 0.0; }}

  //print_matrix( "LAPACK L':", n, n, A, n );
  
  // Stop the timer.
  t = tictoc() - t;

  //show_matrix(a,n,"lower");
  
  // Statistical results.
  /*cout << "Flop: " << n*n*n/3.0 << " operations." << endl;
  cout << "Time: " << t << " seconds. " << endl;
  cout << "Flop Rate: " << 1.0e-6*(n*n*n/3.0)/t << " MFlops" << endl;*/

  // Collect garbage.
  for(i = 0; i < n; i++)
    delete [] a_row[i];
  delete [] a_row;

  for(i = 0; i < n; i++)
    delete [] a_col[i];
  delete [] a_col;
  
  return 0;
}

// Timing routine. Returns the number of seconds elapsed as a double.
double tictoc() 
{
  struct timeval tv;

  gettimeofday(&tv, NULL);

  // Add the number of seconds to the number of (re-scaled) microseconds.
  return (tv.tv_sec + 1.0e-6*tv.tv_usec);
}
