#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include <sys/time.h>

static int max( int a, int b ){
        if (a>b) return(a); else return(b);
}

#ifdef F77_WITH_NO_UNDERSCORE
#define   numroc_      numroc
#define   descinit_    descinit
#define   pdlamch_     pdlamch
#define   pdlange_     pdlange
#define   pdlacpy_     pdlacpy
#define   pdgesv_      pdgesv
#define   pdgemm_      pdgemm
#define   indxg2p_     indxg2p
#endif

extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_get( int context, int request, int* value);
extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void   Cblacs_gridexit( int context);
extern void   Cblacs_exit( int error_code);

extern int    numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void   descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
            int *ictxt, int *lld, int *info);
extern double pdlamch_( int *ictxt , char *cmach);
extern double pdlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);

extern void pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca,
            double *b, int *ib, int *jb, int *descb);
extern void pdgesv_( int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv,
            double *B, int *ib, int *jb, int *descb, int *info);
extern void pdgemm_( char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA,
            double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB,
            double * BETA, double * C, int * IC, int * JC, int * DESCC );
extern int  indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);

extern void pdpotrf_( char *uplo, int *n, double *A, int *ia, int *ja, int *desca, int *info);

int main(int argc, char **argv) {
   int iam, nprocs;
   int myrank_mpi, nprocs_mpi;
   int ictxt, nprow, npcol, myrow, mycol;
   int np, nq, n, nb, nqrhs, nrhs;
   int i, j, k, info, itemp, seed;
   int descA[9], descB[9];
   double *A, *Acpy, *B, *X, *R, eps, *work;
        double AnormF, XnormF, RnormF, BnormF, residF;
   int *ippiv;
   int izero=0,ione=1;
   double mone=(-1.0e0),pone=(1.0e0);
/**/
   double MPIt1, MPIt2, MPIelapsed;
/**/
   MPI_Init( &argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
/**/
   n = 5000; nrhs = 1; nprow = 4; npcol = 4; nb = 64;
/**/
   if (nb>n)
      nb = n;
   if (nprow*npcol>nprocs_mpi){
      if (myrank_mpi==0)
         printf(" **** ERROR : we do not have enough processes available to make a p-by-q process grid ***\n");
         printf(" **** Bye-bye                                                                         ***\n");
      MPI_Finalize(); exit(1);
   }
/**/
   Cblacs_pinfo( &iam, &nprocs ) ;
   Cblacs_get( -1, 0, &ictxt );
   Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
   Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
/**/
/*
 * *
 * *     Work only the process in the process grid
 * *
 * */
   if ((myrow>-1)&(mycol>-1)&(myrow<nprow)&(mycol<npcol)) {
/*
 * *
 * *     Compute the size of the local matrices (thanks to numroc)
 * *
 * */ 
      np    = numroc_( &n   , &nb, &myrow, &izero, &nprow );
      nq    = numroc_( &n   , &nb, &mycol, &izero, &npcol );
      nqrhs = numroc_( &nrhs, &nb, &mycol, &izero, &npcol );
/*
 * *
 * *     Allocate and fill the matrices A and B
 * *
 * */ 

      seed = iam*n*(n+nrhs); srand(seed);
/**/      
      A = (double *)calloc(np*nq,sizeof(double)) ;
      if (A==NULL){ printf("error of memory allocation A on proc %dx%d\n",myrow,mycol); exit(0); }
/**/      
      Acpy = (double *)calloc(np*nq,sizeof(double)) ;
      if (Acpy==NULL){ printf("error of memory allocation Acpy on proc %dx%d\n",myrow,mycol); exit(0); }
/**/      
      B = (double *)calloc(np*nq,sizeof(double)) ;
      if (B==NULL){ printf("error of memory allocation B on proc %dx%d\n",myrow,mycol); exit(0); }
/**/      
      X = (double *)calloc(np*nqrhs,sizeof(double)) ;
      if (X==NULL){ printf("error of memory allocation X on proc %dx%d\n",myrow,mycol); exit(0); }
/**/      
      R = (double *)calloc(np*nqrhs,sizeof(double)) ;
      if (R==NULL){ printf("error of memory allocation R on proc %dx%d\n",myrow,mycol); exit(0); }
/**/      
      ippiv = (int *)calloc(np+nb,sizeof(int)) ;
      if (ippiv==NULL){ printf("error of memory allocation IPIV on proc %dx%d\n",myrow,mycol); exit(0); }
/**/      
      k = 0;
      for (i = 0; i < np; i++) {
         for (j = 0; j < nq; j++) {
				A[k] = abs((double) rand()) / ((double) RAND_MAX)  ;
            k++;   
         }
      }
      k = 0;
      for (i = 0; i < np; i++) {
         for (j = 0; j < nq; j++) {
            B[k] = 0.0;//((double) rand()) / ((double) RAND_MAX) - 0.5 ;
            k++;   
         }
      }
/*
 * *
 * *     Initialize the array descriptor for the matrix A and B
 * *
 * */ 
      itemp = max( 1, np );
      descinit_( descA, &n, &n   , &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
      descinit_( descB, &n, &n   , &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
/*
 * *
 * *     Make a copy of A and the rhs for checking purposes
 * */
      //      pdlacpy_( "All", &n, &n   , A, &ione, &ione, descA, Acpy, &ione, &ione, descA );
      //      pdlacpy_( "All", &n, &n   , B, &ione, &ione, descB, X   , &ione, &ione, descB );
/*
 * **********************************************************************
 * *     Call ScaLAPACK PDGESV routine
 * **********************************************************************
 * */
      if( iam==0 ) {
         printf("                                               \n");
         printf("***********************************************\n");
         printf("  Example of ScaLAPACK routine call: (PDPOTRF)  \n");
         printf("***********************************************\n");
         printf("                                               \n");
         printf("\tn = %d\tnrhs = %d\tprocess grid (%d,%d)\t with blocks: %dx%d\n",n,nrhs,nprow,npcol,nb,nb);
         printf("                                               \n");
      }


      pdgemm_( "N", "T", &n, &n, &n, &pone, A, &ione, &ione, descA, A, &ione, &ione, descA,
            &mone,B , &ione, &ione, descB);
/**/
      MPIt1 = MPI_Wtime();
/**/
      pdpotrf_( "L", &n,  B, &ione, &ione, descB, &info );
/**/
      MPIt2 = MPI_Wtime();
      MPIelapsed=MPIt2-MPIt1;
/**/
	  double flop = (n/1e3*n/1e3*n/1e3/3 ) / (double)(MPIelapsed);
	  printf("GFLOPS %lf \n",flop);

      //if( iam==0 ) {
         printf("\ttime MPI     = %f s\n",MPIelapsed);
      //}

/**/
      if( iam==0 ) {
         printf("                                               \n");
         printf("\tINFO code returned by PDPOTRF = %d              \n",info);
         printf("                                               \n");
      }
/*
 * *     Compute residual ||A * X  - B|| / ( ||X|| * ||A|| * eps * N )
 * *     Froebenius norm
 * */
   /*         pdlacpy_( "All", &n, &nrhs, B, &ione, &ione, descB, R   , &ione, &ione, descB );
            eps = pdlamch_( &ictxt, "Epsilon" );
      pdgemm_( "N", "N", &n, &nrhs, &n, &pone, Acpy, &ione, &ione, descA, X, &ione, &ione, descB,
            &mone, R, &ione, &ione, descB);
      AnormF = pdlange_( "F", &n, &n   , A, &ione, &ione, descA, work);
      BnormF = pdlange_( "F", &n, &nrhs, B, &ione, &ione, descB, work);
      XnormF = pdlange_( "F", &n, &nrhs, X, &ione, &ione, descB, work);
      RnormF = pdlange_( "F", &n, &nrhs, R, &ione, &ione, descB, work);
      residF = RnormF / ( AnormF * XnormF * eps * ((double) n));*/
/**/
   /*   if ( iam==0 ){
         printf("                                                        \n");
         printf("\t||A * X  - B||_F / ( ||X||_F * ||A||_F * eps * N ) = %e \n",residF);
         printf("                                                        \n");
         if (residF<10.0e+0) 
            printf("\tThe answer is correct.                            \n");
         else
            printf("\tThe answer is suspicious.                         \n");
      }*/
/**/
      free(A);
  //    free(Acpy);
      free(B);
  //    free(X);
  //    free(ippiv);
      Cblacs_gridexit( 0 );
   }
/*
 * *     Print ending messages
 * */
   if (iam==0){
      printf("END OF TESTS\n");
      printf("***********************************************\n");
      printf("                                               \n");
   }
   MPI_Finalize();
   exit(0);
}

