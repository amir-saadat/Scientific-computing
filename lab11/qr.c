/*******************************************************************************/
 #include <stdlib.h>
 #include <stdio.h>
 #include <cblas.h>
 #include <math.h>
 
/* Parameters */
 #define M 3000 /*number of rows of A*/
 #define N 32 /*number of cols of A*/
 #define LDA M
 #define LDU M
 #define LDVT N

 #define min(a,b) ((a)>(b)?(b):(a))


/* DGESVD prototype */
 extern void dgesvd_( char* jobu, char* jobvt, int* m, int* n, double* a,
                 int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
                 double* work, int* lwork, int* info );


/* routine for printing a matrix */
 void print_matrix( char* desc, int m, int n, double* a, int lda ) {
         int i, j;
         printf( "\n %s\n", desc );
         for( i = 0; i < m; i++ ) {
                 for( j = 0; j < n; j++ ) printf( " %6.4f", a[i+j*lda] );
                 printf( "\n" );
         }
 }

/* routines for finding the max and min number in array*/
 double maximum(double *a, int n) { 
	    int i; 
	    double Max = a[0]; 
		for(i = 1; i < n; i++) 
			if(a[i] > Max) { 
				Max = a[i]; 
				} 
		return Max; 
 } 
 
 double minimum(double *a, int n) {
	    int i; 
	    double Min = a[0]; 
		for(i = 1; i < n; i++) 
			if(a[i] < Min) { 
				Min = a[i]; 
				} 
		return Min; 
 } 


 void chol_qr_it(double *a){
         /* Locals */
         int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info, info_G,info_SQVT,info_RSQVT; 
		 int lwork, lwork_G, lwork_SQVT;
         double wkopt, wkopt_G;
         double *work, *work_G, *work_SQVT;
         /* Local arrays */
         //double s[min(M,N)], u[LDU*M], vt[LDVT*N];
         
		 /*parameters before while loop*/
		 double *Q, *G, *R;
		 Q=(double *)calloc(LDA*N,sizeof(double));
		 G=(double *)calloc(N*N,sizeof(double));
		 R=(double *)calloc(N*N,sizeof(double));

		 /*parameters for decomposition of G*/
		 int ldu_g = N, LDU_G = N;
		 int ldvt_g = N, LDVT_G = N; 
		 double *s_G , *SQ_s_G , *u_G, *vt_G;

		 s_G =(double *)calloc(N,sizeof(double));
		 SQ_s_G =(double *)calloc(N*N,sizeof(double));
		 u_G = (double *)calloc(LDU_G*N,sizeof(double));
		 vt_G = (double *)calloc(LDVT_G*N,sizeof(double));

		 /*parameters for QR factorization*/
		 double *SQVT, *TAU, *R_SQVT, *RR_SQVT, cond_s_G, *QR;
		 SQVT =(double *)calloc(N*N,sizeof(double));
		 TAU =(double *)calloc(N,sizeof(double));
		 R_SQVT = (double *)calloc(N*N,sizeof(double));
		 RR_SQVT = (double *)calloc(N*N,sizeof(double));
		 QR = (double *)calloc(M*N,sizeof(double));

		 /*parameters involved after inv until end of the function*/
		 double *QQ;
		 QQ=(double *)calloc(LDA*N,sizeof(double));
		 double *I;
		 I=(double *)calloc(N*N,sizeof(double));


		 /* Executable statements */
         //printf( " DGESVD Example Program Results\n" );
         /* Query and allocate the optimal workspace */
         lwork = -1; lwork_G = -1, lwork_SQVT = -1;

		 int i = 0,j,k;
		 double  cn = 200.0;
		 
		 /* Operations before while loop */
		  for (j=0;j<LDA*N;j++)
			 Q[j]=a[j];

		 cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, N, M, 1.0, Q, M, Q, M, 0.0,G,N);
         
		 /* Print G = Q'Q */
		 int ldg = N;
         //print_matrix( "G = Q'Q:", N, N, G,ldg );
		
		 for (j=0;j<N;j++){
			 R[j+j*N]=1.0;
			 I[j+j*N]=1.0;}	
         //print_matrix( "R1:", N, N, R, N );
         //print_matrix( "I:", N, N, R, N );
		 
		 /*allocation of work_G & work_SQVT*/
		 dgesvd_( "All", "All", &n, &n, G, &ldg, s_G, u_G, &ldu_g, vt_G, &ldvt_g, &wkopt_G, &lwork_G,
					&info_G );
		 lwork_G = (int)wkopt_G;
		 lwork_SQVT=lwork_G;
		 work_G = (double*)malloc( lwork_G*sizeof(double) );
		 work_SQVT = (double*)malloc(lwork_SQVT*sizeof(double) );
		

		
		 while (cn > 100) {

			i = i+1;
			printf("iteration step: %d\n",i);
			
			/*SVD of G*/
		//	print_matrix("G",N,N,G,N);
			dgesvd_( "All", "All", &n, &n, G, &ldg, s_G, u_G, &ldu_g, vt_G, &ldvt_g, work_G, &lwork_G,
				&info_G );
			/* Check for convergence */
			if( info_G > 0 ) {
                 printf( "The algorithm computing SVD failed to converge.\n" );
                 exit( 1 );
					}
			for (j=0;j<N;j++)
				SQ_s_G[j+j*N]=sqrt(s_G[j]);
		//	print_matrix( "SQ:", N, N, SQ_s_G, N );

        //    print_matrix( "Left singular vectors (stored columnwise)", n, n, u_G, ldu_g );
		//	  print_matrix( "Right singular vectors (stored rowwise)", n, n, vt_G, ldvt_g );

			
			/*getting matrix sqrt(s)*v'*/
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, SQ_s_G, N, vt_G, N, 0.0,SQVT,N);	
	    //	print_matrix( "SQVT:", N, N, SQVT, N );

			dgeqrf_(&n , &n , SQVT , &n , TAU , work_SQVT , &lwork_SQVT , &info_SQVT);
	    //	print_matrix( "R in QR:", N, N, SQVT, N );
			
			for (j=0;j<N;j++){
				for (k=0;k<N;k++){
					if (j<=k)
				       R_SQVT[j+k*N] = SQVT[j+k*N]; }}
	   //	print_matrix( "R_SQVT:", N, N, R_SQVT, N );

			/* R = r * R */
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, R_SQVT, N, R , N, 0.0,RR_SQVT,N);
	  //	print_matrix( "RR_SQVT:", N, N, RR_SQVT, N );
			
			for (j=0;j<N*N;j++)
				R[j]=RR_SQVT[j];
			
			cond_s_G = maximum(s_G,N)/minimum(s_G,N);
			printf("cond(s): %lf\n",cond_s_G);

			cn = sqrt(cond_s_G);
			
			/* inverse of a triangular matrix*/
			dtrtri_("U" , "N" , &n, R_SQVT ,&n , &info_RSQVT);
			
			if( info_RSQVT > 0 ) {
               printf( "The inverse of r is not possible in iteration %d.\n",i);
               exit( 1 );
				}
	   //  	print_matrix( "inv(r):", N, N, R_SQVT, N );
			
			/* Q = Q * inv(r)*/
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, N, 1.0, Q, M, R_SQVT, N, 0.0,QQ,M);
			
			for (j=0;j<M*N;j++)
				Q[j] = QQ[j];
	  //	print_matrix( "Q new:", M, N, Q, M );

			if ((cn > 100) || (i==1))
				cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, N, M, 1.0, Q, M, Q, M, 0.0,G,N);


		 }/*end of while*/
         
		 //print_matrix("final G: ", N, N, G, N); 
		 /* Print singular values */
         //print_matrix( "Singular values", 1, n, s_G, 1 );
         /* Print left singular vectors */
         //print_matrix( "Left singular vectors (stored columnwise)", n, n, u_G, ldu_g );
         /* Print right singular vectors */
         //print_matrix( "Right singular vectors (stored rowwise)", n, n, vt_G, ldvt_g );
         
		 
		 double Norm1;
		 for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++){
				Norm1 += (G[i*N+j]-I[i*N+j])*(G[i*N+j]-I[i*N+j]);
				}
			}
		 Norm1 = sqrt(Norm1);
		 printf("Norm1: %le\n",Norm1);
		 

		 cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, N, 1.0, Q, M, R , N, 0.0,QR,M);
		 
		 double Norm2;
		 for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++){
				Norm2 += (a[i*N+j]-QR[i*N+j])*(a[i*N+j]-QR[i*N+j]);
				}
			}
		 Norm2 = sqrt(Norm2);
		 printf("Norm2: %le\n",Norm2);
		 
		 
		 /* Free workspace */
         free((double *)work_G);
         exit( 0 );
 }
 

int main() {
		 int i,j;
		 double *a_orig;
		 
		 a_orig = malloc(N*M*sizeof(double));
         
	/*	 double a[] ={
		 
             8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
             9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
             9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
             5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
             3.16,  7.98,  3.01,  5.80,  4.27, -5.31
         };*/



		 for ( i=0; i<M*N; i++ )
		     a_orig[i] = (double)rand() / (double)RAND_MAX;
		 
	
		 chol_qr_it(a_orig);


		 free(a_orig);
		 return 0;
 }
 
