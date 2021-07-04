// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
// includes, project
#include <cutil.h>
#include <cublas.h>
//=============================================================================
extern "C" void  cblas_saxpy(const int, const float, const float *, const int, 
                             const float *, const int);
extern "C" float cblas_snrm2(const int, const float *, const int);
extern "C" float cblas_isamax(const int, const float *, const int);
extern "C" void      sgeqrf_(int*,int*,float*,int*,float*,float*,int*,int*);
extern "C" int strmm_(char*, char *, char*, char *, int *, int *, float *,
                      float *, int *, float *, int *);
extern "C" int sgemm_(char *, char *, int *, int *, int *, float *, float *,
                      int *, float *, int *, float *, float *, int *);

void chol_qr_it(int m, int n, float *A, int lda, float *R);
void chol_qr_it_GPU(int m, int n, float *d_A, int lda, float *d_G, float *R, 
                 float *h_work, int lwork);
//=============================================================================

///////////////////////////////////////////////////////////////////////////////
// Program main
///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char** argv) 
{
    CUT_DEVICE_INIT(argc, argv);

    unsigned int timer = 0;

    /* Matrix size */
    int N, M;                 // NxM would be the size of the matrices 
                              // (M columns) that we would orthogonalize
    float *d_A, *d_G;         // d_A is array for A on the device (GPU)
    float *h_work, *h_tau;    // work space and array tau on the host
    float *h_A, *h_Q1, *h_Q2; // These would be the same NxM matrices 
    float *h_R, *h_G;

    int info[1], lwork, i;

    N  = 131072;
    M  = 128;    

    if (argc != 1)
    for(i = 1; i<argc; i++){	
      if (strcmp("-N", argv[i])==0)
         N = atoi(argv[++i]);
      else if (strcmp("-M", argv[i])==0)
         M = atoi(argv[++i]);
    }
    printf("\nUsage: \n");
    printf("  chol_qr_it -N %d -M %d\n\n", N, M);

    lwork = 2*N*M;

    int n2 = N * M;

    /* Initialize CUBLAS */
    cublasInit();

    /* Allocate host memory for the matrix */
    h_A  = (float*)malloc(n2 * sizeof( h_A[0]));
    h_Q1 = (float*)malloc(n2 * sizeof(h_Q1[0]));
    h_Q2 = (float*)malloc(n2 * sizeof(h_Q2[0]));
   
    h_G = (float*)malloc(M*M * sizeof(h_G[0]));
    h_R = (float*)malloc(M*M * sizeof(h_R[0]));
  
    CUDA_SAFE_CALL( cudaMallocHost( (void**)&h_work, lwork*4) );
  
    h_tau = (float*)malloc(N * sizeof(h_tau[0]));
   
    /* Take a random matrix h_A = h_Q1 = h_Q2 */
    for (i = 0; i < n2; i++) {
        h_A[i] = h_Q1[i] = h_Q2[i] = rand() / (float)RAND_MAX;
    }

    /* Allocate device memory for the matrices */
    cublasAlloc(n2, sizeof(d_A[0]), (void**)&d_A);
    cublasAlloc(M*M, sizeof(d_G[0]), (void**)&d_G);

    // create and start timer
    timer = 0;
    CUT_SAFE_CALL(cutCreateTimer(&timer));
    CUT_SAFE_CALL(cutStartTimer(timer));

    /* =====================================================================
         Performs QR on CPU using LAPACK 
       ===================================================================== */
    sgeqrf_(&N, &M, h_A, &N, h_tau, h_work, &lwork, info);
    if (info[0] < 0)  
       printf("Argument %d of sgeqrf had an illegal value.\n", -info[0]);     

    // stop and destroy timer
    CUT_SAFE_CALL(cutStopTimer(timer));
    printf("CPU Processing time: %f (ms) \n", cutGetTimerValue(timer));
    printf("Speed: %f GFlops \n", 4.*N*M*M/
           (3.*1000000*cutGetTimerValue(timer)));
    CUT_SAFE_CALL(cutDeleteTimer(timer));


    /* Initialize the device matrix with the host matrices */
    cublasSetVector(n2, sizeof(h_Q2[0]), h_Q2, 1, d_A, 1);

    timer = 0;
    CUT_SAFE_CALL(cutCreateTimer(&timer));
    CUT_SAFE_CALL(cutStartTimer(timer));

    /* =====================================================================
         Performs orthogonalization on CPU using chol_qr_it
       ===================================================================== */
    chol_qr_it(N, M, h_Q2, N, h_R);

    // stop and destroy timer
    CUT_SAFE_CALL(cutStopTimer(timer));
    printf("\n\nCPU Processing time: %f (ms) \n", cutGetTimerValue(timer));
    printf("Speed: %f GFlops \n", 4.*N*M*M/
           (3.*1000000*cutGetTimerValue(timer)));
    CUT_SAFE_CALL(cutDeleteTimer(timer));

    float one = 1.f, zero = 0.f;
    sgemm_("t", "n", &M, &M, &N, &one, h_Q2, &N, h_Q2, &N, &zero, h_G, &M);
    for(i=0; i<M*M; i+=(M+1)) h_G[i] -= one;
    printf(" ||I - Q'Q||_F = %e, ||I-Q'Q||_max = %e \n",
            cblas_snrm2(M*M, h_G, 1), cblas_isamax(M*M, h_G, 1));

    strmm_("r", "u", "n", "n", &N, &M, &one, h_R, &M, h_Q2, &N);
    cblas_saxpy(n2, -1.0f, h_Q1, 1, h_Q2, 1);
    printf(" ||A - Q R||_F = %e \n",
            cblas_snrm2(n2, h_Q2, 1));    

    // chol_qr on GPU
    timer = 0;
    CUT_SAFE_CALL(cutCreateTimer(&timer));
    CUT_SAFE_CALL(cutStartTimer(timer));

    /* =====================================================================
         Performs orthogonalization on CPU-GPU using chol_qr_it
       ===================================================================== */
    chol_qr_it_GPU(N, M, d_A, N, d_G, h_R, h_work, lwork);

    // stop and destroy timer
    CUT_SAFE_CALL(cutStopTimer(timer));
    printf("\n\nGPU Processing time: %f (ms) \n", cutGetTimerValue(timer));
    printf("Speed: %f GFlops \n", 4.*N*M*M/
           (3.*1000000*cutGetTimerValue(timer)));
    CUT_SAFE_CALL(cutDeleteTimer(timer));

    /* Read the result back */
    cublasGetVector(n2, sizeof(h_Q2[0]), d_A, 1, h_Q2, 1);

    sgemm_("t", "n", &M, &M, &N, &one, h_Q2, &N, h_Q2, &N, &zero, h_G, &M);
    for(i=0; i<M*M; i+=(M+1)) h_G[i] -= one;
    printf(" ||I - Q'Q||_F = %e, ||I-Q'Q||_max = %e \n",
            cblas_snrm2(M*M, h_G, 1), cblas_isamax(M*M, h_G, 1));

    strmm_("r", "u", "n", "n", &N, &M, &one, h_R, &M, h_Q2, &N);
    cblas_saxpy(n2, -1.0f, h_Q1, 1, h_Q2, 1);
    printf(" ||A - Q R||_F = %e \n",
            cblas_snrm2(n2, h_Q2, 1));

    /* Memory clean up */
    free(h_A);
    free(h_Q1);
    free(h_Q2);
    free(h_R);
    free(h_G);
    CUDA_SAFE_CALL( cudaFreeHost(h_work) );
    free(h_tau);

    cublasFree(d_G);
    cublasFree(d_A);

    /* Shutdown */
    cublasShutdown();

    CUT_EXIT(argc, argv);
}
