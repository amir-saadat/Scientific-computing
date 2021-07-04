#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define dimension 50000
#define Nthreads 5
#define chunk 10

int main (int argc, char *argv[])
{
	int **a,*b,*c,i,j,n;
	n=dimension;

	b = (int*)malloc(n*sizeof(int));
	c = (int*)malloc(n*sizeof(int));
	a = (int**)malloc(n*sizeof(int*));

	srand((unsigned)time(NULL));


	for (i=0;i<n;i++)
		a[i]=(int*)malloc((i+1)*sizeof(int));
	for (i=0;i<n;i++)
		c[i]=0;
	for (i=0;i<n;i++){
		for (j=0;j<=i;j++)
			a[i][j] = (rand() %10);
		b[i] = (rand() %10);
	}

	int tid,MaxNt;	
	MaxNt = omp_get_max_threads();

	if (Nthreads>MaxNt){
		printf("Warning");
		printf("Number of requested threads are greater than max number of available threads\n");
		printf("Requested %d, Maximum %d\n",Nthreads,MaxNt);
	}

	omp_set_num_threads(Nthreads);

	double start, end;
	start = omp_get_wtime();

	#pragma omp parallel shared(n,a,b,c) private(tid,i,j)
	
	{
	tid = omp_get_thread_num();

	
	/*** Matrix Vector Multiplication ***/

	#pragma omp for schedule (dynamic, chunk)
	
	for (i=0; i<n; i++){
		printf("Thread %d did row %d\n",tid,i);
		for (j=0; j<=i; j++)
			c[i]+= a[i][j]*b[j];
	}
	}
	
	/*** End of parallel region ***/

	end = omp_get_wtime();

	printf("Time Passed for doing Parallel Region:%lf\n",(end-start));
	
	char release;
	printf("\n");
	printf("Do you want to see the Matrices a,b and c: (y or n)");
	release = fgetc(stdin);
	if (release=='n'){
		return 0;
	}
	else{
	
	printf("******************************************************\n");
	printf("Matrix a (Matrix is Lower Triangular):\n\n");
	
	for (i=0; i<n; i++){
		for (j=0;j<=i;j++)
			printf("%d ", a[i][j]);
		printf("\n");
	}


	printf("Vector b:\n\n");
	
	for (j=0; j<n; j++)
		printf("%d\n", b[j]);

	
	/*** Resulted Matrix ***/
	printf("Result Vector (Vector c):\n\n");
	
	for (j=0; j<n; j++)
		printf("%d\n", c[j]);

	printf("******************************************************\n");

	
	return 0;
	}
}

