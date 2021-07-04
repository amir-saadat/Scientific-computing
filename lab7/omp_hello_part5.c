/******************************************************************************
* FILE: omp_hello.c
* DESCRIPTION:
*   OpenMP Example - Hello World - C/C++ Version
*   In this simple example, the master thread forks a parallel region.
*   All threads in the team obtain their unique thread number and print it.
*   The master thread only prints the total number of threads.  Two OpenMP
*   library routines are used to obtain the number of threads and each
*   thread's number.
* AUTHOR: Blaise Barney  5/99
* LAST REVISED: 04/06/05
******************************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) 
{
	int nthreads, tid, procs, maxt, inpar, dynamic, nested;

    /* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthreads, tid)
    {

    /* Obtain thread number */
    tid = omp_get_thread_num();


	/* Printing Info from Master thread */
	if (tid == 0) 
		 {
			nthreads = omp_get_num_threads();
			printf("Number of threads = %d\n", nthreads);

			printf("Printing Info From Master Thread (Thread %d)\n", tid);

			/* Get environment information */
			printf("Number of processors available = %d\n", omp_get_num_procs());
			printf("Number of threads being used = %d\n", omp_get_num_threads());
			printf("Max number of threads available = %d\n", omp_get_max_threads());
			if (omp_in_parallel()!=0){
				printf("You are in parallel region\n");
			}
			else{
				printf("You are not in parallel region\n");
			}

			if (omp_get_dynamic()!=0){
				printf("Dynamic threads are enabled\n");
			}
			else{
				printf("Dynamic threads are not enabled\n");
			}

			if (omp_get_nested()!=0){
				printf("Nested parallelism is supported\n");
			}
			else{
				printf("Nested parallelism is not supported\n");
			}


	  }
	
	
/*	printf("Hello World from thread = %d\n", tid);*/

	}

  /* Only master thread does this */

  }  /* All threads join master thread and disband */
