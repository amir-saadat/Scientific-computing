#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <memory.h>
#include <malloc.h>
#include "papi.h"
#include "math.h"

/* multiplication routine*/
void matvec(long int *I, long int *J, double *A, double *v,long int Nc,long int nonz, double *result);


/*main program*/
int main(int argc, char *argv[])
{
	
	/*papi variables*/
	float real_time, proc_time, mflops;
	long long flpins=-1;
	int retval;
	
	FILE *Input,*Output;
    long int Nc, Nr, nonz;
	
	Input=fopen("matrix.output","r");
	fscanf(Input,"%ld %ld %ld",&Nc,&Nr,&nonz);
	
	long int size,s,t=0;
	long int *i;
	long int *CCS,*i_re,*j_re;

	double *v,*v_re;
	double *a,*a_re;

	
	/*reading matrix file*/
	Output=fopen("newmatrix.output","w");


	//printf("%p\n",Output);
	printf("Rows:%ld   Columns:%ld   Non-zeros:%ld \n",Nc,Nr,nonz);
	
	//variable for original matrix
	a  = (double*) malloc (nonz*sizeof(double));
	v= (double*) malloc (Nc*sizeof(double));
	i  = (long int*) malloc (nonz*sizeof(long int));
	
	//variables indicating the reordered index and matrix
	a_re = (double*) malloc (nonz*sizeof(double));
	i_re = (long int*) malloc (nonz*sizeof(long int));
	j_re = (long int*) malloc (nonz*sizeof(long int));
	v_re= (double*) malloc (Nc*sizeof(double));
	
	//variable for Compressed Column Storage
	CCS = (long int*) malloc (nonz*sizeof(long int));
	
	size=Nc;
	s=Nc=t;
	long int j;
	
	/********Reading data from matrix.output and converting it to CCS*******/

	while (!feof(Input))
	{
		fscanf(Input,"%lu %lu %lf \n",&i[s],&j,&a[s]);
	
		j_re[s]=(j-1)/8660+3*((j-1)%8660)+1;
		i_re[s]=(i[s]-1)/8660+3*((i[s]-1)%8660)+1;
		a_re[s]=a[s];
		
		/* CCS format construction*/
		if (j==Nc){
		}
		if (j!=Nc){
				Nc = j;
				CCS[t] = s;
				t++;
			}
		s++;
	}
	
	fclose(Input);

	/*printing out the reordered indeces and Matrix*/
	
	long int k;
	for (k=0; k<nonz; k++){
		fprintf(Output,"%lu\t%lu\t%lf\n",i_re[k],j_re[k],a_re[k]);}
	
	/*vect generation (all elements set to 2 for comparing easily with MATLAB)*/
	int l;
	for (l=0; l<Nc; l++){
		v[l] = 2;
		v_re[l] = 2;
	}

	
	double *result;
	result = (double *) malloc(Nc*sizeof(double));

	retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops);
	
	/*operation to be measured*/
	matvec(i,CCS,a,v,Nc,nonz,result);

	retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops);
	printf("For Original Matrix\n");
	printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t%f\n", real_time, proc_time, flpins, mflops);

	
	/*CCS forming*/
	Nc = 0;int y=0;
	int x;
	for (x=0; x<nonz; x++){
		if (j_re[x] != Nc){
				Nc=j_re[x];
				j_re[y]=x;
				y++;
			}
		}

	double *result_re;
	result_re = (double *) malloc(Nc*sizeof(double));
	
	flpins=-1;PAPI_stop_counters(&flpins, 2);
	retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops);
	
	matvec(i_re,j_re,a_re,v_re,y,nonz,result_re);

	retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops);
	printf("For Reordered Matrix\n");
	printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\%f\n", real_time, proc_time, flpins, mflops);

	/*Norm Calculation*/
	
	double norm = 0.0,norm_re = 0.0;

	for (k=0; k<size; k++)	
		norm += result[k]*result[k];
	for (k=0; k<l; k++)	
		norm_re += result_re[k]*result_re[k];
	
	norm = sqrt(norm);
	norm_re = sqrt(norm_re);
	
	fclose(Output);
	
	printf("\n norm = %lf\n norm_reordered = %lf\n",norm,norm_re);
	PAPI_shutdown();

	/*collecting garbage*/
	free(v);free(v_re);
	free(i);free(CCS);
	free(i_re);free(j_re);
	free(a);free(a_re);

	return 0;
}


void matvec(long int *I,long int *J, double *A, double *v, long int Nc, long int nonz, double *result)
{
	int i,j;
	long int start,end;
	for (i=0;i<(Nc-1);i++)
		{
			start = J[i];
			end = J[i+1];
			for (j=start;j<end;j++)
				{
					result[i] += A[j] * v[I[j]-1];
				}
		}
	for (j=J[Nc-1];j<nonz;j++)
		{
			result[Nc-1] += A[j] * v[I[j]-1];
		}
}

