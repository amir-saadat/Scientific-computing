#include  "mpi.h"
#include  <stdio.h>
#include <stdlib.h>
#include <math.h>


int MemAlloc2d(float ***array, int n, int m) { 
 
	/*the goal is to  allocate the n*m contiguous items */ 
        float *p = (float *)malloc(n*m*sizeof(float)); 
        if (!p) return -1; 
 
      /* allocate the row pointers into the memory */ 
        (*array) = (float **)malloc(n*sizeof(float*)); 
        if (!(*array)) { 
           free(p); 
           return -1; 
        		 } 
 
        /* set up the pointers into the contiguous memory */ 
        int i;
	for (i=0; i<n; i++) 
        (*array)[i] = &(p[i*m]); 
 
        return 0; 
} 
 
int MemFree2d(float ***array) { 
	/* free the memory - because the first element of the array is at the start of whole array*/ 
        free(&((*array)[0][0])); 
 
        /* free the pointers into the memory */ 
        free(*array); 
 
        return 0; 
} 
int main(int  argc,  char  *argv[])
{


	int  i, j, k,  myrank, size;
	MPI_Status status;
	MPI_Init(&argc,  &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,  &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	FILE * infile, * outfile, * data;
	float e;
	
	int n, Pr, Pc, nLocRow, nLocCol,nLocRowal, nLocColal,nLocRowaltr, nLocColaltr;
	/*al is for last row and last column processors which have not same #arguments as the others*/
	data=fopen("data.dat","r");
	fscanf(data,"%i\t%i\t%i\t%i\t%i\t ", &n, &Pr, &Pc, &nLocRow, &nLocCol);
	nLocRowal=nLocRow;
	nLocColal=nLocCol;
	nLocRowaltr=nLocCol;
	nLocColaltr=nLocRow;

	infile=fopen("myfile.dat" , "rb" );
	float **C,**Cp; 
	
	if (size != Pr*Pc) { 
        	printf("program Only works with total number of processes=%i\n",  Pr*Pc); 
        	MPI_Abort(MPI_COMM_WORLD,1); 
			    } 


    /* Array reading and printing */ 
	if (myrank == 0) { 
	
		MemAlloc2d(&C, n, n);
		MemAlloc2d(&Cp, Pr*nLocRow, Pc*nLocCol);
    	
		for (i = 0; i < Pr*nLocRow; i++){ 
    	 	for (j = 0; j < Pc*nLocCol; j++){
				Cp[i][j]=0.0;
			}
		}

    	for (i = 0; i < n; i++) {
    	 	for (j = 0; j < n; j++){
    			fscanf (infile, "%f", &e);
    			C[i][j] = e;
				Cp[i][j]=C[i][j];
    			printf (" %f ",Cp[i][j]);
			}
    			printf ("\n");
   		}

    	if (infile!=NULL) {printf ("Matrix has been successfully read by process %i\n",myrank);}
        } 

	/*Creating the Local Blocks*/

	/*Allocating Memory for local blocks*/

	if (((myrank+1)%Pc==0)&&((myrank+1)!=Pr*Pc))/*row end*/
		nLocColal=n-((myrank%Pc)*nLocCol);
	if (((myrank)>=(Pr*Pc-Pc))&&((myrank+1)!=Pr*Pc))
		nLocRowal=n-((myrank/Pc)*nLocRow);
	if ((myrank+1)==Pr*Pc){
		nLocRowal=n-((myrank/Pc)*nLocRow);
		nLocColal=n-((myrank%Pc)*nLocCol);
	}


	float **blk;
	MemAlloc2d(&blk, nLocRow, nLocCol); 

	
	/* Create a datatype to describe the subarrays of the global array */ 
 
	int Csize[2]    = {Pr*nLocRow, Pc*nLocCol};         
    int blocksize[2] = {nLocRow, nLocCol};       
    MPI_Datatype type; 
	int startsAt[2]= {(myrank/Pc)*nLocRow,(myrank%Pc)*nLocCol};
    MPI_Type_create_subarray(2, Csize, blocksize, startsAt, MPI_ORDER_C, MPI_FLOAT, &type); 
    MPI_Type_commit(&type);

/*    MPI_Datatype type; 
	MPI_Type_vector(nLocRow,nLocCol,n,MPI_FLOAT,&type);
	MPI_Type_commit(&type);*/

	/*Setting the block for zero processor*/
	/*Sending the extracted block to other processors*/
	if (myrank==0) {
		for (i=0;i<nLocRow;i++){
		 	for (j=0;j<nLocCol;j++){
				blk[i][j]  =  Cp[i+(myrank/Pc)*nLocRow][j+(myrank%Pc)*nLocCol];
				}
			}
	for (i=1;i<Pr*Pc;i++)
		MPI_Send(&(Cp[(i/Pc)*nLocRow][(i%Pc)*nLocCol]), 1, type, i, 100, MPI_COMM_WORLD);	
	}

	/*Recieving the corresponding blocks by all processors from root processor(0)*/
	if (myrank!=0){
		MPI_Recv(&(blk[0][0]), nLocRow*nLocCol, MPI_FLOAT , 0, 100, MPI_COMM_WORLD , &status);
	}


    /*Printing out the blocks result*/
    char filename[128];
    sprintf(filename,"block.%i.dat",myrank);
    outfile=fopen(filename, "w");

	for (k=0;k<Pr*Pc;k++){
		if (myrank==k){
			printf("Process %i's block\n",myrank);
			for (i=0;i<nLocRowal;i++){
				for (j=0;j<nLocColal;j++){
					fprintf(outfile,"%f\t",blk[i][j]);
					printf("%f\t",blk[i][j]);
								}
				fprintf(outfile,"\n");
				printf("\n");
							}
			}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	/*partitioning the first result*/
	if (myrank==0)
		printf("\n\n");
	sleep(1);

	/* we don't need "type" anymore*/
	MPI_Type_free(&type);
	
	/*Transposing Matrix*/
	float **blk_tr,**blk_tror;
	MemAlloc2d(&blk_tr, nLocCol, nLocRow); 
	MemAlloc2d(&blk_tror, nLocCol, nLocRow); 
	
	MPI_Datatype col, transpose;
	MPI_Status tr_status,tr_status1;
	MPI_Request req;
	int flag;

	MPI_Type_vector(nLocRow, 1, nLocCol, MPI_FLOAT, &col); 
	MPI_Type_hvector(nLocCol, 1,sizeof(float) , col, &transpose); 
	MPI_Type_commit(&transpose);

	/* just for representing the transpose of each block, I sent it to &blk_tror[0][0]
	 of the rank itself (it is not necessary)*/
	MPI_Sendrecv (&blk[0][0],1,transpose,myrank,101, &blk_tror[0][0],nLocRow*nLocCol,
			MPI_FLOAT,myrank,101,MPI_COMM_WORLD,&tr_status);

	/*now sending and receiving to destination and from source
	 The result would be saved at &blk[0][0]*/
	
	int dest = ((myrank%Pc))*Pr+(myrank/Pc);
	int source = ((myrank%Pr))*Pc+(myrank/Pr);

	MPI_Sendrecv (&blk[0][0],1,transpose,dest,102, &blk_tr[0][0],nLocRow*nLocCol,
			MPI_FLOAT,source,102,MPI_COMM_WORLD,&tr_status);

	
	MPI_Type_free(&transpose);
	MPI_Type_free(&col);
		
	


	/* now all processors print their local transposed block: */ 
	 

	if (((myrank+1)%Pr==0)&&((myrank+1)!=Pr*Pc))/*row end*/
		nLocColaltr=n-((myrank%Pr)*nLocRow);
	if (((myrank)>=(Pr*Pc-Pr))&&((myrank+1)!=Pr*Pc))
		nLocRowaltr=n-((myrank/Pr)*nLocCol);
	if ((myrank+1)==Pr*Pc){
		nLocRowaltr=n-((myrank/Pr)*nLocCol);
		nLocColaltr=n-((myrank%Pr)*nLocRow);
	}


	for (k=0;k<Pr*Pc;k++){
		if (myrank==k){
			printf("Process %i's transposed block\n",myrank);
			for (i=0;i<nLocRowaltr;i++){
				for (j=0;j<nLocColaltr;j++){
					fprintf(outfile,"%f\t",blk_tr[i][j]);
					printf("%f\t",blk_tr[i][j]);
								}
				fprintf(outfile,"\n");
				printf("\n");
							}
			}
		MPI_Barrier(MPI_COMM_WORLD);
	}


	sleep(1);
	if (myrank==0)
		printf("\n\n");

	/*Gathering transposed blocks in root*/
	
	/* create a datatype to describe the subarrays of the global array */ 
	/*I have resized the displacement based on block size, so each Row in 
	  transposed blocks has unit one*/
	
	float **C_tr;
	MemAlloc2d(&C_tr,Pc*nLocCol, Pr*nLocRow);
	
	int sizes[2]    = {Pc*nLocCol, Pr*nLocRow};         /* global size */ 
	int subsizes[2] = {nLocCol,nLocRow};     /* local size */ 
	int starts[2]   = {0,0};                        /* where this one starts */ 
	MPI_Datatype typep, subarrtype; 
	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &typep); 
	MPI_Type_create_resized(typep, 0, nLocRow*sizeof(float), &subarrtype); 
	MPI_Type_commit(&subarrtype); 
								 
	if (myrank == 0){
		float *CTRptr=NULL; 
		CTRptr = &(C_tr[0][0]); 				 
	}

	/* setting the displacements for all processors 
	   based on their position in transposed matrix*/ 
	int sendcounts[Pr*Pc]; 
	int displs[Pr*Pc]; 
												 
	if (myrank == 0) { 
		for (i=0; i<Pr*Pc; i++) sendcounts[i] = 1; 
		int disp = 0; 
		for (i=0; i<Pc; i++) { 
			for (j=0; j<Pr; j++) { 
				displs[i*Pr+j] = disp; 
				disp += 1; 
					          } 
			disp += ((nLocCol)-1)*Pr; 
			       } 
			 } 


	MPI_Gatherv(&(blk_tr[0][0]), nLocRow*nLocCol,  
			MPI_FLOAT, &(C_tr[0][0]), sendcounts, displs, subarrtype, 0, MPI_COMM_WORLD); 

    /*printing out the final transposed matrix*/	
	if (myrank==0){
		printf("Transposed Matrix:\n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++){
				printf (" %f ",C_tr[i][j]);
							}
    		printf ("\n");
					}
		}



	
	MemFree2d(&blk);
	MemFree2d(&blk_tr);MemFree2d(&blk_tror);
	MemFree2d(&C_tr);

	if (myrank==0)
		MemFree2d(&C);
    	 
	
	MPI_Finalize();


    return  0;
}
