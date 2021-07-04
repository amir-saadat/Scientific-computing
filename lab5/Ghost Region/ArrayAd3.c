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


	int  i, j, k,  myrank, size, tag;
	MPI_File  thefile;
	MPI_Status status;
	MPI_Init(&argc,  &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,  &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	FILE * infile, * outfile, * data;
	float e;
	
	int n, Pr, Pc, nLocRow, nLocCol;
	data=fopen("data.dat","r");
	fscanf(data,"%i\t%i\t%i\t%i\t%i\t ", &n, &Pr, &Pc, &nLocRow, &nLocCol);

	infile=fopen("myfile.dat" , "rb" );
	float **C; 
	
	if (size != Pr*Pc) { 
        	printf("program Only works with total number of processes=%i\n",  Pr*Pc); 
        	MPI_Abort(MPI_COMM_WORLD,1); 
			    } 


    /* Arrey reading and printing */ 
	if (myrank == 0) { 
	MemAlloc2d(&C, n, n); 
    	for (i = 0; i < n; i++) {
    	 	for (j = 0; j < n; j++){
    			fscanf (infile, "%f", &e);
    			C[i][j] = e;
    			printf (" %f ",C[i][j]);
			}
    			printf ("\n");
   		}

    	if (infile!=NULL) {printf ("Matrix has been successfully read by process %i\n",myrank);}
        } 

	/*This part could generate matrix C with random numbers
	 for ( i=0; i<n; i++ ){	
		 for ( i=0; i<n; i++ ){
			  C[i][j] = (float)rand() / (float)RAND_MAX;
				}
		}*/
	
	/*Creating the Local Blocks*/

	/*Allocating Memory for local blocks*/
	
    float **blk;
	MemAlloc2d(&blk, nLocRow, nLocCol); 

	
	/* Create a datatype to describe the subarrays of the global array */ 
 
    int Csize[2]    = {n, n};         /* global size */ 
    int blocksize[2] = {nLocRow, nLocCol};     /* local size */  
    MPI_Datatype type; 
	int startsAt[2]= {(myrank/Pc)*nLocRow,(myrank%Pc)*nLocCol};
    MPI_Type_create_subarray(2, Csize, blocksize, startsAt, MPI_ORDER_C, MPI_FLOAT, &type); 
    MPI_Type_commit(&type); 
	
	/*Setting the block for zero processor*/
	/*Sending the extracted block to other processors*/
	if (myrank==0) {
		for (i=0;i<nLocRow;i++){
		 	for (j=0;j<nLocCol;j++){
				blk[i][j]  =  C[i+(myrank/Pc)*nLocRow][j+(myrank%Pc)*nLocCol];
				}
			}
	for (i=1;i<Pr*Pc;i++)
		MPI_Send(&(C[(i/Pc)*nLocRow][(i%Pc)*nLocCol]), 1, type, i, 100, MPI_COMM_WORLD);		
	
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
			for (i=0;i<nLocRow;i++){
				for (j=0;j<nLocCol;j++){
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
	sleep(1);
	if (myrank==0)
		printf("\n\n");

	/* we don't need "type" anymore*/
	MPI_Type_free(&type);
	
	/*Ghost Region*/
	
	/*Allocating memory for possible Ghost regin at Left(LeG), Right(RiG), Up(UpG), Bottom(BoG)*/
	float **LeG;
	float **RiG;
	float **UpG;
	float **BoG;

	MemAlloc2d(&LeG, nLocRow, 1);
	MemAlloc2d(&RiG, nLocRow, 1);
	MemAlloc2d(&UpG, 1, nLocCol);
	MemAlloc2d(&BoG, 1, nLocCol);

	/*required datatypes and requests*/
	MPI_Datatype LeGtype,RiGtype,UpGtype,BoGtype;
	MPI_Request reqRi,reqLe,reqUp,reqBo;
	MPI_Status statusLeG,statusRiG,statusUpG,statusBoG;
	
	/*sending and receiving ghost regions*/
	
	/*sending*/

	/*to right side neighbor*/
	if (((myrank+1)%Pc)!=0){/*so the processor has right neighbor*/
		MPI_Type_vector(nLocRow,1,nLocCol,MPI_FLOAT,&RiGtype);
		MPI_Type_commit(&RiGtype);
		MPI_Isend(&(blk[0][nLocCol-1]), 1, RiGtype, myrank+1, 101, MPI_COMM_WORLD, &reqRi);
		MPI_Type_free(&RiGtype);
					}

	/*to left side neighbor*/
	if (((myrank+1)%Pc)!=1){/*so the processor has left neighbor*/
		MPI_Type_vector(nLocRow,1,nLocCol,MPI_FLOAT,&LeGtype);
		MPI_Type_commit(&LeGtype);
		MPI_Isend(&(blk[0][0]), 1, LeGtype, myrank-1, 102, MPI_COMM_WORLD, &reqLe);
		MPI_Type_free(&LeGtype);
					}
	
	/*to upper side neighbor*/
	if ((myrank-Pc)>=0){/*so the processor has upper neighbor*/
		MPI_Type_vector(nLocCol,1,1,MPI_FLOAT,&UpGtype);
		MPI_Type_commit(&UpGtype);
		MPI_Isend(&(blk[0][0]), 1, UpGtype, myrank-Pc, 103, MPI_COMM_WORLD, &reqUp);
		MPI_Type_free(&UpGtype);
			}

	/*to bottom side neighbor*/
	if ((myrank+Pc)<(Pr*Pc)){/*so the processor has bottom neighbor*/
		MPI_Type_vector(nLocCol,1,1,MPI_FLOAT,&BoGtype);
		MPI_Type_commit(&BoGtype);
	    MPI_Isend(&(blk[nLocRow-1][0]), 1, BoGtype, myrank+Pc, 104, MPI_COMM_WORLD, &reqBo);
		MPI_Type_free(&BoGtype);
			}

	/*receiving*/

	
	/*from the right neighbor*/
	if (((myrank+1)%Pc)!=0){/*so the processor has right neighbor*/
		MPI_Recv(&(RiG[0][0]), nLocRow, MPI_FLOAT, myrank+1, 102, MPI_COMM_WORLD,&statusRiG);
					}

	/*from the left neighbor*/
	if (((myrank+1)%Pc)!=1){/*so the processor has left neighbor*/
		MPI_Recv(&(LeG[0][0]), nLocRow, MPI_FLOAT, myrank-1, 101, MPI_COMM_WORLD,&statusLeG);
					}			
	
	/*from upper neighbor*/
	if ((myrank-Pc)>=0){/*so the processor has upper neighbor*/
		MPI_Recv(&(UpG[0][0]), nLocCol, MPI_FLOAT , myrank-Pc, 104, MPI_COMM_WORLD,&statusUpG);
				    }			
	/*from bottom neighbor*/
	if ((myrank+Pc)<(Pr*Pc)){/*so the processor has bottom neighbor*/
		MPI_Recv(&(BoG[0][0]), nLocCol, MPI_FLOAT, myrank+Pc, 103, MPI_COMM_WORLD,&statusBoG);
					}


	/* now all processors print their local Ghost Region data: */ 
	 

	for (j=0; j<size; j++) {
		if (myrank == j) {
			printf("\n");
			printf("I am block number %i my Ghost Regions are:\n",myrank);
			if (((myrank+1)%Pc)!=0){
				printf("RIGHT GHOST REGION(block.%i):\n",myrank);
				fprintf(outfile,"RIGHT GHOST REGION(block.%i):\n",myrank);
				for (i=0;i<nLocRow;i++){
					printf("%f\n",RiG[i][0]);
					fprintf(outfile,"%f\n",RiG[i][0]);
								} 
				printf("\n");
				fprintf(outfile,"\n");
						} 

			if (((myrank+1)%Pc)!=1){
				printf("LEFT GHOST REGION(block.%i):\n",myrank);
				fprintf(outfile,"LEFT GHOST REGION(block.%i):\n",myrank);
				for (i=0;i<nLocRow;i++){
					printf("%f\n",LeG[i][0]);
					fprintf(outfile,"%f\n",LeG[i][0]);
								} 
				printf("\n"); 
				fprintf(outfile,"\n"); 
						}  

			if ((myrank-Pc)>=0){
				printf("UPPER GHOST REGION(block.%i):\n",myrank);
				fprintf(outfile,"UPPER GHOST REGION(block.%i):\n",myrank);
				for (i=0;i<nLocCol;i++){
					printf("%f\t",UpG[0][i]);
					fprintf(outfile,"%f\t",UpG[0][i]);
								} 
				printf("\n");
				fprintf(outfile,"\n");
						}

			if ((myrank+Pc)<size){
				printf("BOTTOM GHOST REGION(block.%i):\n",myrank);
				fprintf(outfile,"BOTTOM GHOST REGION(block.%i):\n",myrank);
				for (i=0;i<nLocCol;i++){
					printf("%f\t",BoG[0][i]);
					fprintf(outfile,"%f\t",BoG[0][i]);
								} 
				printf("\n");
				fprintf(outfile,"\n");
						}
					
				}
	    MPI_Barrier(MPI_COMM_WORLD); 
	}




	MemFree2d(&blk);
	if (myrank==0)
	MemFree2d(&C);
	MemFree2d(&RiG);MemFree2d(&LeG);MemFree2d(&UpG);MemFree2d(&BoG);
    	 
	
	MPI_Finalize();


    return  0;
}
