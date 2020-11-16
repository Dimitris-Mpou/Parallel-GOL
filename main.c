#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "functions.h"

#define STEPS 100

int main(int argc, char const *argv[]){
  int i, j, size, n, rank, blockSize, stepCount, sqrtP, up, down, left, right, upLeft, upRight, downLeft, downRight, localChanged, totalChanged;
  char *recUp, *recDown, *recLeft, *recRight, *sendLeft, *sendRight, recDownLeft, recDownRight, recUpLeft, recUpRight;  // halo points vars
  double start, end;
  char **cur, **next, **temp;
  MPI_Status status[8];
  MPI_Request sendReq[8], recvReq[8];
  MPI_Datatype LineType;
	MPI_Datatype ColType;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  sqrtP = (int) sqrt(size);
  n = atoi(argv[1]);
  if(n % sqrtP != 0){    // Elegxos N kai p (num procceses)
    if(rank == 0)
      printf("Wrong numbers, blocks can't be equal!\nProgramm ending\n");
    MPI_Finalize();
    return 0;
  }

  blockSize = n / sqrtP;
  recUp = malloc((blockSize+1) * sizeof(char));
  recDown = malloc((blockSize+1) * sizeof(char));
  recLeft = malloc((blockSize+2) * sizeof(char));
  recRight = malloc((blockSize+2) * sizeof(char));
  sendLeft = malloc(blockSize * sizeof(char));
  sendRight = malloc(blockSize * sizeof(char));

  // Arxikopoiisi pinakwn
  cur = malloc((blockSize+2) * sizeof(char *));   // blockSize+2 gia ta halo points
  next = malloc((blockSize+2) * sizeof(char *));
  for(i=0; i<(blockSize+2); i++){
    cur[i] = malloc((blockSize+2) * sizeof(char));
    next[i] = malloc((blockSize+2) * sizeof(char));
  }
  if(cur == NULL || next == NULL){
    printf("Not enough space!!!!\n");
  }

  srand(time(0)+rank);
  for(i=0; i<blockSize+2; i++){
    for(j=0; j<blockSize+2; j++){
      switch (rand()%2) { // 0 -> nekros, 1 -> zontanos
        case 0:
          cur[i][j] = '0';
          break;
        case 1:
          cur[i][j] = '1';
          break;
      }
    }
  }
  findNeighbors(rank, sqrtP, size, &up, &down, &left, &right, &upLeft, &upRight, &downLeft, &downRight);

  MPI_Type_contiguous(blockSize+1, MPI_CHAR, &LineType);   // Create datatypes
  MPI_Type_commit(&LineType);
  MPI_Type_vector(blockSize+2, 1, blockSize+2, MPI_CHAR, &ColType);
  MPI_Type_commit(&ColType);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  for(stepCount=0; stepCount<STEPS; stepCount++){   // Big for loop

      // receive halo points
    MPI_Irecv(&recDownRight, 1, MPI_CHAR, downRight, 0, MPI_COMM_WORLD, &recvReq[3]);	//dr
    MPI_Irecv(&recDownLeft, 1, MPI_CHAR, downLeft, 0, MPI_COMM_WORLD, &recvReq[2]);	//dl
    MPI_Irecv(&recUpRight, 1, MPI_CHAR, upRight, 0, MPI_COMM_WORLD, &recvReq[1]);	//ur
    MPI_Irecv(&recUpLeft, 1, MPI_CHAR, upLeft, 0, MPI_COMM_WORLD, &recvReq[0]);	//ul
    MPI_Irecv(recDown, 1, LineType, down, 0, MPI_COMM_WORLD, &recvReq[5]);	//d
    MPI_Irecv(recUp, 1, LineType, up, 0, MPI_COMM_WORLD, &recvReq[4]);	//u
    MPI_Irecv(recRight, blockSize+2, MPI_CHAR, right, 0, MPI_COMM_WORLD, &recvReq[7]);	//r
    MPI_Irecv(recLeft, blockSize+2, MPI_CHAR, left, 0, MPI_COMM_WORLD, &recvReq[6]);	//l

      // Send halo points
    MPI_Isend(&cur[1][1], 1, MPI_CHAR, upLeft, 0, MPI_COMM_WORLD, &sendReq[0]);	//ul
    MPI_Isend(&cur[1][blockSize], 1, MPI_CHAR, upRight, 0, MPI_COMM_WORLD, &sendReq[1]);	//ur
    MPI_Isend(&cur[blockSize][1], 1, MPI_CHAR, downLeft, 0, MPI_COMM_WORLD, &sendReq[2]);	//dl
    MPI_Isend(&cur[blockSize][blockSize], 1, MPI_CHAR, downRight, 0, MPI_COMM_WORLD, &sendReq[3]);	//dr
    MPI_Isend(cur[1], 1, LineType, up, 0, MPI_COMM_WORLD, &sendReq[4]);	//u
    MPI_Isend(cur[blockSize], 1, LineType, down, 0, MPI_COMM_WORLD, &sendReq[5]);	//d
    MPI_Isend(&cur[0][1], 1, ColType, left, 0, MPI_COMM_WORLD, &sendReq[6]);	//l
    MPI_Isend(&cur[0][blockSize], 1, ColType, right, 0, MPI_COMM_WORLD, &sendReq[7]);	//r

    nextStateIn(cur, next, blockSize);    // Calculate next generation for inners

    MPI_Waitall(8, sendReq, status);
  	MPI_Waitall(8, recvReq, status);

    cur[0][0] = recUpLeft;  // Antigrafw ta halo points ston cur
    cur[0][blockSize+1] = recUpRight;
    cur[blockSize+1][0] = recDownLeft;
    cur[blockSize+1][blockSize+1] = recDownRight;
    for(i=1; i<=blockSize; i++){
      cur[0][i] = recUp[i];
      cur[blockSize+1][i] = recDown[i];
      cur[i][0] = recLeft[i];
      cur[i][blockSize+1] = recRight[i];
    }

    nextStateOut(cur, next, blockSize);    // Calculate next generation for outers
    localChanged = 0;
    if( !((stepCount + 1) % 10) ){      	// Elegxw an to grid den allaxe katse 10 epanalipseis
   		if( sameGrid(cur, next, blockSize) )
     		localChanged = 1;

		MPI_Allreduce(&localChanged, &totalChanged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	 	if(totalChanged == size) {
   			if(rank==0)
     			printf("100%%\nGrid didn't change in %d generation! Programm will continue all the steps\n",stepCount+1);
	  	}
 	  }

    temp = cur;     // swap
    cur = next;
    next = temp;
    if(!rank && stepCount % 25 == 0)
      printf("%d %%\n", 100*stepCount/STEPS);
  }

  free(recUp);
  free(recDown);
  free(recRight);
  free(recLeft);
  free(sendRight);
  free(sendLeft);
  for(i=0; i<blockSize+2; i++){
    free(cur[i]);
    free(next[i]);
  }
  free(cur);
  free(next);

  end = MPI_Wtime();
  if(!rank)
    printf("100%%\nWith %d procceses and %d size of problem time was: %.3f seconds\n", size, n, end-start);

  MPI_Finalize();
  return 0;
}
