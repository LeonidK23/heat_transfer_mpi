#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "src/bmmm.h"

int main(int argc, char *argv[]) {
  int n;
  int i, j;
  int num, rank, block_size;
  // double **A, **B;
  double *block;

  int ret = MPI_Init(&argc, &argv);
  MPI_Status stat;
  if(MPI_SUCCESS != ret)
    MPI_Abort(MPI_COMM_WORLD, ret);

  MPI_Comm_size(MPI_COMM_WORLD, &num) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  n = 4;
  block_size = 2;
  // master rank is the total last rank
  if(rank == num-1){
    int rank_id;
    // A = initialize_matrix(n);
    // B = initialize_matrix(n);
    double A[4][4] = {{1, 2, 3, 4},
                      {5, 6, 7 ,8},
                      {9, 10, 11, 12},
                      {13, 14, 15, 16}};
    double B[4][4] = {{1, 2, 3, 4},
                      {5, 6, 7 ,8},
                      {9, 10, 11, 12},
                      {13, 14, 15, 16}};

    block = (double *)malloc(block_size*block_size*sizeof(double));
    for (i = 0; i < n; i+=block_size)
      for (j = 0; j < n; j+=block_size){
        block = slice_matrix(A, i, j, block_size);
        rank_id = (int)((i*n / block_size + j)/block_size);
        // send blocks to corresponding ranks
        MPI_Send(block, block_size*block_size, MPI_DOUBLE, rank_id, 0, MPI_COMM_WORLD);
      }

    free(block);
    // free(A);
    // free(B);
  }
  else{
    printf("I am rank %i\n", rank);
    block = (double *)malloc(block_size*block_size*sizeof(double));
    MPI_Recv(block, block_size*block_size, MPI_DOUBLE, num-1, 0, MPI_COMM_WORLD, &stat);
    for (i = 0; i < block_size; i++)
      for (j = 0; j < block_size; j++)
        printf("%f\n", block[i*block_size + j]);
    printf("---\n");
  }
  // computing of block matrix on each "slave" rank
  // else {
  //    MPI_Recv(get_arr, n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &stat);
  //    MPI_Send(arr, n, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
  // }




  MPI_Finalize();

  return 0;
}
