#include <stdlib.h>

#include "bmmm.h"

double** initialize_matrix(int n){
  int i, j;
  double **arr;

  arr = (double **)malloc(n * sizeof(double *));
  for (i = 0; i < n; i++)
    arr[i] = (double *)malloc(n * sizeof(double));

   for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
     arr[i][j] = rand()/(double)RAND_MAX;

  return arr;
}

double* slice_matrix(double mat[4][4], int ind_x, int ind_y, int size){
  double *block;
  int i, j;

  block = (double *)malloc(size*size*sizeof(double));
  for (i = 0; i < size; i++){
    for (j = 0; j < size; j++)
      block[i*size + j] = mat[i + ind_x][j + ind_y];
  }

  // for (i = ind_x; i < size + ind_x; i++){
  //   for (j = ind_y; j < size + ind_y; j++){
  //     printf("%f\n", block[i*size+j]);
  //   }
  // }

  return block;
}
