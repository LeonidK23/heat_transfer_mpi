#include <iostream>

#include "auxiliary.hpp"

double* slice_matrix_square(double* grid, int N, int offset, int window_size, int source_temp){
  double* window = new double[window_size*window_size];
  int grid_ind;

  for (int i = 0; i < window_size; i++){
    for (int j = 0; j < window_size; j++){
      grid_ind = offset + i*N + j;
      if (grid_ind < 0)
        window[i*window_size + j] = source_temp;
      else
        window[i*window_size + j] = grid[grid_ind];
      std::cout << grid_ind << ' ';
    }
    // std::cout << '\n';
  }
  // std::cout << "-------------------------------------" << '\n';

  return window;
}


double* slice_matrix_rectangle(double* grid, int N, int rank_id, int block_size,int ghost_size, int source_temp, bool is_border){
  int grid_ind, offset, window_size;
  offset = rank_id*block_size - ghost_size;
  window_size = block_size + 2*ghost_size;

  double* window = new double[window_size*N];

  if (!is_border){
    for (int i = 0; i < N; i++){
      for (int j = 0; j < window_size; j++){
        grid_ind = i*N + j + offset;
        window[i*window_size + j] = grid[grid_ind];
        // std::cout << grid_ind << ' ' << grid[grid_ind] << ' ';
        // std::cout << window[i*window_size + j] << ' ';
      }
      // std::cout << '\n';
    }
    // std::cout << "-------------------------------------" << '\n';
  } else {
    if (rank_id == 0){
      for (int i = 0; i < N; i++){
        for (int j = 0; j < window_size; j++){
          // the leftmost ghost line = source temperature
          if (j < ghost_size){
            grid_ind = i*N + j + offset;
            window[i*window_size + j] = source_temp;
          }
          else {
            grid_ind = i*N + j + offset;
            window[i*window_size + j] = grid[grid_ind];
          }
          // std::cout << grid_ind << ' ' << window[i*window_size + j] << ' ';
          // std::cout << window[i*window_size + j] << ' ';
        }
        // std::cout << '\n';
      }
      // std::cout << "-------------------------------------" << '\n';
    } else {
      std::cout << rank_id << '\n';
      for (int i = 0; i < N; i++){
        for (int j = 0; j < window_size; j++){
          // the rightmost ghost line also = source temperature
          if (j >= window_size - ghost_size){
            grid_ind = i*N + j + offset;
            window[i*window_size + j] = source_temp;
          }
          else {
            grid_ind = i*N + j + offset;
            window[i*window_size + j] = grid[grid_ind];
          }
          // std::cout << grid_ind << ' ' << window[i*window_size + j] << ' ';
          // std::cout << window[i*window_size + j] << ' ';
        }
        // std::cout << '\n';
      }
      // std::cout << "-------------------------------------" << '\n';
    }
  }

  return window;
}
