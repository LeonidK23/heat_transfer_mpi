#include <iostream>
#include <fstream>

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


double* slice_matrix_rectangle(double* grid, int m, int n, int rank_id, int block_size, int ghost_size, double source_temp, bool is_border, int offset){
  int grid_ind, window_size;
  // offset = rank_id*block_size - ghost_size;
  window_size = block_size + 2*ghost_size;

  double* window = new double[window_size*m];

  if (!is_border){
    for (int i = 0; i < m; i++){
      for (int j = 0; j < window_size; j++){
        grid_ind = i*n + j + offset;
        window[i*window_size + j] = grid[grid_ind];
        // std::cout << grid_ind << ' ' << grid[grid_ind] << ' ';
        // std::cout << window[i*window_size + j] << ' ';
      }
      // std::cout << '\n';
    }
    // std::cout << "-------------------------------------" << '\n';
  } else {
    if (rank_id == 0){
      for (int i = 0; i < m; i++){
        for (int j = 0; j < window_size; j++){
          // the leftmost ghost line = source temperature
          if (j < ghost_size){
            grid_ind = i*n + j + offset;
            window[i*window_size + j] = source_temp;
          }
          else {
            grid_ind = i*n + j + offset;
            window[i*window_size + j] = grid[grid_ind];
          }
          // std::cout << grid_ind << ' ' << window[i*window_size + j] << ' ';
          // std::cout << window[i*window_size + j] << ' ';
        }
        // std::cout << '\n';
      }
      // std::cout << "-------------------------------------" << '\n';
    } else {
      for (int i = 0; i < m; i++){
        for (int j = 0; j < window_size; j++){
          // the rightmost ghost line also = source temperature
          if (j >= window_size - ghost_size){
            grid_ind = i*n + j + offset;
            window[i*window_size + j] = source_temp;
          }
          else {
            grid_ind = i*n + j + offset;
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

void insert_block(double* window_matrix, double* ghost_lines, int offset, int m, int n, int gl_m, int gl_n){
  int global_ind;

  for (int i = 0; i < gl_m; i++){
    for (int j = 0; j < gl_n; j++){
      window_matrix[offset + i*n + j] = ghost_lines[i*gl_n + j];
    }
  }
}

double* reshape_grid(double* mat, int N, int block_size){
  double *new_mat = new double[N*N];

  for (int b = 0; b < N/block_size; b++)
    for (int i = 0; i < N; i++){
      for (int j = 0; j < block_size; j++){
        new_mat[b*block_size + i*N + j] = mat[b*N*block_size + i*block_size + j];
      }
    }

  return new_mat;
}

void print_grid(double* grid, int m, int n){
  std::ofstream myfile;
  myfile.open ("results.txt");

  for (int i = 0; i < m; i++){
    for (int j = 0; j < n; j++){
      // std::cout << grid[i*n + j] << ' ';
      myfile << i << " " << j << " " << grid[i*n + j] << '\n';
    }
    // std::cout << '\n';
  }

  myfile.close();
}
