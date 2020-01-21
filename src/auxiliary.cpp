#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "auxiliary.hpp"

double* slice_matrix(double* grid, int n, int rank_id, int block_size_x, int block_size_y, int ghost_size, double source_temp, bool is_border, int offset_x, int offset_y, int n_rank_dim){
  int grid_ind, window_size_x, window_size_y;
  window_size_x = block_size_x + 2*ghost_size;
  window_size_y = block_size_y + 2*ghost_size;

  double* window = new double[window_size_x*window_size_y];

  if (!is_border){
    for (int i = 0; i < window_size_y; i++){
      for (int j = 0; j < window_size_x; j++){
        grid_ind = offset_y*n + offset_x + i*n + j;
        window[i*window_size_x + j] = grid[grid_ind];
      }
      // std::cout << '\n';
    }
    // left border of processors grid
  } else if (rank_id % n_rank_dim == 0){
      for (int i = 0; i < window_size_y; i++){
        for (int j = 0; j < window_size_x; j++){
          grid_ind = offset_y*n + offset_x + i*n + j;
          // if left ghost line, or if left top ghost line or if left bottom ghost line
          if (j < ghost_size || (rank_id == 0 &&  i < ghost_size) || (rank_id == n_rank_dim * (n_rank_dim-1) && i >= window_size_y - ghost_size)){
            window[i*window_size_x + j] = source_temp;
          }
          else {
            window[i*window_size_x + j] = grid[grid_ind];
          }
        }
      }
    } else if ((rank_id+1) % n_rank_dim == 0){
      for (int i = 0; i < window_size_y; i++){
        for (int j = 0; j < window_size_x; j++){
          grid_ind = offset_y*n + offset_x + i*n + j;
          // if right ghost line, or if right top ghost line or if right bottom ghost line
          if (j >= window_size_x - ghost_size || (rank_id == n_rank_dim - 1 && i < ghost_size) || (rank_id == n_rank_dim*n_rank_dim - 1 && i >= window_size_y - ghost_size)){
            window[i*window_size_x + j] = source_temp;
          } else {
            window[i*window_size_x + j] = grid[grid_ind];
          }
        }
      }
    } else {
      // only top or only bottom processor
      for (int i = 0; i < window_size_y; i++){
        for (int j = 0; j < window_size_x; j++){
          grid_ind = offset_y*n + offset_x + i*n + j;
          // if top ghost line, or if bottom ghost line
          if ((rank_id < n_rank_dim && i < ghost_size) || (rank_id >= n_rank_dim * (n_rank_dim-1) && i >= window_size_y - ghost_size)){
            window[i*window_size_x + j] = source_temp;
          } else {
            window[i*window_size_x + j] = grid[grid_ind];
          }
        }
      }
    }

  return window;
}


double* slice_matrix_rectangle(double* grid, int m, int n, int rank_id, int block_size, int ghost_size, double source_temp, bool is_border, int offset){
  int grid_ind, window_size;
  window_size = block_size + 2*ghost_size;

  double* window = new double[window_size*m];

  if (!is_border){
    for (int i = 0; i < m; i++){
      for (int j = 0; j < window_size; j++){
        grid_ind = i*n + j + offset;
        window[i*window_size + j] = grid[grid_ind];
      }
    }
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
        }
      }
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
        }
      }
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


double* reshape_grid_2d(double* mat, int N, int block_size, const int n_blocks){
  double *new_mat = new double[N*N];
  int n_blocks_row = std::sqrt(n_blocks);

  for (int b = 0; b < n_blocks_row; b++)
    for (int b_x = 0; b_x < n_blocks_row; b_x++)
      for (int i = 0; i < block_size; i++)
        for (int j = 0; j < block_size; j++)
          new_mat[b*block_size*N + b_x*block_size + i*N + j] = mat[b*n_blocks_row*block_size*block_size + b_x*block_size*block_size + i*block_size + j];

  return new_mat;
}

void print_grid(double* grid, int m, int n, bool to_file, int iter){
  std::ofstream myfile;
  std::string filename ("data/results_");
  filename += std::to_string(iter)+".txt";

  if (to_file == true)
    myfile.open (filename);

  for (int i = 0; i < m; i++){
    for (int j = 0; j < n; j++){
      if (to_file == true)
        myfile << i << " " << j << " " << grid[i*n + j] << '\n';
      else
        std::cout << grid[i*n + j] << ' ';
    }
    if (to_file == false)
      std::cout << '\n';
  }

  myfile.close();
}
