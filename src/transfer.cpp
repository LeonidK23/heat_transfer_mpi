#include <iostream>
#include <vector>
#include <math.h>

#include "transfer.hpp"

double* heat_transfer_2d(double* grid, int window_size, int ghost_size, const double source_temperature, const double alpha, const double h, const int delta_t){
  int i, j, k, current_index;
  double left, right, bottom, top, current_point;
  double* old_grid;
  int block_size;

  block_size = window_size - 2*ghost_size;
  old_grid = new double[window_size*window_size];
  for (int i = 0; i < window_size*window_size; ++i)
    old_grid[i] = 0;

  // inside the block may be computed ghost_size iterations
  for(k = 0; k < ghost_size; k++){
    std::swap(old_grid, grid);
    for(i = 0; i < window_size; i++)
      // on each iteration reduce size of window because the left- and
      // rightmost ghost columns cannot be updated
      // for(j = k + 1; j <= window_size - (k + 2); j++){
      for(j = 0; j < window_size; j++){
        // if inside of "block to be computed" - compute new value, else copy the old value
        if (j >= k + 1 && j <= window_size - k - 2 && i >= k + 1 && i <= window_size - k - 2){
          // index = row + offset on row + column index
          current_index = i*window_size + j;
          current_point = old_grid[current_index];

          // heat source remains always the same
          if(current_point != source_temperature){

            top = old_grid[current_index - window_size];
            bottom = old_grid[current_index + window_size];
            left = old_grid[current_index - 1];
            right = old_grid[current_index + 1];

            grid[current_index] = current_point + alpha*delta_t*(left+right+top+bottom-4*current_point) / pow(h, 2);
          } else
          grid[current_index] = current_point;
        } else {
          grid[i*window_size + j] = old_grid[i*window_size + j];
        }
      }
    }

    return grid;
}
