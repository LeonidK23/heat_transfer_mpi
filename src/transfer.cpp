#include <iostream>
#include <vector>
#include <math.h>

#include "transfer.hpp"

double* heat_transfer(double* grid, int N, int window_size, int ghost_size, const double source_temperature, const double alpha, const double h){
  int i, j, k, current_index;
  double left, right, bottom, top, current_point;
  double* old_grid;
  int block_size;

  block_size = window_size - 2*ghost_size;
  old_grid = new double[N*window_size];
  for (int i = 0; i < N*window_size; ++i)
    old_grid[i] = 0;

  for(k = 0; k < ghost_size; k++){
    std::swap(old_grid, grid);
    for(i = 0; i < N; i++)
      // on each iteration reduce size of window because the left- and
      // rightmost ghost columns cannot be updated
      for(j = 1; j <= window_size - 2*(k + 1); j++){
        // index = row + offset on row + column index
        current_index = i*window_size + k + j;
        current_point = old_grid[current_index];
        // std::cout << current_index << '\n';
        // std::cout << k << ' ' << i << ' ' << j << '\n';

        if(current_point != source_temperature){
          // perform check on left and left in "for i" loop

          if(i == 0){
            top = source_temperature;
            bottom = old_grid[current_index + window_size];
          } else if(i == N - 1){
            bottom = source_temperature;
            top = old_grid[current_index - window_size];
          } else {
            top = old_grid[current_index - window_size];
            bottom = old_grid[current_index + window_size];
          }
          left = old_grid[current_index - 1];
          right = old_grid[current_index + 1];

          grid[current_index] = current_point + alpha * (left+right+top+bottom-4*current_point) / pow(h, 2);
        } else
          grid[current_index] = current_point;
      }
    }

    return grid;
}
