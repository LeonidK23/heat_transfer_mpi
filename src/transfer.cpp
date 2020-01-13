#include <iostream>
#include <vector>
#include <math.h>

#include "transfer.hpp"

double* heat_transfer(double* grid, int n_iter, int block_size, const int source_temperature, const int alpha, const int h){
  int i, j, k;
  float left, right, bottom, top, current_point;
  double* old_grid;

  old_grid = new double[block_size*block_size];
  for (int i = 0; i < block_size*block_size; ++i)
    old_grid[i] = 0;

  for(k = 0; k < n_iter; k++){
    std::swap(old_grid, grid);
    for(i = 0; i < block_size; i++)
      for(j = 0; j < block_size; j++){
        current_point = old_grid[i*block_size + j];
        // perform check on left and left in "for i" loop
        if(i == 0){
          top = 0;
          bottom = old_grid[(i+1)*block_size + j];
        }
        else if(i == block_size - 1){
          bottom = 0;
          top = old_grid[(i-1)*block_size + j];
        }
        else{
          top = old_grid[(i-1)*block_size + j];
          bottom = old_grid[(i+1)*block_size + j];
        }
        if(j == 0){
          left = 0;
          right = old_grid[i*block_size + j + 1];
        }
        else if(j == block_size - 1){
          right = 0;
          left = old_grid[i*block_size + j - 1];
        }
        else{
          left = old_grid[i*block_size + j - 1];
          right = old_grid[i*block_size + j + 1];
        }

        if(current_point != source_temperature)
          grid[i*block_size + j] = current_point + alpha * (left+right+top+bottom-4*current_point) / pow(h, 2);
        else
          grid[i*block_size + j] = current_point;
      }
    }

    return grid;
}
