#include <iostream>
#include <vector>
#include <math.h>

#include "serial.h"

#define ALPHA 0.2
#define SOURCE 100
#define H 2

std::vector<std::vector<float>> compute_heat_propagation(
      std::vector<std::vector<float>> grid, int n_iter){
  int i, j, k;
  float left, right, bottom, top, current_point;
  int n_rows, n_cols;
  std::vector<std::vector<float>> old_grid;

  n_rows = grid.size();
  n_cols = grid[0].size();
  for(k = 0; k < n_iter; k++){
    old_grid = grid;
    for(i = 0; i < n_rows; i++)
      for(j = 0; j < n_cols; j++){
        current_point = old_grid[i][j];
        // perform check on left and left in "for i" loop
        if(i == 0){
          top = 0;
          bottom = old_grid[i+1][j];
        }
        else if(i == n_rows - 1){
          bottom = 0;
          top = old_grid[i-1][j];
        }
        else{
          top = old_grid[i-1][j];
          bottom = old_grid[i+1][j];
        }
        if(j == 0){
          left = 0;
          right = old_grid[i][j+1];
        }
        else if(j == n_cols - 1){
          right = 0;
          left = old_grid[i][j-1];
        }
        else{
          left = old_grid[i][j-1];
          right = old_grid[i][j+1];
        }

        if(current_point != SOURCE)
          grid[i][j] = current_point + ALPHA * (left+right+top+bottom-4*current_point) / pow(H, 2);
      }
    }

    return grid;
}
