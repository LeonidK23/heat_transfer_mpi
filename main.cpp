#include <iostream>
#include <vector>
#include <random>

#include "src/serial.h"

#define N 8
#define N_SOURCES 5

void print_grid(std::vector<std::vector<float>> grid){
  int i, j;
  int n_rows = grid.size();
  int n_cols = grid[0].size();

  for (i = 0; i < n_rows; i++){
    for(j = 0; j < n_cols; j++)
      std::cout << grid[i][j] << ' ';
    std::cout << '\n';
    }
}

int main() {
  int i, j;
  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<> dis(0, N);
  int rand_row, rand_col;
  std::vector<std::vector<float>> grid(N);

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++){
      grid[i].push_back(0);
    }
  }

  for(i = 0; i < N_SOURCES; i++){
    rand_row = dis(gen);
    rand_col = dis(gen);
    grid[rand_row][rand_col] = SOURCE;
  }

  std::cout << "Start grid:" << '\n';
  print_grid(grid);
  grid = compute_heat_propagation(grid, 1);
  std::cout << "New grid:" << '\n';
  print_grid(grid);
  
  return 0;
}
