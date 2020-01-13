#include <iostream>
#include <random>

#include "src/transfer.hpp"

#define N 8
#define N_SOURCES 5
#define SOURCE_TEMPERATURE 25

void print_grid(double* grid){
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++)
      std::cout << grid[i*N + j] << ' ';
    std::cout << '\n';
  }
}

int main() {
  int i, j;
  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<> dis(0, N);
  int rand_row, rand_col;
  double *grid;

  grid = new double[N*N];
  for (int i = 0; i < N*N; ++i)
        grid[i] = 0;

  for(i = 0; i < N_SOURCES; i++){
    rand_row = dis(gen);
    rand_col = dis(gen);
    grid[rand_row*N + rand_col] = SOURCE_TEMPERATURE;
  }

  std::cout << "Start grid:" << '\n';
  print_grid(grid);
  grid = heat_transfer(grid, 2, N);
  std::cout << "New grid:" << '\n';
  print_grid(grid);

  return 0;
}
