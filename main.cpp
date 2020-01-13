#include <iostream>
#include <random>

#include "src/transfer.hpp"
#include "src/auxiliary.hpp"

void print_grid(double* grid, int N){
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++)
      std::cout << grid[i*N + j] << ' ';
    std::cout << '\n';
  }
}

int main() {
  const int N = 8;
  const int N_SOURCES = 5;
  const int SOURCE_TEMPERATURE = 25;

  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<> dis(0, N);
  int rand_row, rand_col;
  int ghost_zone, window_size;
  double *grid;

  // initialize grid
  grid = new double[N*N];
  for (int i = 0; i < N*N; ++i)
        grid[i] = 0;

  for (int i = 0; i < N_SOURCES; i++){
    rand_row = dis(gen);
    rand_col = dis(gen);
    grid[rand_row*N + rand_col] = SOURCE_TEMPERATURE;
  }

  // split matrix and send to other processors
  int n_proc = 16, n_proc_x = 4, n_proc_y = 4;
  int block_size = 2;
  int offset_without_ghost, offset_with_ghost;
  ghost_zone = 1;
  window_size = block_size + 2*ghost_zone;
  for (int i = 0; i < n_proc_y; i++)
    for (int j = 0; j < n_proc_x; j++){
      // original offset - ghost zone
      std::cout << i*n_proc_x + j << '\n';
      offset_without_ghost = i*N*block_size + j*block_size;
      offset_with_ghost = (i*N*block_size + j*block_size) - ghost_zone*N - ghost_zone;
      grid = slice_matrix(grid, N, offset_with_ghost, window_size);
      // std::cout << i << ' ' << j << ' ' << offset_without_ghost << '\n';
      // std::cout << i << ' ' << j << ' ' << offset_with_ghost << '\n';
      // std::cout << "----------------------------------" << '\n';
    }

  // std::cout << "Start grid:" << '\n';
  // print_grid(grid, N);
  grid = heat_transfer(grid, 2, N, SOURCE_TEMPERATURE);
  // std::cout << "New grid:" << '\n';
  // print_grid(grid, N);

  return 0;
}
