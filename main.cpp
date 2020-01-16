#include <iostream>
#include <random>
#include <mpi.h>

#include "src/transfer.hpp"
#include "src/auxiliary.hpp"

void print_grid(double* grid, int N){
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++)
      std::cout << grid[i*N + j] << ' ';
    std::cout << '\n';
  }
}

int main(int argc, char *argv[]) {
  const int N = 8;
  const int N_SOURCES = 5;
  const int SOURCE_TEMPERATURE = 25;
  const int ALPHA = 0.2;
  const int H = 1;

  int num, rank;
  int ret = MPI_Init(&argc, &argv);
  MPI_Request r, r1, r2;
  MPI_Status stat;

  if(MPI_SUCCESS != ret)
    MPI_Abort(MPI_COMM_WORLD, ret);

  MPI_Comm_size(MPI_COMM_WORLD, &num) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<> dis(0, N);
  int rand_row, rand_col;
  int ghost_zone, window_size;
  double *grid, *window_matrix;

  // int n_proc = 16, n_proc_x = 4, n_proc_y = 4;
  int n_proc = num;
  int block_size = 2;
  ghost_zone = 1;
  window_size = block_size + 2*ghost_zone;
  window_matrix = new double[N*window_size];

  if (rank == num-1){
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
    std::cout << "Start grid:" << '\n';
    print_grid(grid, N);
    // for (int i = 0; i < n_proc_y/2; i++)
    //   for (int j = 0; j < n_proc_x/2; j++){
      for (int rank_id = 0; rank_id < n_proc-1; rank_id++){
        // original offset - ghost zone
        std::cout << "First send to processor " << rank_id << '\n';
        // block_matrix = slice_matrix_rectangle(grid, N, offset_with_ghost, window_size, SOURCE_TEMPERATURE);
        window_matrix = slice_matrix_rectangle(grid, N, rank_id, block_size, ghost_zone, SOURCE_TEMPERATURE,
                                               rank_id == 0 || rank_id == n_proc-1);

        // first send to all working processors
        MPI_Send(window_matrix, N*window_size, MPI_DOUBLE, rank_id, 0, MPI_COMM_WORLD);
      }


      // std::cout << "New grid:" << '\n';
      // print_grid(grid, N);
  } else {
    // first receive from master processor
    MPI_Recv(window_matrix, N*window_size, MPI_DOUBLE, num-1, 0, MPI_COMM_WORLD, &stat);

    // grid = heat_transfer(grid, 2, N, SOURCE_TEMPERATURE, ALPHA, H);
  }

  MPI_Finalize();
  free(window_matrix);
  free(grid);

  return 0;
}
