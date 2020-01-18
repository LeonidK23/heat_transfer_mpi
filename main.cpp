#include <iostream>
#include <random>
#include <mpi.h>

#include "src/transfer.hpp"
#include "src/auxiliary.hpp"

int main(int argc, char *argv[]) {
  const int N = 16;
  const int N_SOURCES = 5;
  const double SOURCE_TEMPERATURE = 25;
  const double ALPHA = 0.2;
  const double H = 1.0f;

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
  int ghost_zone, window_size, offset;
  double *grid, *window_matrix;

  // int n_proc = 16, n_proc_x = 4, n_proc_y = 4;
  int n_proc = num;
  int block_size = N/(n_proc - 1);
  ghost_zone = 4;
  // window = block + ghost_lines
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
    print_grid(grid, N, N);
    // for (int i = 0; i < n_proc_y/2; i++)
    //   for (int j = 0; j < n_proc_x/2; j++){
      for (int rank_id = 0; rank_id < n_proc-1; rank_id++){
        // original offset - ghost zone
        // std::cout << "First send to processor " << rank_id << '\n';
        // block_matrix = slice_matrix_rectangle(grid, N, offset_with_ghost, window_size, SOURCE_TEMPERATURE);
        offset = rank_id*block_size - ghost_zone;
        window_matrix = slice_matrix_rectangle(grid, N, N, rank_id, block_size, ghost_zone, SOURCE_TEMPERATURE,
                                               rank_id == 0 || rank_id == n_proc-2, offset);

        // first send to all working processors
        MPI_Send(window_matrix, N*window_size, MPI_DOUBLE, rank_id, 0, MPI_COMM_WORLD);
      }


      // std::cout << "New grid:" << '\n';
      // print_grid(grid, N);
  } else {
    double *ghost_lines;

    ghost_lines = new double[N*ghost_zone];

    // first receive from master processor
    MPI_Recv(window_matrix, N*window_size, MPI_DOUBLE, num-1, 0, MPI_COMM_WORLD, &stat);

    if (rank == 1){
      std::cout << "---------------------------------" << '\n';
      // print_grid(window_matrix, N, window_size);
      // std::cout << "---------------------------------" << '\n';
      window_matrix = heat_transfer(window_matrix, N, window_size, ghost_zone, SOURCE_TEMPERATURE, ALPHA, H);
      print_grid(window_matrix, N, window_size);
      // slice right ghost lines
      ghost_lines = slice_matrix_rectangle(window_matrix, N, window_size, 0, ghost_zone, 0, SOURCE_TEMPERATURE, false, block_size);
      print_grid(ghost_lines, N, ghost_zone);
    }

  }

  MPI_Finalize();
  // free(window_matrix);
  // free(grid);

  return 0;
}
