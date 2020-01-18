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
  const int N_ITER = 10;

  int num, rank;
  int ret = MPI_Init(&argc, &argv);
  MPI_Request r;
  MPI_Status stat;

  if(MPI_SUCCESS != ret)
    MPI_Abort(MPI_COMM_WORLD, ret);

  MPI_Comm_size(MPI_COMM_WORLD, &num) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<> dis(0, N);
  int rand_row, rand_col;
  int ghost_zone, window_size, offset;
  double *grid, *window_matrix, *block_matrix;

  // int n_proc = 16, n_proc_x = 4, n_proc_y = 4;
  int n_proc = num;
  int block_size = N/(n_proc - 1);
  ghost_zone = 3;
  // window = block + ghost_lines
  window_size = block_size + 2*ghost_zone;
  window_matrix = new double[N*window_size];
  block_matrix = new double[N*block_size];

  double *rbuf;
  if (rank == n_proc-1){
    rbuf = new double[n_proc*N*block_size];
  }

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

    std::cout << "Start grid:" << '\n';
    print_grid(grid, N, N);
    std::cout << "----------------------" << '\n';
    // split matrix and send to other processors
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
  } else {
    double *ghost_lines_left, *ghost_lines_right, *recv_ghostlines_left, *recv_ghostlines_right;

    ghost_lines_right = new double[N*ghost_zone];
    ghost_lines_left = new double[N*ghost_zone];
    recv_ghostlines_left = new double[N*ghost_zone];
    recv_ghostlines_right = new double[N*ghost_zone];

    // first receive from master processor
    MPI_Recv(window_matrix, N*window_size, MPI_DOUBLE, num-1, 0, MPI_COMM_WORLD, &stat);

    for (int k = 0; k < N_ITER; k++){
      window_matrix = heat_transfer(window_matrix, N, window_size, ghost_zone, SOURCE_TEMPERATURE, ALPHA, H);
      // if (rank == 0)
      //   print_grid(window_matrix, N, window_size);
      // slice right ghost lines
      ghost_lines_right = slice_matrix_rectangle(window_matrix, N, window_size, 0, ghost_zone, 0, SOURCE_TEMPERATURE, false, block_size);
      // print_grid(ghost_lines, N, ghost_zone);

      if (ghost_zone < block_size){
        ghost_lines_left = slice_matrix_rectangle(window_matrix, N, window_size, 0, ghost_zone, 0, SOURCE_TEMPERATURE, false, ghost_zone);
        // if (rank == 2)
        //   print_grid(ghost_lines_left, N, ghost_zone);
      } else
        ghost_lines_left = ghost_lines_right;

      // if not the border processors
      if (rank != 0 && rank != n_proc - 2){
        // exchange with the right neighbour
        MPI_Sendrecv(ghost_lines_right, N*ghost_zone, MPI_DOUBLE, rank+1, 0, recv_ghostlines_right,
                     N*ghost_zone, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &stat);

        // exchange with the left neighbour
        MPI_Sendrecv(ghost_lines_left, N*ghost_zone, MPI_DOUBLE, rank-1, 0, recv_ghostlines_left,
                     N*ghost_zone, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &stat);

        // insert reived ghost lines on the left and right side
        insert_block(window_matrix, recv_ghostlines_left, 0, N, window_size, N, ghost_zone);
        insert_block(window_matrix, recv_ghostlines_right, ghost_zone+block_size, N, window_size, N, ghost_zone);

      } else if (rank == 0){
        // exchange only with the right neighbour
        MPI_Sendrecv(ghost_lines_right, N*ghost_zone, MPI_DOUBLE, rank+1, 0, recv_ghostlines_right,
                     N*ghost_zone, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &stat);

        // insert received ghost lines on the right side
        insert_block(window_matrix, recv_ghostlines_right, ghost_zone+block_size, N, window_size, N, ghost_zone);

      } else if (rank == n_proc - 2){
        // exchange only with the left neighbour
        MPI_Sendrecv(ghost_lines_left, N*ghost_zone, MPI_DOUBLE, rank-1, 0, recv_ghostlines_left,
                     N*ghost_zone, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &stat);

        // insert received ghost lines on the left side
        insert_block(window_matrix, recv_ghostlines_left, 0, N, window_size, N, ghost_zone);
      }
    }

    block_matrix = slice_matrix_rectangle(window_matrix, N, window_size, 0, block_size, 0, SOURCE_TEMPERATURE, false, ghost_zone);

    free(ghost_lines_left);
    free(ghost_lines_right);
    free(recv_ghostlines_left);
    free(recv_ghostlines_right);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // gather matrices from all processes
  MPI_Gather(block_matrix, N*block_size, MPI_DOUBLE, rbuf, block_size*N, MPI_DOUBLE, n_proc-1, MPI_COMM_WORLD);
  if (rank == n_proc - 1){
    grid = reshape_grid(rbuf, N, block_size);

    std::cout << "Final grid:" << '\n';
    print_grid(grid, N, N);

    free(rbuf);
    free(grid);
  }

  MPI_Finalize();
  free(window_matrix);
  free(block_matrix);

  return 0;
}
