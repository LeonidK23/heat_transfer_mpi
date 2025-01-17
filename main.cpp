#include <iostream>
#include <random>
#include <mpi.h>
#include <chrono>
#include <cmath>
#include <cstdlib>

#include "src/transfer.hpp"
#include "src/auxiliary.hpp"

using std::chrono::system_clock;
using std::chrono::duration;

int main(int argc, char *argv[]) {
  int N = 216;
  int N_ITER = 400;
  int GHOST_ZONE = 10;
  int DELTA_T = 1;
  double ALPHA = 0.4;
  double H = 2.0f;
  const int SAVE_FREQUENCY = 1;
  const int N_SOURCES = 8;
  const double SOURCE_TEMPERATURE = 25;

  // read command line parameters
  if (argc > 1){
    N = atoi(argv[1]);
    GHOST_ZONE = atoi(argv[2]);
    N_ITER = atoi(argv[3]);
    ALPHA = atof(argv[4]);
    H = atof(argv[5]);
    DELTA_T = atoi(argv[6]);
  }

  int num, rank;
  int ret = MPI_Init(&argc, &argv);
  MPI_Request r;
  MPI_Status stat;

  if(MPI_SUCCESS != ret)
    MPI_Abort(MPI_COMM_WORLD, ret);

  MPI_Comm_size(MPI_COMM_WORLD, &num);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // std::mt19937 gen(std::random_device{}());
  std::mt19937 gen(8);
  std::uniform_real_distribution<> dis(0, N);
  int rand_row, rand_col, offset_x, offset_y;
  int window_size, offset, rank_id;
  double *grid, *window_matrix, *block_matrix;

  int n_proc = num, n_proc_x = std::sqrt(num - 1), n_proc_y = std::sqrt(num - 1);
  int block_size = N/n_proc_x;
  window_size = block_size + 2*GHOST_ZONE;
  window_matrix = new double[N*window_size];
  block_matrix = new double[N*block_size];

  double *rbuf;
  if (rank == n_proc-1){
    // allocate array for gathering of gitter
    rbuf = new double[n_proc*block_size*block_size];
  }

  if (rank == num-1){
    // MASTER PROCESS

    // initialize grid
    grid = new double[N*N];
    for (int i = 0; i < N*N; ++i)
      grid[i] = 0;

    for (int i = 0; i < N_SOURCES; i++){
      rand_row = dis(gen);
      rand_col = dis(gen);
      grid[rand_row*N + rand_col] = SOURCE_TEMPERATURE;
    }

    // measure time for the whole computation
    auto start = system_clock::now();

    // split the initialized matrix and send to "working" processors
    for (int i = 0; i < n_proc_y; i++)
      for (int j = 0; j < n_proc_x; j++){
        rank_id = i*n_proc_x + j;
        offset_x = j*block_size - GHOST_ZONE;
        offset_y = i*block_size - GHOST_ZONE;
        window_matrix = slice_matrix(grid, N, rank_id, block_size, block_size, GHOST_ZONE, SOURCE_TEMPERATURE,
                                               j == 0 || j == n_proc_x - 1 || i == 0 || i == n_proc_y - 1, offset_x, offset_y, n_proc_x);
        MPI_Send(window_matrix, window_size*window_size, MPI_DOUBLE, rank_id, 0, MPI_COMM_WORLD);
      }

    for (int k = 0; k < N_ITER/GHOST_ZONE; k++)
      // gather between results from all processors SAVE_FREQUENCY times
      if (k % SAVE_FREQUENCY == 0){
        MPI_Gather(block_matrix, block_size*block_size, MPI_DOUBLE, rbuf, block_size*block_size, MPI_DOUBLE, n_proc-1, MPI_COMM_WORLD);
        grid = reshape_grid_2d(rbuf, N, block_size, n_proc - 1);
        // save grid to file
        print_grid(grid, N, N, true, k);
      }
    auto end = system_clock::now();
    std::cout << "Computation took: " << duration<double>(end - start).count() << " seconds" << '\n';

  } else {
    // WORKING PROCESS

    MPI_Recv(window_matrix, window_size*window_size, MPI_DOUBLE, num-1, 0, MPI_COMM_WORLD, &stat);

    double *ghost_lines, *recv_ghostlines;

    ghost_lines = new double[block_size*GHOST_ZONE];
    recv_ghostlines = new double[block_size*GHOST_ZONE];

    for (int k = 0; k < N_ITER/GHOST_ZONE; k++){
      // compute heat transfer inside the block (# iterations == GHOST_ZONE)
      window_matrix = heat_transfer_2d(window_matrix, window_size, GHOST_ZONE, SOURCE_TEMPERATURE, ALPHA, H, DELTA_T);

      // if not left border process - crop the data for left neighbour and exchange with him
      if (rank % n_proc_x != 0){
        ghost_lines = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, block_size, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, GHOST_ZONE, 0);
        MPI_Sendrecv(ghost_lines, block_size*GHOST_ZONE, MPI_DOUBLE, rank-1, 0, recv_ghostlines,
                     block_size*GHOST_ZONE, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &stat);
        // insert received data into the block matrix
        insert_block(window_matrix, recv_ghostlines, GHOST_ZONE*window_size, window_size, window_size, block_size, GHOST_ZONE);

         // if not the top border process - exchange with the top left neighbour
         if (rank >= n_proc_x){
           ghost_lines = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, GHOST_ZONE, 0);
           MPI_Sendrecv(ghost_lines, GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x-1, 0, recv_ghostlines,
             GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x-1, 0, MPI_COMM_WORLD, &stat);
           insert_block(window_matrix, recv_ghostlines, 0, window_size, window_size, GHOST_ZONE, GHOST_ZONE);
         }
         // if not the bottom border process - exchange with the bottom left process
         if (rank < n_proc_x*(n_proc_y - 1)){
           ghost_lines = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, block_size, 0);
           MPI_Sendrecv(ghost_lines, GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x-1, 0, recv_ghostlines,
             GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x-1, 0, MPI_COMM_WORLD, &stat);
           insert_block(window_matrix, recv_ghostlines, (GHOST_ZONE+block_size)*window_size, window_size, window_size, GHOST_ZONE, GHOST_ZONE);
         }
      }
      // if not right border processor - exchange with tne right neighbour
      if ((rank + 1)%n_proc_x != 0){
        ghost_lines = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, block_size, 0, SOURCE_TEMPERATURE, false, block_size, GHOST_ZONE, 0);
        MPI_Sendrecv(ghost_lines, block_size*GHOST_ZONE, MPI_DOUBLE, rank+1, 0, recv_ghostlines,
                     block_size*GHOST_ZONE, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &stat);
        insert_block(window_matrix, recv_ghostlines, GHOST_ZONE*window_size+GHOST_ZONE+block_size, window_size, window_size, block_size, GHOST_ZONE);
        // if not top border processor - exchange with tne top right neighbour
        if (rank >= n_proc_x){
          ghost_lines = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, block_size, GHOST_ZONE, 0);
          MPI_Sendrecv(ghost_lines, GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x+1, 0, recv_ghostlines,
             GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x+1, 0, MPI_COMM_WORLD, &stat);
          insert_block(window_matrix, recv_ghostlines, GHOST_ZONE+block_size, window_size, window_size, GHOST_ZONE, GHOST_ZONE);
        }
        // if not bottom border processor - exchange with tne bottom right neighbour
        if (rank < n_proc_x*(n_proc_y - 1)){
          ghost_lines = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, block_size, block_size, 0);
          MPI_Sendrecv(ghost_lines, GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x+1, 0, recv_ghostlines,
             GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x+1, 0, MPI_COMM_WORLD, &stat);
          insert_block(window_matrix, recv_ghostlines, (GHOST_ZONE+block_size)*window_size+GHOST_ZONE+block_size, window_size, window_size, GHOST_ZONE, GHOST_ZONE);
         }
        }
      // if not the top border processor - exchange with the top neighbour
      if (rank >= n_proc_x){
        ghost_lines = slice_matrix(window_matrix, window_size, 0, block_size, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, GHOST_ZONE, 0);
        MPI_Sendrecv(ghost_lines, block_size*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x, 0, recv_ghostlines,
                     block_size*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x, 0, MPI_COMM_WORLD, &stat);
        insert_block(window_matrix, recv_ghostlines, GHOST_ZONE, window_size, window_size, GHOST_ZONE, block_size);
      }
      // if not the bottom border processor - exchange with the bottom neighbour
      if (rank < n_proc_x*(n_proc_y - 1)){
        ghost_lines = slice_matrix(window_matrix, window_size, 0, block_size, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, block_size, 0);

        MPI_Sendrecv(ghost_lines, block_size*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x, 0, recv_ghostlines,
                     block_size*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x, 0, MPI_COMM_WORLD, &stat);
        insert_block(window_matrix, recv_ghostlines, (GHOST_ZONE+block_size)*window_size+GHOST_ZONE, window_size, window_size, GHOST_ZONE, block_size);
      }

      // send the between result to master process
      if (k % SAVE_FREQUENCY == 0){
        // take only the block matrix without ghost lines
        block_matrix = slice_matrix(window_matrix, window_size, 0, block_size, block_size, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, GHOST_ZONE, 0);
        // "send" to master
        MPI_Gather(block_matrix, block_size*block_size, MPI_DOUBLE, rbuf, block_size*block_size, MPI_DOUBLE, n_proc-1, MPI_COMM_WORLD);
      }
    }
    // crop out the final result after all iterations are over
    block_matrix = slice_matrix(window_matrix, window_size, 0, block_size, block_size, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, GHOST_ZONE, 0);

    free(ghost_lines);
    free(recv_ghostlines);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // gather final matrices from all processes
  MPI_Gather(block_matrix, block_size*block_size, MPI_DOUBLE, rbuf, block_size*block_size, MPI_DOUBLE, n_proc-1, MPI_COMM_WORLD);

  if (rank == n_proc - 1){
    grid = reshape_grid_2d(rbuf, N, block_size, n_proc - 1);
    print_grid(grid, N, N, true, 999);
    std::cout << "ready" << '\n';

    free(rbuf);
    free(grid);
  }

  MPI_Finalize();
  free(window_matrix);
  free(block_matrix);

  return 0;
}
