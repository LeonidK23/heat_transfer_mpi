#include <iostream>
#include <random>
#include <mpi.h>
#include <chrono>
#include <cmath>

#include "src/transfer.hpp"
#include "src/auxiliary.hpp"

using std::chrono::system_clock;
using std::chrono::duration;

int main(int argc, char *argv[]) {
  const int N = 9;
  const int N_SOURCES = 8;
  const double SOURCE_TEMPERATURE = 25;
  const double ALPHA = 0.2;
  const double H = 1.0f;
  const int N_ITER = 1;
  const int GHOST_ZONE = 2;

  int num, rank;
  int ret = MPI_Init(&argc, &argv);
  MPI_Request r;
  MPI_Status stat;

  if(MPI_SUCCESS != ret)
    MPI_Abort(MPI_COMM_WORLD, ret);

  MPI_Comm_size(MPI_COMM_WORLD, &num) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // std::mt19937 gen(std::random_device{}());
  std::mt19937 gen(8);
  std::uniform_real_distribution<> dis(0, N);
  int rand_row, rand_col, offset_x, offset_y;
  int window_size, offset, rank_id;
  double *grid, *window_matrix, *block_matrix;

  int n_proc = num, n_proc_x = std::sqrt(num - 1), n_proc_y = std::sqrt(num - 1);
  // int n_proc = num;
  int block_size = N/n_proc_x;
  // window = block + ghost_lines
  window_size = block_size + 2*GHOST_ZONE;
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
    print_grid(grid, N, N, false);
    std::cout << "----------------------" << '\n';
    // split matrix and send to other processors
    for (int i = 0; i < n_proc_y; i++)
      for (int j = 0; j < n_proc_x; j++){
        rank_id = i*n_proc_x + j;
        offset_x = j*block_size - GHOST_ZONE;
        offset_y = i*block_size - GHOST_ZONE;
        window_matrix = slice_matrix(grid, N, rank_id, block_size, block_size, GHOST_ZONE, SOURCE_TEMPERATURE,
                                               j == 0 || j == n_proc_x - 1 || i == 0 || i == n_proc_y - 1, offset_x, offset_y, n_proc_x);
        // first send to all working processors
        MPI_Send(window_matrix, N*window_size, MPI_DOUBLE, rank_id, 0, MPI_COMM_WORLD);
      }
  } else {
    MPI_Recv(window_matrix, N*window_size, MPI_DOUBLE, num-1, 0, MPI_COMM_WORLD, &stat);

    double *ghost_lines_left, *ghost_lines_right, *ghost_lines_top, *ghost_lines_bottom;
    double *ghost_lines_topleft, *ghost_lines_topright, *ghost_lines_bottomleft, *ghost_lines_bottomright;
    double *recv_ghostlines_left, *recv_ghostlines_right, *recv_ghostlines_top, *recv_ghostlines_bottom;
    double *recv_ghostlines_topleft, *recv_ghostlines_topright, *recv_ghostlines_bottomleft, *recv_ghostlines_bottomright;
    // first receive from master processor

    window_matrix = heat_transfer_2d(window_matrix, window_size, GHOST_ZONE, SOURCE_TEMPERATURE, ALPHA, H);

    // slice right block
    ghost_lines_left = new double[block_size*GHOST_ZONE];
    ghost_lines_right = new double[block_size*GHOST_ZONE];
    ghost_lines_top = new double[block_size*GHOST_ZONE];
    ghost_lines_bottom = new double[block_size*GHOST_ZONE];
    ghost_lines_topleft = new double[GHOST_ZONE*GHOST_ZONE];
    ghost_lines_topright = new double[GHOST_ZONE*GHOST_ZONE];
    ghost_lines_bottomleft = new double[GHOST_ZONE*GHOST_ZONE];
    ghost_lines_bottomright = new double[GHOST_ZONE*GHOST_ZONE];
    recv_ghostlines_left = new double[block_size*GHOST_ZONE];
    recv_ghostlines_right = new double[block_size*GHOST_ZONE];
    recv_ghostlines_top = new double[block_size*GHOST_ZONE];
    recv_ghostlines_bottom = new double[block_size*GHOST_ZONE];
    recv_ghostlines_topleft = new double[GHOST_ZONE*GHOST_ZONE];
    recv_ghostlines_topright = new double[GHOST_ZONE*GHOST_ZONE];
    recv_ghostlines_bottomleft = new double[GHOST_ZONE*GHOST_ZONE];
    recv_ghostlines_bottomright = new double[GHOST_ZONE*GHOST_ZONE];

    ghost_lines_right = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, block_size, 0, SOURCE_TEMPERATURE, false, block_size, GHOST_ZONE, 0);
    ghost_lines_left = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, block_size, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, GHOST_ZONE, 0);
    ghost_lines_top = slice_matrix(window_matrix, window_size, 0, block_size, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, GHOST_ZONE, 0);
    ghost_lines_bottom = slice_matrix(window_matrix, window_size, 0, block_size, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, block_size, 0);
    ghost_lines_topleft = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, GHOST_ZONE, 0);
    ghost_lines_topright = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, block_size, GHOST_ZONE, 0);
    ghost_lines_bottomleft = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, GHOST_ZONE, block_size, 0);
    ghost_lines_bottomright = slice_matrix(window_matrix, window_size, 0, GHOST_ZONE, GHOST_ZONE, 0, SOURCE_TEMPERATURE, false, block_size, block_size, 0);

    // if not left border processor - exchange with tne left neighbour
    if (rank % n_proc_x != 0){
      MPI_Sendrecv(ghost_lines_left, block_size*GHOST_ZONE, MPI_DOUBLE, rank-1, 0, recv_ghostlines_left,
                   block_size*GHOST_ZONE, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &stat);
      insert_block(window_matrix, recv_ghostlines_left, GHOST_ZONE*window_size, window_size, window_size, block_size, GHOST_ZONE);
       if (rank >= n_proc_x){
         MPI_Sendrecv(ghost_lines_topleft, GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x-1, 0, recv_ghostlines_topleft,
           GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x-1, 0, MPI_COMM_WORLD, &stat);
         insert_block(window_matrix, recv_ghostlines_topleft, 0, window_size, window_size, GHOST_ZONE, GHOST_ZONE);
       }
       if (rank < n_proc_x*(n_proc_y - 1)){
         MPI_Sendrecv(ghost_lines_bottomleft, GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x-1, 0, recv_ghostlines_bottomleft,
           GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x-1, 0, MPI_COMM_WORLD, &stat);
         insert_block(window_matrix, recv_ghostlines_bottomleft, (GHOST_ZONE+block_size)*window_size, window_size, window_size, GHOST_ZONE, GHOST_ZONE);
       }
    }
    // if not right border processor - exchange with tne right neighbour
    if ((rank + 1)%n_proc_x != 0){
      MPI_Sendrecv(ghost_lines_right, block_size*GHOST_ZONE, MPI_DOUBLE, rank+1, 0, recv_ghostlines_right,
                   block_size*GHOST_ZONE, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &stat);
      insert_block(window_matrix, recv_ghostlines_right, GHOST_ZONE*window_size+GHOST_ZONE+block_size, window_size, window_size, block_size, GHOST_ZONE);
       if (rank >= n_proc_x){
         MPI_Sendrecv(ghost_lines_topright, GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x+1, 0, recv_ghostlines_topright,
           GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x+1, 0, MPI_COMM_WORLD, &stat);
         insert_block(window_matrix, recv_ghostlines_topright, GHOST_ZONE+block_size, window_size, window_size, GHOST_ZONE, GHOST_ZONE);
       }
       if (rank < n_proc_x*(n_proc_y - 1)){
         MPI_Sendrecv(ghost_lines_bottomright, GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x+1, 0, recv_ghostlines_bottomright,
           GHOST_ZONE*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x+1, 0, MPI_COMM_WORLD, &stat);
         insert_block(window_matrix, recv_ghostlines_bottomright, (GHOST_ZONE+block_size)*window_size+GHOST_ZONE+block_size, window_size, window_size, GHOST_ZONE, GHOST_ZONE);
       }
      }
    // if not the top border processor - exchange with the neighbour on the top
    if (rank >= n_proc_x){
      MPI_Sendrecv(ghost_lines_top, block_size*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x, 0, recv_ghostlines_top,
                   block_size*GHOST_ZONE, MPI_DOUBLE, rank-n_proc_x, 0, MPI_COMM_WORLD, &stat);
      insert_block(window_matrix, recv_ghostlines_top, GHOST_ZONE, window_size, window_size, GHOST_ZONE, block_size);
    }
    if (rank < n_proc_x*(n_proc_y - 1)){
      MPI_Sendrecv(ghost_lines_bottom, block_size*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x, 0, recv_ghostlines_bottom,
                   block_size*GHOST_ZONE, MPI_DOUBLE, rank+n_proc_x, 0, MPI_COMM_WORLD, &stat);
      insert_block(window_matrix, recv_ghostlines_bottom, (GHOST_ZONE+block_size)*window_size+GHOST_ZONE, window_size, window_size, GHOST_ZONE, block_size);
    }

    if (rank == 4){
        print_grid(window_matrix, window_size, window_size, false);
        std::cout << "-------------------------" << '\n';
      }

    // free(ghost_lines_left);
    // free(ghost_lines_right);
    // free(ghost_lines_top);
    // free(ghost_lines_bottom);
    // free(ghost_lines_topleft);
    // free(ghost_lines_topright);
    // free(ghost_lines_bottomleft);
    // free(ghost_lines_bottomright);
    // free(recv_ghostlines_left);
    // free(recv_ghostlines_right);
    // free(recv_ghostlines_top);
    // free(recv_ghostlines_bottom);


  //
  //   auto start = system_clock::now();
  //
  //   // for (int k = 0; k < N_ITER/ghost_zone; k++){
  //   for (int k = 0; k < N_ITER; k++){
  //     window_matrix = heat_transfer(window_matrix, N, window_size, ghost_zone, SOURCE_TEMPERATURE, ALPHA, H);
  //     // if (rank == 0)
  //     //   print_grid(window_matrix, N, window_size);
  //     // slice right ghost lines
  //     ghost_lines_right = slice_matrix_rectangle(window_matrix, N, window_size, 0, ghost_zone, 0, SOURCE_TEMPERATURE, false, block_size);
  //     // print_grid(ghost_lines, N, ghost_zone);
  //
  //     if (ghost_zone < block_size){
  //       ghost_lines_left = slice_matrix_rectangle(window_matrix, N, window_size, 0, ghost_zone, 0, SOURCE_TEMPERATURE, false, ghost_zone);
  //       // if (rank == 2)
  //       //   print_grid(ghost_lines_left, N, ghost_zone);
  //     } else
  //       ghost_lines_left = ghost_lines_right;
  //
  //     // if not the border processors
  //     if (rank != 0 && rank != n_proc - 2){
  //       // exchange with the right neighbour
  //       MPI_Sendrecv(ghost_lines_right, N*ghost_zone, MPI_DOUBLE, rank+1, 0, recv_ghostlines_right,
  //                    N*ghost_zone, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &stat);
  //
  //       // exchange with the left neighbour
  //       MPI_Sendrecv(ghost_lines_left, N*ghost_zone, MPI_DOUBLE, rank-1, 0, recv_ghostlines_left,
  //                    N*ghost_zone, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &stat);
  //
  //       // insert reived ghost lines on the left and right side
  //       insert_block(window_matrix, recv_ghostlines_left, 0, N, window_size, N, ghost_zone);
  //       insert_block(window_matrix, recv_ghostlines_right, ghost_zone+block_size, N, window_size, N, ghost_zone);
  //
  //     } else if (rank == 0){
  //       // exchange only with the right neighbour
  //       MPI_Sendrecv(ghost_lines_right, N*ghost_zone, MPI_DOUBLE, rank+1, 0, recv_ghostlines_right,
  //                    N*ghost_zone, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &stat);
  //
  //       // insert received ghost lines on the right side
  //       insert_block(window_matrix, recv_ghostlines_right, ghost_zone+block_size, N, window_size, N, ghost_zone);
  //
  //     } else if (rank == n_proc - 2){
  //       // exchange only with the left neighbour
  //       MPI_Sendrecv(ghost_lines_left, N*ghost_zone, MPI_DOUBLE, rank-1, 0, recv_ghostlines_left,
  //                    N*ghost_zone, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &stat);
  //
  //       // insert received ghost lines on the left side
  //       insert_block(window_matrix, recv_ghostlines_left, 0, N, window_size, N, ghost_zone);
  //     }
  //
  //   }
  //
  //   auto end = system_clock::now();
  //   // std::cout << "Computing time:" << duration<double>(end - start).count() << '\n';
  //
  //   block_matrix = slice_matrix_rectangle(window_matrix, N, window_size, 0, block_size, 0, SOURCE_TEMPERATURE, false, ghost_zone);
  //
  //   free(ghost_lines_left);
  //   free(ghost_lines_right);
  //   free(recv_ghostlines_left);
  //   free(recv_ghostlines_right);
  }
  //
  // MPI_Barrier(MPI_COMM_WORLD);
  //
  // // gather matrices from all processes
  // MPI_Gather(block_matrix, N*block_size, MPI_DOUBLE, rbuf, block_size*N, MPI_DOUBLE, n_proc-1, MPI_COMM_WORLD);
  // if (rank == n_proc - 1){
  //   grid = reshape_grid(rbuf, N, block_size);
  //
  //   // std::cout << "Final grid:" << '\n';
  //   // print_grid(grid, N, N, false);
  //   // std::cout << "ready" << '\n';
  //
  //   free(rbuf);
  //   free(grid);
  // }

  MPI_Finalize();
  free(window_matrix);
  free(block_matrix);

  return 0;
}
