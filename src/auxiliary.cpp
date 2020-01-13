#include "auxiliary.hpp"

#include <iostream>

double* slice_matrix(double* grid, int N, int offset, int window_size){
  double* window = new double[window_size*window_size];

  for (int i = 0; i < window_size; i++){
    for (int j = 0; j < window_size; j++)
      // window[i*window_size + j] = grid[(i+1)*offset + i*window_size + j];
      std::cout << offset + i*N + j << ' ';
    std::cout << '\n';
  }
  std::cout << "-------------------------------------" << '\n';

  return window;
}
