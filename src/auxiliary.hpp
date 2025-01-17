#ifndef AUXILIARY_H
#define AUXILIARY_H

double* slice_matrix(double*, int, int, int, int, int, double, bool, int, int, int);
double* reshape_grid_2d(double*, int, int, const int);
void insert_block(double*, double*, int, int, int, int, int);
void print_grid(double*, int, int, bool, int);


#endif
