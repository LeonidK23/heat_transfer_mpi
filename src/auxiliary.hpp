#ifndef AUXILIARY_H
#define AUXILIARY_H

double* slice_matrix_square(double*, int, int, int, int);
double* slice_matrix_rectangle(double*, int, int, int, int, int, double, bool, int);
double* reshape_grid(double*, int, int);
void insert_block(double*, double*, int, int, int, int, int);
void print_grid(double*, int, int);


#endif
