#ifndef AUXILIARY_H
#define AUXILIARY_H

double* slice_matrix_square(double*, int, int, int, int);
double* slice_matrix_rectangle(double*, int, int, int, int, int, double, bool, int);
void insert_ghost_lines(double*, double*, int, int, int, int, int);
void print_grid(double*, int, int);

#endif
