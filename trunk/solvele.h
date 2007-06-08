#ifndef _SOLVELE_H_
#define _SOLVELE_H_

#include "sparse_matrix.h"

typedef struct Solvele {
    int size;
    SparseMatrix *mat;
    DenseVector *vec;
} Solvele;

extern double eps_sor;
extern double omega_sor;

void solvele_set_matrix(Solvele *self, int i, int j, double val);
void solvele_add_matrix(Solvele *self, int i, int j, double val);
void solvele_set_vector(Solvele *self, int i, double val);
void solvele_add_vector(Solvele *self, int i, double val);
void solvele_print_matrix(Solvele *self);
void solvele_print_vector(Solvele *self);
double *solvele_solve(Solvele *self, int ni, int nj, int nk);
Solvele *solvele_new(int size);

#endif
