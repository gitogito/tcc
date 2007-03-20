#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DenseVector {
    int n;
    double *val;
} DenseVector;

typedef struct SpMatElem {
    double val;
    int idx;
    struct SpMatElem *next;
} SpMatElem;

typedef struct SparseMatrix {
    int n;
    int nnz;
    SpMatElem **elems;
} SparseMatrix;

DenseVector *dvec_new(int n);
void dvec_free(DenseVector *v);
void dvec_set(DenseVector *v, int i, double val);
double dvec_ref(DenseVector *v, int i);
void dvec_print(DenseVector *v);
SparseMatrix *smat_new(int n);
void smat_free(SparseMatrix *m);
void smat_set(SparseMatrix *m, int i, int j, double val);
double smat_ref(SparseMatrix *m, int i, int j);
void smat_print(SparseMatrix *m);
void get_crs(SparseMatrix *a0, DenseVector *b0,
	int **pap, int **pai, double **pax, double **pb);

#ifdef __cplusplus
}
#endif

#endif
