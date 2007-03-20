#include <stdio.h>
#include "sparse_matrix.h"
#include "mem.h"

DenseVector *dvec_new(int n)
{
    DenseVector *p;
    int i;

    p = EALLOC(DenseVector);
    p->n = n;
    p->val = EALLOCN(double, n);
    for (i = 0; i < n; ++i)
	p->val[i] = 0.0;

    return p;
}

void dvec_set(DenseVector *v, int i, double val)
{
    v->val[i] = val;
}

double dvec_ref(DenseVector *v, int i)
{
    return v->val[i];
}

void dvec_print(DenseVector *v)
{
    int i;

    printf("---\n");
    for (i = 0; i < v->n; ++i)
	printf("\t%g", v->val[i]);
    printf("\n");
}

static SpMatElem *elem_new(double val, int idx, SpMatElem *next)
{
    SpMatElem *p;

    p = EALLOC(SpMatElem);
    p->val = val;
    p->idx = idx;
    p->next = next;

    return p;
}

SparseMatrix *smat_new(int n)
{
    SparseMatrix *p;
    int i;

    p = EALLOC(SparseMatrix);
    p->n = n;
    p->nnz = 0;
    p->elems = EALLOCN(SpMatElem *, n);
    for (i = 0; i < n; ++i)
	p->elems[i] = NULL;

    return p;
}

void smat_set(SparseMatrix *m, int i, int j, double val)
{
    SpMatElem *pre;

    if (m->elems[i] == NULL) {
	if (val == 0.0)
	    return;
	m->elems[i] = elem_new(val, j, NULL);
	++(m->nnz);
	return;
    }

    /* first elem */
    if (j < m->elems[i]->idx) {
	if (val != 0.0) {
	    m->elems[i] = elem_new(val, j, m->elems[i]);
	    ++(m->nnz);
	}
	return;
    } else if (j == m->elems[i]->idx) {
	if (val == 0.0) {
	    m->elems[i] = m->elems[i]->next;
	    --(m->nnz);
	} else {
	    m->elems[i]->val = val;
	}
	return;
    }

    /* second and later elements */
    for (pre = m->elems[i]; pre->next != NULL; pre = pre->next) {
	if (j < pre->next->idx) {
	    if (val != 0.0) {
		pre->next = elem_new(val, j, pre->next);
		++(m->nnz);
	    }
	    return;
	} else if (j == pre->next->idx) {
	    if (val == 0.0) {
		pre->next = pre->next->next;
		--(m->nnz);
	    } else {
		pre->next->val = val;
	    }
	    return;
	}
    }

    /* add elem to last */
    if (val != 0.0) {
	pre->next = elem_new(val, j, NULL);
	++(m->nnz);
    }
}

double smat_ref(SparseMatrix *m, int i, int j)
{
    SpMatElem *elem;

    for (elem = m->elems[i]; elem != NULL; elem = elem->next) {
	if (j == elem->idx)
	    return elem->val;
    }
    return 0.0;
}

void smat_print(SparseMatrix *m)
{
    int i, j;

    printf("===\n");
    for (i = 0; i < m->n; ++i) {
	for (j = 0; j < m->n; ++j) {
	    printf("\t%g", smat_ref(m, i, j));
	}
	printf("\n");
    }
}

void get_crs(SparseMatrix *a0, DenseVector *b0,
	int **pap, int **pai, double **pax, double **pb)
{
    int *ap;
    int *ai;
    double *ax;
    double *b;
    int i, k, n;
    SpMatElem *elem;

    *pap = EALLOCN(int, a0->n + 1);
    *pai = EALLOCN(int, a0->nnz);
    *pax = EALLOCN(double, a0->nnz);
    *pb = EALLOCN(double, a0->n);

    ap = *pap;
    ai = *pai;
    ax = *pax;
    b = *pb;

    ap[0] = 0;
    n = 0;
    k = 0;
    for (i = 0; i < a0->n; ++i) {
	for (elem = a0->elems[i]; elem != NULL; elem = elem->next) {
	    ++n;
	    ai[k] = elem->idx;
	    ax[k] = elem->val;
	    ++k;
	}
	ap[i + 1] = n;
	b[i] = b0->val[i];
    }
}
