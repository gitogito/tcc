#include <stdio.h>
#include "sparse_matrix.h"
#include "mem.h"

#define NTH	100

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

void dvec_free(DenseVector *v)
{
    if (v == NULL)
	return;
    FREE(v->val);
    FREE(v);
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
    p->rows = EALLOCN(SpMatRow, n);
    for (i = 0; i < n; ++i) {
	p->rows[i].type = ROW_LIST;
	p->rows[i].urow.list = NULL;
    }

    return p;
}

void smat_free(SparseMatrix *m)
{
    int i;
    SpMatElem *p, *p2;

    if (m == NULL)
	return;
    for (i = 0; i < m->n; ++i) {
	if (m->rows[i].type == ROW_ARY) {
	    p = m->rows[i].urow.list;
	    while (p != NULL) {
		p2 = p;
		p = p->next;
		FREE(p2);
	    }
	} else {
	    FREE(m->rows[i].urow.ary);
	}
    }
    FREE(m->rows);
    FREE(m);
}

void smat_set(SparseMatrix *m, int i, int j, double val)
{
    SpMatElem *pre, *p;

    if (m->rows[i].type == ROW_ARY) {
	if (val != 0.0) {
	    m->rows[i].urow.ary[j] = val;
	    ++(m->nnz);
	}
	return;
    }

    if (m->rows[i].urow.list == NULL) {
	if (val == 0.0)
	    return;
	m->rows[i].urow.list = elem_new(val, j, NULL);
	++(m->nnz);
	return;
    }

    /* first elem */
    if (j < m->rows[i].urow.list->idx) {
	if (val != 0.0) {
	    m->rows[i].urow.list = elem_new(val, j, m->rows[i].urow.list);
	    ++(m->nnz);
	}
	return;
    } else if (j == m->rows[i].urow.list->idx) {
	if (val == 0.0) {
	    p = m->rows[i].urow.list;
	    m->rows[i].urow.list = m->rows[i].urow.list->next;
	    FREE(p);
	    --(m->nnz);
	} else {
	    m->rows[i].urow.list->val = val;
	}
	return;
    }

    /* second and later elements */
    for (pre = m->rows[i].urow.list; pre->next != NULL; pre = pre->next) {
	if (j < pre->next->idx) {
	    if (val != 0.0) {
		pre->next = elem_new(val, j, pre->next);
		++(m->nnz);
	    }
	    return;
	} else if (j == pre->next->idx) {
	    if (val == 0.0) {
		p = pre->next;
		pre->next = pre->next->next;
		FREE(p);
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
    int n;
    double *ary;
    int jj;

    if (m->rows[i].type == ROW_ARY) {
	return m->rows[i].urow.ary[j];
    }

    n = 0;
    for (elem = m->rows[i].urow.list; elem != NULL; elem = elem->next) {
	if (j == elem->idx) {
	    if (n < NTH) {
		return elem->val;
	    }
	    ary = EALLOCN(double, m->n);
	    for (jj = 0; jj < m->n; ++jj)
		ary[jj] = 0.0;
	    for (elem = m->rows[i].urow.list; elem != NULL; elem = elem->next)
		ary[elem->idx] = elem->val;
	    m->rows[i].type = ROW_ARY;
	    m->rows[i].urow.ary = ary;
	    return m->rows[i].urow.ary[j];
	}
	++n;
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
    int i, j, k, n;
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
	if (a0->rows[i].type == ROW_ARY) {
	    for (j = 0; j < a0->n; ++j) {
		if (a0->rows[i].urow.ary[j] != 0.0) {
		    ++n;
		    ai[k] = j;
		    ax[k] = a0->rows[i].urow.ary[j];
		    ++k;
		}
	    }
	} else {
	    for (elem = a0->rows[i].urow.list; elem != NULL; elem = elem->next) {
		++n;
		ai[k] = elem->idx;
		ax[k] = elem->val;
		++k;
	    }
	}
	ap[i + 1] = n;
	b[i] = b0->val[i];
    }
}
