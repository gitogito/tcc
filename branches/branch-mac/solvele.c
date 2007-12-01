#include "config.h"

#ifdef HAVE_UMFPACK_H
#include <umfpack.h>
#endif
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "solvele.h"
#include "sparse_matrix.h"
#include "mem.h"
#include "tc.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define EPS	1.0e-6
#define N	20

double eps_sor;
double omega_sor;

void solvele_set_matrix(Solvele *self, int i, int j, double val)
{
    if (i < 0 || j < 0)
	warn_exit("negative index in solvele_set_matrix");
    if (i >= self->size || j >= self->size)
	warn_exit("too large index in solvele_set_matrix");
    smat_set(self->mat, i, j, val);
}

void solvele_add_matrix(Solvele *self, int i, int j, double val)
{
    solvele_set_matrix(self, i, j, smat_ref(self->mat, i, j) + val);
}

void solvele_set_vector(Solvele *self, int i, double val)
{
    if (i < 0)
	warn_exit("negative index in solvele_set_vector");
    if (i >= self->size)
	warn_exit("too large index in solvele_set_vector");
    dvec_set(self->vec, i, val);
}

void solvele_add_vector(Solvele *self, int i, double val)
{
    solvele_set_vector(self, i, dvec_ref(self->vec, i) + val);
}

void solvele_print_matrix(Solvele *self)
{
    smat_print(self->mat);
}

void solvele_print_vector(Solvele *self)
{
    dvec_print(self->vec);
}

static void sor_without_diag(double *u, int size, int ni, int nj, int nk,
	int *ap, int *ai, double *ax, double *pb)
{
    static char rotate[] = "|/-\\";
    int i, ii, index, pi;
    double eps, omega, new_val, old_val, sum_abs, sum_abs_diff;

    if (opt_v)
        warn("solving using SOR method ...");
    for (i = 0; i < size; ++i) {
        u[i] = 0.0;
    }
    if (opt_e)
        eps = eps_sor;
    else
        eps = EPS;
    if (opt_o) {
        omega = omega_sor;
    } else {
        omega = 2.0 / (1.0 + sqrt(1.0 - 1.0/3.0 * (cos(M_PI / ni) + cos(M_PI / nj) + cos(M_PI / nk))));
        if (omega < 1.9)
            omega = 1.9;
    }
    if (opt_v) {
        warn("SOR epsilon is %g", eps);
        warn("SOR relaxation factor is %g", omega);
    }
    for (ii = 0; /* do nothing */; ++ii) {
        for (index = 0; index < N; ++index) {
            for (i = 0; i < size; ++i) {
                new_val = pb[i];
                for (pi = ap[i]; pi < ap[i + 1]; ++pi) {
		    new_val -= ax[pi] * u[ai[pi]];
                }
                old_val = u[i];
                u[i] += omega * (new_val - old_val);
            }
        }

        sum_abs = sum_abs_diff = 0.0;
        for (i = 0; i < size; ++i) {
            new_val = pb[i];
            for (pi = ap[i]; pi < ap[i + 1]; ++pi) {
		new_val -= ax[pi] * u[ai[pi]];
            }
            old_val = u[i];
            u[i] += omega * (new_val - old_val);
	    sum_abs += fabs(new_val);
	    sum_abs_diff += fabs(new_val - old_val);
	}
        if (sum_abs_diff / sum_abs < eps) {
	    if (opt_v)
		fprintf(stderr, "\rfinished\n");
            break;
        }
	if (opt_v)
	    fprintf(stderr, "\r%c %5.1f%%", rotate[ii % (sizeof(rotate) - 1)],
		    100.0 * eps / (sum_abs_diff / sum_abs));
    }
}

double *solvele_solve(Solvele *self, int ni, int nj, int nk)
{
    int *ap;
    int *ai;
    double *ax;
    double *pb;
    double *u;
#ifdef HAVE_UMFPACK_H
    void *Symbolic, *Numeric;
#endif

    u = EALLOCN(double, self->size);

    if (opt_u) {
#ifdef HAVE_UMFPACK_H
	get_crs(self->mat, self->vec, &ap, &ai, &ax, &pb);
	if (opt_v)
	    warn("solving using UMFPACK ...");
	umfpack_di_symbolic(self->size, self->size, ap, ai, ax, &Symbolic, NULL, NULL);
	umfpack_di_numeric(ap, ai, ax, Symbolic, &Numeric, NULL, NULL);
	umfpack_di_free_symbolic(&Symbolic);
	umfpack_di_solve(UMFPACK_At, ap, ai, ax, u, pb, Numeric, NULL, NULL);
	umfpack_di_free_numeric (&Numeric);
#else
	warn_exit("solver with UMFPACK is not implemented in solvele_solve");
#endif
    } else {
	get_crs_without_diag(self->mat, self->vec, &ap, &ai, &ax, &pb);
        sor_without_diag(u, self->size, ni, nj, nk, ap, ai, ax, pb);
    }
    return u;
}

Solvele *solvele_new(int size)
{
    Solvele *self;

    self = EALLOC(Solvele);
    self->size = size;
    self->mat = smat_new(size);
    self->vec = dvec_new(size);
    return self;
}
