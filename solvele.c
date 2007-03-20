#include "config.h"

#ifdef HAVE_UMFPACK_H
#include "umfpack.h"
#endif
#include <math.h>
#include "solvele.h"
#include "sparse_matrix.h"
#include "mem.h"
#include "tc.h"

#define DBL_EPS	1.0e-15
#define EPS	1.0e-6
#define OMEGA	1.0

void solvele_set_matrix(Solvele *self, int i, int j, double val)
{
    if (i < 0 || j < 0)
	warn_exit("negative index in solvele_set_matrix");
    if (i >= self->size || j >= self->size)
	warn_exit("too large index in solvele_set_matrix");
    smat_set(self->mat, i, j, val);
}

void solvele_set_vector(Solvele *self, int i, double val)
{
    if (i < 0)
	warn_exit("negative index in solvele_set_vector");
    if (i >= self->size)
	warn_exit("too large index in solvele_set_vector");
    dvec_set(self->vec, i, val);
}

void solvele_print_matrix(Solvele *self)
{
    smat_print(self->mat);
}

void solvele_print_vector(Solvele *self)
{
    dvec_print(self->vec);
}

double *solvele_solve(Solvele *self, int use_umfpack_p)
{
    int *ap;
    int *ai;
    double *ax;
    double *pb;
    double *u;
#ifdef HAVE_UMFPACK_H
    void *Symbolic, *Numeric;
#endif
    int size, ok, pi, i, j;
    double v, old_val, new_val, c0;

    get_crs(self->mat, self->vec, &ap, &ai, &ax, &pb);

    u = EALLOCN(double, self->size);

    if (use_umfpack_p) {
#ifdef HAVE_UMFPACK_H
        umfpack_di_symbolic(self->size, self->size, ap, ai, ax, &Symbolic, NULL, NULL);
        umfpack_di_numeric(ap, ai, ax, Symbolic, &Numeric, NULL, NULL);
        umfpack_di_free_symbolic(&Symbolic);
        umfpack_di_solve(UMFPACK_At, ap, ai, ax, u, pb, Numeric, NULL, NULL);
        umfpack_di_free_numeric (&Numeric);
#else
        warn_exit("solver with UMFPACK is not implemented in solvele_solve");
#endif
    } else {
        size = self->size;
        for (i = 0; i < size; ++i) {
            u[i] = 0.0;
        }
        for (;;) {
            ok = 1;
            for (i = 0; i < size; ++i) {
                v = pb[i];
                for (pi = ap[i]; pi < ap[i + 1]; ++pi) {
                    j = ai[pi];
                    if (i == j)
                        c0 = ax[pi];
                    else
                        v -= ax[pi] * u[j];
                }
                new_val = v / c0;
                old_val = u[i];
                u[i] += OMEGA * (new_val - old_val);
                if (ok && fabs(new_val) > DBL_EPS &&
                        fabs(new_val - old_val) > EPS * fabs(new_val))
                {
                    ok = 0;
                }
            }
            if (ok)
                break;
        }
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
