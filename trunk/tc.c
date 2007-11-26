#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <float.h>
#include "tc.h"
#include "sim.h"
#include "solvele.h"
#include "config.h"
#ifdef HAVE_LIBGC
#include "gc.h"
#endif

char *prgname;

int opt_e;
int opt_o;
int opt_u;
int opt_v;
int opt_y;

void bug(char *fmt, ...)
{
    va_list ap;

    fprintf(stderr, "%s: [BUG] ", prgname);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    fflush(stderr);
    abort();
}

void warn(char *fmt, ...)
{
    va_list ap;

    fprintf(stderr, "%s: ", prgname);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    fflush(stderr);
}

void warn_exit(char *fmt, ...)
{
    va_list ap;

    fprintf(stderr, "%s: ", prgname);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    fflush(stderr);
    exit(1);
}

int main(int argc, char **argv)
{
    char *s;
    char *fname;
    double *sol;
    iPoint ipoint;
    IP_TYPE ni, nj, nk;
    double x0, y0, z0;
    double dx, dy, dz;
    IP_TYPE i, j, k;
    double x, y, z;
    double val;
    double min, max;
    int act;

#ifdef HAVE_LIBGC
    GC_init();
#endif

    prgname = argv[0];
    --argc;
    ++argv;

    opt_e = 0;
    opt_o = 0;
    opt_u = 0;
    opt_v = 0;
    opt_y = 0;

    while (argc > 0) {
	if (argv[0][0] != '-')
	    break;
	for (s = argv[0] +1; *s != '\0'; s++) {
	    if (*s == 'e') {
		opt_e = 1;
		if (*(s + 1) != '\0')
		    warn_exit("option 'e' must have an argument");
		--argc;
		++argv;
		eps_sor = atof(argv[0]);
		if (eps_sor <= 0.0)
		    warn_exit("invalid eps %g", eps_sor);
		break;
	    } else if (*s == 'o') {
		opt_o = 1;
		if (*(s + 1) != '\0')
		    warn_exit("option 'o' must have an argument");
		--argc;
		++argv;
		omega_sor = atof(argv[0]);
		if (omega_sor <= 0.0 || omega_sor >= 2.0)
		    warn_exit("invalid omega %g", omega_sor);
		break;
	    } else if (*s == 'u') {
		opt_u = 1;
	    } else if (*s == 'v') {
		opt_v = 1;
	    } else if (*s == 'y') {
		opt_y = 1;
	    } else {
		warn_exit("invalid option: %c", *s);
	    }
	}
	--argc;
	++argv;
    }

    if (argc > 0) {
	fname = *argv;
	--argc;
	++argv;
    } else {
	fname = NULL;
    }
    if (argc > 0)
	warn_exit("too many args");

    sim = sim_new(fname);
    sol = sim_calc();

    ni = world->ni;
    nj = world->nj;
    nk = world->nk;
    x0 = world->x0;
    y0 = world->y0;
    z0 = world->z0;
    dx = world->dx;
    dy = world->dy;
    dz = world->dz;

    max = DBL_MIN;
    min = DBL_MAX;
    for (k = 0; k < nk; ++k) {
	z = z0 + dz * k;
	for (j = 0; j < nj; ++j) {
	    y = y0 + dy * j;
	    for (i = 0; i < ni; ++i) {
		x = x0 + dx * i;
		ipoint.i = i;
		ipoint.j = j;
		ipoint.k = k;
		if (sim_active_p(&ipoint)) {
		    val = sol[world_to_index(&ipoint)];
		    if (val > max)
			max = val;
		    if (val < min)
			min = val;
		}
	    }
	}
    }

    printf("# %d\t%g\t%g\n", ni, x0, x0 + world->xlen);
    printf("# %d\t%g\t%g\n", nj, y0, y0 + world->ylen);
    printf("# %d\t%g\t%g\n", nk, z0, z0 + world->zlen);
    for (k = 0; k < nk; ++k) {
	z = z0 + dz * k;
	for (j = 0; j < nj; ++j) {
	    y = y0 + dy * j;
	    for (i = 0; i < ni; ++i) {
		x = x0 + dx * i;
		ipoint.i = i;
		ipoint.j = j;
		ipoint.k = k;
		if (sim_active_p(&ipoint)) {
		    val = sol[world_to_index(&ipoint)];
		    act = 1;
		} else {
		    val = min - 0.2 * (max - min);
		    act = 0;
		}
		printf("%g\t%g\t%g\t%g\t%d\n", x, y, z, val, act);
	    }
	    putchar('\n');
	}
    }

    return 0;
}
