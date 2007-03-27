#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "sim.h"
#include "solvele.h"

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
    Sim *sim;
    Array3Dd ary;
    int ni, nj, nk;
    double x0, y0, z0;
    double dx, dy, dz;
    int i, j, k;
    double x, y, z;
    double val;
    int act;

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
		argc--;
		argv++;
		eps_sor = atof(argv[0]);
		if (eps_sor <= 0.0)
		    warn_exit("invalid eps %g", eps_sor);
		break;
	    } else if (*s == 'o') {
		opt_o = 1;
		if (*(s + 1) != '\0')
		    warn_exit("option 'o' must have an argument");
		argc--;
		argv++;
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
	argc--;
	argv++;
    }

    sim = sim_new();
    ary = sim_calc(sim);

    ni = sim->world->ni;
    nj = sim->world->nj;
    nk = sim->world->nk;
    x0 = sim->world->x0;
    y0 = sim->world->y0;
    z0 = sim->world->z0;
    dx = sim->world->dx;
    dy = sim->world->dy;
    dz = sim->world->dz;

    printf("# %d\t%g\t%g\n", ni, x0, x0 + sim->world->xlen);
    printf("# %d\t%g\t%g\n", nj, y0, y0 + sim->world->ylen);
    printf("# %d\t%g\t%g\n", nk, z0, z0 + sim->world->zlen);
    for (k = 0; k < nk; ++k) {
        z = z0 + dz * k;
        for (j = 0; j < nj; ++j) {
            y = y0 + dy * j;
            for (i = 0; i < ni; ++i) {
                x = x0 + dx * i;
                if (sim_active_p(sim, get_ipoint(i, j, k))) {
                    val = ary[i][j][k];
                    act = 1;
                } else {
                    val = 0.0;
                    act = 0;
                }
                printf("%g\t%g\t%g\t%g\t%d\n", x, y, z, val, act);
            }
            putchar('\n');
        }
    }

    return 0;
}
