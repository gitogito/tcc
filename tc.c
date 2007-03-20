#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "sim.h"

char *prgname;

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
    Sim *sim;
    Array3Dd ary;
    int ni, nj, nk;
    double dx, dy, dz;
    int i, j, k;
    double x, y, z;
    double val;
    int act;

    prgname = argv[0];
    --argc;
    ++argv;

    sim = sim_new();
    ary = sim_calc(sim);

    ni = sim->world->ni;
    nj = sim->world->nj;
    nk = sim->world->nk;
    dx = sim->world->dx;
    dy = sim->world->dy;
    dz = sim->world->dz;

    printf("# %d\t%g\t%g\n", ni, sim->world->x, sim->world->xlen);
    printf("# %d\t%g\t%g\n", nj, sim->world->y, sim->world->ylen);
    printf("# %d\t%g\t%g\n", nk, sim->world->z, sim->world->zlen);
    for (k = 0; k < nk; ++k) {
        z = dz * k;
        for (j = 0; j < nj; ++j) {
            y = dy * j;
            for (i = 0; i < ni; ++i) {
                x = dx * i;
                if (sim_active_p(sim, get_point(i, j, k))) {
                    val = ary[i][j][k];
                    act = 1;
                } else {
                    val = -1.0;
                    act = 0;
                }
                printf("%g\t%g\t%g\t%g\t%d\n", x, y, z, val, act);
            }
            putchar('\n');
        }
    }

    return 0;
}
