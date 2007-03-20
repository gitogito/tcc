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

    prgname = argv[0];
    --argc;
    ++argv;

    sim = sim_new();
    ary = sim_calc(sim);

    return 0;
}
