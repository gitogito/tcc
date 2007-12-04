#ifndef _MEM_H_
#define _MEM_H_

#include "config.h"

#define EALLOCN(type, n)	((type *) emalloc(sizeof(type) * (n)))
#define EALLOC(type)		(EALLOCN(type, 1))
#define EALLOCN_ATOMIC(type, n)	((type *) emalloc_atomic(sizeof(type) * (n)))
#define EALLOC_ATOMIC(type)	(EALLOCN_ATOMIC(type, 1))
#define ALLOCATE_3D(var, type, ni, nj, nk)	do { \
    int i, j; \
    (var) = EALLOCN(type **, (ni)); \
    (var)[0] = EALLOCN(type *, (ni) * (nj)); \
    (var)[0][0] = EALLOCN(type, (ni) * (nj) * (nk)); \
    for (i = 0; i < (ni); ++i) { \
	(var)[i] = (var)[0] + (nj) * i; \
	(var)[i][0] = (var)[0][0] + ((nj) * (nk)) * i; \
	for (j = 0; j < (nj); ++j) { \
	    (var)[i][j] = (var)[i][0] + (nk) * j; \
	} \
    } \
} while (0)
#define ALLOCATE_3D2(var, type, ni, nj, nk, init)	do { \
    int i, j, k; \
    ALLOCATE_3D((var), type, (ni), (nj), (nk)); \
    for (i = 0; i < (ni); ++i) { \
	for (j = 0; j < (nj); ++j) { \
	    for (k = 0; k < (nk); ++k) { \
		(var)[i][j][k] = init; \
	    } \
	} \
    } \
} while (0)
#define FREE_3D(var)	do { \
    FREE((var)[0][0]); \
    FREE((var)[0]); \
    FREE(var); \
} while (0)

#ifdef HAVE_LIBGC

#include <gc.h>

#define MALLOC(size)		GC_MALLOC(size)
#define MALLOC_ATOMIC(size)	GC_MALLOC_ATOMIC(size)
#define CALLOC(number, size)	GC_MALLOC((number) * (size))
#define REALLOC(ptr, size)	GC_REALLOC((ptr), (size))
#define FREE(ptr)		GC_FREE(ptr)

#else

#include <stdlib.h>

#define MALLOC(size)		malloc(size)
#define MALLOC_ATOMIC(size)	malloc(size)
#define CALLOC(number, size)	calloc((number), (size))
#define REALLOC(ptr, size)	realloc((ptr), (size))
#define FREE(ptr)		free(ptr)

#endif

void *emalloc(size_t size);
void *emalloc_atomic(size_t size);
void *ecalloc(size_t number, size_t size);
void *erealloc(void *ptr, size_t size);
char *estrdup(const char *s);
char *estrndup(const char *s, size_t len);

#endif
