#include <string.h>
#include "mem.h"
#include "tcc.h"

void *emalloc(size_t size)
{
    void *p;

    p = MALLOC(size);
    if (p == NULL)
	warn_exit("can't allocate memory");
    return p;
}

void *emalloc_atomic(size_t size)
{
    void *p;

    p = MALLOC_ATOMIC(size);
    if (p == NULL)
	warn_exit("can't allocate memory");
    return p;
}

void *ecalloc(size_t number, size_t size)
{
    void *p;

    p = CALLOC(number, size);
    if (p == NULL)
	warn_exit("can't allocate memory");
    return p;
}

void *erealloc(void *ptr, size_t size)
{
    void *p;

    p = REALLOC(ptr, size);
    if (p == NULL)
	warn_exit("can't allocate memory");
    return p;
}

char *estrdup(const char *s)
{
    char *p;

    p = EALLOCN_ATOMIC(char, strlen(s) + 1);
    if (p == NULL)
	warn_exit("can't allocate memory");
    strcpy(p, s);
    return p;
}

char *estrndup(const char *s, size_t len)
{
    char *p;

    p = EALLOCN_ATOMIC(char, len + 1);
    if (p == NULL)
	warn_exit("can't allocate memory");
    memcpy(p, s, len + 1);
    p[len] = '\0';
    return p;
}
