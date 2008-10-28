#ifndef BLAS_MALLOC_H
#define BLAS_MALLOC_H 1

#ifdef NO_BLASMALLOC

#include <malloc.h>
#define blas_malloc(s) malloc(s)
#define blas_realloc(p, s) realloc(p,s)
#define blas_free(p) free(p)

#else

/* stddef is needed for size_t */
#include <stddef.h>
void  *blas_malloc(size_t size);
void *blas_realloc(void *p, size_t size);
void blas_free(void *p);


#endif

#endif
