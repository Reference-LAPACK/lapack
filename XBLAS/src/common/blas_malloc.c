#include <stdlib.h>
#include "blas_extended.h"

void *blas_malloc(size_t size)
{
#ifdef BLAS_DEBUG
  void *ptr = malloc(size);
  if (size % sizeof(float) == 0) {
    int n = size / sizeof(float);
    int i;
    for (i = 0; i < n; i++) { 
      ((float *) ptr)[i] = 0.0 / 0.0;
    }
  }
  return ptr;
#else
  return (malloc(size));
#endif
}

void blas_free(void *ptr)
{
  if (ptr != NULL)
    free(ptr);
}

void *blas_realloc(void *ptr, size_t size)
{
  return (realloc(ptr, size));
}
