#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void sprint_vector(const float *x, int n, int inc, const char *name)
{
  const float *x_i = x;
  int i, xi;


  xi = (inc >= 0) ? 0 : (1 - n) * inc;

  if (name) {
    printf("%s = ", name);
  }
  printf("[\n");
  for (i = 0; i < n; i++, xi += inc) {
    printf("  ");
    printf("%16.8e", x_i[xi]);
    printf("\n");
  }
  printf("];\n");

}
void dprint_vector(const double *x, int n, int inc, const char *name)
{
  const double *x_i = x;
  int i, xi;


  xi = (inc >= 0) ? 0 : (1 - n) * inc;

  if (name) {
    printf("%s = ", name);
  }
  printf("[\n");
  for (i = 0; i < n; i++, xi += inc) {
    printf("  ");
    printf("%24.16e", x_i[xi]);
    printf("\n");
  }
  printf("];\n");

}
void cprint_vector(const void *x, int n, int inc, const char *name)
{
  const float *x_i = (float *) x;
  int i, xi;

  inc *= 2;
  xi = (inc >= 0) ? 0 : (1 - n) * inc;

  if (name) {
    printf("%s = ", name);
  }
  printf("[\n");
  for (i = 0; i < n; i++, xi += inc) {
    printf("  ");
    printf("(%16.8e, %16.8e)", x_i[xi], x_i[xi + 1]);
    printf("\n");
  }
  printf("];\n");

}
void zprint_vector(const void *x, int n, int inc, const char *name)
{
  const double *x_i = (double *) x;
  int i, xi;

  inc *= 2;
  xi = (inc >= 0) ? 0 : (1 - n) * inc;

  if (name) {
    printf("%s = ", name);
  }
  printf("[\n");
  for (i = 0; i < n; i++, xi += inc) {
    printf("  ");
    printf("(%24.16e, %24.16e)", x_i[xi], x_i[xi + 1]);
    printf("\n");
  }
  printf("];\n");

}
