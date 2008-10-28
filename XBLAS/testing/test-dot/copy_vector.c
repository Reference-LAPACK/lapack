#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

void scopy_vector(const float *x, int n, int incx, float *y, int incy)
{
  const float *x_i = x;
  float *y_i = y;
  int i, xi, yi;



  xi = (incx >= 0) ? 0 : (1 - n) * incx;
  yi = (incy >= 0) ? 0 : (1 - n) * incy;

  for (i = 0; i < n; i++, xi += incx, yi += incy) {
    y_i[yi] = x_i[xi];
  }

}
void dcopy_vector(const double *x, int n, int incx, double *y, int incy)
{
  const double *x_i = x;
  double *y_i = y;
  int i, xi, yi;



  xi = (incx >= 0) ? 0 : (1 - n) * incx;
  yi = (incy >= 0) ? 0 : (1 - n) * incy;

  for (i = 0; i < n; i++, xi += incx, yi += incy) {
    y_i[yi] = x_i[xi];
  }

}
void ccopy_vector(const void *x, int n, int incx, void *y, int incy)
{
  const float *x_i = (float *) x;
  float *y_i = (float *) y;
  int i, xi, yi;

  incx *= 2;
  incy *= 2;
  xi = (incx >= 0) ? 0 : (1 - n) * incx;
  yi = (incy >= 0) ? 0 : (1 - n) * incy;

  for (i = 0; i < n; i++, xi += incx, yi += incy) {
    y_i[yi] = x_i[xi];
    y_i[yi + 1] = x_i[xi + 1];
  }

}
void zcopy_vector(const void *x, int n, int incx, void *y, int incy)
{
  const double *x_i = (double *) x;
  double *y_i = (double *) y;
  int i, xi, yi;

  incx *= 2;
  incy *= 2;
  xi = (incx >= 0) ? 0 : (1 - n) * incx;
  yi = (incy >= 0) ? 0 : (1 - n) * incy;

  for (i = 0; i < n; i++, xi += incx, yi += incy) {
    y_i[yi] = x_i[xi];
    y_i[yi + 1] = x_i[xi + 1];
  }

}
