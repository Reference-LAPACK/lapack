#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zdot_z_c(enum blas_conj_type conj, int n, const void *alpha,
		   const void *x, int incx, const void *beta,
		   const void *y, int incy, void *r)

/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *     r <- beta * r + alpha * SUM_{i=0, n-1} x[i] * y[i].
 * 
 * Arguments
 * =========
 *  
 * conj   (input) enum blas_conj_type
 *        When x and y are complex vectors, specifies whether vector
 *        components x[i] are used unconjugated or conjugated. 
 * 
 * n      (input) int
 *        The length of vectors x and y.
 * 
 * alpha  (input) const void*
 * 
 * x      (input) const void*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * beta   (input) const void*
 *
 * y      (input) const void*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input/output) void*
 * 
 */
{
  static const char routine_name[] = "BLAS_zdot_z_c";

  int i, ix = 0, iy = 0;
  double *r_i = (double *) r;
  const double *x_i = (double *) x;
  const float *y_i = (float *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double x_ii[2];
  float y_ii[2];
  double r_v[2];
  double prod[2];
  double sum[2];
  double tmp1[2];
  double tmp2[2];


  /* Test the input parameters. */
  if (n < 0)
    BLAS_error(routine_name, -2, n, NULL);
  else if (incx == 0)
    BLAS_error(routine_name, -5, incx, NULL);
  else if (incy == 0)
    BLAS_error(routine_name, -8, incy, NULL);

  /* Immediate return. */
  if (((beta_i[0] == 1.0 && beta_i[1] == 0.0))
      && (n == 0 || (alpha_i[0] == 0.0 && alpha_i[1] == 0.0)))
    return;



  r_v[0] = r_i[0];
  r_v[1] = r_i[0 + 1];
  sum[0] = sum[1] = 0.0;
  incx *= 2;
  incy *= 2;
  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incy;

  if (conj == blas_conj) {
    for (i = 0; i < n; ++i) {
      x_ii[0] = x_i[ix];
      x_ii[1] = x_i[ix + 1];
      y_ii[0] = y_i[iy];
      y_ii[1] = y_i[iy + 1];
      x_ii[1] = -x_ii[1];
      {
	prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
	prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
      }				/* prod = x[i]*y[i] */
      sum[0] = sum[0] + prod[0];
      sum[1] = sum[1] + prod[1];	/* sum = sum+prod */
      ix += incx;
      iy += incy;
    }				/* endfor */
  } else {
    /* do not conjugate */

    for (i = 0; i < n; ++i) {
      x_ii[0] = x_i[ix];
      x_ii[1] = x_i[ix + 1];
      y_ii[0] = y_i[iy];
      y_ii[1] = y_i[iy + 1];

      {
	prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
	prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
      }				/* prod = x[i]*y[i] */
      sum[0] = sum[0] + prod[0];
      sum[1] = sum[1] + prod[1];	/* sum = sum+prod */
      ix += incx;
      iy += incy;
    }				/* endfor */
  }

  {
    tmp1[0] = (double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
    tmp1[1] = (double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
  }				/* tmp1 = sum*alpha */
  {
    tmp2[0] = (double) r_v[0] * beta_i[0] - (double) r_v[1] * beta_i[1];
    tmp2[1] = (double) r_v[0] * beta_i[1] + (double) r_v[1] * beta_i[0];
  }				/* tmp2 = r*beta */
  tmp1[0] = tmp1[0] + tmp2[0];
  tmp1[1] = tmp1[1] + tmp2[1];	/* tmp1 = tmp1+tmp2 */
  ((double *) r)[0] = tmp1[0];
  ((double *) r)[1] = tmp1[1];	/* r = tmp1 */



}
