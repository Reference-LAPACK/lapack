#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_cdot_c_s(enum blas_conj_type conj, int n, const void *alpha,
		   const void *x, int incx, const void *beta,
		   const float *y, int incy, void *r)

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
 * y      (input) const float*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input/output) void*
 * 
 */
{
  static const char routine_name[] = "BLAS_cdot_c_s";

  int i, ix = 0, iy = 0;
  float *r_i = (float *) r;
  const float *x_i = (float *) x;
  const float *y_i = y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float x_ii[2];
  float y_ii;
  float r_v[2];
  float prod[2];
  float sum[2];
  float tmp1[2];
  float tmp2[2];


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

  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incy;

  if (conj == blas_conj) {
    for (i = 0; i < n; ++i) {
      x_ii[0] = x_i[ix];
      x_ii[1] = x_i[ix + 1];
      y_ii = y_i[iy];
      x_ii[1] = -x_ii[1];
      {
	prod[0] = x_ii[0] * y_ii;
	prod[1] = x_ii[1] * y_ii;
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
      y_ii = y_i[iy];

      {
	prod[0] = x_ii[0] * y_ii;
	prod[1] = x_ii[1] * y_ii;
      }				/* prod = x[i]*y[i] */
      sum[0] = sum[0] + prod[0];
      sum[1] = sum[1] + prod[1];	/* sum = sum+prod */
      ix += incx;
      iy += incy;
    }				/* endfor */
  }

  {
    tmp1[0] = sum[0] * alpha_i[0] - sum[1] * alpha_i[1];
    tmp1[1] = sum[0] * alpha_i[1] + sum[1] * alpha_i[0];
  }
  /* tmp1 = sum*alpha */
  {
    tmp2[0] = r_v[0] * beta_i[0] - r_v[1] * beta_i[1];
    tmp2[1] = r_v[0] * beta_i[1] + r_v[1] * beta_i[0];
  }
  /* tmp2 = r*beta */
  tmp1[0] = tmp1[0] + tmp2[0];
  tmp1[1] = tmp1[1] + tmp2[1];	/* tmp1 = tmp1+tmp2 */
  ((float *) r)[0] = tmp1[0];
  ((float *) r)[1] = tmp1[1];	/* r = tmp1 */



}
