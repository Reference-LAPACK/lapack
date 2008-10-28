#include "blas_extended.h"
#include "blas_extended_private.h"

void s_r_truth2(enum blas_conj_type conj, int n, float alpha,
		const float *x, int incx, float beta,
		const float *head_y, const float *tail_y, int incy, float *r,
		double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *   r_true <- beta * r + alpha * SUM_{i=0, n-1} x[i] * (head_y[i] + tail_y[i])
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
 * alpha  (input) float
 * 
 * x      (input) const float*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * beta   (input) float
 *
 * head_y
 * tail_y (input) const float*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input) float*
 *
 * head_r_true
 * tail_r_true (output) float*
 *        The truth.
 * 
 */
{

  int i, ix = 0, iy = 0;
  float *r_i = r;
  const float *x_i = x;
  const float *head_y_i = head_y;
  const float *tail_y_i = tail_y;
  float alpha_i = alpha;
  float beta_i = beta;
  float x_ii;
  float y_ii;
  float r_v;
  double head_prod, tail_prod;
  double head_sum, tail_sum;
  double head_tmp1, tail_tmp1;
  double head_tmp2, tail_tmp2;
  FPU_FIX_DECL;

  /* Immediate return */
  if (n < 0) {
    *head_r_true = *tail_r_true = 0.0;
    return;
  }

  FPU_FIX_START;

  r_v = r_i[0];
  head_sum = tail_sum = 0.0;


  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incy;

  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = head_y_i[iy];

    head_prod = (double) x_ii *y_ii;
    tail_prod = 0.0;		/* prod = x[i] * head_y[i] */
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_sum + head_prod;
      bv = s1 - head_sum;
      s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_sum + tail_prod;
      bv = t1 - tail_sum;
      t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_sum = t1 + t2;
      tail_sum = t2 - (head_sum - t1);
    }				/* sum = sum+prod */
    y_ii = tail_y_i[iy];
    head_prod = (double) x_ii *y_ii;
    tail_prod = 0.0;		/* prod = x[i] * tail_y[i] */
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_sum + head_prod;
      bv = s1 - head_sum;
      s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_sum + tail_prod;
      bv = t1 - tail_sum;
      t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_sum = t1 + t2;
      tail_sum = t2 - (head_sum - t1);
    }				/* sum = sum+prod */
    ix += incx;
    iy += incy;
  }				/* endfor */


  {
    double dt = (double) alpha_i;
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

      con = head_sum * split;
      a11 = con - head_sum;
      a11 = con - a11;
      a21 = head_sum - a11;
      con = dt * split;
      b1 = con - dt;
      b1 = con - b1;
      b2 = dt - b1;

      c11 = head_sum * dt;
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = tail_sum * dt;
      t1 = c11 + c2;
      t2 = (c2 - (t1 - c11)) + c21;

      head_tmp1 = t1 + t2;
      tail_tmp1 = t2 - (head_tmp1 - t1);
    }
  }				/* tmp1 = sum*alpha */
  head_tmp2 = (double) r_v *beta_i;
  tail_tmp2 = 0.0;		/* tmp2 = r*beta */
  {
    /* Compute double-double = double-double + double-double. */
    double bv;
    double s1, s2, t1, t2;

    /* Add two hi words. */
    s1 = head_tmp1 + head_tmp2;
    bv = s1 - head_tmp1;
    s2 = ((head_tmp2 - bv) + (head_tmp1 - (s1 - bv)));

    /* Add two lo words. */
    t1 = tail_tmp1 + tail_tmp2;
    bv = t1 - tail_tmp1;
    t2 = ((tail_tmp2 - bv) + (tail_tmp1 - (t1 - bv)));

    s2 += t1;

    /* Renormalize (s1, s2)  to  (t1, s2) */
    t1 = s1 + s2;
    s2 = s2 - (t1 - s1);

    t2 += s2;

    /* Renormalize (t1, t2)  */
    head_tmp1 = t1 + t2;
    tail_tmp1 = t2 - (head_tmp1 - t1);
  }				/* tmp1 = tmp1+tmp2 */
  *head_r_true = head_tmp1;
  *tail_r_true = tail_tmp1;	/* *r = tmp1 */

  FPU_FIX_STOP;

}
void d_r_truth2(enum blas_conj_type conj, int n, double alpha,
		const double *x, int incx, double beta,
		const double *head_y, const double *tail_y, int incy,
		double *r, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *   r_true <- beta * r + alpha * SUM_{i=0, n-1} x[i] * (head_y[i] + tail_y[i])
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
 * alpha  (input) double
 * 
 * x      (input) const double*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * beta   (input) double
 *
 * head_y
 * tail_y (input) const double*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input) double*
 *
 * head_r_true
 * tail_r_true (output) double*
 *        The truth.
 * 
 */
{

  int i, ix = 0, iy = 0;
  double *r_i = r;
  const double *x_i = x;
  const double *head_y_i = head_y;
  const double *tail_y_i = tail_y;
  double alpha_i = alpha;
  double beta_i = beta;
  double x_ii;
  double y_ii;
  double r_v;
  double head_prod, tail_prod;
  double head_sum, tail_sum;
  double head_tmp1, tail_tmp1;
  double head_tmp2, tail_tmp2;
  FPU_FIX_DECL;

  /* Immediate return */
  if (n < 0) {
    *head_r_true = *tail_r_true = 0.0;
    return;
  }

  FPU_FIX_START;

  r_v = r_i[0];
  head_sum = tail_sum = 0.0;


  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incy;

  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = head_y_i[iy];

    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = x_ii * split;
      a1 = con - x_ii;
      a1 = con - a1;
      a2 = x_ii - a1;
      con = y_ii * split;
      b1 = con - y_ii;
      b1 = con - b1;
      b2 = y_ii - b1;

      head_prod = x_ii * y_ii;
      tail_prod = (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
    }				/* prod = x[i] * head_y[i] */
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_sum + head_prod;
      bv = s1 - head_sum;
      s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_sum + tail_prod;
      bv = t1 - tail_sum;
      t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_sum = t1 + t2;
      tail_sum = t2 - (head_sum - t1);
    }				/* sum = sum+prod */
    y_ii = tail_y_i[iy];
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = x_ii * split;
      a1 = con - x_ii;
      a1 = con - a1;
      a2 = x_ii - a1;
      con = y_ii * split;
      b1 = con - y_ii;
      b1 = con - b1;
      b2 = y_ii - b1;

      head_prod = x_ii * y_ii;
      tail_prod = (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
    }				/* prod = x[i] * tail_y[i] */
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_sum + head_prod;
      bv = s1 - head_sum;
      s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_sum + tail_prod;
      bv = t1 - tail_sum;
      t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_sum = t1 + t2;
      tail_sum = t2 - (head_sum - t1);
    }				/* sum = sum+prod */
    ix += incx;
    iy += incy;
  }				/* endfor */


  {
    /* Compute double-double = double-double * double. */
    double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

    con = head_sum * split;
    a11 = con - head_sum;
    a11 = con - a11;
    a21 = head_sum - a11;
    con = alpha_i * split;
    b1 = con - alpha_i;
    b1 = con - b1;
    b2 = alpha_i - b1;

    c11 = head_sum * alpha_i;
    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

    c2 = tail_sum * alpha_i;
    t1 = c11 + c2;
    t2 = (c2 - (t1 - c11)) + c21;

    head_tmp1 = t1 + t2;
    tail_tmp1 = t2 - (head_tmp1 - t1);
  }				/* tmp1 = sum*alpha */
  {
    /* Compute double_double = double * double. */
    double a1, a2, b1, b2, con;

    con = r_v * split;
    a1 = con - r_v;
    a1 = con - a1;
    a2 = r_v - a1;
    con = beta_i * split;
    b1 = con - beta_i;
    b1 = con - b1;
    b2 = beta_i - b1;

    head_tmp2 = r_v * beta_i;
    tail_tmp2 = (((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
  }				/* tmp2 = r*beta */
  {
    /* Compute double-double = double-double + double-double. */
    double bv;
    double s1, s2, t1, t2;

    /* Add two hi words. */
    s1 = head_tmp1 + head_tmp2;
    bv = s1 - head_tmp1;
    s2 = ((head_tmp2 - bv) + (head_tmp1 - (s1 - bv)));

    /* Add two lo words. */
    t1 = tail_tmp1 + tail_tmp2;
    bv = t1 - tail_tmp1;
    t2 = ((tail_tmp2 - bv) + (tail_tmp1 - (t1 - bv)));

    s2 += t1;

    /* Renormalize (s1, s2)  to  (t1, s2) */
    t1 = s1 + s2;
    s2 = s2 - (t1 - s1);

    t2 += s2;

    /* Renormalize (t1, t2)  */
    head_tmp1 = t1 + t2;
    tail_tmp1 = t2 - (head_tmp1 - t1);
  }				/* tmp1 = tmp1+tmp2 */
  *head_r_true = head_tmp1;
  *tail_r_true = tail_tmp1;	/* *r = tmp1 */

  FPU_FIX_STOP;

}
void c_r_truth2(enum blas_conj_type conj, int n, const void *alpha,
		const void *x, int incx, const void *beta,
		const void *head_y, const void *tail_y, int incy, void *r,
		double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *   r_true <- beta * r + alpha * SUM_{i=0, n-1} x[i] * (head_y[i] + tail_y[i])
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
 * head_y
 * tail_y (input) const void*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input) void*
 *
 * head_r_true
 * tail_r_true (output) void*
 *        The truth.
 * 
 */
{

  int i, ix = 0, iy = 0;
  float *r_i = (float *) r;
  const float *x_i = (float *) x;
  const float *head_y_i = (float *) head_y;
  const float *tail_y_i = (float *) tail_y;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float x_ii[2];
  float y_ii[2];
  float r_v[2];
  double head_prod[2], tail_prod[2];
  double head_sum[2], tail_sum[2];
  double head_tmp1[2], tail_tmp1[2];
  double head_tmp2[2], tail_tmp2[2];
  FPU_FIX_DECL;

  /* Immediate return */
  if (n < 0) {
    head_r_true[0] = tail_r_true[0] = head_r_true[1] = tail_r_true[1] = 0.0;
    return;
  }

  FPU_FIX_START;

  r_v[0] = r_i[0];
  r_v[1] = r_i[0 + 1];
  head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;
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
      y_ii[0] = head_y_i[iy];
      y_ii[1] = head_y_i[iy + 1];
      x_ii[1] = -x_ii[1];
      {
	double head_e1, tail_e1;
	double d1;
	double d2;
	/* Real part */
	d1 = (double) x_ii[0] * y_ii[0];
	d2 = (double) -x_ii[1] * y_ii[1];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_prod[0] = head_e1;
	tail_prod[0] = tail_e1;
	/* imaginary part */
	d1 = (double) x_ii[0] * y_ii[1];
	d2 = (double) x_ii[1] * y_ii[0];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_prod[1] = head_e1;
	tail_prod[1] = tail_e1;
      }				/* prod = x[i] * head_y[i] */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_sum[0];
	tail_a = tail_sum[0];
	head_b = head_prod[0];
	tail_b = tail_prod[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[0] = head_t;
	tail_sum[0] = tail_t;
	/* Imaginary part */
	head_a = head_sum[1];
	tail_a = tail_sum[1];
	head_b = head_prod[1];
	tail_b = tail_prod[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[1] = head_t;
	tail_sum[1] = tail_t;
      }				/* sum = sum+prod */
      y_ii[0] = tail_y_i[iy];
      y_ii[1] = tail_y_i[iy + 1];
      {
	double head_e1, tail_e1;
	double d1;
	double d2;
	/* Real part */
	d1 = (double) x_ii[0] * y_ii[0];
	d2 = (double) -x_ii[1] * y_ii[1];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_prod[0] = head_e1;
	tail_prod[0] = tail_e1;
	/* imaginary part */
	d1 = (double) x_ii[0] * y_ii[1];
	d2 = (double) x_ii[1] * y_ii[0];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_prod[1] = head_e1;
	tail_prod[1] = tail_e1;
      }				/* prod = x[i] * tail_y[i] */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_sum[0];
	tail_a = tail_sum[0];
	head_b = head_prod[0];
	tail_b = tail_prod[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[0] = head_t;
	tail_sum[0] = tail_t;
	/* Imaginary part */
	head_a = head_sum[1];
	tail_a = tail_sum[1];
	head_b = head_prod[1];
	tail_b = tail_prod[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[1] = head_t;
	tail_sum[1] = tail_t;
      }				/* sum = sum+prod */
      ix += incx;
      iy += incy;
    }				/* endfor */
  } else {
    /* do not conjugate */

    for (i = 0; i < n; ++i) {
      x_ii[0] = x_i[ix];
      x_ii[1] = x_i[ix + 1];
      y_ii[0] = head_y_i[iy];
      y_ii[1] = head_y_i[iy + 1];

      {
	double head_e1, tail_e1;
	double d1;
	double d2;
	/* Real part */
	d1 = (double) x_ii[0] * y_ii[0];
	d2 = (double) -x_ii[1] * y_ii[1];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_prod[0] = head_e1;
	tail_prod[0] = tail_e1;
	/* imaginary part */
	d1 = (double) x_ii[0] * y_ii[1];
	d2 = (double) x_ii[1] * y_ii[0];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_prod[1] = head_e1;
	tail_prod[1] = tail_e1;
      }				/* prod = x[i] * head_y[i] */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_sum[0];
	tail_a = tail_sum[0];
	head_b = head_prod[0];
	tail_b = tail_prod[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[0] = head_t;
	tail_sum[0] = tail_t;
	/* Imaginary part */
	head_a = head_sum[1];
	tail_a = tail_sum[1];
	head_b = head_prod[1];
	tail_b = tail_prod[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[1] = head_t;
	tail_sum[1] = tail_t;
      }				/* sum = sum+prod */
      y_ii[0] = tail_y_i[iy];
      y_ii[1] = tail_y_i[iy + 1];
      {
	double head_e1, tail_e1;
	double d1;
	double d2;
	/* Real part */
	d1 = (double) x_ii[0] * y_ii[0];
	d2 = (double) -x_ii[1] * y_ii[1];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_prod[0] = head_e1;
	tail_prod[0] = tail_e1;
	/* imaginary part */
	d1 = (double) x_ii[0] * y_ii[1];
	d2 = (double) x_ii[1] * y_ii[0];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_prod[1] = head_e1;
	tail_prod[1] = tail_e1;
      }				/* prod = x[i] * tail_y[i] */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_sum[0];
	tail_a = tail_sum[0];
	head_b = head_prod[0];
	tail_b = tail_prod[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[0] = head_t;
	tail_sum[0] = tail_t;
	/* Imaginary part */
	head_a = head_sum[1];
	tail_a = tail_sum[1];
	head_b = head_prod[1];
	tail_b = tail_prod[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[1] = head_t;
	tail_sum[1] = tail_t;
      }				/* sum = sum+prod */
      ix += incx;
      iy += incy;
    }				/* endfor */
  }

  {
    double cd[2];
    cd[0] = (double) alpha_i[0];
    cd[1] = (double) alpha_i[1];
    {
      /* Compute complex-extra = complex-extra * complex-double. */
      double head_a0, tail_a0;
      double head_a1, tail_a1;
      double head_t1, tail_t1;
      double head_t2, tail_t2;
      head_a0 = head_sum[0];
      tail_a0 = tail_sum[0];
      head_a1 = head_sum[1];
      tail_a1 = tail_sum[1];
      /* real part */
      {
	/* Compute double-double = double-double * double. */
	double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	con = head_a0 * split;
	a11 = con - head_a0;
	a11 = con - a11;
	a21 = head_a0 - a11;
	con = cd[0] * split;
	b1 = con - cd[0];
	b1 = con - b1;
	b2 = cd[0] - b1;

	c11 = head_a0 * cd[0];
	c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	c2 = tail_a0 * cd[0];
	t1 = c11 + c2;
	t2 = (c2 - (t1 - c11)) + c21;

	head_t1 = t1 + t2;
	tail_t1 = t2 - (head_t1 - t1);
      }
      {
	/* Compute double-double = double-double * double. */
	double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	con = head_a1 * split;
	a11 = con - head_a1;
	a11 = con - a11;
	a21 = head_a1 - a11;
	con = cd[1] * split;
	b1 = con - cd[1];
	b1 = con - b1;
	b2 = cd[1] - b1;

	c11 = head_a1 * cd[1];
	c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	c2 = tail_a1 * cd[1];
	t1 = c11 + c2;
	t2 = (c2 - (t1 - c11)) + c21;

	head_t2 = t1 + t2;
	tail_t2 = t2 - (head_t2 - t1);
      }
      head_t2 = -head_t2;
      tail_t2 = -tail_t2;
      {
	/* Compute double-double = double-double + double-double. */
	double bv;
	double s1, s2, t1, t2;

	/* Add two hi words. */
	s1 = head_t1 + head_t2;
	bv = s1 - head_t1;
	s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	/* Add two lo words. */
	t1 = tail_t1 + tail_t2;
	bv = t1 - tail_t1;
	t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	s2 += t1;

	/* Renormalize (s1, s2)  to  (t1, s2) */
	t1 = s1 + s2;
	s2 = s2 - (t1 - s1);

	t2 += s2;

	/* Renormalize (t1, t2)  */
	head_t1 = t1 + t2;
	tail_t1 = t2 - (head_t1 - t1);
      }
      head_tmp1[0] = head_t1;
      tail_tmp1[0] = tail_t1;
      /* imaginary part */
      {
	/* Compute double-double = double-double * double. */
	double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	con = head_a1 * split;
	a11 = con - head_a1;
	a11 = con - a11;
	a21 = head_a1 - a11;
	con = cd[0] * split;
	b1 = con - cd[0];
	b1 = con - b1;
	b2 = cd[0] - b1;

	c11 = head_a1 * cd[0];
	c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	c2 = tail_a1 * cd[0];
	t1 = c11 + c2;
	t2 = (c2 - (t1 - c11)) + c21;

	head_t1 = t1 + t2;
	tail_t1 = t2 - (head_t1 - t1);
      }
      {
	/* Compute double-double = double-double * double. */
	double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	con = head_a0 * split;
	a11 = con - head_a0;
	a11 = con - a11;
	a21 = head_a0 - a11;
	con = cd[1] * split;
	b1 = con - cd[1];
	b1 = con - b1;
	b2 = cd[1] - b1;

	c11 = head_a0 * cd[1];
	c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	c2 = tail_a0 * cd[1];
	t1 = c11 + c2;
	t2 = (c2 - (t1 - c11)) + c21;

	head_t2 = t1 + t2;
	tail_t2 = t2 - (head_t2 - t1);
      }
      {
	/* Compute double-double = double-double + double-double. */
	double bv;
	double s1, s2, t1, t2;

	/* Add two hi words. */
	s1 = head_t1 + head_t2;
	bv = s1 - head_t1;
	s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	/* Add two lo words. */
	t1 = tail_t1 + tail_t2;
	bv = t1 - tail_t1;
	t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	s2 += t1;

	/* Renormalize (s1, s2)  to  (t1, s2) */
	t1 = s1 + s2;
	s2 = s2 - (t1 - s1);

	t2 += s2;

	/* Renormalize (t1, t2)  */
	head_t1 = t1 + t2;
	tail_t1 = t2 - (head_t1 - t1);
      }
      head_tmp1[1] = head_t1;
      tail_tmp1[1] = tail_t1;
    }

  }				/* tmp1 = sum*alpha */
  {
    double head_e1, tail_e1;
    double d1;
    double d2;
    /* Real part */
    d1 = (double) r_v[0] * beta_i[0];
    d2 = (double) -r_v[1] * beta_i[1];
    {
      /* Compute double-double = double + double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = d1 + d2;
      e = t1 - d1;
      t2 = ((d2 - e) + (d1 - (t1 - e)));

      /* The result is t1 + t2, after normalization. */
      head_e1 = t1 + t2;
      tail_e1 = t2 - (head_e1 - t1);
    }
    head_tmp2[0] = head_e1;
    tail_tmp2[0] = tail_e1;
    /* imaginary part */
    d1 = (double) r_v[0] * beta_i[1];
    d2 = (double) r_v[1] * beta_i[0];
    {
      /* Compute double-double = double + double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = d1 + d2;
      e = t1 - d1;
      t2 = ((d2 - e) + (d1 - (t1 - e)));

      /* The result is t1 + t2, after normalization. */
      head_e1 = t1 + t2;
      tail_e1 = t2 - (head_e1 - t1);
    }
    head_tmp2[1] = head_e1;
    tail_tmp2[1] = tail_e1;
  }				/* tmp2 = r*beta */
  {
    double head_t, tail_t;
    double head_a, tail_a;
    double head_b, tail_b;
    /* Real part */
    head_a = head_tmp1[0];
    tail_a = tail_tmp1[0];
    head_b = head_tmp2[0];
    tail_b = tail_tmp2[0];
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_a + head_b;
      bv = s1 - head_a;
      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_a + tail_b;
      bv = t1 - tail_a;
      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_t = t1 + t2;
      tail_t = t2 - (head_t - t1);
    }
    head_tmp1[0] = head_t;
    tail_tmp1[0] = tail_t;
    /* Imaginary part */
    head_a = head_tmp1[1];
    tail_a = tail_tmp1[1];
    head_b = head_tmp2[1];
    tail_b = tail_tmp2[1];
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_a + head_b;
      bv = s1 - head_a;
      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_a + tail_b;
      bv = t1 - tail_a;
      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_t = t1 + t2;
      tail_t = t2 - (head_t - t1);
    }
    head_tmp1[1] = head_t;
    tail_tmp1[1] = tail_t;
  }				/* tmp1 = tmp1+tmp2 */
  head_r_true[0] = head_tmp1[0];
  tail_r_true[0] = tail_tmp1[0];
  head_r_true[1] = head_tmp1[1];
  tail_r_true[1] = tail_tmp1[1];	/* *r = tmp1 */

  FPU_FIX_STOP;

}
void z_r_truth2(enum blas_conj_type conj, int n, const void *alpha,
		const void *x, int incx, const void *beta,
		const void *head_y, const void *tail_y, int incy, void *r,
		double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *   r_true <- beta * r + alpha * SUM_{i=0, n-1} x[i] * (head_y[i] + tail_y[i])
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
 * head_y
 * tail_y (input) const void*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input) void*
 *
 * head_r_true
 * tail_r_true (output) void*
 *        The truth.
 * 
 */
{

  int i, ix = 0, iy = 0;
  double *r_i = (double *) r;
  const double *x_i = (double *) x;
  const double *head_y_i = (double *) head_y;
  const double *tail_y_i = (double *) tail_y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double x_ii[2];
  double y_ii[2];
  double r_v[2];
  double head_prod[2], tail_prod[2];
  double head_sum[2], tail_sum[2];
  double head_tmp1[2], tail_tmp1[2];
  double head_tmp2[2], tail_tmp2[2];
  FPU_FIX_DECL;

  /* Immediate return */
  if (n < 0) {
    head_r_true[0] = tail_r_true[0] = head_r_true[1] = tail_r_true[1] = 0.0;
    return;
  }

  FPU_FIX_START;

  r_v[0] = r_i[0];
  r_v[1] = r_i[0 + 1];
  head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;
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
      y_ii[0] = head_y_i[iy];
      y_ii[1] = head_y_i[iy + 1];
      x_ii[1] = -x_ii[1];
      {
	/* Compute complex-extra = complex-double * complex-double. */
	double head_t1, tail_t1;
	double head_t2, tail_t2;
	/* Real part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[0] * split;
	  a1 = con - x_ii[0];
	  a1 = con - a1;
	  a2 = x_ii[0] - a1;
	  con = y_ii[0] * split;
	  b1 = con - y_ii[0];
	  b1 = con - b1;
	  b2 = y_ii[0] - b1;

	  head_t1 = x_ii[0] * y_ii[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[1] * split;
	  a1 = con - x_ii[1];
	  a1 = con - a1;
	  a2 = x_ii[1] - a1;
	  con = y_ii[1] * split;
	  b1 = con - y_ii[1];
	  b1 = con - b1;
	  b2 = y_ii[1] - b1;

	  head_t2 = x_ii[1] * y_ii[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	head_t2 = -head_t2;
	tail_t2 = -tail_t2;
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_prod[0] = head_t1;
	tail_prod[0] = tail_t1;
	/* Imaginary part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[1] * split;
	  a1 = con - x_ii[1];
	  a1 = con - a1;
	  a2 = x_ii[1] - a1;
	  con = y_ii[0] * split;
	  b1 = con - y_ii[0];
	  b1 = con - b1;
	  b2 = y_ii[0] - b1;

	  head_t1 = x_ii[1] * y_ii[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[0] * split;
	  a1 = con - x_ii[0];
	  a1 = con - a1;
	  a2 = x_ii[0] - a1;
	  con = y_ii[1] * split;
	  b1 = con - y_ii[1];
	  b1 = con - b1;
	  b2 = y_ii[1] - b1;

	  head_t2 = x_ii[0] * y_ii[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_prod[1] = head_t1;
	tail_prod[1] = tail_t1;
      }				/* prod = x[i] * head_y[i] */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_sum[0];
	tail_a = tail_sum[0];
	head_b = head_prod[0];
	tail_b = tail_prod[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[0] = head_t;
	tail_sum[0] = tail_t;
	/* Imaginary part */
	head_a = head_sum[1];
	tail_a = tail_sum[1];
	head_b = head_prod[1];
	tail_b = tail_prod[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[1] = head_t;
	tail_sum[1] = tail_t;
      }				/* sum = sum+prod */
      y_ii[0] = tail_y_i[iy];
      y_ii[1] = tail_y_i[iy + 1];
      {
	/* Compute complex-extra = complex-double * complex-double. */
	double head_t1, tail_t1;
	double head_t2, tail_t2;
	/* Real part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[0] * split;
	  a1 = con - x_ii[0];
	  a1 = con - a1;
	  a2 = x_ii[0] - a1;
	  con = y_ii[0] * split;
	  b1 = con - y_ii[0];
	  b1 = con - b1;
	  b2 = y_ii[0] - b1;

	  head_t1 = x_ii[0] * y_ii[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[1] * split;
	  a1 = con - x_ii[1];
	  a1 = con - a1;
	  a2 = x_ii[1] - a1;
	  con = y_ii[1] * split;
	  b1 = con - y_ii[1];
	  b1 = con - b1;
	  b2 = y_ii[1] - b1;

	  head_t2 = x_ii[1] * y_ii[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	head_t2 = -head_t2;
	tail_t2 = -tail_t2;
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_prod[0] = head_t1;
	tail_prod[0] = tail_t1;
	/* Imaginary part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[1] * split;
	  a1 = con - x_ii[1];
	  a1 = con - a1;
	  a2 = x_ii[1] - a1;
	  con = y_ii[0] * split;
	  b1 = con - y_ii[0];
	  b1 = con - b1;
	  b2 = y_ii[0] - b1;

	  head_t1 = x_ii[1] * y_ii[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[0] * split;
	  a1 = con - x_ii[0];
	  a1 = con - a1;
	  a2 = x_ii[0] - a1;
	  con = y_ii[1] * split;
	  b1 = con - y_ii[1];
	  b1 = con - b1;
	  b2 = y_ii[1] - b1;

	  head_t2 = x_ii[0] * y_ii[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_prod[1] = head_t1;
	tail_prod[1] = tail_t1;
      }				/* prod = x[i] * tail_y[i] */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_sum[0];
	tail_a = tail_sum[0];
	head_b = head_prod[0];
	tail_b = tail_prod[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[0] = head_t;
	tail_sum[0] = tail_t;
	/* Imaginary part */
	head_a = head_sum[1];
	tail_a = tail_sum[1];
	head_b = head_prod[1];
	tail_b = tail_prod[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[1] = head_t;
	tail_sum[1] = tail_t;
      }				/* sum = sum+prod */
      ix += incx;
      iy += incy;
    }				/* endfor */
  } else {
    /* do not conjugate */

    for (i = 0; i < n; ++i) {
      x_ii[0] = x_i[ix];
      x_ii[1] = x_i[ix + 1];
      y_ii[0] = head_y_i[iy];
      y_ii[1] = head_y_i[iy + 1];

      {
	/* Compute complex-extra = complex-double * complex-double. */
	double head_t1, tail_t1;
	double head_t2, tail_t2;
	/* Real part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[0] * split;
	  a1 = con - x_ii[0];
	  a1 = con - a1;
	  a2 = x_ii[0] - a1;
	  con = y_ii[0] * split;
	  b1 = con - y_ii[0];
	  b1 = con - b1;
	  b2 = y_ii[0] - b1;

	  head_t1 = x_ii[0] * y_ii[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[1] * split;
	  a1 = con - x_ii[1];
	  a1 = con - a1;
	  a2 = x_ii[1] - a1;
	  con = y_ii[1] * split;
	  b1 = con - y_ii[1];
	  b1 = con - b1;
	  b2 = y_ii[1] - b1;

	  head_t2 = x_ii[1] * y_ii[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	head_t2 = -head_t2;
	tail_t2 = -tail_t2;
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_prod[0] = head_t1;
	tail_prod[0] = tail_t1;
	/* Imaginary part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[1] * split;
	  a1 = con - x_ii[1];
	  a1 = con - a1;
	  a2 = x_ii[1] - a1;
	  con = y_ii[0] * split;
	  b1 = con - y_ii[0];
	  b1 = con - b1;
	  b2 = y_ii[0] - b1;

	  head_t1 = x_ii[1] * y_ii[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[0] * split;
	  a1 = con - x_ii[0];
	  a1 = con - a1;
	  a2 = x_ii[0] - a1;
	  con = y_ii[1] * split;
	  b1 = con - y_ii[1];
	  b1 = con - b1;
	  b2 = y_ii[1] - b1;

	  head_t2 = x_ii[0] * y_ii[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_prod[1] = head_t1;
	tail_prod[1] = tail_t1;
      }				/* prod = x[i] * head_y[i] */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_sum[0];
	tail_a = tail_sum[0];
	head_b = head_prod[0];
	tail_b = tail_prod[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[0] = head_t;
	tail_sum[0] = tail_t;
	/* Imaginary part */
	head_a = head_sum[1];
	tail_a = tail_sum[1];
	head_b = head_prod[1];
	tail_b = tail_prod[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[1] = head_t;
	tail_sum[1] = tail_t;
      }				/* sum = sum+prod */
      y_ii[0] = tail_y_i[iy];
      y_ii[1] = tail_y_i[iy + 1];
      {
	/* Compute complex-extra = complex-double * complex-double. */
	double head_t1, tail_t1;
	double head_t2, tail_t2;
	/* Real part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[0] * split;
	  a1 = con - x_ii[0];
	  a1 = con - a1;
	  a2 = x_ii[0] - a1;
	  con = y_ii[0] * split;
	  b1 = con - y_ii[0];
	  b1 = con - b1;
	  b2 = y_ii[0] - b1;

	  head_t1 = x_ii[0] * y_ii[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[1] * split;
	  a1 = con - x_ii[1];
	  a1 = con - a1;
	  a2 = x_ii[1] - a1;
	  con = y_ii[1] * split;
	  b1 = con - y_ii[1];
	  b1 = con - b1;
	  b2 = y_ii[1] - b1;

	  head_t2 = x_ii[1] * y_ii[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	head_t2 = -head_t2;
	tail_t2 = -tail_t2;
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_prod[0] = head_t1;
	tail_prod[0] = tail_t1;
	/* Imaginary part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[1] * split;
	  a1 = con - x_ii[1];
	  a1 = con - a1;
	  a2 = x_ii[1] - a1;
	  con = y_ii[0] * split;
	  b1 = con - y_ii[0];
	  b1 = con - b1;
	  b2 = y_ii[0] - b1;

	  head_t1 = x_ii[1] * y_ii[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii[0] * split;
	  a1 = con - x_ii[0];
	  a1 = con - a1;
	  a2 = x_ii[0] - a1;
	  con = y_ii[1] * split;
	  b1 = con - y_ii[1];
	  b1 = con - b1;
	  b2 = y_ii[1] - b1;

	  head_t2 = x_ii[0] * y_ii[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_prod[1] = head_t1;
	tail_prod[1] = tail_t1;
      }				/* prod = x[i] * tail_y[i] */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_sum[0];
	tail_a = tail_sum[0];
	head_b = head_prod[0];
	tail_b = tail_prod[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[0] = head_t;
	tail_sum[0] = tail_t;
	/* Imaginary part */
	head_a = head_sum[1];
	tail_a = tail_sum[1];
	head_b = head_prod[1];
	tail_b = tail_prod[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_sum[1] = head_t;
	tail_sum[1] = tail_t;
      }				/* sum = sum+prod */
      ix += incx;
      iy += incy;
    }				/* endfor */
  }

  {
    /* Compute complex-extra = complex-extra * complex-double. */
    double head_a0, tail_a0;
    double head_a1, tail_a1;
    double head_t1, tail_t1;
    double head_t2, tail_t2;
    head_a0 = head_sum[0];
    tail_a0 = tail_sum[0];
    head_a1 = head_sum[1];
    tail_a1 = tail_sum[1];
    /* real part */
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

      con = head_a0 * split;
      a11 = con - head_a0;
      a11 = con - a11;
      a21 = head_a0 - a11;
      con = alpha_i[0] * split;
      b1 = con - alpha_i[0];
      b1 = con - b1;
      b2 = alpha_i[0] - b1;

      c11 = head_a0 * alpha_i[0];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = tail_a0 * alpha_i[0];
      t1 = c11 + c2;
      t2 = (c2 - (t1 - c11)) + c21;

      head_t1 = t1 + t2;
      tail_t1 = t2 - (head_t1 - t1);
    }
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

      con = head_a1 * split;
      a11 = con - head_a1;
      a11 = con - a11;
      a21 = head_a1 - a11;
      con = alpha_i[1] * split;
      b1 = con - alpha_i[1];
      b1 = con - b1;
      b2 = alpha_i[1] - b1;

      c11 = head_a1 * alpha_i[1];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = tail_a1 * alpha_i[1];
      t1 = c11 + c2;
      t2 = (c2 - (t1 - c11)) + c21;

      head_t2 = t1 + t2;
      tail_t2 = t2 - (head_t2 - t1);
    }
    head_t2 = -head_t2;
    tail_t2 = -tail_t2;
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_t1 + head_t2;
      bv = s1 - head_t1;
      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_t1 + tail_t2;
      bv = t1 - tail_t1;
      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_t1 = t1 + t2;
      tail_t1 = t2 - (head_t1 - t1);
    }
    head_tmp1[0] = head_t1;
    tail_tmp1[0] = tail_t1;
    /* imaginary part */
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

      con = head_a1 * split;
      a11 = con - head_a1;
      a11 = con - a11;
      a21 = head_a1 - a11;
      con = alpha_i[0] * split;
      b1 = con - alpha_i[0];
      b1 = con - b1;
      b2 = alpha_i[0] - b1;

      c11 = head_a1 * alpha_i[0];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = tail_a1 * alpha_i[0];
      t1 = c11 + c2;
      t2 = (c2 - (t1 - c11)) + c21;

      head_t1 = t1 + t2;
      tail_t1 = t2 - (head_t1 - t1);
    }
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

      con = head_a0 * split;
      a11 = con - head_a0;
      a11 = con - a11;
      a21 = head_a0 - a11;
      con = alpha_i[1] * split;
      b1 = con - alpha_i[1];
      b1 = con - b1;
      b2 = alpha_i[1] - b1;

      c11 = head_a0 * alpha_i[1];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = tail_a0 * alpha_i[1];
      t1 = c11 + c2;
      t2 = (c2 - (t1 - c11)) + c21;

      head_t2 = t1 + t2;
      tail_t2 = t2 - (head_t2 - t1);
    }
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_t1 + head_t2;
      bv = s1 - head_t1;
      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_t1 + tail_t2;
      bv = t1 - tail_t1;
      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_t1 = t1 + t2;
      tail_t1 = t2 - (head_t1 - t1);
    }
    head_tmp1[1] = head_t1;
    tail_tmp1[1] = tail_t1;
  }
  /* tmp1 = sum*alpha */
  {
    /* Compute complex-extra = complex-double * complex-double. */
    double head_t1, tail_t1;
    double head_t2, tail_t2;
    /* Real part */
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = r_v[0] * split;
      a1 = con - r_v[0];
      a1 = con - a1;
      a2 = r_v[0] - a1;
      con = beta_i[0] * split;
      b1 = con - beta_i[0];
      b1 = con - b1;
      b2 = beta_i[0] - b1;

      head_t1 = r_v[0] * beta_i[0];
      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = r_v[1] * split;
      a1 = con - r_v[1];
      a1 = con - a1;
      a2 = r_v[1] - a1;
      con = beta_i[1] * split;
      b1 = con - beta_i[1];
      b1 = con - b1;
      b2 = beta_i[1] - b1;

      head_t2 = r_v[1] * beta_i[1];
      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    head_t2 = -head_t2;
    tail_t2 = -tail_t2;
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_t1 + head_t2;
      bv = s1 - head_t1;
      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_t1 + tail_t2;
      bv = t1 - tail_t1;
      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_t1 = t1 + t2;
      tail_t1 = t2 - (head_t1 - t1);
    }
    head_tmp2[0] = head_t1;
    tail_tmp2[0] = tail_t1;
    /* Imaginary part */
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = r_v[1] * split;
      a1 = con - r_v[1];
      a1 = con - a1;
      a2 = r_v[1] - a1;
      con = beta_i[0] * split;
      b1 = con - beta_i[0];
      b1 = con - b1;
      b2 = beta_i[0] - b1;

      head_t1 = r_v[1] * beta_i[0];
      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = r_v[0] * split;
      a1 = con - r_v[0];
      a1 = con - a1;
      a2 = r_v[0] - a1;
      con = beta_i[1] * split;
      b1 = con - beta_i[1];
      b1 = con - b1;
      b2 = beta_i[1] - b1;

      head_t2 = r_v[0] * beta_i[1];
      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_t1 + head_t2;
      bv = s1 - head_t1;
      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_t1 + tail_t2;
      bv = t1 - tail_t1;
      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_t1 = t1 + t2;
      tail_t1 = t2 - (head_t1 - t1);
    }
    head_tmp2[1] = head_t1;
    tail_tmp2[1] = tail_t1;
  }				/* tmp2 = r*beta */
  {
    double head_t, tail_t;
    double head_a, tail_a;
    double head_b, tail_b;
    /* Real part */
    head_a = head_tmp1[0];
    tail_a = tail_tmp1[0];
    head_b = head_tmp2[0];
    tail_b = tail_tmp2[0];
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_a + head_b;
      bv = s1 - head_a;
      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_a + tail_b;
      bv = t1 - tail_a;
      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_t = t1 + t2;
      tail_t = t2 - (head_t - t1);
    }
    head_tmp1[0] = head_t;
    tail_tmp1[0] = tail_t;
    /* Imaginary part */
    head_a = head_tmp1[1];
    tail_a = tail_tmp1[1];
    head_b = head_tmp2[1];
    tail_b = tail_tmp2[1];
    {
      /* Compute double-double = double-double + double-double. */
      double bv;
      double s1, s2, t1, t2;

      /* Add two hi words. */
      s1 = head_a + head_b;
      bv = s1 - head_a;
      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

      /* Add two lo words. */
      t1 = tail_a + tail_b;
      bv = t1 - tail_a;
      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

      s2 += t1;

      /* Renormalize (s1, s2)  to  (t1, s2) */
      t1 = s1 + s2;
      s2 = s2 - (t1 - s1);

      t2 += s2;

      /* Renormalize (t1, t2)  */
      head_t = t1 + t2;
      tail_t = t2 - (head_t - t1);
    }
    head_tmp1[1] = head_t;
    tail_tmp1[1] = tail_t;
  }				/* tmp1 = tmp1+tmp2 */
  head_r_true[0] = head_tmp1[0];
  tail_r_true[0] = tail_tmp1[0];
  head_r_true[1] = head_tmp1[1];
  tail_r_true[1] = tail_tmp1[1];	/* *r = tmp1 */

  FPU_FIX_STOP;

}
