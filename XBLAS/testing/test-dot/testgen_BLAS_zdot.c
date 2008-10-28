#include <stdlib.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"


/* Complex-Complex Multiplication */
void z_mul(double a[], double b[], double c[])
{
  double cr, ci;
  cr = a[0] * b[0] - a[1] * b[1];
  ci = a[1] * b[0] + a[0] * b[1];
  c[0] = cr;
  c[1] = ci;
}

/* Complex Division c = a/b */
void z_div(double a[], double b[], double c[])
{
  double ratio, den;
  double abr, abi, cr, ci;

  if ((abr = b[0]) < 0.)
    abr = -abr;
  if ((abi = b[1]) < 0.)
    abi = -abi;
  if (abr <= abi) {
    if (abi == 0) {
      BLAS_error("z_div: division by zero", 0, 0, NULL);
    }
    ratio = b[0] / b[1];
    den = b[1] * (1 + ratio * ratio);
    cr = (a[0] * ratio + a[1]) / den;
    ci = (a[1] * ratio - a[0]) / den;
  } else {
    ratio = b[1] / b[0];
    den = b[0] * (1 + ratio * ratio);
    cr = (a[0] + a[1] * ratio) / den;
    ci = (a[1] - a[0] * ratio) / den;
  }
  c[0] = cr;
  c[1] = ci;
}

static double ulp(double a)
/*
 * Purpose 
 * =======
 * 
 * Compute the unit last place of a double precision number.
 */
{
  double f;
  int e;
  f = frexp(a, &e);
  return power(2, e - BITS_D);
}


static double rand_half_1(int l_bits, int *seed)
/*
 * Purpose
 * =======
 * 
 * Generate random number in the interval [0.5, 1). l_bits specifies that only
 * the leading l_bits are nonzero.
 * 
 */
{
  double a = xrand(seed);       /* [0,1] */
  a /= 2.;
  a += 0.5;
  if (l_bits < BITS_D) {
    double s = power(2, l_bits);
    double t = a / s;           /* shift right l_bits */
    t = (t + a) - a;            /* cancel trailing bits */
    a = t * s;                  /* shift back */
  }
  return a;
}


static void r_truth(enum blas_conj_type conj, int n, void *alpha, const void *x, int incx, void *beta, const void *y, int incy, void *r,        /* input */
                    double *r_true_l, double *r_true_t)
{
  int i, ix = 0, iy = 0;
  double *r_i = (double *) r;
  const double *x_i = (double *) x;
  const double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double x_ii[2];
  double y_ii[2];
  double r_v[2];
  double prod_l[2], prod_t[2];
  double sum_l[2], sum_t[2];
  double tmp1_l[2], tmp1_t[2];
  double tmp2_l[2], tmp2_t[2];

  /* Immediate return */
  if (n < 0) {
    r_true_l[0] = r_true_l[1] = r_true_t[0] = r_true_t[1] = 0.0;
    return;
  }

  r_v[0] = r_i[0];
  r_v[1] = r_i[0 + 1];
  sum_l[0] = sum_l[1] = sum_t[0] = sum_t[1] = 0.0;      /* sum = 0 */
  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    if (conj == blas_conj)
      x_ii[1] = -x_ii[1];
    y_ii[0] = y_i[iy];
    y_ii[1] = y_i[iy + 1];
    {
      /*
       * Compute complex-extra = complex-double * complex-double.
       */
      double t1_l, t1_t;
      double t2_l, t2_t;
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

        t1_l = x_ii[0] * y_ii[0];
        t1_t = (((a1 * b1 - t1_l) + a1 * b2) + a2 * b1) + a2 * b2;
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

        t2_l = x_ii[1] * y_ii[1];
        t2_t = (((a1 * b1 - t2_l) + a1 * b2) + a2 * b1) + a2 * b2;
      }
      t2_l = -t2_l;
      t2_t = -t2_t;
      {
        /*
         * Compute double-double = double-double + double-double.
         */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = t1_l + t2_l;
        e = t1 - t1_l;
        t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

        /* The result is t1 + t2, after normalization. */
        t1_l = t1 + t2;
        t1_t = t2 - (t1_l - t1);
      }
      prod_l[0] = t1_l;
      prod_t[0] = t1_t;
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

        t1_l = x_ii[1] * y_ii[0];
        t1_t = (((a1 * b1 - t1_l) + a1 * b2) + a2 * b1) + a2 * b2;
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

        t2_l = x_ii[0] * y_ii[1];
        t2_t = (((a1 * b1 - t2_l) + a1 * b2) + a2 * b1) + a2 * b2;
      }
      {
        /*
         * Compute double-double = double-double + double-double.
         */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = t1_l + t2_l;
        e = t1 - t1_l;
        t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

        /* The result is t1 + t2, after normalization. */
        t1_l = t1 + t2;
        t1_t = t2 - (t1_l - t1);
      }
      prod_l[1] = t1_l;
      prod_t[1] = t1_t;
    }                           /* prod = x[i]*y[i] */
    {
      double t_l, t_t;
      double a_l, a_t;
      double b_l, b_t;
      /* Real part */
      a_l = sum_l[0];
      a_t = sum_t[0];
      b_l = prod_l[0];
      b_t = prod_t[0];
      {
        /*
         * Compute double-double = double-double + double-double.
         */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = a_l + b_l;
        e = t1 - a_l;
        t2 = ((b_l - e) + (a_l - (t1 - e))) + a_t + b_t;

        /* The result is t1 + t2, after normalization. */
        t_l = t1 + t2;
        t_t = t2 - (t_l - t1);
      }
      sum_l[0] = t_l;
      sum_t[0] = t_t;
      /* Imaginary part */
      a_l = sum_l[1];
      a_t = sum_t[1];
      b_l = prod_l[1];
      b_t = prod_t[1];
      {
        /*
         * Compute double-double = double-double + double-double.
         */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = a_l + b_l;
        e = t1 - a_l;
        t2 = ((b_l - e) + (a_l - (t1 - e))) + a_t + b_t;

        /* The result is t1 + t2, after normalization. */
        t_l = t1 + t2;
        t_t = t2 - (t_l - t1);
      }
      sum_l[1] = t_l;
      sum_t[1] = t_t;
    }                           /* sum = sum+prod */
    ix += 2;
    iy += 2;
  }                             /* endfor */
  {
    /* Compute complex-extra = complex-extra * complex-double. */
    double a0_l, a0_t;
    double a1_l, a1_t;
    double t1_l, t1_t;
    double t2_l, t2_t;
    a0_l = sum_l[0];
    a0_t = sum_t[0];
    a1_l = sum_l[1];
    a1_t = sum_t[1];
    /* real part */
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

      con = a0_l * split;
      a11 = con - a0_l;
      a11 = con - a11;
      a21 = a0_l - a11;
      con = alpha_i[0] * split;
      b1 = con - alpha_i[0];
      b1 = con - b1;
      b2 = alpha_i[0] - b1;

      c11 = a0_l * alpha_i[0];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = a0_t * alpha_i[0];
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

      con = a1_l * split;
      a11 = con - a1_l;
      a11 = con - a11;
      a21 = a1_l - a11;
      con = alpha_i[1] * split;
      b1 = con - alpha_i[1];
      b1 = con - b1;
      b2 = alpha_i[1] - b1;

      c11 = a1_l * alpha_i[1];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = a1_t * alpha_i[1];
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      t2_l = t1 + t2;
      t2_t = t2 - (t2_l - t1);
    }
    t2_l = -t2_l;
    t2_t = -t2_t;
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = t1_l + t2_l;
      e = t1 - t1_l;
      t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

      /* The result is t1 + t2, after normalization. */
      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    tmp1_l[0] = t1_l;
    tmp1_t[0] = t1_t;
    /* imaginary part */
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

      con = a1_l * split;
      a11 = con - a1_l;
      a11 = con - a11;
      a21 = a1_l - a11;
      con = alpha_i[0] * split;
      b1 = con - alpha_i[0];
      b1 = con - b1;
      b2 = alpha_i[0] - b1;

      c11 = a1_l * alpha_i[0];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = a1_t * alpha_i[0];
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

      con = a0_l * split;
      a11 = con - a0_l;
      a11 = con - a11;
      a21 = a0_l - a11;
      con = alpha_i[1] * split;
      b1 = con - alpha_i[1];
      b1 = con - b1;
      b2 = alpha_i[1] - b1;

      c11 = a0_l * alpha_i[1];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = a0_t * alpha_i[1];
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      t2_l = t1 + t2;
      t2_t = t2 - (t2_l - t1);
    }
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = t1_l + t2_l;
      e = t1 - t1_l;
      t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

      /* The result is t1 + t2, after normalization. */
      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    tmp1_l[1] = t1_l;
    tmp1_t[1] = t1_t;
  }
  /* tmp1 = sum*alpha */
  {
    /* Compute complex-extra = complex-double * complex-double. */
    double t1_l, t1_t;
    double t2_l, t2_t;
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

      t1_l = r_v[0] * beta_i[0];
      t1_t = (((a1 * b1 - t1_l) + a1 * b2) + a2 * b1) + a2 * b2;
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

      t2_l = r_v[1] * beta_i[1];
      t2_t = (((a1 * b1 - t2_l) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    t2_l = -t2_l;
    t2_t = -t2_t;
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = t1_l + t2_l;
      e = t1 - t1_l;
      t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

      /* The result is t1 + t2, after normalization. */
      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    tmp2_l[0] = t1_l;
    tmp2_t[0] = t1_t;
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

      t1_l = r_v[1] * beta_i[0];
      t1_t = (((a1 * b1 - t1_l) + a1 * b2) + a2 * b1) + a2 * b2;
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

      t2_l = r_v[0] * beta_i[1];
      t2_t = (((a1 * b1 - t2_l) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = t1_l + t2_l;
      e = t1 - t1_l;
      t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

      /* The result is t1 + t2, after normalization. */
      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    tmp2_l[1] = t1_l;
    tmp2_t[1] = t1_t;
  }                             /* tmp2 = r*beta */
  {
    double t_l, t_t;
    double a_l, a_t;
    double b_l, b_t;
    /* Real part */
    a_l = tmp1_l[0];
    a_t = tmp1_t[0];
    b_l = tmp2_l[0];
    b_t = tmp2_t[0];
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = a_l + b_l;
      e = t1 - a_l;
      t2 = ((b_l - e) + (a_l - (t1 - e))) + a_t + b_t;

      /* The result is t1 + t2, after normalization. */
      t_l = t1 + t2;
      t_t = t2 - (t_l - t1);
    }
    tmp1_l[0] = t_l;
    tmp1_t[0] = t_t;
    /* Imaginary part */
    a_l = tmp1_l[1];
    a_t = tmp1_t[1];
    b_l = tmp2_l[1];
    b_t = tmp2_t[1];
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = a_l + b_l;
      e = t1 - a_l;
      t2 = ((b_l - e) + (a_l - (t1 - e))) + a_t + b_t;

      /* The result is t1 + t2, after normalization. */
      t_l = t1 + t2;
      t_t = t2 - (t_l - t1);
    }
    tmp1_l[1] = t_l;
    tmp1_t[1] = t_t;
  }                             /* tmp1 = tmp1+tmp2 */

  /* Return r_truth = tmp1 */
  r_true_l[0] = tmp1_l[0];
  r_true_l[1] = tmp1_l[1];
  r_true_t[0] = tmp1_t[0];
  r_true_t[1] = tmp1_t[1];
}                               /* end r_truth */


static void
gen_y_to_cancel(int k, int n, enum blas_conj_type conj,
                void *alpha, void *x, void *y)
/*
 * Purpose
 * =======
 * 
 * Generate Y(i)'s from k to n-1 to cancel as much as possible.
 * 
 */
{
  int i, ii;
  double zero[2] = { 0.0, 0.0 };
  double r[2] = { 0.0, 0.0 };
  double tmp[2], tmp_l[2], tmp_t[2];
  double *x_i = x, *y_i = y;
  double r_true_l[2], r_true_t[2];

  for (i = k; i < n; ++i) {
    /* y[i] = -rtmp / (alpha * x[i]); */
    r_truth(conj, i, alpha, x, 1, zero, y, 1, r, r_true_l, r_true_t);
    ii = 2 * i;
    if (conj == blas_conj) {
      tmp[0] = x_i[ii];
      tmp[1] = -x_i[ii + 1];
      z_mul((double *) alpha, tmp, tmp);
    } else {
      z_mul((double *) alpha, &x_i[ii], tmp);
    }
    if (tmp[0] == 0. && tmp[1] == 0.)
      y_i[ii] = y_i[ii + 1] = 0.;
    else {
      z_dddivd(r_true_l, r_true_t, tmp, tmp_l, tmp_t);
      y_i[ii] = -tmp_l[0];
      y_i[ii + 1] = -tmp_l[1];
    }
  }
}


static void
gen_r_to_cancel(int n, enum blas_conj_type conj,
                void *alpha, void *beta, void *x, void *y, void *r, int *seed)
/*
 * Purpose
 * =======
 * 
 * Generate r to cancel as much as possible.
 * 
 */
{
  double zero[2] = { 0.0, 0.0 };
  double rtmp[2] = { 0.0, 0.0 };
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  double r_true_l[2], r_true_t[2];

  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
    r_i[0] = xrand(seed);
    r_i[1] = xrand(seed);
  } else {
    r_truth(conj, n, alpha, x, 1, zero, y, 1, rtmp, r_true_l, r_true_t);
    z_dddivd(r_true_l, r_true_t, beta_i, r_true_l, r_true_t);
    r_i[0] = -r_true_l[0];
    r_i[1] = -r_true_l[1];
  }
}



void
testgen_BLAS_zdot(int n, int n_fix2, int n_mix, int norm,
                  enum blas_conj_type conj,
                  void *alpha, int alpha_flag, void *beta, int beta_flag,
                  void *x, void *y, int *seed,
                  void *r, double r_true_l[], double r_true_t[])
/*
 * Purpose
 * =======
 *
 * This routine generates the test vectors X and Y for C_ZDOT.
 *
 * Arguments
 * =========
 * 
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int* 
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double[]
 *         The leading (real,imaginary) parts of the truth in double-double.
 *
 * r_true_t (output) double[]
 *         The trailing (real,imaginary) parts of the truth in double-double.
 *
 */
{
  int B, frees, y_free, i, ii, k, s;
  double zero[2] = { 0.0, 0.0 };
  double a, b, eps, eps_out;
  double f[2], rtmp[2];
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  double *x_i = (double *) x, *y_i = (double *) y;

  if (alpha_flag == 0) {
    alpha_i[0] = xrand(seed);
    alpha_i[1] = xrand(seed);
  }
  if (beta_flag == 0) {
    beta_i[0] = xrand(seed);
    beta_i[1] = xrand(seed);
  }

  y_free = n - n_fix2;
  k = 2 * n_fix2;
  eps_out = power(2, -BITS_D);

  /*
   * Compute the number of bits in the prefix sum:
   *     alpha * SUM_{i=0,n_fix2-1}(x[i] * y[i])
   */
  r_i[0] = r_i[1] = 0.0;
  r_truth(conj, n_fix2, alpha, x, 1, zero, y, 1, r, r_true_l, r_true_t);
  B = FixedBits(r_true_l[0], r_true_t[0]);      /* real */
  B = MAX(B, FixedBits(r_true_l[1], r_true_t[1]));      /* imag */

  /* Pick r at random */
  r_i[0] = xrand(seed);
  r_i[1] = xrand(seed);

  /* Pick the free X(i)'s at random. */
  for (i = n_fix2 + n_mix; i < n; ++i) {
    ii = 2 * i;
    x_i[ii] = xrand(seed);
    x_i[ii + 1] = xrand(seed);
  }

  if (alpha_flag == 1 && alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
    /* Pick the free Y(i)'s at random. */
    for (i = n_fix2; i < n; ++i) {
      ii = 2 * i;
      y_i[ii] = xrand(seed);
      y_i[ii + 1] = xrand(seed);
    }
    /* Compute r_truth in double-double */
    r_truth(conj, n, alpha, x, 1, beta, y, 1, r, r_true_l, r_true_t);
    return;
  }

  if (beta_flag == 1 && beta_i[0] == 0.0 && beta_i[1] == 0.0) {
    if (B == 0) {               /* Assume alpha is not very big . */
      switch (y_free) {
      case 0:
        break;
      case 1:
        y_i[k] = xrand(seed);
        y_i[k + 1] = xrand(seed);
        break;
      case 2:
        /*
         * Make SUM_{i=0,1}(x[k+i] * y[k+i]) small ... alpha * eps^2
         */
        if (n_mix == 0) {       /* Both x[k] and x[k+1] free. */
          /* x[k]*y[k] + x[k+1]*y[k+1] = eps_out^2, small.
           * For complex, each real number is multiplied by (i+1),
           * the result is 2i * eps_out^2.
           */
          a = rand_half_1(26, seed);    /* strictly < 1 */
          x_i[k] = a;           /* real */
          x_i[k + 1] = a;       /* imag */
          y_i[k] = a;
          y_i[k + 1] = a;
          x_i[k + 2] = a + eps_out;     /* exact */
          x_i[k + 3] = a + eps_out;     /* exact */
          y_i[k + 2] = -a + eps_out;    /* exact */
          y_i[k + 3] = -a + eps_out;    /* exact */
        } else if (n_mix == 1) {        /* x[k] fixed, x[k+1] free. */
          /* x[k]*y[k] + x[k+1]*y[k+1] = (eps1 + eps2*i)^2 */
          a = x_i[k];           /* real */
          b = x_i[k + 1];       /* imag */
          if (conj == blas_conj)
            b = -b;
          y_i[k] = a;
          y_i[k + 1] = b;
          eps = ulp(a);
          x_i[k + 2] = a + eps; /* exact */
          y_i[k + 2] = -a + eps;        /* exact */
          eps = ulp(b);
          x_i[k + 3] = b + eps; /* exact */
          if (conj == blas_conj)
            x_i[k + 3] = -x_i[k + 3];
          y_i[k + 3] = -b + eps;        /* exact */
        } else {                /* Both x[k] and x[k+1] fixed; cancel 53 bits. */
          y_i[k] = xrand(seed);
          y_i[k + 1] = xrand(seed);
          gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, x, y);
        }
        break;
      case 3:
        /*
         * Make SUM_{i=0,2}(x[k+i] * y[k+i]) small
         * ... x[k]*y[k] = -x[k+2]*y[k+2] exact, x[k+1]*y[k+1] small
         */
        y_i[k] = -x_i[k + 4];
        y_i[k + 1] = -x_i[k + 5];
        y_i[k + 4] = x_i[k];
        y_i[k + 5] = x_i[k + 1];
        rtmp[0] = x_i[k];
        if (conj == blas_conj) {
          rtmp[1] = -x_i[k + 1];
          y_i[k + 1] = -y_i[k + 1];
          y_i[k + 5] = -y_i[k + 5];
        } else {
          rtmp[1] = x_i[k + 1];
        }
        z_mul(rtmp, &y_i[k], f);
        s = 100;                /* Should be random between 40 and 100. */
        b = power(2, -s);
        f[0] *= b;
        f[1] *= b;
        rtmp[0] = x_i[k + 2];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 3];
        else
          rtmp[1] = x_i[k + 3];
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[k + 2] = y_i[k + 3] = 0.;
        else
          z_div(f, rtmp, &y_i[k + 2]);
        break;
      case 4:
        /*
         * Make SUM_{i=0,3}(x[k+i] * y[k+i]) small
         * ... x[k]*y[k] = -x[k+3]*y[k+3] exact, x[k+1]*y[k+1] small,
         *     x[k+2]*y[k+2] = 0.
         */
        y_i[k] = -x_i[k + 6];
        y_i[k + 1] = -x_i[k + 7];
        y_i[k + 6] = x_i[k];
        y_i[k + 7] = x_i[k + 1];
        rtmp[0] = x_i[k];
        if (conj == blas_conj) {
          rtmp[1] = -x_i[k + 1];
          y_i[k + 1] = -y_i[k + 1];
          y_i[k + 7] = -y_i[k + 7];
        } else {
          rtmp[1] = x_i[k + 1];
        }
        z_mul(rtmp, &y_i[k], f);
        s = 100;                /* Should be random between 40 and 100. */
        b = power(2, -s);
        f[0] *= b;
        f[1] *= b;
        rtmp[0] = x_i[k + 2];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 3];
        else
          rtmp[1] = x_i[k + 3];
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[k + 2] = y_i[k + 3] = 0.;
        else
          z_div(f, rtmp, &y_i[k + 2]);
        y_i[k + 4] = 0.0;
        y_i[k + 5] = 0.0;
        break;
      default:                 /* y_free >= 5 */
        /*
         * Make SUM_{i=0,n-1}(x[k+i] * y[k+i]) small.
         * Use last 2 to cancel bits, and leading ones to add bits.
         * ... Cancel >= 106 bits.
         */
        y_i[k] = xrand(seed);
        y_i[k + 1] = xrand(seed);
        rtmp[0] = x_i[k];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 1];
        else
          rtmp[1] = x_i[k + 1];
        z_mul(rtmp, &y_i[k], rtmp);
        z_mul(alpha_i, rtmp, rtmp);
        s = 30;
        b = power(2, -s);
        for (i = n_fix2 + 1; i < n - 2; ++i) {
          rtmp[0] *= b;
          rtmp[1] *= b;
          ii = 2 * i;
          if (conj == blas_conj) {
            f[0] = x_i[ii];
            f[1] = -x_i[ii + 1];
          } else {
            f[0] = x_i[ii];
            f[1] = x_i[ii + 1];
          }
          z_mul(alpha_i, f, f);
          if (f[0] == 0. && f[1] == 0.)
            y_i[ii] = y_i[ii + 1] = 0.;
          else
            z_div(rtmp, f, &y_i[ii]);
        }
        gen_y_to_cancel(n - 2, n, conj, alpha, x, y);
      }                         /* end switch */
    } else {                    /* B > 0 */
      if (B >= BITS_E) {        /* Choose Y(i)'s to cancel. */
        gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
      } else {                  /* At least y[n-1] is free. */
        if (y_free == 1) {
          /* Cancel min(B,53) bits. */
          gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
        } else {                /* >= 2 frees. */
          /*
           * There are 2 possibilities:
           * (1) use 1 to add bits, and y_free-1 to cancel 53*(y_free-1) bits
           * (2) use all to cancel min(B, 53*y_free) bits
           * Goal is to maximize the # of bits cancelled. By equating (1)
           * and (2), we find the crossover point is y_free = B/53 + 1.
           */
          if (y_free > B / (double) BITS_D + 1) {       /* Use scheme (1) */
            f[0] = f[1] = 0.0;
            r_truth(conj, n_fix2, alpha, x, 1, zero, y, 1, f,
                    r_true_l, r_true_t);
            f[0] = r_true_l[0];
            f[1] = r_true_l[1];
            s = 100;            /* Should be random between 40 and 100. */
            b = power(2, -s);
            f[0] *= b;
            f[1] *= b;
            rtmp[0] = x_i[k];
            if (conj == blas_conj)
              rtmp[1] = -x_i[k + 1];
            else
              rtmp[1] = x_i[k + 1];
            z_mul(alpha_i, rtmp, rtmp);
            if (rtmp[0] == 0. && rtmp[1] == 0.)
              y_i[k] = y_i[k + 1] = 0.;
            else
              z_div(f, rtmp, &y_i[k]);
            gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, x, y);
          } else {              /* Use scheme (2) */
            gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
          }
        }
      }                         /* end else B < 106 */
    }                           /* end else B > 0 */

    /* Compute r_truth in double-double */
    r_truth(conj, n, alpha, x, 1, zero, y, 1, r, r_true_l, r_true_t);
    return;
  }

  /* Now, beta is non-zero. */
  if (B == 0) {
    switch (y_free) {
    case 0:
      break;
    case 1:
      /* Make alpha*x[k]*y[k] + beta*r small. */
      /* Count number of frees in alpha, x[k], and beta. */
      frees = 0;
      if (alpha_flag == 0)
        ++frees;
      if (beta_flag == 0)
        ++frees;
      if (n_mix == 0)
        ++frees;
      if (frees >= 2) {
        /* alpha*x[k]*y[k] + beta*r = -alpha * eps_out^2
         * For complex, each real number is multiplied by (i+1) to
         * yield final result -2i * alpha * eps_out^2.
         */
        a = rand_half_1(26, seed);      /* strictly < 1, only leading 26 bits */
        r_i[0] = 0.0;           /* real */
        r_i[1] = -a * a * 2;    /* imag, exact */
        if (beta_flag == 1) {   /* beta fixed, alpha must be free */
          alpha_i[0] = beta_i[0];
          alpha_i[1] = beta_i[1];
          x_i[k] = a + eps_out; /* real, exact */
          x_i[k + 1] = a + eps_out;     /* imag, exact */
          if (conj == blas_conj)
            x_i[k + 1] = -x_i[k + 1];
          y_i[k] = a - eps_out; /* exact */
          y_i[k + 1] = a - eps_out;     /* exact */
        } else if (n_mix == 1) {        /* x[k] fixed */
          beta_i[0] = x_i[k];
          beta_i[1] = x_i[k + 1];
          if (conj == blas_conj)
            beta_i[1] = -beta_i[1];
          alpha_i[0] = a + eps_out;     /* real, exact */
          alpha_i[1] = a + eps_out;     /* imag, exact */
          y_i[k] = a - eps_out; /* exact */
          y_i[k + 1] = a - eps_out;     /* exact */
        } else {                /* alpha fixed or free, x[k] and beta free */
          beta_i[0] = alpha_i[0];
          beta_i[1] = alpha_i[1];
          x_i[k] = a + eps_out; /* real, exact */
          x_i[k + 1] = a + eps_out;     /* imag, exact */
          if (conj == blas_conj)
            x_i[k + 1] = -x_i[k + 1];
          y_i[k] = a - eps_out; /* exact */
          y_i[k + 1] = a - eps_out;     /* exact */
        }
      } else {                  /* Cancel 53 bits. */
        y_i[k] = xrand(seed);
        y_i[k + 1] = xrand(seed);
        gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
      }
      break;
    case 2:                    /* Actual frees = 3 */
      /* Make SUM_{i=0,1}(alpha * x[k+i] * y[k+i]) + beta*r small. */
      /* Count number of frees in alpha, x[k], and beta. */
      frees = 0;
      if (alpha_flag == 0)
        ++frees;
      if (beta_flag == 0)
        ++frees;
      if (n_mix == 0)
        ++frees;
      if (frees > 0) {
        /*
         * Make alpha*x[k]*y[k] = -beta*r exact, alpha*x[k+1]*y[k+1] small.
         */
        y_i[k] = -1.0;
        y_i[k + 1] = 0.0;
        if (alpha_flag == 0) {  /* alpha free */
          alpha_i[0] = beta_i[0];
          alpha_i[1] = beta_i[1];
          r_i[0] = x_i[k];
          r_i[1] = x_i[k + 1];
          if (conj == blas_conj)
            r_i[1] = -r_i[1];
        } else if (beta_flag == 0) {    /* beta free */
          beta_i[0] = alpha_i[0];
          beta_i[1] = alpha_i[1];
          r_i[0] = x_i[k];
          r_i[1] = x_i[k + 1];
          if (conj == blas_conj)
            r_i[1] = -r_i[1];
        } else {                /* x[k] free */
          x_i[k] = beta_i[0];
          x_i[k + 1] = beta_i[1];
          if (conj == blas_conj)
            x_i[k + 1] = -x_i[k + 1];
          r_i[0] = alpha_i[0];
          r_i[1] = alpha_i[1];
        }
        rtmp[0] = x_i[k];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 1];
        else
          rtmp[1] = x_i[k + 1];
        z_mul(rtmp, &y_i[k], f);
        z_mul(alpha_i, f, f);
        s = 100;                /* Should be random between 40 and 100. */
        b = power(2, -s);
        f[0] *= b;
        f[1] *= b;
        rtmp[0] = x_i[k + 2];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 3];
        else
          rtmp[1] = x_i[k + 3];
        z_mul(alpha_i, rtmp, rtmp);
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[k + 2] = y_i[k + 3] = 0.;
        else
          z_div(f, rtmp, &y_i[k + 2]);
      } else {                  /* Cancel 53 bits. */
        y_i[k] = xrand(seed);
        y_i[k + 1] = xrand(seed);
        gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, x, y);
        gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
      }
      break;
    default:                   /* Actual frees >= 4 */
      /*
       * Make SUM_{i=0,n-1}(alpha * x[k+i] * y[k+i]) + beta*r small.
       * Use last 2 ( Y(n) and r) to cancel bits, and leading ones to
       * add bits ... Cancel >= 106 bits.
       */
      y_i[k] = xrand(seed);
      y_i[k + 1] = xrand(seed);
      rtmp[0] = x_i[k];
      if (conj == blas_conj)
        rtmp[1] = -x_i[k + 1];
      else
        rtmp[1] = x_i[k + 1];
      z_mul(rtmp, &y_i[k], f);
      z_mul(alpha_i, f, f);
      s = 30;
      b = power(2, -s);
      for (i = n_fix2 + 1; i < n - 1; ++i) {
        f[0] *= b;
        f[1] *= b;
        ii = 2 * i;
        rtmp[0] = x_i[ii];
        if (conj == blas_conj)
          rtmp[1] = -x_i[ii + 1];
        else
          rtmp[1] = x_i[ii + 1];
        z_mul(alpha_i, rtmp, rtmp);
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[ii] = y_i[ii + 1] = 0.;
        else
          z_div(f, rtmp, &y_i[ii]);
      }
      gen_y_to_cancel(n - 1, n, conj, alpha, x, y);
      gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
      break;
    }                           /* end switch */
  } else {                      /* B > 0 */
    if (B >= BITS_E) {          /* Choose Y(i)'s and r to cancel. */
      gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
      gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
    } else {                    /* >= 2 frees. Use y[k] to add bits. */
      frees = y_free + 1;
      /*
       * There are 2 possibilities:
       * (1) use 1 to add bits, and y_free-1 to cancel 53*(y_free-1) bits
       * (2) use all to cancel min(B, 53*y_free) bits
       * Goal is to maximize the # of bits cancelled. By equating (1)
       * and (2), we find the crossover point is y_free = B/53 + 1.
       */
      if (frees > B / (double) BITS_D + 1) {    /* Use scheme (1) */
        f[0] = f[1] = 0.0;
        r_truth(conj, n_fix2, alpha, x, 1, zero, y, 1, f, r_true_l, r_true_t);
        f[0] = r_true_l[0];
        f[1] = r_true_l[1];
        s = 100;                /* Should be random between 40 and 100. */
        b = power(2, -s);
        f[0] *= b;
        f[1] *= b;
        rtmp[0] = x_i[k];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 1];
        else
          rtmp[1] = x_i[k + 1];
        z_mul(alpha_i, rtmp, rtmp);
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[k] = y_i[k + 1] = 0.;
        else
          z_div(f, rtmp, &y_i[k]);
        gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, x, y);
      } else {                  /* Use scheme (2) */
        gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
      }
      gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
    }
  }

  /* Compute r_truth in double-double */
  r_truth(conj, n, alpha, x, 1, beta, y, 1, r, r_true_l, r_true_t);
}                               /* testgen_BLAS_zdot */
