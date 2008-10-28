#include <stdlib.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

extern double xrand(int *);
extern int FixedBits(double, double);

void r_truth(enum blas_conj_type conj, int n, float alpha, const double *x_l, const double *x_t, float beta, const float *y, float *r,  /* input */
             double *r_true_l, double *r_true_t)
{
  int i, ix = 0, iy = 0;
  float *r_i = r;
  float alpha_i = alpha;
  float beta_i = beta;
  double dt = (double) alpha_i;
  double x_ii_l, x_ii_t;
  float y_ii;
  float r_v;
  double prod_l, prod_t;
  double sum_l, sum_t;
  double tmp1_l, tmp1_t;
  double tmp2_l, tmp2_t;

  /*printf("r_truth0> r=%e\n", *r); */

  /* Immediate return */
  if (n < 0) {
    *r_true_l = *r_true_t = 0.0;
    return;
  }
  r_v = r_i[0];
  sum_l = sum_t = 0.0;          /* sum = 0 */


  for (i = 0; i < n; ++i) {
    x_ii_l = x_l[ix];
    x_ii_t = x_t[ix];
    y_ii = (double) y[iy];
    /* prod = x[i]*y[i] */
    ddmuld(x_ii_l, x_ii_t, y_ii, &prod_l, &prod_t);
    /* sum = sum+prod */
    ddadd(sum_l, sum_t, prod_l, prod_t, &sum_l, &sum_t);
    ++ix;
    ++iy;
  }                             /* endfor */
  dt = (double) alpha_i;
  ddmuld(sum_l, sum_t, dt, &tmp1_l, &tmp1_t);   /* tmp1 = sum*alpha */

  tmp2_l = (double) r_v *beta_i;
  tmp2_t = 0.0;                 /* tmp2 = r*beta */
  /* tmp1 = tmp1+tmp2 */
  ddadd(tmp1_l, tmp1_t, tmp2_l, tmp2_t, &tmp1_l, &tmp1_t);

  /*printf("r_truth> r=%e\n", *r); */

  /* Return r_truth */
  *r_true_l = tmp1_l;
  *r_true_t = tmp1_t;

}                               /* end r_truth */

void
gen_y_to_cancel(int k, int n, enum blas_conj_type conj,
                float alpha, double *x_l, double *x_t, float *y)
/*
 * Purpose
 * =======
 *
 * Generate Y(i)'s from k to n-1 to cancel as much as possible.
 *
 */
{
  int i;
  float rtmp = 0.0;
  double x_i_l, x_i_t;
  double r_true_l, r_true_t;
  double tmp_l, tmp_t;

  for (i = k; i < n; ++i) {
    r_truth(conj, i, alpha, x_l, x_t, 0.0, y, &rtmp, &r_true_l, &r_true_t);
    x_i_l = x_l[i];
    x_i_t = x_t[i];
    ddmuld(x_i_l, x_i_t, alpha, &tmp_l, &tmp_t);
    if (tmp_l == 0. && tmp_t == 0.)
      y[i] = 0.;
    else {
      dddiv(r_true_l, r_true_t, tmp_l, tmp_t, &x_i_l, &x_i_t);
      y[i] = (float) -x_i_l;
    }
  }
}                               /* end gen_y_to_cancel */


void
gen_r_to_cancel(int n, enum blas_conj_type conj,
                float alpha, float beta,
                double *x_l, double *x_t, float *y, float *r, int *seed)
/*
 * Purpose
 * =======
 *
 * Generate r to cancel as much as possible.
 *
 */
{
  float rtmp;
  double r_true_l, r_true_t;

  if (beta == 0.0) {
    *r = xrand(seed);
  } else {
    rtmp = 0.0;
    r_truth(conj, n, alpha, x_l, x_t, 0.0, y, &rtmp, &r_true_l, &r_true_t);
    *r = -r_true_l / beta;
  }
}                               /* end gen_r_to_cancel */

void
testgen_BLAS_sdot_x(int n, int n_fix2, int n_mix, int norm,
                    enum blas_conj_type conj,
                    float *alpha, int alpha_flag, float *beta, int beta_flag,
                    double *x_l, double *x_t, float *y, int *seed,
                    float *r, double *r_true_l, double *r_true_t)
/*
 * Purpose
 * =======
 *
 * This routine generates the test vectors X and Y for DOT.
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
 * alpha   (input/output) float*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) float*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x_l     (input) double*
 *         The leading part of x vector.
 *
 * x_t     (input) double*
 *         The trailing part of x vector.
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 *
 * r       (output) float*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int B, frees, y_free, i, k, s;
  float rtmps;
  double rtmpd, f;
  double tmp_l, tmp_t;

  if (n_fix2 + n_mix < n)
    BLAS_error("testgen_BLAS_sdot_x", 0, 0, NULL);

  if (alpha_flag == 0)
    *alpha = xrand(seed);
  if (beta_flag == 0)
    *beta = xrand(seed);

  y_free = n - n_fix2;
  k = n_fix2;

  /*
   * Compute the number of bits in the prefix sum: alpha *
   * SUM_{i=0,n_fix2-1}(x[i] * y[i])
   */
  *r = 0.0;
  r_truth(conj, n_fix2, *alpha, x_l, x_t, 0.0, y, r, r_true_l, r_true_t);
  B = FixedBits(*r_true_l, *r_true_t);

  /* Pick r at random */
  *r = xrand(seed);

  /*printf("testgen_s_dot_x> r=%e\n", *r); */

  if (alpha_flag == 1 && *alpha == 0.0) {
    /* Pick the free Y(i)'s at random. */
    for (i = n_fix2; i < n; ++i)
      y[i] = xrand(seed);
    /* Compute r_truth in double-double */
    r_truth(conj, n, *alpha, x_l, x_t, *beta, y, r, r_true_l, r_true_t);
    return;
  }
  if (beta_flag == 1 && *beta == 0.0) {
    if (B == 0) {               /* Assume alpha is not very big . */
      switch (y_free) {
      case 0:
        break;
      case 1:
        y[k] = xrand(seed);
        break;
      case 2:
        /* Both x[k] and x[k+1] fixed; cancel 24 bits. */
        y[k] = xrand(seed);
        gen_y_to_cancel(k + 1, n, conj, *alpha, x_l, x_t, y);
        break;
      case 3:
        /*
         * Make SUM_{i=0,2}(x[k+i] * y[k+i]) small ... x[k]*y[k] =
         * -x[k+2]*y[k+2] exact, x[k+1]*y[k+1] small
         */
        y[k] = -x_l[k + 2];
        y[k + 2] = x_l[k];
        tmp_l = x_t[k] * y[k] + x_t[k + 2] * y[k + 2];
        if (tmp_l == 0.0) {
          ddmuld(x_l[k], x_t[k], y[k], &f, &tmp_t);
          s = 100;              /* Should be random between 40 and 100. */
          tmp_l = f * power(2, -s);
          tmp_t = 0.0;
          if (x_l[k + 1] == 0. && x_t[k + 1] == 0.)
            y[k + 1] = 0.;
          else {
            dddiv(tmp_l, tmp_t, x_l[k + 1], x_t[k + 1], &tmp_l, &tmp_t);
            y[k + 1] = (float) -tmp_l;
          }
        } else {
          y[k] = xrand(seed);   /* add bits */
          gen_y_to_cancel(k + 1, n, conj, *alpha, x_l, x_t, y);
        }
        break;
      case 4:
        /*
         * Make SUM_{i=0,3}(x[k+i] * y[k+i]) small ... x[k]*y[k] =
         * -x[k+3]*y[k+3] exact, x[k+1]*y[k+1] small, x[k+2]*y[k+2] =
         * 0.
         */
        y[k] = -x_l[k + 3];
        y[k + 3] = x_l[k];
        tmp_l = x_t[k] * y[k] + x_t[k + 3] * y[k + 3];
        if (tmp_l == 0.0) {
          ddmuld(x_l[k], x_t[k], y[k], &f, &tmp_t);
          s = 100;              /* Should be random between 40 and 100. */
          tmp_l = f * power(2, -s);
          tmp_t = 0.0;
          if (x_l[k + 1] == 0. && x_t[k + 1] == 0.)
            y[k + 1] = 0.;
          else {
            dddiv(tmp_l, tmp_t, x_l[k + 1], x_t[k + 1], &tmp_l, &tmp_t);
            y[k + 1] = (float) -tmp_l;
          }
          y[k + 2] = 0.0;
        } else {
          y[k] = xrand(seed);   /* add bits */
          gen_y_to_cancel(k + 1, n, conj, *alpha, x_l, x_t, y);
        }
        break;
      default:                 /* y_free >= 5 */
        /*
         * Make alpha * SUM_{i=0,n-1}(x[k+i] * y[k+i]) small ...
         * Cancel >= 72 bits.
         */
        if (y_free <= 7) {
          /* Use first 2 to add bits, rest to cancal bits. */
          y[k] = xrand(seed);
          f = *alpha * x_l[k] * y[k];
          s = 100;              /* Should be random between 40 and 100. */
          if (*alpha * x_l[k + 1] == 0.)
            y[k + 1] = 0.;
          else
            y[k + 1] = f * power(2, -s) / (*alpha * x_l[k + 1]);
          gen_y_to_cancel(k + 2, n, conj, *alpha, x_l, x_t, y);
        } else {
          /*
           * Use last 5 to cancel bits, and leading ones to add
           * bits.
           */
          y[k] = xrand(seed);
          rtmpd = *alpha * x_l[k] * y[k];
          s = 30;
          for (i = k + 1; i < n - 5; ++i) {
            rtmpd *= power(2, -s);
            if (*alpha * x_l[i] == 0.)
              y[i] = 0.;
            else
              y[i] = rtmpd / (*alpha * x_l[i]);
          }
          gen_y_to_cancel(n - 5, n, conj, *alpha, x_l, x_t, y);
        }
        break;
      }                         /* end switch */
    } else {                    /* B > 0 */
      if (B >= BITS_E) {        /* Choose Y(i)'s to cancel. */
        gen_y_to_cancel(k, n, conj, *alpha, x_l, x_t, y);
      } else {                  /* At least y[n-1] is free. */
        if (y_free == 1) {
          /* Cancel min(B,24) bits. */
          gen_y_to_cancel(k, n, conj, *alpha, x_l, x_t, y);
        } else {                /* >= 2 frees. */
          /*
           * There are 2 possibilities: (1) use 1 to add bits, and
           * y_free-1 to cancel 24*(y_free-1) bits (2) use all to
           * cancel min(B, 24*y_free) bits Goal is to maximize the
           * # of bits cancelled. By equating (1) and (2), we find
           * the crossover point is y_free = B/24 + 1.
           */
          if (y_free > B / 24.0 + 1) {  /* Use scheme (1) */
            rtmps = 0.0;
            r_truth(conj, k, *alpha, x_l, x_t, 0.0, y,
                    &rtmps, r_true_l, r_true_t);
            s = 100;            /* Should be random between 40 and 100. */
            if (*alpha * x_l[k] == 0.0)
              y[k] = 0.;
            else
              y[k] = *r_true_l * power(2, -s) / (*alpha * x_l[k]);
            gen_y_to_cancel(k + 1, n, conj, *alpha, x_l, x_t, y);
          } else {              /* Use scheme (2) */
            gen_y_to_cancel(k, n, conj, *alpha, x_l, x_t, y);
          }
        }
      }                         /* end else B < 106 */
    }                           /* end else B > 0 */

    /* Compute r_truth in double-double */
    r_truth(conj, n, *alpha, x_l, x_t, *beta, y, r, r_true_l, r_true_t);
    return;
  }
  /* if beta == 0 */
  /* Now, beta is non-zero. */
  if (B == 0) {
    switch (y_free) {
    case 0:
      break;
    case 1:
      /* Make alpha*x[k]*y[k] + beta*r small. */
      /* Count number of frees in alpha and beta. */
      y[k] = xrand(seed);
      gen_r_to_cancel(n, conj, *alpha, *beta, x_l, x_t, y, r, seed);
      break;
    case 2:
      /* Make alpha * SUM_{i=0,1}(x[k+i] * y[k+i]) + beta*r small. */
      /* Count number of frees in alpha and beta. */
      frees = 0;
      if (alpha_flag == 0)
        ++frees;
      if (beta_flag == 0)
        ++frees;
      if (frees > 0 && x_t[k] == 0.0) {
        /*
         * Make alpha*x[k]*y[k] = -beta*r exact, alpha*x[k+1]*y[k+1]
         * small.
         */
        y[k] = -1.0;
        if (alpha_flag == 0)
          *alpha = *beta;
        else if (beta_flag == 0)
          *beta = *alpha;
        *r = x_l[k];
        /*printf("testgen_s_dot_x2> r=%e\n", *r); */
        f = *alpha * x_l[k] * y[k];
        s = 100;                /* Should be random between 40 and 100. */
        if (*alpha * x_l[k + 1] == 0.)
          y[k + 1] = 0.;
        else
          y[k + 1] = f * power(2, -s) / (*alpha * x_l[k + 1]);
      } else {                  /* Cancel 48 bits. */
        y[k] = xrand(seed);
        gen_y_to_cancel(k + 1, n, conj, *alpha, x_l, x_t, y);
        gen_r_to_cancel(n, conj, *alpha, *beta, x_l, x_t, y, r, seed);
      }
      break;
    case 3:
      /*
       * Make alpha * SUM_{i=0,2}(x[k+i] * y[k+i]) + beta*r small ...
       * x[k]*y[k] = -x[k+2]*y[k+2] exact, x[k+1]*y[k+1] small
       */
      y[k] = -x_l[k + 2];
      y[k + 2] = x_l[k];
      tmp_l = x_t[k] * y[k] + x_t[k + 2] * y[k + 2];
      if (tmp_l == 0.0) {
        ddmuld(x_l[k], x_t[k], y[k], &f, &tmp_t);
        s = 100;                /* Should be random between 40 and 100. */
        tmp_l = f * power(2, -s);
        tmp_t = 0.0;
        if (x_l[k + 1] == 0. && x_t[k + 1] == 0.)
          y[k + 1] = 0.;
        else {
          dddiv(tmp_l, tmp_t, x_l[k + 1], x_t[k + 1], &tmp_l, &tmp_t);
          y[k + 1] = (float) -tmp_l;
        }
        *r = 0.0;
        /*printf("testgen_s_dot_x3> r=%e\n", *r); */
      } else {
        y[k] = xrand(seed);     /* add bits */
        gen_y_to_cancel(k + 1, n, conj, *alpha, x_l, x_t, y);
        gen_r_to_cancel(n, conj, *alpha, *beta, x_l, x_t, y, r, seed);
      }
      break;
    default:                   /* Actual frees >= 5 */
      /*
       * Make SUM_{i=0,n-1}(alpha * x[k+i] * y[k+i]) + beta*r small.
       * ... Cancel >= 72 bits.
       */
      if (y_free <= 6) {
        /* Use 2 to add bits, rest to cancel bits. */
        y[k] = xrand(seed);
        f = *alpha * x_l[k] * y[k];
        s = 100;                /* Should be random between 40 and 100. */
        if (*alpha * x_l[k + 1] == 0.)
          y[k + 1] = 0.;
        else
          y[k + 1] = f * power(2, -s) / (*alpha * x_l[k + 1]);
        gen_y_to_cancel(k + 2, n, conj, *alpha, x_l, x_t, y);
        gen_r_to_cancel(n, conj, *alpha, *beta, x_l, x_t, y, r, seed);
      } else {
        /*
         * Use last 5 (4 Y(i)'s and r) to cancel bits, and leading
         * ones to add bits.
         */
        y[k] = xrand(seed);
        rtmpd = *alpha * x_l[k] * y[k];
        s = 30;
        for (i = k + 1; i < n - 4; ++i) {
          rtmpd *= power(2, -s);
          if (*alpha * x_l[i] == 0.)
            y[i] = 0.;
          else
            y[i] = rtmpd / (*alpha * x_l[i]);
        }
        gen_y_to_cancel(n - 4, n, conj, *alpha, x_l, x_t, y);
        gen_r_to_cancel(n, conj, *alpha, *beta, x_l, x_t, y, r, seed);
      }
      break;
    }                           /* end switch */
  } else {                      /* B > 0 */
    if (B >= BITS_E) {          /* Choose Y(i)'s and r to cancel. */
      gen_y_to_cancel(k, n, conj, *alpha, x_l, x_t, y);
      gen_r_to_cancel(n, conj, *alpha, *beta, x_l, x_t, y, r, seed);
    } else {                    /* >= 2 frees. Use y[k] to add bits. */
      frees = y_free + 1;
      /*
       * There are 2 possibilities: (1) use 1 to add bits, and frees-1
       * to cancel 24*(frees-1) bits (2) use all to cancel min(B,
       * 24*frees) bits Goal is to maximize the # of bits cancelled. By
       * equating (1) and (2), we find the crossover point is frees =
       * B/24 + 1.
       */
      if (frees > B / 24.0 + 1) {       /* Use scheme (1) */
        rtmps = 0.0;
        r_truth(conj, k, *alpha, x_l, x_t, 0.0, y,
                &rtmps, r_true_l, r_true_t);
        s = 100;                /* Should be random between 40 and 100. */
        if (*alpha * x_l[k] == 0.)
          y[k] = 0.;
        else
          y[k] = *r_true_l * power(2, -s) / (*alpha * x_l[k]);
        gen_y_to_cancel(k + 1, n, conj, *alpha, x_l, x_t, y);
      } else {                  /* Use scheme (2) */
        gen_y_to_cancel(k, n, conj, *alpha, x_l, x_t, y);
      }
      gen_r_to_cancel(n, conj, *alpha, *beta, x_l, x_t, y, r, seed);
    }
  }

  /* Compute r_truth in double-double */
  r_truth(conj, n, *alpha, x_l, x_t, *beta, y, r, r_true_l, r_true_t);

}                               /* end testgen_BLAS_sdot_x */
