#include <stdlib.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"


static float rand_half_1(int, int *);
static double ulp(float);
static void r_truth(enum blas_conj_type, int, float, const float *, int,
                    float, const float *, int, float *, double *, double *);
static void gen_y_to_cancel(int, int, enum blas_conj_type, float,
                            float *, float *);
static void gen_r_to_cancel(int, enum blas_conj_type, float, float,
                            float *, float *, float *, int *);

void
testgen_BLAS_sdot(int n, int n_fix2, int n_mix, int norm,
                  enum blas_conj_type conj,
                  float *alpha, int alpha_flag, float *beta, int beta_flag,
                  float *x, float *y, int *seed,
                  float *r, double *r_true_l, double *r_true_t)
/*
 * Purpose
 * =======
 *
 * This routine generates the test vectors X and Y for C_SDOT.
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
 * x       (input/output) float*
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
  int B, frees, y_free, i, k, m, s;
  float a, rtmps;
  double rtmpd, eps, eps_out;
  double f;

  if (alpha_flag == 0)
    *alpha = xrand(seed);
  if (beta_flag == 0)
    *beta = xrand(seed);

  y_free = n - n_fix2;
  k = n_fix2;
  eps_out = power(2, -BITS_S);

  /* Compute the number of bits in the prefix sum:
   *     alpha * SUM_{i=0,n_fix2-1}(x[i] * y[i])
   */
  *r = 0.0;
  r_truth(conj, n_fix2, *alpha, x, 1, 0.0, y, 1, r, r_true_l, r_true_t);
  B = FixedBits(*r_true_l, *r_true_t);

  /* Pick r at random */
  *r = xrand(seed);

  /* Pick the free X(i)'s at random. */
  for (i = n_fix2 + n_mix; i < n; ++i)
    x[i] = xrand(seed);

  if (alpha_flag == 1 && *alpha == 0.0) {
    /* Pick the free Y(i)'s at random. */
    for (i = n_fix2; i < n; ++i)
      y[i] = xrand(seed);
    /* Compute r_truth in double-double */
    r_truth(conj, n, *alpha, x, 1, *beta, y, 1, r, r_true_l, r_true_t);
    return;
  }

  if (beta_flag == 1 && *beta == 0.0) {
    if (B == 0) {               /* Assume alpha is not very big . */
      switch (y_free) {
      case 0:
        break;
      case 1:
        y[n_fix2] = xrand(seed);
        break;
      case 2:
        /* Make SUM_{i=0,1}(x[k+i] * y[k+i]) small ...
           ... = alpha * eps^2                           */
        if (n_mix == 0) {       /* Both x[k] and x[k+1] free. */
          /* x[k]*y[k] + x[k+1]*y[k+1] = eps_out^2, small */
          a = rand_half_1(BITS_S, seed);        /* strictly < 1 */
          x[k] = a;
          y[k] = a;
          x[k + 1] = a + eps_out;       /* exact */
          y[k + 1] = -a + eps_out;      /* exact */
        } else if (n_mix == 1) {        /* x[k] fixed, x[k+1] free. */
          /* x[k]*y[k] + x[k+1]*y[k+1] = eps_out^2 */
          a = x[k];
          y[k] = a;
          eps = ulp(a);
          x[k + 1] = a + eps;   /* exact */
          y[k + 1] = -a + eps;  /* exact */
        } else {                /* Both x[k] and x[k+1] fixed; cancel 24 bits. */
          y[k] = xrand(seed);
          gen_y_to_cancel(k + 1, n, conj, *alpha, x, y);
        }
        break;
      case 3:
        /* Make SUM_{i=0,2}(x[k+i] * y[k+i]) small
           ... x[k]*y[k] = -x[k+2]*y[k+2] exact, x[k+1]*y[k+1] small */
        y[k] = -x[k + 2];
        y[k + 2] = x[k];
        f = x[k] * y[k];
        s = 100;                /* Should be random between 40 and 100. */
        if (x[k + 1] == 0.0)
          y[k + 1] = f;
        else
          y[k + 1] = f * power(2, -s) / x[k + 1];
        break;
      case 4:
        /* Make SUM_{i=0,3}(x[k+i] * y[k+i]) small
           ... x[k]*y[k] = -x[k+3]*y[k+3] exact, x[k+1]*y[k+1] small,
           x[k+2]*y[k+2] = 0. */
        y[k] = -x[k + 3];
        y[k + 3] = x[k];
        f = x[k] * y[k];
        s = 100;                /* Should be random between 40 and 100. */
        if (x[k + 1] == 0.0)
          y[k + 1] = f;
        else
          y[k + 1] = f * power(2, -s) / x[k + 1];
        y[k + 2] = 0.0;
        break;
      default:                 /* y_free >= 5 */
        /* Make SUM_{i=0,n-1}(x[k+i] * y[k+i]) small
           ... Cancel >= 72 bits. */
        a = xrand(seed);        /* Flip a coin */
        if (a < 0.5) {
          if (y_free <= 7) {
            /* Use first 2 to add bits, rest to cancal bits. */
            y[k] = xrand(seed);
            f = *alpha * x[k] * y[k];
            s = 100;            /* Should be random between 40 and 100. */
            if (*alpha * x[k + 1] == 0.0)
              y[k + 1] = 0.0;
            else
              y[k + 1] = f * power(2, -s) / (*alpha * x[k + 1]);
            gen_y_to_cancel(k + 2, n, conj, *alpha, x, y);
          } else {
            /* Use last 5 to cancel bits, and leading ones to add bits. */
            y[k] = xrand(seed);
            rtmpd = *alpha * x[k] * y[k];
            s = 30;
            for (i = k + 1; i < n - 5; ++i) {
              rtmpd *= power(2, -s);
              if (*alpha * x[i] == 0.0)
                y[i] = 0.0;
              else
                y[i] = rtmpd / (*alpha * x[i]);
            }
            gen_y_to_cancel(n - 5, n, conj, *alpha, x, y);
          }
        } else {                /* Use the scheme similar to case 3. */
          m = (y_free - 1) / 2;
          rtmpd = 0.0;
          for (i = 0; i < m; ++i) {
            y[k + i] = -x[n - i - 1];
            rtmpd = MAX(rtmpd, fabs(x[k + i] * y[k + i]));
          }
          for (i = 0; i < m; ++i)
            y[n - i - 1] = x[k + i];
          s = 100;              /* Should be random between 40 and 100. */
          if (x[k + m] == 0.0)
            y[k + m] = 0.0;
          else
            y[k + m] = rtmpd * power(2, -s) / x[k + m];
          if (y_free % 2 == 0)
            y[k + m + 1] = 0.0;
        }
        break;
      }                         /* end switch */
    } else {                    /* B > 0 */
      if (B >= BITS_E) {        /* Choose Y(i)'s to cancel. */
        gen_y_to_cancel(k, n, conj, *alpha, x, y);
      } else {                  /* At least y[n-1] is free. */
        if (y_free == 1) {
          /* Cancel min(B,24) bits. */
          gen_y_to_cancel(k, n, conj, *alpha, x, y);
        } else {                /* >= 2 frees. */
          /* There are 2 possibilities:
           * (1) use 1 to add bits, and y_free-1 to cancel 24*(y_free-1) bits
           * (2) use all to cancel min(B, 24*y_free) bits
           * Goal is to maximize the # of bits cancelled. By equating (1)
           * and (2), we find the crossover point is y_free = B/24 + 1.
           */
          if (y_free > B / 24.0 + 1) {  /* Use scheme (1) */
            rtmps = 0.0;
            BLAS_sdot_x(conj, k, *alpha, x, 1, 0.0, y, 1,
                        &rtmps, blas_prec_extra);
            rtmpd = rtmps;
            s = 100;            /* Should be random between 40 and 100. */
            y[k] = rtmpd * power(2, -s) / (*alpha * x[k]);
            gen_y_to_cancel(k + 1, n, conj, *alpha, x, y);
          } else {              /* Use scheme (2) */
            gen_y_to_cancel(k, n, conj, *alpha, x, y);
          }
        }
      }                         /* end else B < 106 */
    }                           /* end else B > 0 */

    /* Compute r_truth in double-double */
    r_truth(conj, n, *alpha, x, 1, *beta, y, 1, r, r_true_l, r_true_t);
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
      /* Count number of frees in alpha, x[k], and beta. */
      frees = 0;
      if (alpha_flag == 0)
        ++frees;
      if (beta_flag == 0)
        ++frees;
      if (n_mix == 0)
        ++frees;
      if (frees >= 2) {
        /* alpha*x[k]*y[k] + beta*r = -alpha * eps_out^2 */
        a = rand_half_1(12, seed);      /* strictly < 1, only leading 12 bits */
        if (alpha_flag == 1) {  /* alpha fixed */
          *beta = *alpha;
          x[k] = a + eps_out;   /* exact */
          y[k] = a - eps_out;   /* exact */
        } else if (beta_flag == 1) {    /* beta fixed */
          *alpha = *beta;
          x[k] = a + eps_out;   /* exact */
          y[k] = a - eps_out;   /* exact */
        } else if (n_mix == 1) {        /* x[k] fixed */
          *beta = x[k];
          *alpha = a + eps_out; /* exact */
          y[k] = a - eps_out;   /* exact */
        } else {                /* alpha, x[k], and beta all free */
          *beta = *alpha;
          x[k] = a + eps_out;   /* exact */
          y[k] = a - eps_out;   /* exact */
        }
        *r = -a * a;            /* exact */
      } else {                  /* Cancel 24 bits. */
        y[k] = xrand(seed);
        gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
      }
      break;
    case 2:
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
        /* Make alpha*x[k]*y[k] = -beta*r exact, alpha*x[k+1]*y[k+1] small. */
        y[k] = -1.0;
        if (alpha_flag == 0) {  /* alpha free */
          *alpha = *beta;
          *r = x[k];
        } else if (beta_flag == 0) {    /* beta free */
          *beta = *alpha;
          *r = x[k];
        } else {                /* x[k] free */
          x[k] = *beta;
          *r = *alpha;
        }
        f = *alpha * x[k] * y[k];
        s = 100;                /* Should be random between 40 and 100. */
        y[k + 1] = f * power(2, -s) / (*alpha * x[k + 1]);
      } else {                  /* Cancel 48 bits. */
        y[k] = xrand(seed);
        gen_y_to_cancel(k + 1, n, conj, *alpha, x, y);
        gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
      }
      break;
    case 3:
      /* Make SUM_{i=0,2}(alpha * x[k+i] * y[k+i]) + beta*r small
         ... x[k]*y[k] = -x[k+2]*y[k+2] exact, x[k+1]*y[k+1] small */
      y[k] = -x[k + 2];
      y[k + 2] = x[k];
      f = x[k] * y[k];
      s = 100;                  /* Should be random between 40 and 100. */
      if (x[k + 1] == 0.0)
        y[k + 1] = f;
      else
        y[k + 1] = f * power(2, -s) / x[k + 1];
      *r = 0.0;
      break;
    default:                   /* Actual frees >= 5 */
      /* Make SUM_{i=0,n-1}(alpha * x[k+i] * y[k+i]) + beta*r small.
         ... Cancel >= 72 bits. */
      a = xrand(seed);          /* Flip a coin */
      if (a < 0.5) {
        if (y_free <= 6) {
          /* Use 2 to add bits, rest to cancel bits. */
          y[k] = xrand(seed);
          f = *alpha * x[k] * y[k];
          s = 100;              /* Should be random between 40 and 100. */
          y[k + 1] = f * power(2, -s) / (*alpha * x[k + 1]);
          gen_y_to_cancel(k + 2, n, conj, *alpha, x, y);
          gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
        } else {
          /* Use last 5 (4 Y(i)'s and r) to cancel bits, and leading ones
             to add bits. */
          y[k] = xrand(seed);
          rtmpd = *alpha * x[k] * y[k];
          s = 30;
          for (i = k + 1; i < n - 4; ++i) {
            rtmpd *= power(2, -s);
            if (*alpha * x[i] == 0.0)
              y[i] = 0.0;
            else
              y[i] = rtmpd / (*alpha * x[i]);
          }
          gen_y_to_cancel(n - 4, n, conj, *alpha, x, y);
          gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
        }
      } else {                  /* Use the scheme similar to case 3. */
        /* Make the middle one tiny, and the rest cancel each other
           from both ends. */
        *r = 0.0;
        m = (y_free - 1) / 2;
        rtmpd = 0.0;
        for (i = 0; i < m; ++i) {
          y[k + i] = -x[n - i - 1];
          rtmpd = MAX(rtmpd, fabs(x[k + i] * y[k + i]));
        }
        for (i = 0; i < m; ++i)
          y[n - i - 1] = x[k + i];
        s = 100;                /* Should be random between 40 and 100. */
        if (x[k + m] == 0.0)
          y[k + m] = 0.;
        else
          y[k + m] = rtmpd * power(2, -s) / x[k + m];
        if (y_free % 2 == 0)
          y[k + m + 1] = 0.0;
      }
      break;
    }                           /* end switch */
  } else {                      /* B > 0 */
    if (B >= BITS_E) {          /* Choose Y(i)'s and r to cancel. */
      gen_y_to_cancel(k, n, conj, *alpha, x, y);
      gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
    } else {                    /* >= 2 frees. Use y[k] to add bits. */
      frees = y_free + 1;
      /* There are 2 possibilities:
       * (1) use 1 to add bits, and frees-1 to cancel 24*(frees-1) bits
       * (2) use all to cancel min(B, 24*frees) bits
       * Goal is to maximize the # of bits cancelled. By equating (1)
       * and (2), we find the crossover point is frees = B/24 + 1.
       */
      if (frees > B / 24.0 + 1) {       /* Use scheme (1) */
        rtmps = 0.0;
        BLAS_sdot_x(conj, k, *alpha, x, 1, 0.0, y, 1,
                    &rtmps, blas_prec_extra);
        s = 100;                /* Should be random between 40 and 100. */
        if (*alpha * x[k] == 0.0)
          y[k] = 0.;
        else
          y[k] = rtmps * power(2, -s) / (*alpha * x[k]);
        gen_y_to_cancel(k + 1, n, conj, *alpha, x, y);
      } else {                  /* Use scheme (2) */
        gen_y_to_cancel(k, n, conj, *alpha, x, y);
      }
      gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
    }
  }

  /* Compute r_truth in double-double */
  r_truth(conj, n, *alpha, x, 1, *beta, y, 1, r, r_true_l, r_true_t);
}

static double ulp(float a)
/*
 * Purpose
 * =======
 *
 * Compute the unit last place of a single precision number.
 */
{
  double f;
  int e;
  f = frexp(a, &e);
  return power(2, e - BITS_S);
}

static float rand_half_1(int l_bits, int *seed)
/*
 * Purpose
 * =======
 *
 * Generate random number in the interval [0.5, 1).
 * l_bits specifies that only the leading l_bits are nonzero.
 *
 */
{
  float a = xrand(seed);        /* [0,1] */
  a /= 2.;
  a += 0.5;
  if (l_bits < BITS_S) {
    float s = power(2, l_bits);
    float t = a / s;            /* shift right l_bits */
    t = (t + a) - a;            /* cancel trailing bits */
    a = t * s;                  /* shift back */
  }
  return a;
}

static void gen_y_to_cancel(int k, int n, enum blas_conj_type conj,
                            float alpha, float *x, float *y)
/*
 * Purpose
 * =======
 *
 * Generate Y(i)'s from k to n-1 to cancel as much as possible.
 *
 */
{
  int i;
  float rtmp;

  for (i = k; i < n; ++i) {
    rtmp = 0.0;
    BLAS_sdot_x(conj, i, alpha, x, 1, 0.0, y, 1, &rtmp, blas_prec_extra);
    if (alpha * x[i] == 0.0)
      y[i] = 0.;
    else
      y[i] = -rtmp / (alpha * x[i]);
  }
}


static void
gen_r_to_cancel(int n, enum blas_conj_type conj, float alpha,
                float beta, float *x, float *y, float *r, int *seed)
/*
 * Purpose
 * =======
 *
 * Generate r to cancel as much as possible.
 *
 */
{
  float rtmp;

  if (beta == 0.0)
    *r = xrand(seed);
  else {
    rtmp = 0.0;
    BLAS_sdot_x(conj, n, alpha, x, 1, 0.0, y, 1, &rtmp, blas_prec_extra);
    *r = -rtmp / beta;
  }
}


static void r_truth(enum blas_conj_type conj, int n, float alpha, const float *x, int incx, float beta, const float *y, int incy, float *r,     /* input */
                    double *r_true_l, double *r_true_t)
{
  int i, ix = 0, iy = 0;
  float *r_i = r;
  const float *x_i = x;
  const float *y_i = y;
  float alpha_i = alpha;
  float beta_i = beta;
  float x_ii;
  float y_ii;
  float r_v;
  double prod_l, prod_t;
  double sum_l, sum_t;
  double tmp1_l, tmp1_t;
  double tmp2_l, tmp2_t;

  /* Immediate return */
  if (n < 0) {
    *r_true_l = *r_true_t = 0.0;
    return;
  }
  r_v = r_i[0];
  sum_l = sum_t = 0.0;          /* sum = 0 */


  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incx;
  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = y_i[iy];
    prod_l = (double) x_ii *y_ii;
    prod_t = 0.0;               /* prod = x[i]*y[i] */
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = sum_l + prod_l;
      e = t1 - sum_l;
      t2 = ((prod_l - e) + (sum_l - (t1 - e))) + sum_t + prod_t;

      /* The result is t1 + t2, after normalization. */
      sum_l = t1 + t2;
      sum_t = t2 - (sum_l - t1);
    }                           /* sum = sum+prod */
    ix += incx;
    iy += incy;
  }                             /* endfor */
  {
    double dt = (double) alpha_i;
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

      con = sum_l * split;
      a11 = con - sum_l;
      a11 = con - a11;
      a21 = sum_l - a11;
      con = dt * split;
      b1 = con - dt;
      b1 = con - b1;
      b2 = dt - b1;

      c11 = sum_l * dt;
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = sum_t * dt;
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      tmp1_l = t1 + t2;
      tmp1_t = t2 - (tmp1_l - t1);
    }
  }                             /* tmp1 = sum*alpha */
  tmp2_l = (double) r_v *beta_i;
  tmp2_t = 0.0;                 /* tmp2 = r*beta */
  {
    /* Compute double-double = double-double + double-double. */
    double e, t1, t2;

    /* Knuth trick. */
    t1 = tmp1_l + tmp2_l;
    e = t1 - tmp1_l;
    t2 = ((tmp2_l - e) + (tmp1_l - (t1 - e))) + tmp1_t + tmp2_t;

    /* The result is t1 + t2, after normalization. */
    tmp1_l = t1 + t2;
    tmp1_t = t2 - (tmp1_l - t1);
  }                             /* tmp1 = tmp1+tmp2 */

  /* Return r_truth */
  *r_true_l = tmp1_l;
  *r_true_t = tmp1_t;
}                               /* end r_truth */
