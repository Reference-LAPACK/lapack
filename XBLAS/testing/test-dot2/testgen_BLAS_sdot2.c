#include <stdlib.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"


static float rand_half_1(int, int *);
static void gen_y_to_cancel(int, int, enum blas_conj_type, float,
                            float *, float *, float *);
static void gen_r_to_cancel(int, enum blas_conj_type, float, float,
                            float *, float *, float *, float *, int *);

void
testgen_BLAS_sdot2(int n, int n_fix2, int n_mix, int norm,
		   enum blas_conj_type conj,
		   float *alpha, int alpha_flag, float *beta, int beta_flag,
		   float *head_x, float *tail_x, float *y, int *seed,
		   float *r, double *r_true_l, double *r_true_t)
/*
 * Purpose
 * =======
 *
 * This routine generates the test vectors (head_X, tail_X) and Y
 * for GEMV2 testing. X contains leading part (head_X) and trainling
 * part (tail_X).
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
 * head_x
 * tail_x  (input/output) double*
 *         Leading and trailing parts of X.
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
  float a, rtmps;
  double rtmpd, eps_out;
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
  s_r_truth2(conj, n_fix2, *alpha, y, 1, 0.0, head_x, tail_x, 1, r,
	     r_true_l, r_true_t);
  B = FixedBits(*r_true_l, *r_true_t);

  /* Pick r at random */
  *r = xrand(seed);

  /* Pick the free X(i)'s at random. */
  for (i = n_fix2 + n_mix; i < n; ++i) {
      head_x[i] = xrand(seed);
      tail_x[i] = xrand(seed) * power(2, -BITS_S);
  }

  if (alpha_flag == 1 && *alpha == 0.0) {
      /* Pick the free Y(i)'s at random. */
      for (i = n_fix2; i < n; ++i) y[i] = xrand(seed);
      /* Compute r_truth in double-double */
      s_r_truth2(conj, n, *alpha, y, 1, *beta, head_x, tail_x, 1, r,
		 r_true_l, r_true_t);
      return;
  }

  if (beta_flag == 1 && *beta == 0.0) { /* alpha != 0 */
      if (B == 0) {     /* prefix sum is zero */
	  switch (y_free) {
	  case 0: break;
	  case 1:
	      y[n_fix2] = xrand(seed);
	      break;
	  case 2:
	      /*
	       * Make SUM_{i=0,1}(x[k+i] * y[k+i]) small ...
	       * head_x[k]*[k] = -head_x[k+1]*y[k+1] exact,
	       * tail_x[k]*y[k] + tail_x[k+1]*y[k+1] small
	       */
	      if (n_mix == 0) {       /* Both x[k] and x[k+1] free. */
		  /* head_x[k]*y[k] + head_x[k+1]*y[k+1] = eps_out^2, small,
		     tail_x[k]*y[k] + tail_x[k+1]*y[k+1] = 0, exact */
		  a = rand_half_1(26, seed);    /* strictly < 1 */
		  head_x[k] = a;
		  y[k] = a;
		  head_x[k + 1] = a + eps_out;  /* exact */
		  y[k + 1] = -a + eps_out;      /* exact */
		  f = power(2, -BITS_S);
		  tail_x[k] = y[k + 1] * f;
		  tail_x[k + 1] = -y[k] * f;
	      } else {  /* x[k] fixed, x[k+1] free. */
		  y[k] = xrand(seed);
		  gen_y_to_cancel(k + 1, n, conj, *alpha, head_x, tail_x, y);
	      }
	      break;
	  default:               /* y_free >= 3 */
	      /*
	       * Make SUM_{i=0,n-1}(x[k+i] * y[k+i]) small
	       * ... Cancel >= 48 bits.
	       */
	      /* Use last 2 to cancel bits, and leading ones to add bits. */
	      y[k] = xrand(seed);
	      rtmpd = *alpha * head_x[k] * y[k] + *alpha * tail_x[k] * y[k];
	      s = 50;
	      for (i = k + 1; i < n - 2; ++i) {
		  rtmpd *= power(2, -s);
		  f = *alpha * head_x[i] + *alpha * tail_x[i];
		  if (f == 0.) y[i] = 0.;
		  else y[i] = rtmpd / f;
	      }
	      gen_y_to_cancel(n - 2, n, conj, *alpha, head_x, tail_x, y);
	      break;
	  }   /* end switch */
      } else {   /* B > 0 */
	  if (B >= BITS_E) {        /* Choose Y(i)'s to cancel. */
	      gen_y_to_cancel(k, n, conj, *alpha, head_x, tail_x, y);
	  } else {                  /* At least y[n-1] is free. */
	      if (y_free == 1) {
		  /* Cancel min(B,24) bits. */
		  gen_y_to_cancel(k, n, conj, *alpha, head_x, tail_x, y);
	      } else {                /* >= 2 frees. */
		  /* There are 2 possibilities:
		   * (1) use 1 to add bits, and y_free-1 to cancel 
		   *     24*(y_free-1) bits
		   * (2) use all to cancel min(B, 24*y_free) bits
		   * Goal is to maximize the # of bits cancelled. By equating
		   * (1) and (2), we find the crossover point is
		   * y_free = B/24 + 1.
		   */
		  if (y_free > B / 24.0 + 1) {  /* Use scheme (1) */
		      rtmps = 0.0;
		      BLAS_sdot2_x(conj, k, *alpha, y, 1, 0.0, head_x, tail_x,
				   1, &rtmps, blas_prec_extra);
		      rtmpd = rtmps;
		      s = 60;     /* Should be random between 40 and 100. */
		      f = *alpha * head_x[k] + *alpha *tail_x[k];
		      if (f == 0.) y[k] = 0.;
		      else y[k] = rtmpd * power(2, -s) / f;
		      gen_y_to_cancel(k+1, n, conj, *alpha, head_x, tail_x, y);
		  } else {              /* Use scheme (2) */
		      gen_y_to_cancel(k, n, conj, *alpha, head_x, tail_x, y);
		  }
	      }
	  }   /* end else B < 106 */
      }   /* end else B > 0 */

    /* Compute r_truth in double-double */
    s_r_truth2(conj, n, *alpha, y, 1, *beta, head_x, tail_x, 1, r,
	       r_true_l, r_true_t);
    return;
  }

  /* Now, beta is non-zero. */
  if (B == 0) { /* prefix sum is zero. */
      switch (y_free) {
      case 0: break;
      case 1:
	  /* Make alpha*x[k]*y[k] + beta*r small. */
	  /* Count number of frees in alpha, x[k], and beta. */
	  frees = 0;
	  if (alpha_flag == 0) ++frees;
	  if (beta_flag == 0)  ++frees;
	  if (n_mix == 0) ++frees;
	  if (frees >= 2 && n_mix == 0) {  /* x[k] must be free */
	      /* Make alpha*head_x[k]*y[k] = -beta*r exact,
		 alpha*tail_x[k]*y[k] small.         */
	      a = rand_half_1(26, seed); /* strictly<1, only leading 26 bits */
	      *r = -a * a;            /* exact */
	      if (alpha_flag == 1) {  /* alpha fixed */
		  *beta = *alpha;
	      } else if (beta_flag == 1) { /* beta fixed */
		  *alpha = *beta;
	      }
	      head_x[k] = a;
	      y[k] = a;
	      f = *alpha * head_x[k] * y[k];
	      s = 60;      /* Should be random between 40 and 100. */
	      if (*alpha * y[k] == 0.) tail_x[k] = 0.;
	      else tail_x[k] = f * power(2, -s) / (*alpha * y[k]);
	  } else {                  /* Cancel 24 bits. */
	      y[k] = xrand(seed);
	      gen_r_to_cancel(n, conj, *alpha, *beta, head_x, tail_x, y,
			      r, seed);
	  }
	  break;
      default:                   /* Actual frees >= 3 */
	  /*
	   * Make SUM_{i=0,n-1}(alpha * x[k+i] * y[k+i]) + beta*r small.
	   * ... Cancel >= 72 bits.
	   * Use last 3 ( Y(n-1), Y(n) and r) to cancel bits, and
	   * leading ones to add bits.
	   */
	  y[k] = xrand(seed);
	  rtmpd = *alpha * head_x[k] * y[k] + *alpha * tail_x[k] * y[k];
	  s = 30;
	  for (i = k + 1; i < n - 2; ++i) {
	      rtmpd *= power(2, -s);
	      f = *alpha * head_x[i] + *alpha * tail_x[i];
	      if (f == 0.) y[i] = 0.;
	      else y[i] = rtmpd / f;
	  }
	  gen_y_to_cancel(n - 2, n, conj, *alpha, head_x, tail_x, y);
	  gen_r_to_cancel(n, conj, *alpha, *beta, head_x, tail_x, y, r, seed);
	  break;
      }  /* end switch */
  } else {  /* B > 0 */
      if (B >= BITS_E) {          /* Choose Y(i)'s and r to cancel. */
	  gen_y_to_cancel(k, n, conj, *alpha, head_x, tail_x, y);
	  gen_r_to_cancel(n, conj, *alpha, *beta, head_x, tail_x, y, r, seed);
      } else {                    /* >= 2 frees. Use y[k] to add bits. */
	  frees = y_free + 1;
	  /* 
	   * There are 2 possibilities:
	   * (1) use 1 to add bits, and frees-1 to cancel 24*(frees-1) bits
	   * (2) use all to cancel min(B, 24*frees) bits
	   * Goal is to maximize the # of bits cancelled. By equating
	   * (1) and (2), we find the crossover point is frees = B/24 + 1.
	   */
	  if (frees > B / (double) BITS_S + 1) {       /* Use scheme (1) */
	      rtmps = 0.0;
	      BLAS_sdot2_x(conj, k, *alpha, y, 1, 0.0, head_x, tail_x, 1,
			   &rtmps, blas_prec_extra);
	      s = 60;           /* Should be random between 40 and 100. */
	      f = *alpha * head_x[k] + *alpha * tail_x[k];
	      if (f == 0.) y[k] = 0.;
	      else y[k] = rtmps * power(2, -s) / f;
	      gen_y_to_cancel(k + 1, n, conj, *alpha, head_x, tail_x, y);
	  } else {               /* Use scheme (2) */
	      gen_y_to_cancel(k, n, conj, *alpha, head_x, tail_x, y);
	  }
	  gen_r_to_cancel(n, conj, *alpha, *beta, head_x, tail_x, y, r, seed);
      }
  }

  /* Compute r_truth in double-double */
  s_r_truth2(conj, n, *alpha, y, 1, *beta, head_x, tail_x, 1, r,
	     r_true_l, r_true_t);
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
			    float alpha, float *head_x, float *tail_x, float *y)
/*
 * Purpose
 * =======
 *
 * Generate Y(i)'s from k to n-1 to cancel as much as possible.
 *
 */
{
  int i;
  float rtmp, f;

  for (i = k; i < n; ++i) {
      rtmp = 0.0;
      BLAS_sdot2_x(conj, i, alpha, y, 1, 0.0, head_x, tail_x, 1,
		   &rtmp, blas_prec_extra);
      f = alpha * head_x[i] + alpha * tail_x[i];
      if (f == 0.) y[i] = 0.;
      else y[i] = -rtmp / f;
  }
}

static void
gen_r_to_cancel(int n, enum blas_conj_type conj, float alpha,
                float beta, float *head_x, float *tail_x, float *y,
		float *r, int *seed)
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
      BLAS_sdot2_x(conj, n, alpha, y, 1, 0.0, head_x, tail_x, 1, &rtmp,
		   blas_prec_extra);
      *r = -rtmp / beta;
  }
}

