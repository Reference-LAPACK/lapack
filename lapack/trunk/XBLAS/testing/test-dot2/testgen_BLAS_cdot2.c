#include <stdlib.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"


/* Complex-Complex Multiplication */
void c_mul(float a[], float b[], float c[])
{
  float cr, ci;
  cr = a[0] * b[0] - a[1] * b[1];
  ci = a[1] * b[0] + a[0] * b[1];
  c[0] = cr;
  c[1] = ci;
}

/* Complex Division c = a/b */
void c_div(float a[], float b[], float c[])
{
  float ratio, den;
  float abr, abi, cr, ci;

  if ((abr = b[0]) < 0.)
    abr = -abr;
  if ((abi = b[1]) < 0.)
    abi = -abi;
  if (abr <= abi) {
    if (abi == 0) {
      BLAS_error("c_div: division by zero", 0, 0, NULL);
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

static void
gen_y_to_cancel(int k, int n, enum blas_conj_type conj,
                void *alpha, void *head_x, void *tail_x, void *y)
/*
 * Purpose
 * =======
 * 
 * Generate Y(i)'s from k to n-1 to cancel as much as possible.
 * 
 */
{
  int i, ii;
  float zero[2] = { 0.0, 0.0 };
  float r[2] = { 0.0, 0.0 };
  float tmp_l[2], tmp_t[2];
  double tmpd[2], tmpd_l[2], tmpd_t[2]; 
  float *head_x_i = head_x, *tail_x_i = tail_x, *y_i = y;
  double r_true_l[2], r_true_t[2];

  for (i = k; i < n; ++i) {
      /* y[i] = -rtmp / (alpha * x[i]); */
      c_r_truth2(conj, i, alpha, y, 1, zero, head_x, tail_x, 1, r,
		 r_true_l, r_true_t);
      ii = 2 * i;
      if (conj == blas_conj) {
	  tmp_l[0] = head_x_i[ii];
	  tmp_l[1] = -head_x_i[ii + 1];
	  c_mul((float *) alpha, tmp_l, tmp_l);
	  tmp_t[0] = tail_x_i[ii];
	  tmp_t[1] = -tail_x_i[ii + 1];
	  c_mul((float *) alpha, tmp_t, tmp_t);
      } else {
	  c_mul((float *) alpha, &head_x_i[ii], tmp_l);
	  c_mul((float *) alpha, &tail_x_i[ii], tmp_t);
      }
      tmpd[0] = tmp_l[0] + tmp_t[0];
      tmpd[1] = tmp_l[1] + tmp_t[1];

      if (tmpd[0] == 0. && tmpd[1] == 0.)
	  y_i[ii] = y_i[ii + 1] = 0.;
      else {
	  z_dddivd(r_true_l, r_true_t, tmpd, tmpd_l, tmpd_t);
	  y_i[ii] = -tmpd_l[0];
	  y_i[ii + 1] = -tmpd_l[1];
      }
  }
}


static void
gen_r_to_cancel(int n, enum blas_conj_type conj,
                void *alpha, void *beta, void *head_x, void *tail_x, 
		void *y, void *r, int *seed)
/*
 * Purpose
 * =======
 * 
 * Generate r to cancel as much as possible.
 * 
 */
{
  float zero[2] = { 0.0, 0.0 };
  float rtmp[2] = { 0.0, 0.0 };
  float *beta_i = (float *) beta;
  double beta_d[2];
  float *r_i = (float *) r;
  double r_true_l[2], r_true_t[2];

  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
      r_i[0] = xrand(seed);
      r_i[1] = xrand(seed);
  } else {
      c_r_truth2(conj, n, alpha, y, 1, zero, head_x, tail_x, 1, rtmp,
		 r_true_l, r_true_t);
      beta_d[0] = beta_i[0];
      beta_d[1] = beta_i[1];
      z_dddivd(r_true_l, r_true_t, beta_d, r_true_l, r_true_t);
      r_i[0] = -r_true_l[0];
      r_i[1] = -r_true_l[1];
  }
}


void
testgen_BLAS_cdot2(int n, int n_fix2, int n_mix, int norm,
		   enum blas_conj_type conj,
		   void *alpha, int alpha_flag, void *beta, int beta_flag,
		   void *head_x, void *tail_x, void *y, int *seed,
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
 * head_x
 * tail_x  (input/output) void*
 *         Leading and trailing parts of X.
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
  float zero[2] = { 0.0, 0.0 };
  float a, b;
  double eps_out;
  float f[2], rtmp[2], rtmp_t[2];
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;
  float *r_i = (float *) r;
  float *head_x_i = (float *) head_x, *tail_x_i = (float *) tail_x;
  float *y_i = (float *) y;

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
  c_r_truth2(conj, n_fix2, alpha, y, 1, zero, head_x, tail_x, 1, r,
	     r_true_l, r_true_t);
  B = FixedBits(r_true_l[0], r_true_t[0]);         /* real */
  B = MAX(B, FixedBits(r_true_l[1], r_true_t[1])); /* imag */

  /* Pick r at random */
  r_i[0] = xrand(seed);
  r_i[1] = xrand(seed);

  /* Pick the free X(i)'s at random. */
  for (i = n_fix2 + n_mix; i < n; ++i) {
      ii = 2 * i;
      head_x_i[ii] = xrand(seed);
      head_x_i[ii + 1] = xrand(seed);
      a = power(2, -BITS_S);
      tail_x_i[ii] = xrand(seed) * a;
      tail_x_i[ii + 1] = xrand(seed) * a;
  }

  if (alpha_flag == 1 && alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
      /* Pick the free Y(i)'s at random. */
      for (i = n_fix2; i < n; ++i) {
	  ii = 2 * i;
	  y_i[ii] = xrand(seed);
	  y_i[ii + 1] = xrand(seed);
      }
      /* Compute r_truth in double-double */
      c_r_truth2(conj, n, alpha, y, 1, beta, head_x, tail_x, 1, r,
		 r_true_l, r_true_t);
      return;
  }

  if (beta_flag==1 && beta_i[0] == 0.0 && beta_i[1] == 0.0) { /* alpha != 0 */
      if (B == 0) {   /* prefix sum is zero. */
	  switch (y_free) {
	  case 0: break;
	  case 1:
	      y_i[k] = xrand(seed);
	      y_i[k + 1] = xrand(seed);
	      break;
	  case 2:
	      /*
	       * Make SUM_{i=0,1}(x[k+i] * y[k+i]) small ... alpha * eps^2
	       * head_x[k]*y[k] = -head_x[k+1]*y[k+1] exact,
	       * tail_x[k]*y[k] + tail_x[k+1]*y[k+1] small
	       * For complex, each real number is multiplied by (i+1),
	       * the result is 2i * eps_out^2.
	       */
	      if (n_mix == 0) {       /* Both x[k] and x[k+1] free. */
		  /* head_x[k]*y[k] + head_x[k+1]*y[k+1] = eps_out^2, small,
		   * tail_x[k]*y[k] + tail_x[k+1]*y[k+1] = 0, exact.
		   */
		  a = rand_half_1(BITS_S, seed);    /* strictly < 1 */
		  head_x_i[k] = a;           /* real */
		  head_x_i[k + 1] = a;       /* imag */
		  y_i[k] = a;
		  y_i[k + 1] = a;
		  head_x_i[k + 2] = a + eps_out;     /* exact */
		  head_x_i[k + 3] = a + eps_out;     /* exact */
		  y_i[k + 2] = -a + eps_out;    /* exact */
		  y_i[k + 3] = -a + eps_out;    /* exact */
		  b = power(2, -BITS_S);        /* shift right */
		  tail_x_i[k] = y_i[k+2] * b;
		  tail_x_i[k+1] = tail_x_i[k];
		  tail_x_i[k+2] = -y_i[k] * b;
		  tail_x_i[k+3] = tail_x_i[k+2];
	      } else {  /* x[k] fixed, x[k+1] fixed or free. */
		  y_i[k] = xrand(seed);
		  y_i[k + 1] = xrand(seed);
		  gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, head_x, tail_x, y);
	      }
	      break;
	  default:                 /* y_free >= 3 */
	      /*
	       * Make SUM_{i=0,n-1}(x[k+i] * y[k+i]) small.
	       * Use last 2 to cancel bits, and leading ones to add bits.
	       * ... Cancel >= 48 bits.
	       */
	      y_i[k] = xrand(seed);
	      y_i[k + 1] = xrand(seed);
	      rtmp[0] = head_x_i[k];	  /* leading part */
	      if (conj == blas_conj) rtmp[1] = -head_x_i[k + 1];
	      else rtmp[1] = head_x_i[k + 1];
	      c_mul(rtmp, &y_i[k], rtmp);
	      c_mul(alpha_i, rtmp, rtmp);
	      rtmp_t[0] = tail_x_i[k];    /* imaginary part */
	      if (conj == blas_conj) rtmp_t[1] = -tail_x_i[k + 1];
	      else rtmp_t[1] = tail_x_i[k + 1];
	      c_mul(rtmp_t, &y_i[k], rtmp_t);
	      c_mul(alpha_i, rtmp_t, rtmp_t);
	      rtmp[0] += rtmp_t[0];
	      rtmp[1] += rtmp_t[1];
	      s = 40;
	      b = power(2, -s);
	      for (i = n_fix2 + 1; i < n - 2; ++i) {
		  ii = 2 * i;
		  rtmp[0] *= b;
		  rtmp[1] *= b;
		  f[0] = head_x_i[ii];
		  if (conj == blas_conj) f[1] = -head_x_i[ii + 1];
		  else f[1] = head_x_i[ii + 1];
		  c_mul(alpha_i, f, f);
		  rtmp_t[0] = tail_x_i[ii];
		  if (conj == blas_conj) rtmp_t[1] = -tail_x_i[ii + 1];
		  else f[1] = tail_x_i[ii + 1];
		  c_mul(alpha_i, rtmp_t, rtmp_t);
		  f[0] += rtmp_t[0];
		  f[1] += rtmp_t[1];
		  if (f[0] == 0. && f[1] == 0.)
		      y_i[ii] = y_i[ii + 1] = 0.;
		  else
		      c_div(rtmp, f, &y_i[ii]);
	      }
	      gen_y_to_cancel(n - 2, n, conj, alpha, head_x, tail_x, y);
	      break;
	  }                         /* end switch */
      } else {      /* B > 0 */
	  if (B >= BITS_E) {        /* Choose Y(i)'s to cancel. */
	      gen_y_to_cancel(n_fix2, n, conj, alpha, head_x, tail_x, y);
	  } else {                  /* At least y[n-1] is free. */
	      if (y_free == 1) {
		  /* Cancel min(B,53) bits. */
		  gen_y_to_cancel(n_fix2, n, conj, alpha, head_x, tail_x, y);
	      } else {              /* >= 2 frees in Y. */
		  /*
		   * There are 2 possibilities:
		   * (1) use 1 to add bits, and y_free-1 to cancel
		   *     53*(y_free-1) bits
		   * (2) use all to cancel min(B, 53*y_free) bits
		   * Goal is to maximize the # of bits cancelled. By equating
		   * (1) and (2), we find the crossover point is
		   * y_free = B/53 + 1.
		   */
		  if (y_free > B / (double) BITS_D + 1) { /* Use scheme (1) */
		      f[0] = f[1] = 0.0;
		      c_r_truth2(conj, n_fix2, alpha, y, 1, zero,
				 head_x, tail_x, 1, f, r_true_l, r_true_t);
		      f[0] = r_true_l[0];
		      f[1] = r_true_l[1];
		      s = 100;  /* Should be random between 40 and 100. */
		      b = power(2, -s);
		      f[0] *= b;
		      f[1] *= b;
		      rtmp[0] = head_x_i[k];  /* leading part */
		      if (conj == blas_conj) rtmp[1] = -head_x_i[k + 1];
		      else rtmp[1] = head_x_i[k + 1];
		      c_mul(alpha_i, rtmp, rtmp);
		      rtmp_t[0] = tail_x_i[k];  /* trailing part */
		      if (conj == blas_conj) rtmp_t[1] = -tail_x_i[k + 1];
		      else rtmp_t[1] = tail_x_i[k + 1];
		      c_mul(alpha_i, rtmp_t, rtmp_t);
		      rtmp[0] += rtmp_t[0];
		      rtmp[1] += rtmp_t[1];
		      if (rtmp[0] == 0. && rtmp[1] == 0.)
			  y_i[k] = y_i[k + 1] = 0.;
		      else
			  c_div(f, rtmp, &y_i[k]);
		      gen_y_to_cancel(n_fix2 + 1, n, conj, alpha,
				      head_x, tail_x, y);
		  } else {              /* Use scheme (2) */
		      gen_y_to_cancel(n_fix2, n, conj, alpha,
				      head_x, tail_x, y);
		  }
	      }
	  }     /* end else B < 106 */
      }         /* end else B > 0 */

      /* Compute r_truth in double-double */
      c_r_truth2(conj, n, alpha, y, 1, zero, head_x, tail_x, 1, r,
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
	  if (beta_flag == 0) ++frees;
	  if (n_mix == 0) ++frees;
	  if (frees >= 2 && n_mix == 0) { /* x[k] must be free */
	      /* Make alpha*head_x[k]*y[k] = -beta*r exact,
	       * alpha*tail_x[k]*y[k] small.
	       * For complex, each real number is multiplied by (i+1).
	       */
	      a = rand_half_1(26, seed); /* strictly<1, only leading 26 bits */
	      r_i[0] = -a * a;    /* real, exact */
	      r_i[1] = r_i[0];    /* imag */
	      if (alpha_flag == 1) {  /* alpha fixed, beta must be free */
		  beta_i[0] = alpha_i[0];
		  beta_i[1] = alpha_i[1];
	      } else if (beta_flag == 1) { /* beta fixed, alpha must be free */
		  alpha_i[0] = beta_i[0];
		  alpha_i[1] = beta_i[1];
	      }
	      head_x_i[k] = a;   /* real, exact */
	      head_x_i[k+1] = a; /* imag, exact */
	      if (conj == blas_conj) head_x_i[k + 1] = -head_x_i[k + 1];
	      y_i[k] = a;     /* exact */
	      y_i[k + 1] = a; /* exact */
	      /* f = *alpha * head_x[k] * y[k]; */
	      c_mul(alpha_i, &head_x_i[k], f);
	      c_mul(f, &y_i[k], f);
	      c_mul(alpha_i, &y_i[k], rtmp);
	      s = power(2, -100); /* Should be random between 40 and 100. */
	      if (rtmp[0] == 0. && rtmp[1] == 0.) {
		  tail_x_i[k] = 0.;
		  tail_x_i[k+1] = 0.;
	      } else {
		  f[0] *= s;
		  f[1] *= s;
		  c_div(f, rtmp, &tail_x_i[k]);
	      }
	  } else {     /* Cancel 53 bits. */
	      y_i[k] = xrand(seed);
	      y_i[k+1] = xrand(seed);
	      gen_r_to_cancel(n, conj, alpha, beta, head_x, tail_x, y,
			      r, seed);
	  }
	  break;
      default:                 /* Actual frees >= 3 */
	  /*
	   * Make SUM_{i=0,n-1}(alpha * x[k+i] * y[k+i]) + beta*r small.
	   * Use last 2 ( Y(n) and r) to cancel bits, and leading ones to
	   * add bits ... Cancel >= 106 bits.
	   */
	  y_i[k] = xrand(seed);
	  y_i[k + 1] = xrand(seed);
	  rtmp[0] = head_x_i[k];         /* leading part */
	  if (conj == blas_conj) rtmp[1] = -head_x_i[k + 1];
	  else rtmp[1] = head_x_i[k + 1];
	  c_mul(rtmp, &y_i[k], f);
	  c_mul(alpha_i, f, rtmp);
	  rtmp_t[0] = tail_x_i[k];       /* trailing part */
	  if (conj == blas_conj) rtmp_t[1] = -tail_x_i[k + 1];
	  else rtmp_t[1] = tail_x_i[k + 1];
	  c_mul(rtmp_t, &y_i[k], f);
	  c_mul(alpha_i, f, rtmp_t);
	  rtmp[0] += rtmp_t[0];
	  rtmp[1] += rtmp_t[1];
	  s = 30;
	  b = power(2, -s);
	  for (i = n_fix2 + 1; i < n - 1; ++i) {
	      ii = 2 * i;
	      rtmp[0] *= b;
	      rtmp[1] *= b;
	      f[0] = head_x_i[ii];      /* leading part */
	      if (conj == blas_conj) f[1] = -head_x_i[ii + 1];
	      else f[1] = head_x_i[ii + 1];
	      c_mul(alpha_i, f, f);
	      rtmp_t[0] = tail_x_i[ii]; /* trailing part */
	      if (conj == blas_conj) rtmp_t[1] = -tail_x_i[ii + 1];
	      else rtmp_t[1] = tail_x_i[ii + 1];
	      c_mul(alpha_i, rtmp_t, rtmp_t);
	      f[0] += rtmp_t[0];
	      f[1] += rtmp_t[1];
	      if (f[0] == 0. && f[1] == 0.)
		  y_i[ii] = y_i[ii + 1] = 0.;
	      else
		  c_div(rtmp, f, &y_i[ii]);
	  }
	  gen_y_to_cancel(n - 1, n, conj, alpha, head_x, tail_x, y);
	  gen_r_to_cancel(n, conj, alpha, beta, head_x, tail_x, y, r, seed);
	  break;
      }  /* end switch */
  } else {    /* B > 0 */
      if (B >= BITS_E) {        /* Choose Y(i)'s and r to cancel. */
	  gen_y_to_cancel(n_fix2, n, conj, alpha, head_x, tail_x, y);
	  gen_r_to_cancel(n, conj, alpha, beta, head_x, tail_x, y, r, seed);
      } else {                  /* >= 2 frees. Use y[k] to add bits. */
	  frees = y_free + 1;
	  /*
	   * There are 2 possibilities:
	   * (1) use 1 to add bits, and y_free-1 to cancel 24*(y_free-1) bits
	   * (2) use all to cancel min(B, 24*y_free) bits
	   * Goal is to maximize the # of bits cancelled. By equating (1)
	   * and (2), we find the crossover point is y_free = B/24 + 1.
	   */
	  if (frees > B / (float) BITS_S + 1) {    /* Use scheme (1) */
	      f[0] = f[1] = 0.0;
	      c_r_truth2(conj, n_fix2, alpha, y, 1, zero, head_x, tail_x, 1, f,
			 r_true_l, r_true_t);
	      f[0] = r_true_l[0];
	      f[1] = r_true_l[1];
	      s = 60;    /* Should be random between 40 and 100. */
	      b = power(2, -s);
	      f[0] *= b;
	      f[1] *= b;
	      rtmp[0] = head_x_i[k];    /* leading part */
	      if (conj == blas_conj) rtmp[1] = -head_x_i[k + 1];
	      else rtmp[1] = -head_x_i[k + 1];
	      c_mul(alpha_i, rtmp, rtmp);
	      rtmp_t[0] = tail_x_i[k];  /* trailing part */
	      if (conj == blas_conj) rtmp_t[1] = -tail_x_i[k + 1];
	      else rtmp_t[1] = -tail_x_i[k + 1];
	      c_mul(alpha_i, rtmp_t, rtmp_t);
	      rtmp[0] += rtmp_t[0];
	      rtmp[1] += rtmp_t[1];
	      if (rtmp[0] == 0. && rtmp[1] == 0.) {
		  y_i[k] = y_i[k + 1] = 0.;
	      } else {
		  c_div(f, rtmp, &y_i[k]);
	      }
	      gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, head_x, tail_x, y);
	  } else {                  /* Use scheme (2) */
	      gen_y_to_cancel(n_fix2, n, conj, alpha, head_x, tail_x, y);
	  }
	  gen_r_to_cancel(n, conj, alpha, beta, head_x, tail_x, y, r, seed);
      }
  }

  /* Compute r_truth in double-double */
  c_r_truth2(conj, n, alpha, y, 1, beta, head_x, tail_x, 1, r,
	     r_true_l, r_true_t);
}                               /* testgen_BLAS_zdot */
