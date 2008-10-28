#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

double power(int i1, int i2)
{
  int i, j;
  double r = 1.0;

  if (i2 < 0)
    j = -i2;
  else
    j = i2;
  for (i = 0; i < j; ++i)
    r *= i1;
  if (i2 < 0)
    r = 1. / r;

  return r;
}

double xrand(int *is)
/*
 * XRAND returns a uniformly distributed pseudorandom number in (0, 1).
 * IS is the seed, and is changed with each call.  The period of this
 * linear congruential generator is 2^26, according to Knuth vol. 2.
 */
{
  double s1, s2, ret_val;
#define f7    78125.0           /* 5.d0 ** 7 */
#define r26   1.4901161193847656e-8     /* 2^(-26) */
#define r28   3.7252902984619141e-9     /* 2^(-28) */
#define t28   268435456.0       /* 2^28 */

  s1 = *is;
  s2 = fmod(f7 * s1, t28);
  ret_val = (s1 + r26 * s2) * r28;
  *is = s2;

  return ret_val;
}

int FixedBits(double r_true_l, double r_true_t)
/*
 * Purpose
 * =======
 *
 * Compute the number of fixed bits in r_true.
 * (r_true_l, r_true_t) is double-double representation of r_true.
 *
 */
{

  int b;                        /* Number of fixed bits in r_true */
  int i, k;
  double tmp_l, tmp_t;
  double res[5], t, temp;

  b = k = 0;
  res[0] = r_true_l;
  while (res[k] != 0.0) {       /* Each time cancel 53 bits */
    tmp_l = r_true_l;
    tmp_t = r_true_t;
    for (i = 0; i <= k; ++i) {  /* tmp = tmp - res[i] */
      t = -res[i];
      {
        /* Compute double-double = double-double + double. */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = tmp_l + t;
        e = t1 - tmp_l;
        t2 = ((t - e) + (tmp_l - (t1 - e))) + tmp_t;

        /* The result is t1 + t2, after normalization. */
        tmp_l = t1 + t2;
        tmp_t = t2 - (tmp_l - t1);
      }
    }
    b += BITS_D;
    ++k;
    res[k] = tmp_l;
  }

  if (k > 0) {
    int e;
    double f;

    /* Use Kahan's trick for the last nonzero residual. */
    b -= BITS_D;
    --k;
    t = fabs(res[k]);
    f = frexp(t, &e);           /* 1/2 <= f < 1 */
    /* g = f * 2;              1 <= g < 2 */
    for (i = 1;; ++i) {         /* Compute number of bits in g */
      t = f * power(2, i);      /* Shift left i bits */
      temp = floor(t);
      t -= temp;
      if (t == 0.0)
        break;
      ++b;
    }
  }
  return b;
}                               /* FixedBits */


void ddadd(double dda_l, double dda_t, double ddb_l, double ddb_t,
           double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 *
 * This subroutine computes ddc = dda + ddb.
 *
 * Taken from D. H. Bailey's ddfun90.f.
 *
 */
  double e, t1, t2;

  /* Compute dda + ddb using Knuth's trick. */
  t1 = dda_l + ddb_l;
  e = t1 - dda_l;
  t2 = ((ddb_l - e) + (dda_l - (t1 - e))) + dda_t + ddb_t;

  /* The result is t1 + t2, after normalization. */
  *ddc_l = t1 + t2;
  *ddc_t = t2 - (*ddc_l - t1);

}                               /* end ddadd */

void ddmuld(double dda_l, double dda_t, double db,
            double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 *
 * This routine multiplies the DD number DDA by the DP number DB to yield
 * the DD product DDC.
 *
 * Taken from D. H. Bailey's ddfun90.f.
 *
 */
  double a1, a2, b1, b2, cona, conb, c11, c21, c2, e, t1, t2;

  /* This splits dda(1) and db into high-order and low-order words. */
  cona = dda_l * split;
  conb = db * split;
  a1 = cona - (cona - dda_l);
  b1 = conb - (conb - db);
  a2 = dda_l - a1;
  b2 = db - b1;

  /* Multilply dda(1) * db using Dekker's method. */
  c11 = dda_l * db;
  c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

  /* Compute dda(2) * db (only high-order word is needed). */
  c2 = dda_t * db;

  /* Compute (c11, c21) + c2 using Knuth's trick. */
  t1 = c11 + c2;
  e = t1 - c11;
  t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

  /* The result is t1 + t2, after normalization. */
  *ddc_l = t1 + t2;
  *ddc_t = t2 - (*ddc_l - t1);

}                               /* end ddmuld */

void dddiv(double dda_l, double dda_t,
           double ddb_l, double ddb_t, double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 * 
 * This divides the DD number DDA by the DD number DDB to yield the DD
 * quotient DDC.
 * 
 * Taken from D. H. Bailey's ddfun90.f.
 *
 */
  double a1, a2, b1, b2, cona, conb, c11, c2, c21, e,
    s1, s2, t1, t2, t11, t12, t21, t22;

  /* Compute a DP approximation to the quotient. */
  s1 = dda_l / ddb_l;

  /* This splits s1 and ddb(1) into high-order and low-order words. */
  cona = s1 * split;
  conb = ddb_l * split;
  a1 = cona - (cona - s1);
  b1 = conb - (conb - ddb_l);
  a2 = s1 - a1;
  b2 = ddb_l - b1;

  /* Multiply s1 * ddb(1) using Dekker's method. */
  c11 = s1 * ddb_l;
  c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

  /* Compute s1 * ddb(2) (only high-order word is needed). */
  c2 = s1 * ddb_t;

  /* Compute (c11, c21) + c2 using Knuth's trick. */
  t1 = c11 + c2;
  e = t1 - c11;
  t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

  /* The result is t1 + t2, after normalization. */
  t12 = t1 + t2;
  t22 = t2 - (t12 - t1);

  /* Compute dda - (t12, t22) using Knuth's trick. */
  t11 = dda_l - t12;
  e = t11 - dda_l;
  t21 = ((-t12 - e) + (dda_l - (t11 - e))) + dda_t - t22;

  /* Compute high-order word of (t11, t21) and divide by ddb(1). */
  s2 = (t11 + t21) / ddb_l;

  /* The result is s1 + s2, after normalization. */
  *ddc_l = s1 + s2;
  *ddc_t = s2 - (*ddc_l - s1);

}                               /* end dddiv */


void z_ddmuld(double *dda_l, double *dda_t, double *db,
              double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 * 
 * This multiplies the complex DD number DDA by the complex D number DB to
 * yield the complex DD product DDC.
 * 
 */
  double t_l, t_t, t1_l, t1_t, t2_l, t2_t;
  /* real part */
  ddmuld(dda_l[0], dda_t[0], db[0], &t1_l, &t1_t);
  ddmuld(dda_l[1], dda_t[1], -db[1], &t2_l, &t2_t);
  ddadd(t1_l, t1_t, t2_l, t2_t, &t_l, &t_t);
  ddc_l[0] = t_l;
  ddc_t[0] = t_t;
  /* imaginary part */
  ddmuld(dda_l[1], dda_t[1], db[0], &t1_l, &t1_t);
  ddmuld(dda_l[0], dda_t[0], db[1], &t2_l, &t2_t);
  ddadd(t1_l, t1_t, t2_l, t2_t, &t_l, &t_t);
  ddc_l[1] = t_l;
  ddc_t[1] = t_t;
}

void z_dddivd(double *dda_l, double *dda_t, double *db,
              double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 * 
 * This divides the complex DD number DDA by the complex DP number DB to
 * yield the complex DD quotient DDC.
 * 
 */
  double db_conj[2];
  double t_l[2], t_t[2];
  double d_l, d_t;

  /* b_r^2 + b_i^2 in double-double */
  d_l = db[0];
  d_t = 0.0;
  ddmuld(d_l, d_t, d_l, &t_l[0], &t_t[0]);
  d_l = db[1];
  d_t = 0.0;
  ddmuld(d_l, d_t, d_l, &t_l[1], &t_t[1]);
  ddadd(t_l[0], t_t[0], t_l[1], t_t[1], &d_l, &d_t);

  db_conj[0] = db[0];
  db_conj[1] = -db[1];
  z_ddmuld(dda_l, dda_t, db_conj, t_l, t_t);
  /*printf("\tz_dddivd() -> a * conj(b) = %10.8e + %10.8ei\n",
     t_l[0],t_l[1]); */
  dddiv(t_l[0], t_t[0], d_l, d_t, &ddc_l[0], &ddc_t[0]);
  dddiv(t_l[1], t_t[1], d_l, d_t, &ddc_l[1], &ddc_t[1]);
}
