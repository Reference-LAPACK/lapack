#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
#define NORM_START -1
#define NORM_END    1
#define INC_START  -2
#define INC_END     2
#define PREC_START  0
#define PREC_END    2



void do_test_ssum_x(int n, int ntests, int *seed, double thresh,
		    int debug, float test_prob, double *min_ratio,
		    double *max_ratio, int *num_bad_ratio, int *num_tests)

/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on sum  
 *
 * Arguments
 * =========
 *
 * n         (input) int
 *           The size of vector being tested
 *
 * ntests    (input) int
 *           The number of tests to run for each set of attributes.
 *
 * seed      (input/output) int         
 *           The seed for the random number generator used in testgen().
 *
 * thresh    (input) double
 *           When the ratio returned from test() exceeds the specified
 *           threshold, the current size, r_true, r_comp, and ratio will be
 *           printed.  (Since ratio is supposed to be O(1), we can set thresh
 *           to ~10.)
 *
 * debug     (input) int
 *           If debug=3, print summary 
 *           If debug=2, print summary only if the number of bad ratios > 0
 *           If debug=1, print complete info if tests fail
 *           If debug=0, return max ratio
 *
 * min_ratio (output) double
 *           The minimum ratio
 * 
 * num_bad_ratio (output) int
 *               The number of tests fail; they are above the threshold.
 *
 * num_tests (output) int
 *           The number of tests is being performed.
 *
 */
{

  /* function name */
  const char fname[] = "BLAS_ssum_x";

  int i, j;
  int norm;
  int xi, incx_val, incx;
  double ratio_max, ratio_min;
  double ratio;
  int saved_seed;
  double eps_int;
  double un_int;

  int test_count;
  int bad_ratio_count;

  float *x;
  float *x_gen;
  float x_elem;

  double head_sum_true, tail_sum_true;

  float sum;

  int prec_val;
  enum blas_prec_type prec;
  int x_gen_i, incx_gen;

  FPU_FIX_DECL;


  /* test for bad arguments */
  if (n < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* if there is nothing to test, return all zero */
  if (n == 0 || ntests == 0) {
    *min_ratio = 0.0;
    *max_ratio = 0.0;
    *num_bad_ratio = 0;
    *num_tests = 0;
    return;
  }

  FPU_FIX_START;

  incx_gen = 1;


  /* get space for calculation */
  x = (float *) blas_malloc(n * 2 * sizeof(float));
  if (n * 2 > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_gen = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* initialization */
  ratio_min = 1e308;
  ratio_max = 0.0;
  test_count = 0;
  bad_ratio_count = 0;


  /* varying extra precs */
  for (prec_val = PREC_START; prec_val <= PREC_END; prec_val++) {
    switch (prec_val) {
    case 0:
      eps_int = power(2, -BITS_S);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_single),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_single));
      prec = blas_prec_single;
      break;
    case 1:
      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;
      break;
    case 2:
    default:
      eps_int = power(2, -BITS_E);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
      prec = blas_prec_extra;
      break;
    }

    for (norm = NORM_START; norm <= NORM_END; norm++) {

      for (incx_val = INC_START; incx_val <= INC_END; incx_val++) {

	if (incx_val == 0)
	  continue;

	for (i = 0; i < ntests; i++) {

	  /* For the sake of speed, we throw out this case at random */
	  if (xrand(seed) >= test_prob)
	    continue;

	  saved_seed = *seed;

	  BLAS_ssum_testgen(n, norm, x_gen, seed,
			    &head_sum_true, &tail_sum_true);

	  test_count++;

	  incx = incx_val;


	  xi = 0;
	  if (incx < 0)
	    xi = -(n - 1) * incx;

	  /* copy x_gen to x */
	  incx_gen = 1;

	  for (j = 0, x_gen_i = 0; j < n;
	       j++, x_gen_i += incx_gen, xi += incx) {
	    x_elem = x_gen[x_gen_i];
	    x[xi] = x_elem;
	  }

	  FPU_FIX_STOP;
	  BLAS_ssum_x(n, x, incx_val, &sum, prec);
	  FPU_FIX_START;

	  test_BLAS_ssum(n, sum, head_sum_true, tail_sum_true, x, incx_val,
			 eps_int, un_int, &ratio);

	  /* The !<= below causes NaN error to be detected.
	     Note that (NaN > thresh) is always false */
	  if (!(ratio <= thresh)) {
	    bad_ratio_count++;

	    if (debug == 3) {
	      printf("Seed = %d\n", saved_seed);
	      printf("n = %d\n", n);
	      printf("norm = %d\n", norm);
	      switch (prec) {
	      case blas_prec_single:
		printf("single ");
		break;
	      case blas_prec_double:
		printf("double ");
		break;
	      case blas_prec_indigenous:
		printf("indigenous ");
		break;
	      case blas_prec_extra:
		printf("extra ");
		break;
	      }

	      printf("incx = %d\n", incx);

	      /* print out the vector */
	      xi = (incx < 0) ? -(n - 1) * incx : 0;
	      printf(" [ ");
	      for (j = 0; j < n; j++, xi += incx) {
		printf("%16.8e", x[xi]);
	      }
	      printf("]\n");

	      printf("sum = ");
	      printf("%16.8e", sum);
	      printf("\n");
	      printf("ratio = %.4e\n", ratio);
	      printf("sum_true = ");
	      printf("[%24.16e %24.16e]", head_sum_true, tail_sum_true);

	    }			/* end of if (debug == 3) */
	  }
	  /* end of if (ratio > thresh) */
	  if (ratio > ratio_max)
	    ratio_max = ratio;
	  if (ratio != 0.0 && ratio < ratio_min)
	    ratio_min = ratio;

	}			/* end of ntests loop */
      }				/* end of incx loop */
    }				/* end of norm loop */

  }				/* end of prec loop */

  FPU_FIX_STOP;

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

  blas_free(x);
  blas_free(x_gen);


}
void do_test_dsum_x(int n, int ntests, int *seed, double thresh,
		    int debug, float test_prob, double *min_ratio,
		    double *max_ratio, int *num_bad_ratio, int *num_tests)

/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on sum  
 *
 * Arguments
 * =========
 *
 * n         (input) int
 *           The size of vector being tested
 *
 * ntests    (input) int
 *           The number of tests to run for each set of attributes.
 *
 * seed      (input/output) int         
 *           The seed for the random number generator used in testgen().
 *
 * thresh    (input) double
 *           When the ratio returned from test() exceeds the specified
 *           threshold, the current size, r_true, r_comp, and ratio will be
 *           printed.  (Since ratio is supposed to be O(1), we can set thresh
 *           to ~10.)
 *
 * debug     (input) int
 *           If debug=3, print summary 
 *           If debug=2, print summary only if the number of bad ratios > 0
 *           If debug=1, print complete info if tests fail
 *           If debug=0, return max ratio
 *
 * min_ratio (output) double
 *           The minimum ratio
 * 
 * num_bad_ratio (output) int
 *               The number of tests fail; they are above the threshold.
 *
 * num_tests (output) int
 *           The number of tests is being performed.
 *
 */
{

  /* function name */
  const char fname[] = "BLAS_dsum_x";

  int i, j;
  int norm;
  int xi, incx_val, incx;
  double ratio_max, ratio_min;
  double ratio;
  int saved_seed;
  double eps_int;
  double un_int;

  int test_count;
  int bad_ratio_count;

  double *x;
  double *x_gen;
  double x_elem;

  double head_sum_true, tail_sum_true;

  double sum;

  int prec_val;
  enum blas_prec_type prec;
  int x_gen_i, incx_gen;

  FPU_FIX_DECL;


  /* test for bad arguments */
  if (n < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* if there is nothing to test, return all zero */
  if (n == 0 || ntests == 0) {
    *min_ratio = 0.0;
    *max_ratio = 0.0;
    *num_bad_ratio = 0;
    *num_tests = 0;
    return;
  }

  FPU_FIX_START;

  incx_gen = 1;


  /* get space for calculation */
  x = (double *) blas_malloc(n * 2 * sizeof(double));
  if (n * 2 > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_gen = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* initialization */
  ratio_min = 1e308;
  ratio_max = 0.0;
  test_count = 0;
  bad_ratio_count = 0;


  /* varying extra precs */
  for (prec_val = PREC_START; prec_val <= PREC_END; prec_val++) {
    switch (prec_val) {
    case 0:
      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;
      break;
    case 1:
      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;
      break;
    case 2:
    default:
      eps_int = power(2, -BITS_E);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
      prec = blas_prec_extra;
      break;
    }

    for (norm = NORM_START; norm <= NORM_END; norm++) {

      for (incx_val = INC_START; incx_val <= INC_END; incx_val++) {

	if (incx_val == 0)
	  continue;

	for (i = 0; i < ntests; i++) {

	  /* For the sake of speed, we throw out this case at random */
	  if (xrand(seed) >= test_prob)
	    continue;

	  saved_seed = *seed;

	  BLAS_dsum_testgen(n, norm, x_gen, seed,
			    &head_sum_true, &tail_sum_true);

	  test_count++;

	  incx = incx_val;


	  xi = 0;
	  if (incx < 0)
	    xi = -(n - 1) * incx;

	  /* copy x_gen to x */
	  incx_gen = 1;

	  for (j = 0, x_gen_i = 0; j < n;
	       j++, x_gen_i += incx_gen, xi += incx) {
	    x_elem = x_gen[x_gen_i];
	    x[xi] = x_elem;
	  }

	  FPU_FIX_STOP;
	  BLAS_dsum_x(n, x, incx_val, &sum, prec);
	  FPU_FIX_START;

	  test_BLAS_dsum(n, sum, head_sum_true, tail_sum_true, x, incx_val,
			 eps_int, un_int, &ratio);

	  /* The !<= below causes NaN error to be detected.
	     Note that (NaN > thresh) is always false */
	  if (!(ratio <= thresh)) {
	    bad_ratio_count++;

	    if (debug == 3) {
	      printf("Seed = %d\n", saved_seed);
	      printf("n = %d\n", n);
	      printf("norm = %d\n", norm);
	      switch (prec) {
	      case blas_prec_single:
		printf("single ");
		break;
	      case blas_prec_double:
		printf("double ");
		break;
	      case blas_prec_indigenous:
		printf("indigenous ");
		break;
	      case blas_prec_extra:
		printf("extra ");
		break;
	      }

	      printf("incx = %d\n", incx);

	      /* print out the vector */
	      xi = (incx < 0) ? -(n - 1) * incx : 0;
	      printf(" [ ");
	      for (j = 0; j < n; j++, xi += incx) {
		printf("%24.16e", x[xi]);
	      }
	      printf("]\n");

	      printf("sum = ");
	      printf("%24.16e", sum);
	      printf("\n");
	      printf("ratio = %.4e\n", ratio);
	      printf("sum_true = ");
	      printf("[%24.16e %24.16e]", head_sum_true, tail_sum_true);

	    }			/* end of if (debug == 3) */
	  }
	  /* end of if (ratio > thresh) */
	  if (ratio > ratio_max)
	    ratio_max = ratio;
	  if (ratio != 0.0 && ratio < ratio_min)
	    ratio_min = ratio;

	}			/* end of ntests loop */
      }				/* end of incx loop */
    }				/* end of norm loop */

  }				/* end of prec loop */

  FPU_FIX_STOP;

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

  blas_free(x);
  blas_free(x_gen);


}
void do_test_csum_x(int n, int ntests, int *seed, double thresh,
		    int debug, float test_prob, double *min_ratio,
		    double *max_ratio, int *num_bad_ratio, int *num_tests)

/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on sum  
 *
 * Arguments
 * =========
 *
 * n         (input) int
 *           The size of vector being tested
 *
 * ntests    (input) int
 *           The number of tests to run for each set of attributes.
 *
 * seed      (input/output) int         
 *           The seed for the random number generator used in testgen().
 *
 * thresh    (input) double
 *           When the ratio returned from test() exceeds the specified
 *           threshold, the current size, r_true, r_comp, and ratio will be
 *           printed.  (Since ratio is supposed to be O(1), we can set thresh
 *           to ~10.)
 *
 * debug     (input) int
 *           If debug=3, print summary 
 *           If debug=2, print summary only if the number of bad ratios > 0
 *           If debug=1, print complete info if tests fail
 *           If debug=0, return max ratio
 *
 * min_ratio (output) double
 *           The minimum ratio
 * 
 * num_bad_ratio (output) int
 *               The number of tests fail; they are above the threshold.
 *
 * num_tests (output) int
 *           The number of tests is being performed.
 *
 */
{

  /* function name */
  const char fname[] = "BLAS_csum_x";

  int i, j;
  int norm;
  int xi, incx_val, incx;
  double ratio_max, ratio_min;
  double ratio;
  int saved_seed;
  double eps_int;
  double un_int;

  int test_count;
  int bad_ratio_count;

  float *x;
  float *x_gen;
  float x_elem[2];

  double head_sum_true[2], tail_sum_true[2];

  float sum[2];

  int prec_val;
  enum blas_prec_type prec;
  int x_gen_i, incx_gen;

  FPU_FIX_DECL;


  /* test for bad arguments */
  if (n < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* if there is nothing to test, return all zero */
  if (n == 0 || ntests == 0) {
    *min_ratio = 0.0;
    *max_ratio = 0.0;
    *num_bad_ratio = 0;
    *num_tests = 0;
    return;
  }

  FPU_FIX_START;

  incx_gen = 1;
  incx_gen *= 2;

  /* get space for calculation */
  x = (float *) blas_malloc(n * 2 * sizeof(float) * 2);
  if (n * 2 > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* initialization */
  ratio_min = 1e308;
  ratio_max = 0.0;
  test_count = 0;
  bad_ratio_count = 0;


  /* varying extra precs */
  for (prec_val = PREC_START; prec_val <= PREC_END; prec_val++) {
    switch (prec_val) {
    case 0:
      eps_int = power(2, -BITS_S);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_single),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_single));
      prec = blas_prec_single;
      break;
    case 1:
      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;
      break;
    case 2:
    default:
      eps_int = power(2, -BITS_E);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
      prec = blas_prec_extra;
      break;
    }

    for (norm = NORM_START; norm <= NORM_END; norm++) {

      for (incx_val = INC_START; incx_val <= INC_END; incx_val++) {

	if (incx_val == 0)
	  continue;

	for (i = 0; i < ntests; i++) {

	  /* For the sake of speed, we throw out this case at random */
	  if (xrand(seed) >= test_prob)
	    continue;

	  saved_seed = *seed;

	  BLAS_csum_testgen(n, norm, x_gen, seed,
			    head_sum_true, tail_sum_true);

	  test_count++;

	  incx = incx_val;
	  incx *= 2;

	  xi = 0;
	  if (incx < 0)
	    xi = -(n - 1) * incx;

	  /* copy x_gen to x */
	  incx_gen = 1;
	  incx_gen *= 2;
	  for (j = 0, x_gen_i = 0; j < n;
	       j++, x_gen_i += incx_gen, xi += incx) {
	    x_elem[0] = x_gen[x_gen_i];
	    x_elem[1] = x_gen[x_gen_i + 1];
	    x[xi] = x_elem[0];
	    x[xi + 1] = x_elem[1];
	  }

	  FPU_FIX_STOP;
	  BLAS_csum_x(n, x, incx_val, sum, prec);
	  FPU_FIX_START;

	  test_BLAS_csum(n, sum, head_sum_true, tail_sum_true, x, incx_val,
			 eps_int, un_int, &ratio);

	  /* The !<= below causes NaN error to be detected.
	     Note that (NaN > thresh) is always false */
	  if (!(ratio <= thresh)) {
	    bad_ratio_count++;

	    if (debug == 3) {
	      printf("Seed = %d\n", saved_seed);
	      printf("n = %d\n", n);
	      printf("norm = %d\n", norm);
	      switch (prec) {
	      case blas_prec_single:
		printf("single ");
		break;
	      case blas_prec_double:
		printf("double ");
		break;
	      case blas_prec_indigenous:
		printf("indigenous ");
		break;
	      case blas_prec_extra:
		printf("extra ");
		break;
	      }

	      printf("incx = %d\n", incx);

	      /* print out the vector */
	      xi = (incx < 0) ? -(n - 1) * incx : 0;
	      printf(" [ ");
	      for (j = 0; j < n; j++, xi += incx) {
		printf("(%16.8e, %16.8e)", x[xi], x[xi + 1]);
	      }
	      printf("]\n");

	      printf("sum = ");
	      printf("(%16.8e, %16.8e)", sum[0], sum[1]);
	      printf("\n");
	      printf("ratio = %.4e\n", ratio);
	      printf("sum_true = ");
	      printf("([%24.16e  %24.16e], [%24.16e %24.16e])",
		     head_sum_true[0], tail_sum_true[0], head_sum_true[1],
		     tail_sum_true[1]);

	    }			/* end of if (debug == 3) */
	  }
	  /* end of if (ratio > thresh) */
	  if (ratio > ratio_max)
	    ratio_max = ratio;
	  if (ratio != 0.0 && ratio < ratio_min)
	    ratio_min = ratio;

	}			/* end of ntests loop */
      }				/* end of incx loop */
    }				/* end of norm loop */

  }				/* end of prec loop */

  FPU_FIX_STOP;

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

  blas_free(x);
  blas_free(x_gen);


}
void do_test_zsum_x(int n, int ntests, int *seed, double thresh,
		    int debug, float test_prob, double *min_ratio,
		    double *max_ratio, int *num_bad_ratio, int *num_tests)

/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on sum  
 *
 * Arguments
 * =========
 *
 * n         (input) int
 *           The size of vector being tested
 *
 * ntests    (input) int
 *           The number of tests to run for each set of attributes.
 *
 * seed      (input/output) int         
 *           The seed for the random number generator used in testgen().
 *
 * thresh    (input) double
 *           When the ratio returned from test() exceeds the specified
 *           threshold, the current size, r_true, r_comp, and ratio will be
 *           printed.  (Since ratio is supposed to be O(1), we can set thresh
 *           to ~10.)
 *
 * debug     (input) int
 *           If debug=3, print summary 
 *           If debug=2, print summary only if the number of bad ratios > 0
 *           If debug=1, print complete info if tests fail
 *           If debug=0, return max ratio
 *
 * min_ratio (output) double
 *           The minimum ratio
 * 
 * num_bad_ratio (output) int
 *               The number of tests fail; they are above the threshold.
 *
 * num_tests (output) int
 *           The number of tests is being performed.
 *
 */
{

  /* function name */
  const char fname[] = "BLAS_zsum_x";

  int i, j;
  int norm;
  int xi, incx_val, incx;
  double ratio_max, ratio_min;
  double ratio;
  int saved_seed;
  double eps_int;
  double un_int;

  int test_count;
  int bad_ratio_count;

  double *x;
  double *x_gen;
  double x_elem[2];

  double head_sum_true[2], tail_sum_true[2];

  double sum[2];

  int prec_val;
  enum blas_prec_type prec;
  int x_gen_i, incx_gen;

  FPU_FIX_DECL;


  /* test for bad arguments */
  if (n < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* if there is nothing to test, return all zero */
  if (n == 0 || ntests == 0) {
    *min_ratio = 0.0;
    *max_ratio = 0.0;
    *num_bad_ratio = 0;
    *num_tests = 0;
    return;
  }

  FPU_FIX_START;

  incx_gen = 1;
  incx_gen *= 2;

  /* get space for calculation */
  x = (double *) blas_malloc(n * 2 * sizeof(double) * 2);
  if (n * 2 > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* initialization */
  ratio_min = 1e308;
  ratio_max = 0.0;
  test_count = 0;
  bad_ratio_count = 0;


  /* varying extra precs */
  for (prec_val = PREC_START; prec_val <= PREC_END; prec_val++) {
    switch (prec_val) {
    case 0:
      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;
      break;
    case 1:
      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;
      break;
    case 2:
    default:
      eps_int = power(2, -BITS_E);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_extra),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_extra));
      prec = blas_prec_extra;
      break;
    }

    for (norm = NORM_START; norm <= NORM_END; norm++) {

      for (incx_val = INC_START; incx_val <= INC_END; incx_val++) {

	if (incx_val == 0)
	  continue;

	for (i = 0; i < ntests; i++) {

	  /* For the sake of speed, we throw out this case at random */
	  if (xrand(seed) >= test_prob)
	    continue;

	  saved_seed = *seed;

	  BLAS_zsum_testgen(n, norm, x_gen, seed,
			    head_sum_true, tail_sum_true);

	  test_count++;

	  incx = incx_val;
	  incx *= 2;

	  xi = 0;
	  if (incx < 0)
	    xi = -(n - 1) * incx;

	  /* copy x_gen to x */
	  incx_gen = 1;
	  incx_gen *= 2;
	  for (j = 0, x_gen_i = 0; j < n;
	       j++, x_gen_i += incx_gen, xi += incx) {
	    x_elem[0] = x_gen[x_gen_i];
	    x_elem[1] = x_gen[x_gen_i + 1];
	    x[xi] = x_elem[0];
	    x[xi + 1] = x_elem[1];
	  }

	  FPU_FIX_STOP;
	  BLAS_zsum_x(n, x, incx_val, sum, prec);
	  FPU_FIX_START;

	  test_BLAS_zsum(n, sum, head_sum_true, tail_sum_true, x, incx_val,
			 eps_int, un_int, &ratio);

	  /* The !<= below causes NaN error to be detected.
	     Note that (NaN > thresh) is always false */
	  if (!(ratio <= thresh)) {
	    bad_ratio_count++;

	    if (debug == 3) {
	      printf("Seed = %d\n", saved_seed);
	      printf("n = %d\n", n);
	      printf("norm = %d\n", norm);
	      switch (prec) {
	      case blas_prec_single:
		printf("single ");
		break;
	      case blas_prec_double:
		printf("double ");
		break;
	      case blas_prec_indigenous:
		printf("indigenous ");
		break;
	      case blas_prec_extra:
		printf("extra ");
		break;
	      }

	      printf("incx = %d\n", incx);

	      /* print out the vector */
	      xi = (incx < 0) ? -(n - 1) * incx : 0;
	      printf(" [ ");
	      for (j = 0; j < n; j++, xi += incx) {
		printf("(%24.16e, %24.16e)", x[xi], x[xi + 1]);
	      }
	      printf("]\n");

	      printf("sum = ");
	      printf("(%24.16e, %24.16e)", sum[0], sum[1]);
	      printf("\n");
	      printf("ratio = %.4e\n", ratio);
	      printf("sum_true = ");
	      printf("([%24.16e  %24.16e], [%24.16e %24.16e])",
		     head_sum_true[0], tail_sum_true[0], head_sum_true[1],
		     tail_sum_true[1]);

	    }			/* end of if (debug == 3) */
	  }
	  /* end of if (ratio > thresh) */
	  if (ratio > ratio_max)
	    ratio_max = ratio;
	  if (ratio != 0.0 && ratio < ratio_min)
	    ratio_min = ratio;

	}			/* end of ntests loop */
      }				/* end of incx loop */
    }				/* end of norm loop */

  }				/* end of prec loop */

  FPU_FIX_STOP;

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

  blas_free(x);
  blas_free(x_gen);


}

int main(int argc, char **argv)
{
  int nsizes, ntests, debug;
  double thresh, test_prob;
  double total_min_ratio, total_max_ratio;
  int total_bad_ratios;
  int seed, num_bad_ratio, num_tests;
  int total_tests, nr_failed_routines = 0, nr_routines = 0;
  double min_ratio, max_ratio;
  const char *base_routine = "sum";
  char *fname;
  int n;


  if (argc != 6) {
    printf("Usage:\n");
    printf("do_test_sum <nsizes> <ntests> <thresh> <debug> <test_prob>\n");
    printf("   <nsizes>: number of sizes to be run.\n");
    printf
      ("   <ntests>: the number of tests performed for each set of attributes\n");
    printf
      ("   <thresh>: to catch bad ratios if it is greater than <thresh>\n");
    printf("    <debug>: 0, 1, 2, or 3; \n");
    printf("        if 0, no printing \n");
    printf("        if 1, print error summary only if tests fail\n");
    printf("        if 2, print error summary for each n\n");
    printf("        if 3, print complete info each test fails \n");
    printf("<test_prob>: probability of preforming a given \n");
    printf("           test case: 0.0 does no tests, 1.0 does all tests\n");
    return -1;
  } else {
    nsizes = atoi(argv[1]);
    ntests = atoi(argv[2]);
    thresh = atof(argv[3]);
    debug = atoi(argv[4]);
    test_prob = atof(argv[5]);
  }

  seed = 1999;

  if (nsizes < 0 || ntests < 0 || debug < 0 || debug > 3)
    BLAS_error("Testing sum", 0, 0, NULL);

  printf("Testing %s...\n", base_routine);
  printf("INPUT: nsizes = %d, ntests = %d, thresh = %4.2f, debug = %d\n\n",
	 nsizes, ntests, thresh, debug);





  fname = "BLAS_ssum_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;

  for (n = 0; n <= nsizes; n++) {
    do_test_ssum_x(n, ntests, &seed, thresh, debug,
		   test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		   &num_tests);
    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("  n = %d:  ", n);
      printf("bad/total = %d/%d, min_ratio = %g, max_ratio = %g\n",
	     num_bad_ratio, num_tests, min_ratio, max_ratio);
    }
    total_tests += num_tests;
    total_bad_ratios += num_bad_ratio;
    if (total_min_ratio > min_ratio)
      total_min_ratio = min_ratio;
    if (total_max_ratio < max_ratio)
      total_max_ratio = max_ratio;
  }

  nr_routines++;
  if (total_bad_ratios == 0)
    printf("PASS> ");
  else {
    printf("FAIL> ");
    nr_failed_routines++;
  }
  printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n",
	 fname, total_bad_ratios, total_tests, max_ratio);
  printf("\n");

  fname = "BLAS_dsum_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;

  for (n = 0; n <= nsizes; n++) {
    do_test_dsum_x(n, ntests, &seed, thresh, debug,
		   test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		   &num_tests);
    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("  n = %d:  ", n);
      printf("bad/total = %d/%d, min_ratio = %g, max_ratio = %g\n",
	     num_bad_ratio, num_tests, min_ratio, max_ratio);
    }
    total_tests += num_tests;
    total_bad_ratios += num_bad_ratio;
    if (total_min_ratio > min_ratio)
      total_min_ratio = min_ratio;
    if (total_max_ratio < max_ratio)
      total_max_ratio = max_ratio;
  }

  nr_routines++;
  if (total_bad_ratios == 0)
    printf("PASS> ");
  else {
    printf("FAIL> ");
    nr_failed_routines++;
  }
  printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n",
	 fname, total_bad_ratios, total_tests, max_ratio);
  printf("\n");

  fname = "BLAS_csum_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;

  for (n = 0; n <= nsizes; n++) {
    do_test_csum_x(n, ntests, &seed, thresh, debug,
		   test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		   &num_tests);
    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("  n = %d:  ", n);
      printf("bad/total = %d/%d, min_ratio = %g, max_ratio = %g\n",
	     num_bad_ratio, num_tests, min_ratio, max_ratio);
    }
    total_tests += num_tests;
    total_bad_ratios += num_bad_ratio;
    if (total_min_ratio > min_ratio)
      total_min_ratio = min_ratio;
    if (total_max_ratio < max_ratio)
      total_max_ratio = max_ratio;
  }

  nr_routines++;
  if (total_bad_ratios == 0)
    printf("PASS> ");
  else {
    printf("FAIL> ");
    nr_failed_routines++;
  }
  printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n",
	 fname, total_bad_ratios, total_tests, max_ratio);
  printf("\n");

  fname = "BLAS_zsum_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;

  for (n = 0; n <= nsizes; n++) {
    do_test_zsum_x(n, ntests, &seed, thresh, debug,
		   test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		   &num_tests);
    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("  n = %d:  ", n);
      printf("bad/total = %d/%d, min_ratio = %g, max_ratio = %g\n",
	     num_bad_ratio, num_tests, min_ratio, max_ratio);
    }
    total_tests += num_tests;
    total_bad_ratios += num_bad_ratio;
    if (total_min_ratio > min_ratio)
      total_min_ratio = min_ratio;
    if (total_max_ratio < max_ratio)
      total_max_ratio = max_ratio;
  }

  nr_routines++;
  if (total_bad_ratios == 0)
    printf("PASS> ");
  else {
    printf("FAIL> ");
    nr_failed_routines++;
  }
  printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n",
	 fname, total_bad_ratios, total_tests, max_ratio);
  printf("\n");



  printf("\n");
  if (nr_failed_routines)
    printf("FAILED ");
  else
    printf("PASSED ");
  printf("%-10s: FAIL/TOTAL = %d/%d\n",
	 base_routine, nr_failed_routines, nr_routines);

  return 0;
}
