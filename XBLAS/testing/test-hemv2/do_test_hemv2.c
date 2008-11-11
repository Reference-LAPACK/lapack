#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

/* 0 -- 1 */
#define UPLO_START 0
#define UPLO_END   1

/* 0 -- 1 */
#define ORDER_START  0
#define ORDER_END    1

/* 0 -- 2 */
#define ALPHA_START  0
#define ALPHA_END    2

/* 0 -- 2 */
#define BETA_START   0
#define BETA_END     2

/* -1 -- 1 */
#define NORM_START   -1
#define NORM_END     1

/* 0 -- 2 */
#define LDA_START    0
#define LDA_END      2

/* 0 -- 2 */
#define PREC_START   0
#define PREC_END     2

/* 0 -- 1 */
#define RANDOMIZE_START 0
#define RANDOMIZE_END   1

/* -2 -- 2 (Stride) */
#define INCX_START -2
#define INCX_END 2

/* -2 -- 2 (Stride) */
#define INCY_START -2
#define INCY_END 2

#define NUM_DATA 7










void do_test_zhemv2_z_c
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_zhemv2_z_c";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;



  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  float *head_x;
  float *tail_x;
  double *y;
  double *a_vec;
  double *y_gen;
  float *head_x_gen;
  float *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n * n * sizeof(double) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;

      /* vary norm -- underflow, approx 1, overflow */
      for (norm = NORM_START; norm <= NORM_END; norm++) {

	/* number of tests */
	for (test_no = 0; test_no < ntests; test_no++) {

	  /* vary storage format */
	  for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary lda = n, n+1, 2*n */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		saved_seed = *seed;
		/* For the sake of speed, we throw out this case at random */
		if (xrand(seed) >= test_prob)
		  continue;

		alpha_flag = 0;
		switch (alpha_val) {
		case 0:
		  alpha[0] = alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		case 1:
		  alpha[0] = 1.0;
		  alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		}
		beta_flag = 0;
		switch (beta_val) {
		case 0:
		  beta[0] = beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		case 1:
		  beta[0] = 1.0;
		  beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		}

		/* finally we are here to generate the test case */
		BLAS_zhemv2_z_c_testgen(norm, order_type,
					uplo_type, n, &alpha, alpha_flag,
					&beta, beta_flag, a, lda, head_x_gen,
					tail_x_gen, y_gen, seed, head_r_true,
					tail_r_true);
		test_count++;

		/* vary incx = -2, -1, 1, 2 */
		for (incx_val = INCX_START; incx_val <= INCX_END; incx_val++) {

		  incx = incx_val;
		  if (0 == incx)
		    continue;

		  /* vary incy = -2, -1, 1, 2 */
		  for (incy_val = INCY_START; incy_val <= INCY_END;
		       incy_val++) {

		    incy = incy_val;
		    if (0 == incy)
		      continue;

		    /* copy generated vector with appropriate incs. */
		    zcopy_vector(y_gen, n, 1, y, incy);
		    ccopy_vector(head_x_gen, n, 1, head_x, incx);
		    ccopy_vector(tail_x_gen, n, 1, tail_x, incx);

		    /* call hemv2 routines to be tested */
		    FPU_FIX_STOP;
		    BLAS_zhemv2_z_c(order_type,
				    uplo_type, n, alpha, a, lda, head_x,
				    tail_x, incx, beta, y, incy);
		    FPU_FIX_START;

		    /* now compute the ratio using test_BLAS_xdot */
		    /* copy a row from A, use x, run dot test */

		    incyi = incy;
		    incyi *= 2;
		    yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		    for (i = 0, yi = yi0, ri = 0;
			 i < n; i++, yi += incyi, ri += incri) {
		      zhe_copy_row(order_type, uplo_type, blas_left_side,
				   n, a, lda, a_vec, i);

		      /* just use the x vector - it was unchanged (in theory) */
		      rin[0] = y_gen[i];
		      rin[1] = y_gen[i + 1];
		      rout[0] = y[yi];
		      rout[1] = y[yi + 1];
		      head_r_true_elem[0] = head_r_true[ri];
		      head_r_true_elem[1] = head_r_true[ri + 1];
		      tail_r_true_elem[0] = tail_r_true[ri];
		      tail_r_true_elem[1] = tail_r_true[ri + 1];

		      test_BLAS_zdot2_z_c(n, blas_no_conj, alpha, beta,
					  rin, rout, head_r_true_elem,
					  tail_r_true_elem, a_vec, 1, head_x,
					  tail_x, incx, eps_int, un_int,
					  &ratios[i]);

		      /* take the max ratio */
		      if (i == 0) {
			ratio = ratios[0];

			/* The !<= below causes NaN errors to be included.
			 * Note that (NaN > 0) is false */
		      } else if (!(ratios[i] <= ratio)) {
			ratio = ratios[i];
		      }

		    }		/* end of dot-test loop */


		    /* The !<= below causes NaN errors to be included.
		     * Note that (NaN > 0) is false */
		    if (!(ratio <= thresh)) {

		      if (debug == 3) {
			printf("\n\t\tTest # %d\n", test_count);
			printf("y type : z, a type : z, x type : c\n");
			printf("Seed = %d\t", saved_seed);
			printf("n %d\n", n);
			printf("LDA %d  INCX %d  INCY %d\n", lda, incx, incx);

			if (order_type == blas_rowmajor)
			  printf("row ");
			else
			  printf("col ");

			if (uplo_type == blas_upper)
			  printf("upper ");
			else
			  printf("lower ");

			printf("NORM %d, ALPHA %d, BETA %d\n",
			       norm, alpha_val, beta_val);

			/* print out info */
			printf("alpha = ");
			printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			printf("   ");
			printf("beta = ");
			printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			printf("\n");

			printf("a\n");
			zhe_print_matrix(a, n, lda, order_type, uplo_type);
			cprint_vector(head_x, n, incx, "head_x");
			cprint_vector(tail_x, n, incx, "tail_x");
			zprint_vector(y_gen, n, incy, "y_gen");
			zprint_vector(y, n, incy, "y");
			zprint_vector(head_r_true, n, 1, "head_r_true");
			dprint_vector(ratios, n, 1, "ratios");
			printf("ratio = %g\n", ratio);
		      }
		      bad_ratio_count++;
		      if (bad_ratio_count >= MAX_BAD_TESTS) {
			printf("\ntoo many failures, exiting....");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		      if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			printf("\nFlagrant ratio error, exiting...");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		    }

		    if (!(ratio <= ratio_max))
		      ratio_max = ratio;
		    if (ratio != 0.0 && !(ratio >= ratio_min))
		      ratio_min = ratio;

		  }		/* end of incy loop */

		}		/* end of incx loop */

	      }			/* end of lda loop */

	    }			/* end of uplo loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zhemv2_c_z
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_zhemv2_c_z";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;



  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  double alpha[2];
  double beta[2];
  float *a;
  double *head_x;
  double *tail_x;
  double *y;
  float *a_vec;
  double *y_gen;
  double *head_x_gen;
  double *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n * n * sizeof(float) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;

      /* vary norm -- underflow, approx 1, overflow */
      for (norm = NORM_START; norm <= NORM_END; norm++) {

	/* number of tests */
	for (test_no = 0; test_no < ntests; test_no++) {

	  /* vary storage format */
	  for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary lda = n, n+1, 2*n */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		saved_seed = *seed;
		/* For the sake of speed, we throw out this case at random */
		if (xrand(seed) >= test_prob)
		  continue;

		alpha_flag = 0;
		switch (alpha_val) {
		case 0:
		  alpha[0] = alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		case 1:
		  alpha[0] = 1.0;
		  alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		}
		beta_flag = 0;
		switch (beta_val) {
		case 0:
		  beta[0] = beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		case 1:
		  beta[0] = 1.0;
		  beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		}

		/* finally we are here to generate the test case */
		BLAS_zhemv2_c_z_testgen(norm, order_type,
					uplo_type, n, &alpha, alpha_flag,
					&beta, beta_flag, a, lda, head_x_gen,
					tail_x_gen, y_gen, seed, head_r_true,
					tail_r_true);
		test_count++;

		/* vary incx = -2, -1, 1, 2 */
		for (incx_val = INCX_START; incx_val <= INCX_END; incx_val++) {

		  incx = incx_val;
		  if (0 == incx)
		    continue;

		  /* vary incy = -2, -1, 1, 2 */
		  for (incy_val = INCY_START; incy_val <= INCY_END;
		       incy_val++) {

		    incy = incy_val;
		    if (0 == incy)
		      continue;

		    /* copy generated vector with appropriate incs. */
		    zcopy_vector(y_gen, n, 1, y, incy);
		    zcopy_vector(head_x_gen, n, 1, head_x, incx);
		    zcopy_vector(tail_x_gen, n, 1, tail_x, incx);

		    /* call hemv2 routines to be tested */
		    FPU_FIX_STOP;
		    BLAS_zhemv2_c_z(order_type,
				    uplo_type, n, alpha, a, lda, head_x,
				    tail_x, incx, beta, y, incy);
		    FPU_FIX_START;

		    /* now compute the ratio using test_BLAS_xdot */
		    /* copy a row from A, use x, run dot test */

		    incyi = incy;
		    incyi *= 2;
		    yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		    for (i = 0, yi = yi0, ri = 0;
			 i < n; i++, yi += incyi, ri += incri) {
		      che_copy_row(order_type, uplo_type, blas_left_side,
				   n, a, lda, a_vec, i);

		      /* just use the x vector - it was unchanged (in theory) */
		      rin[0] = y_gen[i];
		      rin[1] = y_gen[i + 1];
		      rout[0] = y[yi];
		      rout[1] = y[yi + 1];
		      head_r_true_elem[0] = head_r_true[ri];
		      head_r_true_elem[1] = head_r_true[ri + 1];
		      tail_r_true_elem[0] = tail_r_true[ri];
		      tail_r_true_elem[1] = tail_r_true[ri + 1];

		      test_BLAS_zdot2_c_z(n, blas_no_conj, alpha, beta,
					  rin, rout, head_r_true_elem,
					  tail_r_true_elem, a_vec, 1, head_x,
					  tail_x, incx, eps_int, un_int,
					  &ratios[i]);

		      /* take the max ratio */
		      if (i == 0) {
			ratio = ratios[0];

			/* The !<= below causes NaN errors to be included.
			 * Note that (NaN > 0) is false */
		      } else if (!(ratios[i] <= ratio)) {
			ratio = ratios[i];
		      }

		    }		/* end of dot-test loop */


		    /* The !<= below causes NaN errors to be included.
		     * Note that (NaN > 0) is false */
		    if (!(ratio <= thresh)) {

		      if (debug == 3) {
			printf("\n\t\tTest # %d\n", test_count);
			printf("y type : z, a type : c, x type : z\n");
			printf("Seed = %d\t", saved_seed);
			printf("n %d\n", n);
			printf("LDA %d  INCX %d  INCY %d\n", lda, incx, incx);

			if (order_type == blas_rowmajor)
			  printf("row ");
			else
			  printf("col ");

			if (uplo_type == blas_upper)
			  printf("upper ");
			else
			  printf("lower ");

			printf("NORM %d, ALPHA %d, BETA %d\n",
			       norm, alpha_val, beta_val);

			/* print out info */
			printf("alpha = ");
			printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			printf("   ");
			printf("beta = ");
			printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			printf("\n");

			printf("a\n");
			che_print_matrix(a, n, lda, order_type, uplo_type);
			zprint_vector(head_x, n, incx, "head_x");
			zprint_vector(tail_x, n, incx, "tail_x");
			zprint_vector(y_gen, n, incy, "y_gen");
			zprint_vector(y, n, incy, "y");
			zprint_vector(head_r_true, n, 1, "head_r_true");
			dprint_vector(ratios, n, 1, "ratios");
			printf("ratio = %g\n", ratio);
		      }
		      bad_ratio_count++;
		      if (bad_ratio_count >= MAX_BAD_TESTS) {
			printf("\ntoo many failures, exiting....");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		      if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			printf("\nFlagrant ratio error, exiting...");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		    }

		    if (!(ratio <= ratio_max))
		      ratio_max = ratio;
		    if (ratio != 0.0 && !(ratio >= ratio_min))
		      ratio_min = ratio;

		  }		/* end of incy loop */

		}		/* end of incx loop */

	      }			/* end of lda loop */

	    }			/* end of uplo loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zhemv2_c_c
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_zhemv2_c_c";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;



  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  double alpha[2];
  double beta[2];
  float *a;
  float *head_x;
  float *tail_x;
  double *y;
  float *a_vec;
  double *y_gen;
  float *head_x_gen;
  float *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n * n * sizeof(float) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;

      /* vary norm -- underflow, approx 1, overflow */
      for (norm = NORM_START; norm <= NORM_END; norm++) {

	/* number of tests */
	for (test_no = 0; test_no < ntests; test_no++) {

	  /* vary storage format */
	  for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary lda = n, n+1, 2*n */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		saved_seed = *seed;
		/* For the sake of speed, we throw out this case at random */
		if (xrand(seed) >= test_prob)
		  continue;

		alpha_flag = 0;
		switch (alpha_val) {
		case 0:
		  alpha[0] = alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		case 1:
		  alpha[0] = 1.0;
		  alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		}
		beta_flag = 0;
		switch (beta_val) {
		case 0:
		  beta[0] = beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		case 1:
		  beta[0] = 1.0;
		  beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		}

		/* finally we are here to generate the test case */
		BLAS_zhemv2_c_c_testgen(norm, order_type,
					uplo_type, n, &alpha, alpha_flag,
					&beta, beta_flag, a, lda, head_x_gen,
					tail_x_gen, y_gen, seed, head_r_true,
					tail_r_true);
		test_count++;

		/* vary incx = -2, -1, 1, 2 */
		for (incx_val = INCX_START; incx_val <= INCX_END; incx_val++) {

		  incx = incx_val;
		  if (0 == incx)
		    continue;

		  /* vary incy = -2, -1, 1, 2 */
		  for (incy_val = INCY_START; incy_val <= INCY_END;
		       incy_val++) {

		    incy = incy_val;
		    if (0 == incy)
		      continue;

		    /* copy generated vector with appropriate incs. */
		    zcopy_vector(y_gen, n, 1, y, incy);
		    ccopy_vector(head_x_gen, n, 1, head_x, incx);
		    ccopy_vector(tail_x_gen, n, 1, tail_x, incx);

		    /* call hemv2 routines to be tested */
		    FPU_FIX_STOP;
		    BLAS_zhemv2_c_c(order_type,
				    uplo_type, n, alpha, a, lda, head_x,
				    tail_x, incx, beta, y, incy);
		    FPU_FIX_START;

		    /* now compute the ratio using test_BLAS_xdot */
		    /* copy a row from A, use x, run dot test */

		    incyi = incy;
		    incyi *= 2;
		    yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		    for (i = 0, yi = yi0, ri = 0;
			 i < n; i++, yi += incyi, ri += incri) {
		      che_copy_row(order_type, uplo_type, blas_left_side,
				   n, a, lda, a_vec, i);

		      /* just use the x vector - it was unchanged (in theory) */
		      rin[0] = y_gen[i];
		      rin[1] = y_gen[i + 1];
		      rout[0] = y[yi];
		      rout[1] = y[yi + 1];
		      head_r_true_elem[0] = head_r_true[ri];
		      head_r_true_elem[1] = head_r_true[ri + 1];
		      tail_r_true_elem[0] = tail_r_true[ri];
		      tail_r_true_elem[1] = tail_r_true[ri + 1];

		      test_BLAS_zdot2_c_c(n, blas_no_conj, alpha, beta,
					  rin, rout, head_r_true_elem,
					  tail_r_true_elem, a_vec, 1, head_x,
					  tail_x, incx, eps_int, un_int,
					  &ratios[i]);

		      /* take the max ratio */
		      if (i == 0) {
			ratio = ratios[0];

			/* The !<= below causes NaN errors to be included.
			 * Note that (NaN > 0) is false */
		      } else if (!(ratios[i] <= ratio)) {
			ratio = ratios[i];
		      }

		    }		/* end of dot-test loop */


		    /* The !<= below causes NaN errors to be included.
		     * Note that (NaN > 0) is false */
		    if (!(ratio <= thresh)) {

		      if (debug == 3) {
			printf("\n\t\tTest # %d\n", test_count);
			printf("y type : z, a type : c, x type : c\n");
			printf("Seed = %d\t", saved_seed);
			printf("n %d\n", n);
			printf("LDA %d  INCX %d  INCY %d\n", lda, incx, incx);

			if (order_type == blas_rowmajor)
			  printf("row ");
			else
			  printf("col ");

			if (uplo_type == blas_upper)
			  printf("upper ");
			else
			  printf("lower ");

			printf("NORM %d, ALPHA %d, BETA %d\n",
			       norm, alpha_val, beta_val);

			/* print out info */
			printf("alpha = ");
			printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			printf("   ");
			printf("beta = ");
			printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			printf("\n");

			printf("a\n");
			che_print_matrix(a, n, lda, order_type, uplo_type);
			cprint_vector(head_x, n, incx, "head_x");
			cprint_vector(tail_x, n, incx, "tail_x");
			zprint_vector(y_gen, n, incy, "y_gen");
			zprint_vector(y, n, incy, "y");
			zprint_vector(head_r_true, n, 1, "head_r_true");
			dprint_vector(ratios, n, 1, "ratios");
			printf("ratio = %g\n", ratio);
		      }
		      bad_ratio_count++;
		      if (bad_ratio_count >= MAX_BAD_TESTS) {
			printf("\ntoo many failures, exiting....");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		      if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			printf("\nFlagrant ratio error, exiting...");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		    }

		    if (!(ratio <= ratio_max))
		      ratio_max = ratio;
		    if (ratio != 0.0 && !(ratio >= ratio_min))
		      ratio_min = ratio;

		  }		/* end of incy loop */

		}		/* end of incx loop */

	      }			/* end of lda loop */

	    }			/* end of uplo loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_chemv2_c_s
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_chemv2_c_s";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;



  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *head_x;
  float *tail_x;
  float *y;
  float *a_vec;
  float *y_gen;
  float *head_x_gen;
  float *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n * n * sizeof(float) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


      eps_int = power(2, -BITS_S);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_single),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_single));
      prec = blas_prec_single;

      /* vary norm -- underflow, approx 1, overflow */
      for (norm = NORM_START; norm <= NORM_END; norm++) {

	/* number of tests */
	for (test_no = 0; test_no < ntests; test_no++) {

	  /* vary storage format */
	  for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary lda = n, n+1, 2*n */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		saved_seed = *seed;
		/* For the sake of speed, we throw out this case at random */
		if (xrand(seed) >= test_prob)
		  continue;

		alpha_flag = 0;
		switch (alpha_val) {
		case 0:
		  alpha[0] = alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		case 1:
		  alpha[0] = 1.0;
		  alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		}
		beta_flag = 0;
		switch (beta_val) {
		case 0:
		  beta[0] = beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		case 1:
		  beta[0] = 1.0;
		  beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		}

		/* finally we are here to generate the test case */
		BLAS_chemv2_c_s_testgen(norm, order_type,
					uplo_type, n, &alpha, alpha_flag,
					&beta, beta_flag, a, lda, head_x_gen,
					tail_x_gen, y_gen, seed, head_r_true,
					tail_r_true);
		test_count++;

		/* vary incx = -2, -1, 1, 2 */
		for (incx_val = INCX_START; incx_val <= INCX_END; incx_val++) {

		  incx = incx_val;
		  if (0 == incx)
		    continue;

		  /* vary incy = -2, -1, 1, 2 */
		  for (incy_val = INCY_START; incy_val <= INCY_END;
		       incy_val++) {

		    incy = incy_val;
		    if (0 == incy)
		      continue;

		    /* copy generated vector with appropriate incs. */
		    ccopy_vector(y_gen, n, 1, y, incy);
		    scopy_vector(head_x_gen, n, 1, head_x, incx);
		    scopy_vector(tail_x_gen, n, 1, tail_x, incx);

		    /* call hemv2 routines to be tested */
		    FPU_FIX_STOP;
		    BLAS_chemv2_c_s(order_type,
				    uplo_type, n, alpha, a, lda, head_x,
				    tail_x, incx, beta, y, incy);
		    FPU_FIX_START;

		    /* now compute the ratio using test_BLAS_xdot */
		    /* copy a row from A, use x, run dot test */

		    incyi = incy;
		    incyi *= 2;
		    yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		    for (i = 0, yi = yi0, ri = 0;
			 i < n; i++, yi += incyi, ri += incri) {
		      che_copy_row(order_type, uplo_type, blas_left_side,
				   n, a, lda, a_vec, i);

		      /* just use the x vector - it was unchanged (in theory) */
		      rin[0] = y_gen[i];
		      rin[1] = y_gen[i + 1];
		      rout[0] = y[yi];
		      rout[1] = y[yi + 1];
		      head_r_true_elem[0] = head_r_true[ri];
		      head_r_true_elem[1] = head_r_true[ri + 1];
		      tail_r_true_elem[0] = tail_r_true[ri];
		      tail_r_true_elem[1] = tail_r_true[ri + 1];

		      test_BLAS_cdot2_c_s(n, blas_no_conj, alpha, beta,
					  rin, rout, head_r_true_elem,
					  tail_r_true_elem, a_vec, 1, head_x,
					  tail_x, incx, eps_int, un_int,
					  &ratios[i]);

		      /* take the max ratio */
		      if (i == 0) {
			ratio = ratios[0];

			/* The !<= below causes NaN errors to be included.
			 * Note that (NaN > 0) is false */
		      } else if (!(ratios[i] <= ratio)) {
			ratio = ratios[i];
		      }

		    }		/* end of dot-test loop */


		    /* The !<= below causes NaN errors to be included.
		     * Note that (NaN > 0) is false */
		    if (!(ratio <= thresh)) {

		      if (debug == 3) {
			printf("\n\t\tTest # %d\n", test_count);
			printf("y type : c, a type : c, x type : s\n");
			printf("Seed = %d\t", saved_seed);
			printf("n %d\n", n);
			printf("LDA %d  INCX %d  INCY %d\n", lda, incx, incx);

			if (order_type == blas_rowmajor)
			  printf("row ");
			else
			  printf("col ");

			if (uplo_type == blas_upper)
			  printf("upper ");
			else
			  printf("lower ");

			printf("NORM %d, ALPHA %d, BETA %d\n",
			       norm, alpha_val, beta_val);

			/* print out info */
			printf("alpha = ");
			printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			printf("   ");
			printf("beta = ");
			printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			printf("\n");

			printf("a\n");
			che_print_matrix(a, n, lda, order_type, uplo_type);
			sprint_vector(head_x, n, incx, "head_x");
			sprint_vector(tail_x, n, incx, "tail_x");
			cprint_vector(y_gen, n, incy, "y_gen");
			cprint_vector(y, n, incy, "y");
			zprint_vector(head_r_true, n, 1, "head_r_true");
			dprint_vector(ratios, n, 1, "ratios");
			printf("ratio = %g\n", ratio);
		      }
		      bad_ratio_count++;
		      if (bad_ratio_count >= MAX_BAD_TESTS) {
			printf("\ntoo many failures, exiting....");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		      if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			printf("\nFlagrant ratio error, exiting...");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		    }

		    if (!(ratio <= ratio_max))
		      ratio_max = ratio;
		    if (ratio != 0.0 && !(ratio >= ratio_min))
		      ratio_min = ratio;

		  }		/* end of incy loop */

		}		/* end of incx loop */

	      }			/* end of lda loop */

	    }			/* end of uplo loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zhemv2_z_d
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_zhemv2_z_d";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;



  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *head_x;
  double *tail_x;
  double *y;
  double *a_vec;
  double *y_gen;
  double *head_x_gen;
  double *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n * n * sizeof(double) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (double *) blas_malloc(2 * n * sizeof(double));
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (double *) blas_malloc(2 * n * sizeof(double));
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


      eps_int = power(2, -BITS_D);
      un_int = pow((double) BLAS_fpinfo_x(blas_base, blas_prec_double),
		   (double) BLAS_fpinfo_x(blas_emin, blas_prec_double));
      prec = blas_prec_double;

      /* vary norm -- underflow, approx 1, overflow */
      for (norm = NORM_START; norm <= NORM_END; norm++) {

	/* number of tests */
	for (test_no = 0; test_no < ntests; test_no++) {

	  /* vary storage format */
	  for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary lda = n, n+1, 2*n */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		saved_seed = *seed;
		/* For the sake of speed, we throw out this case at random */
		if (xrand(seed) >= test_prob)
		  continue;

		alpha_flag = 0;
		switch (alpha_val) {
		case 0:
		  alpha[0] = alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		case 1:
		  alpha[0] = 1.0;
		  alpha[1] = 0.0;
		  alpha_flag = 1;
		  break;
		}
		beta_flag = 0;
		switch (beta_val) {
		case 0:
		  beta[0] = beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		case 1:
		  beta[0] = 1.0;
		  beta[1] = 0.0;
		  beta_flag = 1;
		  break;
		}

		/* finally we are here to generate the test case */
		BLAS_zhemv2_z_d_testgen(norm, order_type,
					uplo_type, n, &alpha, alpha_flag,
					&beta, beta_flag, a, lda, head_x_gen,
					tail_x_gen, y_gen, seed, head_r_true,
					tail_r_true);
		test_count++;

		/* vary incx = -2, -1, 1, 2 */
		for (incx_val = INCX_START; incx_val <= INCX_END; incx_val++) {

		  incx = incx_val;
		  if (0 == incx)
		    continue;

		  /* vary incy = -2, -1, 1, 2 */
		  for (incy_val = INCY_START; incy_val <= INCY_END;
		       incy_val++) {

		    incy = incy_val;
		    if (0 == incy)
		      continue;

		    /* copy generated vector with appropriate incs. */
		    zcopy_vector(y_gen, n, 1, y, incy);
		    dcopy_vector(head_x_gen, n, 1, head_x, incx);
		    dcopy_vector(tail_x_gen, n, 1, tail_x, incx);

		    /* call hemv2 routines to be tested */
		    FPU_FIX_STOP;
		    BLAS_zhemv2_z_d(order_type,
				    uplo_type, n, alpha, a, lda, head_x,
				    tail_x, incx, beta, y, incy);
		    FPU_FIX_START;

		    /* now compute the ratio using test_BLAS_xdot */
		    /* copy a row from A, use x, run dot test */

		    incyi = incy;
		    incyi *= 2;
		    yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		    for (i = 0, yi = yi0, ri = 0;
			 i < n; i++, yi += incyi, ri += incri) {
		      zhe_copy_row(order_type, uplo_type, blas_left_side,
				   n, a, lda, a_vec, i);

		      /* just use the x vector - it was unchanged (in theory) */
		      rin[0] = y_gen[i];
		      rin[1] = y_gen[i + 1];
		      rout[0] = y[yi];
		      rout[1] = y[yi + 1];
		      head_r_true_elem[0] = head_r_true[ri];
		      head_r_true_elem[1] = head_r_true[ri + 1];
		      tail_r_true_elem[0] = tail_r_true[ri];
		      tail_r_true_elem[1] = tail_r_true[ri + 1];

		      test_BLAS_zdot2_z_d(n, blas_no_conj, alpha, beta,
					  rin, rout, head_r_true_elem,
					  tail_r_true_elem, a_vec, 1, head_x,
					  tail_x, incx, eps_int, un_int,
					  &ratios[i]);

		      /* take the max ratio */
		      if (i == 0) {
			ratio = ratios[0];

			/* The !<= below causes NaN errors to be included.
			 * Note that (NaN > 0) is false */
		      } else if (!(ratios[i] <= ratio)) {
			ratio = ratios[i];
		      }

		    }		/* end of dot-test loop */


		    /* The !<= below causes NaN errors to be included.
		     * Note that (NaN > 0) is false */
		    if (!(ratio <= thresh)) {

		      if (debug == 3) {
			printf("\n\t\tTest # %d\n", test_count);
			printf("y type : z, a type : z, x type : d\n");
			printf("Seed = %d\t", saved_seed);
			printf("n %d\n", n);
			printf("LDA %d  INCX %d  INCY %d\n", lda, incx, incx);

			if (order_type == blas_rowmajor)
			  printf("row ");
			else
			  printf("col ");

			if (uplo_type == blas_upper)
			  printf("upper ");
			else
			  printf("lower ");

			printf("NORM %d, ALPHA %d, BETA %d\n",
			       norm, alpha_val, beta_val);

			/* print out info */
			printf("alpha = ");
			printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			printf("   ");
			printf("beta = ");
			printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			printf("\n");

			printf("a\n");
			zhe_print_matrix(a, n, lda, order_type, uplo_type);
			dprint_vector(head_x, n, incx, "head_x");
			dprint_vector(tail_x, n, incx, "tail_x");
			zprint_vector(y_gen, n, incy, "y_gen");
			zprint_vector(y, n, incy, "y");
			zprint_vector(head_r_true, n, 1, "head_r_true");
			dprint_vector(ratios, n, 1, "ratios");
			printf("ratio = %g\n", ratio);
		      }
		      bad_ratio_count++;
		      if (bad_ratio_count >= MAX_BAD_TESTS) {
			printf("\ntoo many failures, exiting....");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		      if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			printf("\nFlagrant ratio error, exiting...");
			printf("\nTesting and compilation");
			printf(" are incomplete\n\n");
			goto end;
		      }
		    }

		    if (!(ratio <= ratio_max))
		      ratio_max = ratio;
		    if (ratio != 0.0 && !(ratio >= ratio_min))
		      ratio_min = ratio;

		  }		/* end of incy loop */

		}		/* end of incx loop */

	      }			/* end of lda loop */

	    }			/* end of uplo loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_chemv2_x
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_chemv2_x";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;

  int prec_val;

  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *head_x;
  float *tail_x;
  float *y;
  float *a_vec;
  float *y_gen;
  float *head_x_gen;
  float *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n * n * sizeof(float) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


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

	/* vary norm -- underflow, approx 1, overflow */
	for (norm = NORM_START; norm <= NORM_END; norm++) {

	  /* number of tests */
	  for (test_no = 0; test_no < ntests; test_no++) {

	    /* vary storage format */
	    for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary lda = n, n+1, 2*n */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		  saved_seed = *seed;
		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  alpha_flag = 0;
		  switch (alpha_val) {
		  case 0:
		    alpha[0] = alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  case 1:
		    alpha[0] = 1.0;
		    alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  }
		  beta_flag = 0;
		  switch (beta_val) {
		  case 0:
		    beta[0] = beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  case 1:
		    beta[0] = 1.0;
		    beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  }

		  /* finally we are here to generate the test case */
		  BLAS_chemv2_testgen(norm, order_type,
				      uplo_type, n, &alpha, alpha_flag, &beta,
				      beta_flag, a, lda, head_x_gen,
				      tail_x_gen, y_gen, seed, head_r_true,
				      tail_r_true);
		  test_count++;

		  /* vary incx = -2, -1, 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    /* vary incy = -2, -1, 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      /* copy generated vector with appropriate incs. */
		      ccopy_vector(y_gen, n, 1, y, incy);
		      ccopy_vector(head_x_gen, n, 1, head_x, incx);
		      ccopy_vector(tail_x_gen, n, 1, tail_x, incx);

		      /* call hemv2 routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_chemv2_x(order_type,
				    uplo_type, n, alpha, a, lda, head_x,
				    tail_x, incx, beta, y, incy, prec);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run dot test */

		      incyi = incy;
		      incyi *= 2;
		      yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		      for (i = 0, yi = yi0, ri = 0;
			   i < n; i++, yi += incyi, ri += incri) {
			che_copy_row(order_type, uplo_type, blas_left_side,
				     n, a, lda, a_vec, i);

			/* just use the x vector - it was unchanged (in theory) */
			rin[0] = y_gen[i];
			rin[1] = y_gen[i + 1];
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_cdot2(n, blas_no_conj, alpha, beta,
					rin, rout, head_r_true_elem,
					tail_r_true_elem, a_vec, 1, head_x,
					tail_x, incx, eps_int, un_int,
					&ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];

			  /* The !<= below causes NaN errors to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}

		      }		/* end of dot-test loop */


		      /* The !<= below causes NaN errors to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : c, a type : c, x type : c\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d\n", n);
			  printf("LDA %d  INCX %d  INCY %d\n", lda, incx,
				 incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  if (uplo_type == blas_upper)
			    printf("upper ");
			  else
			    printf("lower ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			  printf("\n");

			  printf("a\n");
			  che_print_matrix(a, n, lda, order_type, uplo_type);
			  cprint_vector(head_x, n, incx, "head_x");
			  cprint_vector(tail_x, n, incx, "tail_x");
			  cprint_vector(y_gen, n, incy, "y_gen");
			  cprint_vector(y, n, incy, "y");
			  zprint_vector(head_r_true, n, 1, "head_r_true");
			  dprint_vector(ratios, n, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			}
			bad_ratio_count++;
			if (bad_ratio_count >= MAX_BAD_TESTS) {
			  printf("\ntoo many failures, exiting....");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
			if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			  printf("\nFlagrant ratio error, exiting...");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
		      }

		      if (!(ratio <= ratio_max))
			ratio_max = ratio;
		      if (ratio != 0.0 && !(ratio >= ratio_min))
			ratio_min = ratio;

		    }		/* end of incy loop */

		  }		/* end of incx loop */

		}		/* end of lda loop */

	      }			/* end of uplo loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zhemv2_x
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_zhemv2_x";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;

  int prec_val;

  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *head_x;
  double *tail_x;
  double *y;
  double *a_vec;
  double *y_gen;
  double *head_x_gen;
  double *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n * n * sizeof(double) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


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

	/* vary norm -- underflow, approx 1, overflow */
	for (norm = NORM_START; norm <= NORM_END; norm++) {

	  /* number of tests */
	  for (test_no = 0; test_no < ntests; test_no++) {

	    /* vary storage format */
	    for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary lda = n, n+1, 2*n */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		  saved_seed = *seed;
		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  alpha_flag = 0;
		  switch (alpha_val) {
		  case 0:
		    alpha[0] = alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  case 1:
		    alpha[0] = 1.0;
		    alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  }
		  beta_flag = 0;
		  switch (beta_val) {
		  case 0:
		    beta[0] = beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  case 1:
		    beta[0] = 1.0;
		    beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  }

		  /* finally we are here to generate the test case */
		  BLAS_zhemv2_testgen(norm, order_type,
				      uplo_type, n, &alpha, alpha_flag, &beta,
				      beta_flag, a, lda, head_x_gen,
				      tail_x_gen, y_gen, seed, head_r_true,
				      tail_r_true);
		  test_count++;

		  /* vary incx = -2, -1, 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    /* vary incy = -2, -1, 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      /* copy generated vector with appropriate incs. */
		      zcopy_vector(y_gen, n, 1, y, incy);
		      zcopy_vector(head_x_gen, n, 1, head_x, incx);
		      zcopy_vector(tail_x_gen, n, 1, tail_x, incx);

		      /* call hemv2 routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zhemv2_x(order_type,
				    uplo_type, n, alpha, a, lda, head_x,
				    tail_x, incx, beta, y, incy, prec);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run dot test */

		      incyi = incy;
		      incyi *= 2;
		      yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		      for (i = 0, yi = yi0, ri = 0;
			   i < n; i++, yi += incyi, ri += incri) {
			zhe_copy_row(order_type, uplo_type, blas_left_side,
				     n, a, lda, a_vec, i);

			/* just use the x vector - it was unchanged (in theory) */
			rin[0] = y_gen[i];
			rin[1] = y_gen[i + 1];
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot2(n, blas_no_conj, alpha, beta,
					rin, rout, head_r_true_elem,
					tail_r_true_elem, a_vec, 1, head_x,
					tail_x, incx, eps_int, un_int,
					&ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];

			  /* The !<= below causes NaN errors to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}

		      }		/* end of dot-test loop */


		      /* The !<= below causes NaN errors to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : z, x type : z\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d\n", n);
			  printf("LDA %d  INCX %d  INCY %d\n", lda, incx,
				 incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  if (uplo_type == blas_upper)
			    printf("upper ");
			  else
			    printf("lower ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");

			  printf("a\n");
			  zhe_print_matrix(a, n, lda, order_type, uplo_type);
			  zprint_vector(head_x, n, incx, "head_x");
			  zprint_vector(tail_x, n, incx, "tail_x");
			  zprint_vector(y_gen, n, incy, "y_gen");
			  zprint_vector(y, n, incy, "y");
			  zprint_vector(head_r_true, n, 1, "head_r_true");
			  dprint_vector(ratios, n, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			}
			bad_ratio_count++;
			if (bad_ratio_count >= MAX_BAD_TESTS) {
			  printf("\ntoo many failures, exiting....");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
			if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			  printf("\nFlagrant ratio error, exiting...");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
		      }

		      if (!(ratio <= ratio_max))
			ratio_max = ratio;
		      if (ratio != 0.0 && !(ratio >= ratio_min))
			ratio_min = ratio;

		    }		/* end of incy loop */

		  }		/* end of incx loop */

		}		/* end of lda loop */

	      }			/* end of uplo loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zhemv2_z_c_x
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_zhemv2_z_c_x";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;

  int prec_val;

  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  float *head_x;
  float *tail_x;
  double *y;
  double *a_vec;
  double *y_gen;
  float *head_x_gen;
  float *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n * n * sizeof(double) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


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

	/* vary norm -- underflow, approx 1, overflow */
	for (norm = NORM_START; norm <= NORM_END; norm++) {

	  /* number of tests */
	  for (test_no = 0; test_no < ntests; test_no++) {

	    /* vary storage format */
	    for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary lda = n, n+1, 2*n */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		  saved_seed = *seed;
		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  alpha_flag = 0;
		  switch (alpha_val) {
		  case 0:
		    alpha[0] = alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  case 1:
		    alpha[0] = 1.0;
		    alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  }
		  beta_flag = 0;
		  switch (beta_val) {
		  case 0:
		    beta[0] = beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  case 1:
		    beta[0] = 1.0;
		    beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  }

		  /* finally we are here to generate the test case */
		  BLAS_zhemv2_z_c_testgen(norm, order_type,
					  uplo_type, n, &alpha, alpha_flag,
					  &beta, beta_flag, a, lda,
					  head_x_gen, tail_x_gen, y_gen, seed,
					  head_r_true, tail_r_true);
		  test_count++;

		  /* vary incx = -2, -1, 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    /* vary incy = -2, -1, 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      /* copy generated vector with appropriate incs. */
		      zcopy_vector(y_gen, n, 1, y, incy);
		      ccopy_vector(head_x_gen, n, 1, head_x, incx);
		      ccopy_vector(tail_x_gen, n, 1, tail_x, incx);

		      /* call hemv2 routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zhemv2_z_c_x(order_type,
					uplo_type, n, alpha, a, lda, head_x,
					tail_x, incx, beta, y, incy, prec);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run dot test */

		      incyi = incy;
		      incyi *= 2;
		      yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		      for (i = 0, yi = yi0, ri = 0;
			   i < n; i++, yi += incyi, ri += incri) {
			zhe_copy_row(order_type, uplo_type, blas_left_side,
				     n, a, lda, a_vec, i);

			/* just use the x vector - it was unchanged (in theory) */
			rin[0] = y_gen[i];
			rin[1] = y_gen[i + 1];
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot2_z_c(n, blas_no_conj, alpha, beta,
					    rin, rout, head_r_true_elem,
					    tail_r_true_elem, a_vec, 1,
					    head_x, tail_x, incx, eps_int,
					    un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];

			  /* The !<= below causes NaN errors to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}

		      }		/* end of dot-test loop */


		      /* The !<= below causes NaN errors to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : z, x type : c\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d\n", n);
			  printf("LDA %d  INCX %d  INCY %d\n", lda, incx,
				 incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  if (uplo_type == blas_upper)
			    printf("upper ");
			  else
			    printf("lower ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");

			  printf("a\n");
			  zhe_print_matrix(a, n, lda, order_type, uplo_type);
			  cprint_vector(head_x, n, incx, "head_x");
			  cprint_vector(tail_x, n, incx, "tail_x");
			  zprint_vector(y_gen, n, incy, "y_gen");
			  zprint_vector(y, n, incy, "y");
			  zprint_vector(head_r_true, n, 1, "head_r_true");
			  dprint_vector(ratios, n, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			}
			bad_ratio_count++;
			if (bad_ratio_count >= MAX_BAD_TESTS) {
			  printf("\ntoo many failures, exiting....");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
			if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			  printf("\nFlagrant ratio error, exiting...");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
		      }

		      if (!(ratio <= ratio_max))
			ratio_max = ratio;
		      if (ratio != 0.0 && !(ratio >= ratio_min))
			ratio_min = ratio;

		    }		/* end of incy loop */

		  }		/* end of incx loop */

		}		/* end of lda loop */

	      }			/* end of uplo loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zhemv2_c_z_x
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_zhemv2_c_z_x";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;

  int prec_val;

  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  double alpha[2];
  double beta[2];
  float *a;
  double *head_x;
  double *tail_x;
  double *y;
  float *a_vec;
  double *y_gen;
  double *head_x_gen;
  double *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n * n * sizeof(float) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


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

	/* vary norm -- underflow, approx 1, overflow */
	for (norm = NORM_START; norm <= NORM_END; norm++) {

	  /* number of tests */
	  for (test_no = 0; test_no < ntests; test_no++) {

	    /* vary storage format */
	    for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary lda = n, n+1, 2*n */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		  saved_seed = *seed;
		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  alpha_flag = 0;
		  switch (alpha_val) {
		  case 0:
		    alpha[0] = alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  case 1:
		    alpha[0] = 1.0;
		    alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  }
		  beta_flag = 0;
		  switch (beta_val) {
		  case 0:
		    beta[0] = beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  case 1:
		    beta[0] = 1.0;
		    beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  }

		  /* finally we are here to generate the test case */
		  BLAS_zhemv2_c_z_testgen(norm, order_type,
					  uplo_type, n, &alpha, alpha_flag,
					  &beta, beta_flag, a, lda,
					  head_x_gen, tail_x_gen, y_gen, seed,
					  head_r_true, tail_r_true);
		  test_count++;

		  /* vary incx = -2, -1, 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    /* vary incy = -2, -1, 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      /* copy generated vector with appropriate incs. */
		      zcopy_vector(y_gen, n, 1, y, incy);
		      zcopy_vector(head_x_gen, n, 1, head_x, incx);
		      zcopy_vector(tail_x_gen, n, 1, tail_x, incx);

		      /* call hemv2 routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zhemv2_c_z_x(order_type,
					uplo_type, n, alpha, a, lda, head_x,
					tail_x, incx, beta, y, incy, prec);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run dot test */

		      incyi = incy;
		      incyi *= 2;
		      yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		      for (i = 0, yi = yi0, ri = 0;
			   i < n; i++, yi += incyi, ri += incri) {
			che_copy_row(order_type, uplo_type, blas_left_side,
				     n, a, lda, a_vec, i);

			/* just use the x vector - it was unchanged (in theory) */
			rin[0] = y_gen[i];
			rin[1] = y_gen[i + 1];
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot2_c_z(n, blas_no_conj, alpha, beta,
					    rin, rout, head_r_true_elem,
					    tail_r_true_elem, a_vec, 1,
					    head_x, tail_x, incx, eps_int,
					    un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];

			  /* The !<= below causes NaN errors to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}

		      }		/* end of dot-test loop */


		      /* The !<= below causes NaN errors to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : c, x type : z\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d\n", n);
			  printf("LDA %d  INCX %d  INCY %d\n", lda, incx,
				 incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  if (uplo_type == blas_upper)
			    printf("upper ");
			  else
			    printf("lower ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");

			  printf("a\n");
			  che_print_matrix(a, n, lda, order_type, uplo_type);
			  zprint_vector(head_x, n, incx, "head_x");
			  zprint_vector(tail_x, n, incx, "tail_x");
			  zprint_vector(y_gen, n, incy, "y_gen");
			  zprint_vector(y, n, incy, "y");
			  zprint_vector(head_r_true, n, 1, "head_r_true");
			  dprint_vector(ratios, n, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			}
			bad_ratio_count++;
			if (bad_ratio_count >= MAX_BAD_TESTS) {
			  printf("\ntoo many failures, exiting....");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
			if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			  printf("\nFlagrant ratio error, exiting...");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
		      }

		      if (!(ratio <= ratio_max))
			ratio_max = ratio;
		      if (ratio != 0.0 && !(ratio >= ratio_min))
			ratio_min = ratio;

		    }		/* end of incy loop */

		  }		/* end of incx loop */

		}		/* end of lda loop */

	      }			/* end of uplo loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zhemv2_c_c_x
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_zhemv2_c_c_x";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;

  int prec_val;

  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  double alpha[2];
  double beta[2];
  float *a;
  float *head_x;
  float *tail_x;
  double *y;
  float *a_vec;
  double *y_gen;
  float *head_x_gen;
  float *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n * n * sizeof(float) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


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

	/* vary norm -- underflow, approx 1, overflow */
	for (norm = NORM_START; norm <= NORM_END; norm++) {

	  /* number of tests */
	  for (test_no = 0; test_no < ntests; test_no++) {

	    /* vary storage format */
	    for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary lda = n, n+1, 2*n */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		  saved_seed = *seed;
		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  alpha_flag = 0;
		  switch (alpha_val) {
		  case 0:
		    alpha[0] = alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  case 1:
		    alpha[0] = 1.0;
		    alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  }
		  beta_flag = 0;
		  switch (beta_val) {
		  case 0:
		    beta[0] = beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  case 1:
		    beta[0] = 1.0;
		    beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  }

		  /* finally we are here to generate the test case */
		  BLAS_zhemv2_c_c_testgen(norm, order_type,
					  uplo_type, n, &alpha, alpha_flag,
					  &beta, beta_flag, a, lda,
					  head_x_gen, tail_x_gen, y_gen, seed,
					  head_r_true, tail_r_true);
		  test_count++;

		  /* vary incx = -2, -1, 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    /* vary incy = -2, -1, 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      /* copy generated vector with appropriate incs. */
		      zcopy_vector(y_gen, n, 1, y, incy);
		      ccopy_vector(head_x_gen, n, 1, head_x, incx);
		      ccopy_vector(tail_x_gen, n, 1, tail_x, incx);

		      /* call hemv2 routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zhemv2_c_c_x(order_type,
					uplo_type, n, alpha, a, lda, head_x,
					tail_x, incx, beta, y, incy, prec);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run dot test */

		      incyi = incy;
		      incyi *= 2;
		      yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		      for (i = 0, yi = yi0, ri = 0;
			   i < n; i++, yi += incyi, ri += incri) {
			che_copy_row(order_type, uplo_type, blas_left_side,
				     n, a, lda, a_vec, i);

			/* just use the x vector - it was unchanged (in theory) */
			rin[0] = y_gen[i];
			rin[1] = y_gen[i + 1];
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot2_c_c(n, blas_no_conj, alpha, beta,
					    rin, rout, head_r_true_elem,
					    tail_r_true_elem, a_vec, 1,
					    head_x, tail_x, incx, eps_int,
					    un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];

			  /* The !<= below causes NaN errors to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}

		      }		/* end of dot-test loop */


		      /* The !<= below causes NaN errors to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : c, x type : c\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d\n", n);
			  printf("LDA %d  INCX %d  INCY %d\n", lda, incx,
				 incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  if (uplo_type == blas_upper)
			    printf("upper ");
			  else
			    printf("lower ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");

			  printf("a\n");
			  che_print_matrix(a, n, lda, order_type, uplo_type);
			  cprint_vector(head_x, n, incx, "head_x");
			  cprint_vector(tail_x, n, incx, "tail_x");
			  zprint_vector(y_gen, n, incy, "y_gen");
			  zprint_vector(y, n, incy, "y");
			  zprint_vector(head_r_true, n, 1, "head_r_true");
			  dprint_vector(ratios, n, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			}
			bad_ratio_count++;
			if (bad_ratio_count >= MAX_BAD_TESTS) {
			  printf("\ntoo many failures, exiting....");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
			if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			  printf("\nFlagrant ratio error, exiting...");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
		      }

		      if (!(ratio <= ratio_max))
			ratio_max = ratio;
		      if (ratio != 0.0 && !(ratio >= ratio_min))
			ratio_min = ratio;

		    }		/* end of incy loop */

		  }		/* end of incx loop */

		}		/* end of lda loop */

	      }			/* end of uplo loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_chemv2_c_s_x
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_chemv2_c_s_x";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;

  int prec_val;

  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *head_x;
  float *tail_x;
  float *y;
  float *a_vec;
  float *y_gen;
  float *head_x_gen;
  float *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n * n * sizeof(float) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


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

	/* vary norm -- underflow, approx 1, overflow */
	for (norm = NORM_START; norm <= NORM_END; norm++) {

	  /* number of tests */
	  for (test_no = 0; test_no < ntests; test_no++) {

	    /* vary storage format */
	    for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary lda = n, n+1, 2*n */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		  saved_seed = *seed;
		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  alpha_flag = 0;
		  switch (alpha_val) {
		  case 0:
		    alpha[0] = alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  case 1:
		    alpha[0] = 1.0;
		    alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  }
		  beta_flag = 0;
		  switch (beta_val) {
		  case 0:
		    beta[0] = beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  case 1:
		    beta[0] = 1.0;
		    beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  }

		  /* finally we are here to generate the test case */
		  BLAS_chemv2_c_s_testgen(norm, order_type,
					  uplo_type, n, &alpha, alpha_flag,
					  &beta, beta_flag, a, lda,
					  head_x_gen, tail_x_gen, y_gen, seed,
					  head_r_true, tail_r_true);
		  test_count++;

		  /* vary incx = -2, -1, 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    /* vary incy = -2, -1, 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      /* copy generated vector with appropriate incs. */
		      ccopy_vector(y_gen, n, 1, y, incy);
		      scopy_vector(head_x_gen, n, 1, head_x, incx);
		      scopy_vector(tail_x_gen, n, 1, tail_x, incx);

		      /* call hemv2 routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_chemv2_c_s_x(order_type,
					uplo_type, n, alpha, a, lda, head_x,
					tail_x, incx, beta, y, incy, prec);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run dot test */

		      incyi = incy;
		      incyi *= 2;
		      yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		      for (i = 0, yi = yi0, ri = 0;
			   i < n; i++, yi += incyi, ri += incri) {
			che_copy_row(order_type, uplo_type, blas_left_side,
				     n, a, lda, a_vec, i);

			/* just use the x vector - it was unchanged (in theory) */
			rin[0] = y_gen[i];
			rin[1] = y_gen[i + 1];
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_cdot2_c_s(n, blas_no_conj, alpha, beta,
					    rin, rout, head_r_true_elem,
					    tail_r_true_elem, a_vec, 1,
					    head_x, tail_x, incx, eps_int,
					    un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];

			  /* The !<= below causes NaN errors to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}

		      }		/* end of dot-test loop */


		      /* The !<= below causes NaN errors to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : c, a type : c, x type : s\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d\n", n);
			  printf("LDA %d  INCX %d  INCY %d\n", lda, incx,
				 incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  if (uplo_type == blas_upper)
			    printf("upper ");
			  else
			    printf("lower ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			  printf("\n");

			  printf("a\n");
			  che_print_matrix(a, n, lda, order_type, uplo_type);
			  sprint_vector(head_x, n, incx, "head_x");
			  sprint_vector(tail_x, n, incx, "tail_x");
			  cprint_vector(y_gen, n, incy, "y_gen");
			  cprint_vector(y, n, incy, "y");
			  zprint_vector(head_r_true, n, 1, "head_r_true");
			  dprint_vector(ratios, n, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			}
			bad_ratio_count++;
			if (bad_ratio_count >= MAX_BAD_TESTS) {
			  printf("\ntoo many failures, exiting....");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
			if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			  printf("\nFlagrant ratio error, exiting...");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
		      }

		      if (!(ratio <= ratio_max))
			ratio_max = ratio;
		      if (ratio != 0.0 && !(ratio >= ratio_min))
			ratio_min = ratio;

		    }		/* end of incy loop */

		  }		/* end of incx loop */

		}		/* end of lda loop */

	      }			/* end of uplo loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zhemv2_z_d_x
  (int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "do_test_zhemv2_z_d_x";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;

  int prec_val;

  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *head_x;
  double *tail_x;
  double *y;
  double *a_vec;
  double *y_gen;
  double *head_x_gen;
  double *tail_x_gen;
  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0)
    BLAS_error(fname, -1, n, NULL);
  if (ntests < 0)
    BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;
  incri *= 2;

  FPU_FIX_START;

  y = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  y_gen = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && y_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x_gen = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && head_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x_gen = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && tail_x_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n * n * sizeof(double) * 2);
  if (2 * n * n > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_x = (double *) blas_malloc(2 * n * sizeof(double));
  if (2 * n > 0 && head_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  tail_x = (double *) blas_malloc(2 * n * sizeof(double));
  if (2 * n > 0 && tail_x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha[0] = alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha[0] = 1.0;
      alpha[1] = 0.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta[0] = beta[1] = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta[0] = 1.0;
	beta[1] = 0.0;
	beta_flag = 1;
	break;
      }


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

	/* vary norm -- underflow, approx 1, overflow */
	for (norm = NORM_START; norm <= NORM_END; norm++) {

	  /* number of tests */
	  for (test_no = 0; test_no < ntests; test_no++) {

	    /* vary storage format */
	    for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary lda = n, n+1, 2*n */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? n : (lda_val == 1) ? n + 1 : 2 * n;

		  saved_seed = *seed;
		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  alpha_flag = 0;
		  switch (alpha_val) {
		  case 0:
		    alpha[0] = alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  case 1:
		    alpha[0] = 1.0;
		    alpha[1] = 0.0;
		    alpha_flag = 1;
		    break;
		  }
		  beta_flag = 0;
		  switch (beta_val) {
		  case 0:
		    beta[0] = beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  case 1:
		    beta[0] = 1.0;
		    beta[1] = 0.0;
		    beta_flag = 1;
		    break;
		  }

		  /* finally we are here to generate the test case */
		  BLAS_zhemv2_z_d_testgen(norm, order_type,
					  uplo_type, n, &alpha, alpha_flag,
					  &beta, beta_flag, a, lda,
					  head_x_gen, tail_x_gen, y_gen, seed,
					  head_r_true, tail_r_true);
		  test_count++;

		  /* vary incx = -2, -1, 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    /* vary incy = -2, -1, 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      /* copy generated vector with appropriate incs. */
		      zcopy_vector(y_gen, n, 1, y, incy);
		      dcopy_vector(head_x_gen, n, 1, head_x, incx);
		      dcopy_vector(tail_x_gen, n, 1, tail_x, incx);

		      /* call hemv2 routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zhemv2_z_d_x(order_type,
					uplo_type, n, alpha, a, lda, head_x,
					tail_x, incx, beta, y, incy, prec);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run dot test */

		      incyi = incy;
		      incyi *= 2;
		      yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

		      for (i = 0, yi = yi0, ri = 0;
			   i < n; i++, yi += incyi, ri += incri) {
			zhe_copy_row(order_type, uplo_type, blas_left_side,
				     n, a, lda, a_vec, i);

			/* just use the x vector - it was unchanged (in theory) */
			rin[0] = y_gen[i];
			rin[1] = y_gen[i + 1];
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot2_z_d(n, blas_no_conj, alpha, beta,
					    rin, rout, head_r_true_elem,
					    tail_r_true_elem, a_vec, 1,
					    head_x, tail_x, incx, eps_int,
					    un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];

			  /* The !<= below causes NaN errors to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}

		      }		/* end of dot-test loop */


		      /* The !<= below causes NaN errors to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : z, x type : d\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d\n", n);
			  printf("LDA %d  INCX %d  INCY %d\n", lda, incx,
				 incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  if (uplo_type == blas_upper)
			    printf("upper ");
			  else
			    printf("lower ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");

			  printf("a\n");
			  zhe_print_matrix(a, n, lda, order_type, uplo_type);
			  dprint_vector(head_x, n, incx, "head_x");
			  dprint_vector(tail_x, n, incx, "tail_x");
			  zprint_vector(y_gen, n, incy, "y_gen");
			  zprint_vector(y, n, incy, "y");
			  zprint_vector(head_r_true, n, 1, "head_r_true");
			  dprint_vector(ratios, n, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			}
			bad_ratio_count++;
			if (bad_ratio_count >= MAX_BAD_TESTS) {
			  printf("\ntoo many failures, exiting....");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
			if (!(ratio <= TOTAL_FAILURE_THRESHOLD)) {
			  printf("\nFlagrant ratio error, exiting...");
			  printf("\nTesting and compilation");
			  printf(" are incomplete\n\n");
			  goto end;
			}
		      }

		      if (!(ratio <= ratio_max))
			ratio_max = ratio;
		      if (ratio != 0.0 && !(ratio >= ratio_min))
			ratio_min = ratio;

		    }		/* end of incy loop */

		  }		/* end of incx loop */

		}		/* end of lda loop */

	      }			/* end of uplo loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(y);
  blas_free(a);
  blas_free(y_gen);
  blas_free(head_x);
  blas_free(tail_x);
  blas_free(head_x_gen);
  blas_free(tail_x_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

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
  const char *base_routine = "hemv2";
  char *fname;
  int n;

  int i;
  int n_data[NUM_DATA][1] = { {4}, {2}, {3}, {8}, {10}, {1}, {7} };

  if (argc != 6) {
    printf("Usage:\n");
    printf("do_test_hemv2 <nsizes> <ntests> <thresh> <debug> <test_prob>\n");
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
    BLAS_error("Testing hemv2", 0, 0, NULL);

  printf("Testing %s...\n", base_routine);
  printf("INPUT: nsizes = %d, ntests = %d, thresh = %4.2f, debug = %d\n\n",
	 nsizes, ntests, thresh, debug);



  if (nsizes < 0 || nsizes > NUM_DATA)
    BLAS_error("do_test_hemv2", -1, nsizes, NULL);

  fname = "BLAS_zhemv2_z_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_zhemv2_z_c(n, ntests, &seed, thresh, debug,
		       test_prob,
		       &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_zhemv2_c_z";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_zhemv2_c_z(n, ntests, &seed, thresh, debug,
		       test_prob,
		       &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_zhemv2_c_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_zhemv2_c_c(n, ntests, &seed, thresh, debug,
		       test_prob,
		       &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_chemv2_c_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_chemv2_c_s(n, ntests, &seed, thresh, debug,
		       test_prob,
		       &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_zhemv2_z_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_zhemv2_z_d(n, ntests, &seed, thresh, debug,
		       test_prob,
		       &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_chemv2_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_chemv2_x(n, ntests, &seed, thresh, debug,
		     test_prob,
		     &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_zhemv2_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_zhemv2_x(n, ntests, &seed, thresh, debug,
		     test_prob,
		     &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_zhemv2_z_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_zhemv2_z_c_x(n, ntests, &seed, thresh, debug,
			 test_prob,
			 &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_zhemv2_c_z_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_zhemv2_c_z_x(n, ntests, &seed, thresh, debug,
			 test_prob,
			 &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_zhemv2_c_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_zhemv2_c_c_x(n, ntests, &seed, thresh, debug,
			 test_prob,
			 &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_chemv2_c_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_chemv2_c_s_x(n, ntests, &seed, thresh, debug,
			 test_prob,
			 &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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

  fname = "BLAS_zhemv2_z_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    do_test_zhemv2_z_d_x(n, ntests, &seed, thresh, debug,
			 test_prob,
			 &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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
  if (nr_failed_routines)
    printf("FAILED ");
  else
    printf("PASSED ");
  printf("%-10s: FAIL/TOTAL = %d/%d\n",
	 base_routine, nr_failed_routines, nr_routines);

  return 0;
}
