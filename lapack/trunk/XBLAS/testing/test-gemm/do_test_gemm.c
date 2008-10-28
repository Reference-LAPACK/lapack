#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"







/* 0 -- 2 */
#define TRANSA_START    0
#define TRANSA_END      2

/* 0 -- 2 */
#define TRANSB_START    0
#define TRANSB_END      2

/* 0 -- 1 */
#define ORDER_START     0
#define ORDER_END       1

/* 0 -- 2 */
#define ALPHA_START     0
#define ALPHA_END       2

/* 0 -- 2 */
#define BETA_START      0
#define BETA_END        2

/* -1 -- 1 */
#define NORM_START     -1
#define NORM_END        1

/* 0 -- 2 */
#define LDA_START       0
#define LDA_END         2

/* 0 -- 2 */
#define LDB_START       0
#define LDB_END         2

/* 0 -- 2 */
#define LDC_START       0
#define LDC_END         2

/* 0 -- 2 */
#define PREC_START      0
#define PREC_END        2

/* 0 -- 1 */
#define RANDOMIZE_START 0
#define RANDOMIZE_END   1

#define NUM_DATA 9


void do_test_dgemm_d_s(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_dgemm_d_s";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha;
  double beta;
  double *a;
  float *b;
  double *c;
  double *a_vec;
  float *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;




  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha = 1.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta = 1.0;
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_dgemm_d_s_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			dge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_dgemm_d_s(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;




			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  dge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    sge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_ddot_d_s(k, blas_no_conj,
					       alpha, beta, c_gen[cij],
					       c[cij],
					       head_r_true[cij],
					       tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");

			    dge_print_matrix(a, m, k, lda, order, "A");
			    sge_print_matrix(b, k, n, ldb, order, "B");
			    dge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    dge_print_matrix(c, m, n, ldc, order, "C");
			    dge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dgemm_s_d(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_dgemm_s_d";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha;
  double beta;
  float *a;
  double *b;
  double *c;
  float *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;




  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha = 1.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta = 1.0;
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_dgemm_s_d_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			dge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_dgemm_s_d(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;




			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  sge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    dge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_ddot_s_d(k, blas_no_conj,
					       alpha, beta, c_gen[cij],
					       c[cij],
					       head_r_true[cij],
					       tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");

			    sge_print_matrix(a, m, k, lda, order, "A");
			    dge_print_matrix(b, k, n, ldb, order, "B");
			    dge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    dge_print_matrix(c, m, n, ldc, order, "C");
			    dge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dgemm_s_s(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_dgemm_s_s";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha;
  double beta;
  float *a;
  float *b;
  double *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;




  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha = 1.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta = 1.0;
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_dgemm_s_s_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			dge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_dgemm_s_s(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;




			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  sge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    sge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_ddot_s_s(k, blas_no_conj,
					       alpha, beta, c_gen[cij],
					       c[cij],
					       head_r_true[cij],
					       tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");

			    sge_print_matrix(a, m, k, lda, order, "A");
			    sge_print_matrix(b, k, n, ldb, order, "B");
			    dge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    dge_print_matrix(c, m, n, ldc, order, "C");
			    dge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_z_c(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_z_c";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  float *b;
  double *c;
  double *a_vec;
  float *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;
  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_zgemm_z_c_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_zgemm_z_c(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;

			incci *= 2;
			inccij *= 2;

			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  zge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    cge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_zdot_z_c(k, blas_no_conj,
					       alpha, beta, &c_gen[cij],
					       &c[cij],
					       &head_r_true[cij],
					       &tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");

			    zge_print_matrix(a, m, k, lda, order, "A");
			    cge_print_matrix(b, k, n, ldb, order, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order, "C");
			    zge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_c_z(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_c_z";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  float *a;
  double *b;
  double *c;
  float *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;
  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_zgemm_c_z_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_zgemm_c_z(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;

			incci *= 2;
			inccij *= 2;

			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  cge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    zge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_zdot_c_z(k, blas_no_conj,
					       alpha, beta, &c_gen[cij],
					       &c[cij],
					       &head_r_true[cij],
					       &tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");

			    cge_print_matrix(a, m, k, lda, order, "A");
			    zge_print_matrix(b, k, n, ldb, order, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order, "C");
			    zge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_c_c(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_c_c";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  float *a;
  float *b;
  double *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;
  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_zgemm_c_c_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_zgemm_c_c(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;

			incci *= 2;
			inccij *= 2;

			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  cge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    cge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_zdot_c_c(k, blas_no_conj,
					       alpha, beta, &c_gen[cij],
					       &c[cij],
					       &head_r_true[cij],
					       &tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");

			    cge_print_matrix(a, m, k, lda, order, "A");
			    cge_print_matrix(b, k, n, ldb, order, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order, "C");
			    zge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cgemm_c_s(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_cgemm_c_s";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *b;
  float *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  float *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;

  incc *= 2;

  /* allocate memory for arrays */
  c = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_cgemm_c_s_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			cge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_cgemm_c_s(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;

			incci *= 2;
			inccij *= 2;

			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  cge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    sge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_cdot_c_s(k, blas_no_conj,
					       alpha, beta, &c_gen[cij],
					       &c[cij],
					       &head_r_true[cij],
					       &tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			    printf("\n");

			    cge_print_matrix(a, m, k, lda, order, "A");
			    sge_print_matrix(b, k, n, ldb, order, "B");
			    cge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    cge_print_matrix(c, m, n, ldc, order, "C");
			    zge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cgemm_s_c(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_cgemm_s_c";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *b;
  float *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  float *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;

  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_cgemm_s_c_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			cge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_cgemm_s_c(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;

			incci *= 2;
			inccij *= 2;

			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  sge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    cge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_cdot_s_c(k, blas_no_conj,
					       alpha, beta, &c_gen[cij],
					       &c[cij],
					       &head_r_true[cij],
					       &tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			    printf("\n");

			    sge_print_matrix(a, m, k, lda, order, "A");
			    cge_print_matrix(b, k, n, ldb, order, "B");
			    cge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    cge_print_matrix(c, m, n, ldc, order, "C");
			    zge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cgemm_s_s(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_cgemm_s_s";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *b;
  float *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  float *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;


  incc *= 2;

  /* allocate memory for arrays */
  c = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_cgemm_s_s_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			cge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_cgemm_s_s(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;

			incci *= 2;
			inccij *= 2;

			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  sge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    sge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_cdot_s_s(k, blas_no_conj,
					       alpha, beta, &c_gen[cij],
					       &c[cij],
					       &head_r_true[cij],
					       &tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			    printf("\n");

			    sge_print_matrix(a, m, k, lda, order, "A");
			    sge_print_matrix(b, k, n, ldb, order, "B");
			    cge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    cge_print_matrix(c, m, n, ldc, order, "C");
			    zge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_z_d(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_z_d";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *b;
  double *c;
  double *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;

  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_zgemm_z_d_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_zgemm_z_d(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;

			incci *= 2;
			inccij *= 2;

			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  zge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    dge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_zdot_z_d(k, blas_no_conj,
					       alpha, beta, &c_gen[cij],
					       &c[cij],
					       &head_r_true[cij],
					       &tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");

			    zge_print_matrix(a, m, k, lda, order, "A");
			    dge_print_matrix(b, k, n, ldb, order, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order, "C");
			    zge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_d_z(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_d_z";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *b;
  double *c;
  double *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;

  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_zgemm_d_z_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_zgemm_d_z(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;

			incci *= 2;
			inccij *= 2;

			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  dge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    zge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_zdot_d_z(k, blas_no_conj,
					       alpha, beta, &c_gen[cij],
					       &c[cij],
					       &head_r_true[cij],
					       &tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");

			    dge_print_matrix(a, m, k, lda, order, "A");
			    zge_print_matrix(b, k, n, ldb, order, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order, "C");
			    zge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_d_d(int m, int n, int k, int ntests, int *seed,
		       double thresh, int debug, float test_prob,
		       double *min_ratio, double *max_ratio,
		       int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_d_d";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *b;
  double *c;
  double *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;


  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	    if (order_val == 0) {
	      /* row major storage */
	      lda_0 = k;
	      ldb_0 = n;
	      ldc_1 = n;
	      tda_0 = m;
	      tdb_0 = k;
	      order = blas_rowmajor;
	    } else {
	      /* column major storage */
	      lda_0 = m;
	      ldb_0 = k;
	      ldc_1 = m;
	      tda_0 = k;
	      tdb_0 = n;
	      order = blas_colmajor;
	    }

	    /* vary transpositions of A */
	    for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		 transa_val++) {

	      transa = (transa_val == 0) ? blas_no_trans :
		(transa_val == 1) ? blas_trans : blas_conj_trans;

	      if (transa == blas_no_trans) {
		lda_1 = lda_0;
	      } else {
		lda_1 = tda_0;
	      }

	      /* vary transpositions of B */
	      for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		   transb_val++) {

		transb = (transb_val == 0) ? blas_no_trans :
		  (transb_val == 1) ? blas_trans : blas_conj_trans;

		if (transb == blas_no_trans) {
		  ldb_1 = ldb_0;
		} else {
		  ldb_1 = tdb_0;
		}

		/* vary lda = k, k+1, 2*k */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? lda_1 :
		    (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		  /* vary ldb = n, n+1, 2*n */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    ldb = (ldb_val == 0) ? ldb_1 :
		      (ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      ldc = (ldc_val == 0) ? ldc_1 :
			(ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			/* finally we are here to generate the test case */
			BLAS_zgemm_d_d_testgen(norm, order,
					       transa, transb, m, n, k,
					       randomize_val, &alpha,
					       alpha_flag, a, lda, &beta,
					       beta_flag, b, ldb, c, ldc,
					       seed, head_r_true,
					       tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			/* call GEMM routines to be tested */
			FPU_FIX_STOP;
			BLAS_zgemm_d_d(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if (order == blas_colmajor) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;

			incci *= 2;
			inccij *= 2;

			for (i = 0, ci = 0, ri = 0; i < m;
			     i++, ci += incci, ri += incri) {
			  dge_copy_row(order, transa, m, k, a, lda, a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    dge_copy_col(order, transb,
					 k, n, b, ldb, b_vec, j);

			    test_BLAS_zdot_d_d(k, blas_no_conj,
					       alpha, beta, &c_gen[cij],
					       &c[cij],
					       &head_r_true[cij],
					       &tail_r_true[cij], a_vec, 1,
					       b_vec, 1, eps_int, un_int,
					       &ratios[rij]);

			    /* take the max ratio */
			    if (rij == 0) {
			      ratio = ratios[0];
			      /* The !<= below causes NaN error to be detected.
			         Note that (NaN > thresh) is always false. */
			    } else if (!(ratios[rij] <= ratio)) {
			      ratio = ratios[rij];
			    }

			  }
			}	/* end of dot-test loop */

			/* Increase the number of bad ratio, if the ratio
			   is bigger than the threshold.
			   The !<= below causes NaN error to be detected.
			   Note that (NaN > thresh) is always false. */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {

			    printf("\nm %d   n %d   k %d\n", m, n, k);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    printf("A:");
			    switch (transa) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }
			    printf("B:");
			    switch (transb) {
			    case blas_no_trans:
			      printf("no_trans ");
			      break;
			    case blas_trans:
			      printf("trans ");
			      break;
			    case blas_conj_trans:
			      printf("conj_trans ");
			      break;
			    }

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    switch (order) {
			    case blas_rowmajor:
			      printf("row_major ");
			      break;
			    case blas_colmajor:
			      printf("col_major ");
			      break;
			    }


			    if (randomize_val == 0)
			      printf("Not randomized\n");
			    else
			      printf("Randomized\n");

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");

			    dge_print_matrix(a, m, k, lda, order, "A");
			    dge_print_matrix(b, k, n, ldb, order, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order, "C");
			    zge_print_matrix(head_r_true, m, n, ldc, order,
					     "truth");

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
			    printf("\nFlagrant ratio %e, exiting...", ratio);
			    printf("\nTesting and compilation");
			    printf(" are incomplete\n\n");
			    goto end;
			  }
			}

			if (ratio > ratio_max)
			  ratio_max = ratio;

			if (ratio != 0.0 && ratio < ratio_min)
			  ratio_min = ratio;

		      }		/* end of randomize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of transb loop */

	    }			/* end of transa loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_sgemm_x(int m, int n, int k, int ntests, int *seed,
		     double thresh, int debug, float test_prob,
		     double *min_ratio, double *max_ratio, int *num_bad_ratio,
		     int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_sgemm_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  float alpha;
  float beta;
  float *a;
  float *b;
  float *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  float *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;




  /* allocate memory for arrays */
  c = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha = 1.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta = 1.0;
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_sgemm_testgen(norm, order,
					     transa, transb, m, n, k,
					     randomize_val, &alpha,
					     alpha_flag, a, lda, &beta,
					     beta_flag, b, ldb, c, ldc, seed,
					     head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  sge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_sgemm_x(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;




			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    sge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      sge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_sdot(k, blas_no_conj,
					     alpha, beta, c_gen[cij],
					     c[cij],
					     head_r_true[cij],
					     tail_r_true[cij], a_vec, 1,
					     b_vec, 1, eps_int, un_int,
					     &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("%16.8e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%16.8e", beta);;
			      printf("\n");

			      sge_print_matrix(a, m, k, lda, order, "A");
			      sge_print_matrix(b, k, n, ldb, order, "B");
			      sge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      sge_print_matrix(c, m, n, ldc, order, "C");
			      dge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dgemm_x(int m, int n, int k, int ntests, int *seed,
		     double thresh, int debug, float test_prob,
		     double *min_ratio, double *max_ratio, int *num_bad_ratio,
		     int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_dgemm_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha;
  double beta;
  double *a;
  double *b;
  double *c;
  double *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;




  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha = 1.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta = 1.0;
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_dgemm_testgen(norm, order,
					     transa, transb, m, n, k,
					     randomize_val, &alpha,
					     alpha_flag, a, lda, &beta,
					     beta_flag, b, ldb, c, ldc, seed,
					     head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  dge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_dgemm_x(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;




			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    dge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      dge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_ddot(k, blas_no_conj,
					     alpha, beta, c_gen[cij],
					     c[cij],
					     head_r_true[cij],
					     tail_r_true[cij], a_vec, 1,
					     b_vec, 1, eps_int, un_int,
					     &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("%24.16e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%24.16e", beta);;
			      printf("\n");

			      dge_print_matrix(a, m, k, lda, order, "A");
			      dge_print_matrix(b, k, n, ldb, order, "B");
			      dge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      dge_print_matrix(c, m, n, ldc, order, "C");
			      dge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cgemm_x(int m, int n, int k, int ntests, int *seed,
		     double thresh, int debug, float test_prob,
		     double *min_ratio, double *max_ratio, int *num_bad_ratio,
		     int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_cgemm_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *b;
  float *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  float *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;
  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_cgemm_testgen(norm, order,
					     transa, transb, m, n, k,
					     randomize_val, &alpha,
					     alpha_flag, a, lda, &beta,
					     beta_flag, b, ldb, c, ldc, seed,
					     head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  cge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_cgemm_x(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    cge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      cge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_cdot(k, blas_no_conj,
					     alpha, beta, &c_gen[cij],
					     &c[cij],
					     &head_r_true[cij],
					     &tail_r_true[cij], a_vec, 1,
					     b_vec, 1, eps_int, un_int,
					     &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			      printf("\n");

			      cge_print_matrix(a, m, k, lda, order, "A");
			      cge_print_matrix(b, k, n, ldb, order, "B");
			      cge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      cge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_x(int m, int n, int k, int ntests, int *seed,
		     double thresh, int debug, float test_prob,
		     double *min_ratio, double *max_ratio, int *num_bad_ratio,
		     int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *b;
  double *c;
  double *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;
  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_zgemm_testgen(norm, order,
					     transa, transb, m, n, k,
					     randomize_val, &alpha,
					     alpha_flag, a, lda, &beta,
					     beta_flag, b, ldb, c, ldc, seed,
					     head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zgemm_x(order, transa,
				       transb, m, n, k, alpha, a, lda, b, ldb,
				       beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    zge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      zge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_zdot(k, blas_no_conj,
					     alpha, beta, &c_gen[cij],
					     &c[cij],
					     &head_r_true[cij],
					     &tail_r_true[cij], a_vec, 1,
					     b_vec, 1, eps_int, un_int,
					     &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      zge_print_matrix(a, m, k, lda, order, "A");
			      zge_print_matrix(b, k, n, ldb, order, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dgemm_d_s_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_dgemm_d_s_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha;
  double beta;
  double *a;
  float *b;
  double *c;
  double *a_vec;
  float *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;




  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha = 1.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta = 1.0;
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_dgemm_d_s_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  dge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_dgemm_d_s_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;




			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    dge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      sge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_ddot_d_s(k, blas_no_conj,
						 alpha, beta, c_gen[cij],
						 c[cij],
						 head_r_true[cij],
						 tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("%24.16e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%24.16e", beta);;
			      printf("\n");

			      dge_print_matrix(a, m, k, lda, order, "A");
			      sge_print_matrix(b, k, n, ldb, order, "B");
			      dge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      dge_print_matrix(c, m, n, ldc, order, "C");
			      dge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dgemm_s_d_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_dgemm_s_d_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha;
  double beta;
  float *a;
  double *b;
  double *c;
  float *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;




  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha = 1.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta = 1.0;
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_dgemm_s_d_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  dge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_dgemm_s_d_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;




			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    sge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      dge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_ddot_s_d(k, blas_no_conj,
						 alpha, beta, c_gen[cij],
						 c[cij],
						 head_r_true[cij],
						 tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("%24.16e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%24.16e", beta);;
			      printf("\n");

			      sge_print_matrix(a, m, k, lda, order, "A");
			      dge_print_matrix(b, k, n, ldb, order, "B");
			      dge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      dge_print_matrix(c, m, n, ldc, order, "C");
			      dge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dgemm_s_s_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_dgemm_s_s_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha;
  double beta;
  float *a;
  float *b;
  double *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;




  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {

    alpha_flag = 0;
    switch (alpha_val) {
    case 0:
      alpha = 0.0;
      alpha_flag = 1;
      break;
    case 1:
      alpha = 1.0;
      alpha_flag = 1;
      break;
    }

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      beta_flag = 0;
      switch (beta_val) {
      case 0:
	beta = 0.0;
	beta_flag = 1;
	break;
      case 1:
	beta = 1.0;
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_dgemm_s_s_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  dge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_dgemm_s_s_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;




			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    sge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      sge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_ddot_s_s(k, blas_no_conj,
						 alpha, beta, c_gen[cij],
						 c[cij],
						 head_r_true[cij],
						 tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("%24.16e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%24.16e", beta);;
			      printf("\n");

			      sge_print_matrix(a, m, k, lda, order, "A");
			      sge_print_matrix(b, k, n, ldb, order, "B");
			      dge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      dge_print_matrix(c, m, n, ldc, order, "C");
			      dge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_z_c_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_z_c_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  float *b;
  double *c;
  double *a_vec;
  float *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;
  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_zgemm_z_c_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zgemm_z_c_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    zge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      cge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_zdot_z_c(k, blas_no_conj,
						 alpha, beta, &c_gen[cij],
						 &c[cij],
						 &head_r_true[cij],
						 &tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      zge_print_matrix(a, m, k, lda, order, "A");
			      cge_print_matrix(b, k, n, ldb, order, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_c_z_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_c_z_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  float *a;
  double *b;
  double *c;
  float *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;
  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_zgemm_c_z_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zgemm_c_z_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    cge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      zge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_zdot_c_z(k, blas_no_conj,
						 alpha, beta, &c_gen[cij],
						 &c[cij],
						 &head_r_true[cij],
						 &tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      cge_print_matrix(a, m, k, lda, order, "A");
			      zge_print_matrix(b, k, n, ldb, order, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_c_c_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_c_c_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  float *a;
  float *b;
  double *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;
  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_zgemm_c_c_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zgemm_c_c_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    cge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      cge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_zdot_c_c(k, blas_no_conj,
						 alpha, beta, &c_gen[cij],
						 &c[cij],
						 &head_r_true[cij],
						 &tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      cge_print_matrix(a, m, k, lda, order, "A");
			      cge_print_matrix(b, k, n, ldb, order, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cgemm_c_s_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_cgemm_c_s_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *b;
  float *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  float *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;

  incc *= 2;

  /* allocate memory for arrays */
  c = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_cgemm_c_s_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  cge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_cgemm_c_s_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    cge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      sge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_cdot_c_s(k, blas_no_conj,
						 alpha, beta, &c_gen[cij],
						 &c[cij],
						 &head_r_true[cij],
						 &tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			      printf("\n");

			      cge_print_matrix(a, m, k, lda, order, "A");
			      sge_print_matrix(b, k, n, ldb, order, "B");
			      cge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      cge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cgemm_s_c_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_cgemm_s_c_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *b;
  float *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  float *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;

  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_cgemm_s_c_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  cge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_cgemm_s_c_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    sge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      cge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_cdot_s_c(k, blas_no_conj,
						 alpha, beta, &c_gen[cij],
						 &c[cij],
						 &head_r_true[cij],
						 &tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			      printf("\n");

			      sge_print_matrix(a, m, k, lda, order, "A");
			      cge_print_matrix(b, k, n, ldb, order, "B");
			      cge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      cge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cgemm_s_s_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_cgemm_s_s_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  float alpha[2];
  float beta[2];
  float *a;
  float *b;
  float *c;
  float *a_vec;
  float *b_vec;

  /* generated test values for c */
  float *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;


  incc *= 2;

  /* allocate memory for arrays */
  c = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * m * k * sizeof(float));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * k * n * sizeof(float));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_cgemm_s_s_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  cge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_cgemm_s_s_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    sge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      sge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_cdot_s_s(k, blas_no_conj,
						 alpha, beta, &c_gen[cij],
						 &c[cij],
						 &head_r_true[cij],
						 &tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			      printf("\n");

			      sge_print_matrix(a, m, k, lda, order, "A");
			      sge_print_matrix(b, k, n, ldb, order, "B");
			      cge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      cge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_z_d_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_z_d_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *b;
  double *c;
  double *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;
  inca *= 2;

  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double) * 2);
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_zgemm_z_d_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zgemm_z_d_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    zge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      dge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_zdot_z_d(k, blas_no_conj,
						 alpha, beta, &c_gen[cij],
						 &c[cij],
						 &head_r_true[cij],
						 &tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      zge_print_matrix(a, m, k, lda, order, "A");
			      dge_print_matrix(b, k, n, ldb, order, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_d_z_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_d_z_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *b;
  double *c;
  double *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;

  incb *= 2;
  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double) * 2);
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_zgemm_d_z_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zgemm_d_z_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    dge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      zge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_zdot_d_z(k, blas_no_conj,
						 alpha, beta, &c_gen[cij],
						 &c[cij],
						 &head_r_true[cij],
						 &tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      dge_print_matrix(a, m, k, lda, order, "A");
			      zge_print_matrix(b, k, n, ldb, order, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zgemm_d_d_x(int m, int n, int k, int ntests, int *seed,
			 double thresh, int debug, float test_prob,
			 double *min_ratio, double *max_ratio,
			 int *num_bad_ratio, int *num_tests)
{

  /* Function name */
  const char fname[] = "BLAS_zgemm_d_d_x";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;		/* number of tests done so far   */
  int bad_ratio_count;		/* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;

  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;

  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

  double alpha[2];
  double beta[2];
  double *a;
  double *b;
  double *c;
  double *a_vec;
  double *b_vec;

  /* generated test values for c */
  double *c_gen;

  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;

  FPU_FIX_START;

  inca = incb = incc = 1;


  incc *= 2;

  /* allocate memory for arrays */
  c = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  c_gen = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && c_gen == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * m * k * sizeof(double));
  if (2 * m * k > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * k * n * sizeof(double));
  if (2 * k * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && ratios == NULL) {
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

	      if (order_val == 0) {
		/* row major storage */
		lda_0 = k;
		ldb_0 = n;
		ldc_1 = n;
		tda_0 = m;
		tdb_0 = k;
		order = blas_rowmajor;
	      } else {
		/* column major storage */
		lda_0 = m;
		ldb_0 = k;
		ldc_1 = m;
		tda_0 = k;
		tdb_0 = n;
		order = blas_colmajor;
	      }

	      /* vary transpositions of A */
	      for (transa_val = TRANSA_START; transa_val <= TRANSA_END;
		   transa_val++) {

		transa = (transa_val == 0) ? blas_no_trans :
		  (transa_val == 1) ? blas_trans : blas_conj_trans;

		if (transa == blas_no_trans) {
		  lda_1 = lda_0;
		} else {
		  lda_1 = tda_0;
		}

		/* vary transpositions of B */
		for (transb_val = TRANSB_START; transb_val <= TRANSB_END;
		     transb_val++) {

		  transb = (transb_val == 0) ? blas_no_trans :
		    (transb_val == 1) ? blas_trans : blas_conj_trans;

		  if (transb == blas_no_trans) {
		    ldb_1 = ldb_0;
		  } else {
		    ldb_1 = tdb_0;
		  }

		  /* vary lda = k, k+1, 2*k */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? lda_1 :
		      (lda_val == 1) ? lda_1 + 1 : 2 * lda_1;

		    /* vary ldb = n, n+1, 2*n */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      ldb = (ldb_val == 0) ? ldb_1 :
			(ldb_val == 1) ? ldb_1 + 1 : 2 * ldb_1;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			ldc = (ldc_val == 0) ? ldc_1 :
			  (ldc_val == 1) ? ldc_1 + 1 : 2 * ldc_1;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  /* finally we are here to generate the test case */
			  BLAS_zgemm_d_d_testgen(norm, order,
						 transa, transb, m, n, k,
						 randomize_val, &alpha,
						 alpha_flag, a, lda, &beta,
						 beta_flag, b, ldb, c, ldc,
						 seed, head_r_true,
						 tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

			  /* call GEMM routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zgemm_d_d_x(order, transa,
					   transb, m, n, k, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if (order == blas_colmajor) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;

			  incci *= 2;
			  inccij *= 2;

			  for (i = 0, ci = 0, ri = 0; i < m;
			       i++, ci += incci, ri += incri) {
			    dge_copy_row(order, transa,
					 m, k, a, lda, a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      dge_copy_col(order, transb,
					   k, n, b, ldb, b_vec, j);

			      test_BLAS_zdot_d_d(k, blas_no_conj,
						 alpha, beta, &c_gen[cij],
						 &c[cij],
						 &head_r_true[cij],
						 &tail_r_true[cij], a_vec, 1,
						 b_vec, 1, eps_int, un_int,
						 &ratios[rij]);

			      /* take the max ratio */
			      if (rij == 0) {
				ratio = ratios[0];
				/* The !<= below causes NaN error to be detected.
				   Note that (NaN > thresh) is always false. */
			      } else if (!(ratios[rij] <= ratio)) {
				ratio = ratios[rij];
			      }

			    }
			  }	/* end of dot-test loop */

			  /* Increase the number of bad ratio, if the ratio
			     is bigger than the threshold.
			     The !<= below causes NaN error to be detected.
			     Note that (NaN > thresh) is always false. */
			  if (!(ratio <= thresh)) {

			    if (debug == 3) {

			      printf("\nm %d   n %d   k %d\n", m, n, k);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      printf("A:");
			      switch (transa) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }
			      printf("B:");
			      switch (transb) {
			      case blas_no_trans:
				printf("no_trans ");
				break;
			      case blas_trans:
				printf("trans ");
				break;
			      case blas_conj_trans:
				printf("conj_trans ");
				break;
			      }

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      switch (order) {
			      case blas_rowmajor:
				printf("row_major ");
				break;
			      case blas_colmajor:
				printf("col_major ");
				break;
			      }
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

			      if (randomize_val == 0)
				printf("Not randomized\n");
			      else
				printf("Randomized\n");

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      dge_print_matrix(a, m, k, lda, order, "A");
			      dge_print_matrix(b, k, n, ldb, order, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order, "C");
			      zge_print_matrix(head_r_true, m, n, ldc, order,
					       "truth");

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
			      printf("\nFlagrant ratio %e, exiting...",
				     ratio);
			      printf("\nTesting and compilation");
			      printf(" are incomplete\n\n");
			      goto end;
			    }
			  }

			  if (ratio > ratio_max)
			    ratio_max = ratio;

			  if (ratio != 0.0 && ratio < ratio_min)
			    ratio_min = ratio;

			}	/* end of randomize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of transb loop */

	      }			/* end of transa loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

end:
  FPU_FIX_STOP;

  blas_free(c);
  blas_free(a);
  blas_free(b);
  blas_free(c_gen);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(b_vec);

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
  const char *base_routine = "gemm";
  char *fname;
  int n;

  int m, k, i;
  int mnk_data[NUM_DATA][3] = { {2, 3, 4}, {3, 6, 4}, {5, 1, 7},
  {6, 15, 4}, {5, 2, 3}, {8, 4, 1},
  {1, 3, 1}, {8, 8, 8}, {1, 1, 1}
  };

  if (argc != 6) {
    printf("Usage:\n");
    printf("do_test_gemm <nsizes> <ntests> <thresh> <debug> <test_prob>\n");
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
    BLAS_error("Testing gemm", 0, 0, NULL);

  printf("Testing %s...\n", base_routine);
  printf("INPUT: nsizes = %d, ntests = %d, thresh = %4.2f, debug = %d\n\n",
	 nsizes, ntests, thresh, debug);





  fname = "BLAS_dgemm_d_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_dgemm_d_s(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_dgemm_s_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_dgemm_s_d(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_dgemm_s_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_dgemm_s_s(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_z_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_z_c(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_c_z";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_c_z(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_c_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_c_c(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_cgemm_c_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_cgemm_c_s(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_cgemm_s_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_cgemm_s_c(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_cgemm_s_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_cgemm_s_s(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_z_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_z_d(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_d_z";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_d_z(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_d_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_d_d(m, n, k, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_sgemm_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_sgemm_x(m, n, k, ntests, &seed, thresh, debug,
		    test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		    &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_dgemm_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_dgemm_x(m, n, k, ntests, &seed, thresh, debug,
		    test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		    &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_cgemm_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_cgemm_x(m, n, k, ntests, &seed, thresh, debug,
		    test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		    &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_x(m, n, k, ntests, &seed, thresh, debug,
		    test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		    &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_dgemm_d_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_dgemm_d_s_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_dgemm_s_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_dgemm_s_d_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_dgemm_s_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_dgemm_s_s_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_z_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_z_c_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_c_z_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_c_z_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_c_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_c_c_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_cgemm_c_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_cgemm_c_s_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_cgemm_s_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_cgemm_s_c_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_cgemm_s_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_cgemm_s_s_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_z_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_z_d_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_d_z_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_d_z_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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

  fname = "BLAS_zgemm_d_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    do_test_zgemm_d_d_x(m, n, k, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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
