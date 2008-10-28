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
#define SIDE_START 0
#define SIDE_END   1

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
#define LDB_START    0
#define LDB_END      2

/* 0 -- 2 */
#define LDC_START    0
#define LDC_END      2

/* 0 -- 2 */
#define PREC_START   0
#define PREC_END     2

/* 0 -- 1 */
#define RANDOMIZE_START 0
#define RANDOMIZE_END   1

#define NUM_DATA 7










void do_test_dsymm_d_s
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dsymm_d_s";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_dsymm_d_s_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			dge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_dsymm_d_s(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;



			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  dsy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      sge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      sge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin = c_gen[cij];
			    rout = c[cij];
			    head_r_true_elem = head_r_true[cij];
			    tail_r_true_elem = tail_r_true[cij];

			    test_BLAS_ddot_d_s(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");

			    printf("a\n");
			    dsy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    sge_print_matrix(b, m, n, ldb, order_type, "B");
			    dge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    dge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_dsymm_s_d
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dsymm_s_d";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_dsymm_s_d_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			dge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_dsymm_s_d(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;



			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  ssy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      dge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      dge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin = c_gen[cij];
			    rout = c[cij];
			    head_r_true_elem = head_r_true[cij];
			    tail_r_true_elem = tail_r_true[cij];

			    test_BLAS_ddot_s_d(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");

			    printf("a\n");
			    ssy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    dge_print_matrix(b, m, n, ldb, order_type, "B");
			    dge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    dge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_dsymm_s_s
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dsymm_s_s";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_dsymm_s_s_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			dge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_dsymm_s_s(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
			  incci = 1;
			  inccij = ldc;
			} else {
			  incci = ldc;
			  inccij = 1;
			}

			incri = incci;
			incrij = inccij;



			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  ssy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      sge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      sge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin = c_gen[cij];
			    rout = c[cij];
			    head_r_true_elem = head_r_true[cij];
			    tail_r_true_elem = tail_r_true[cij];

			    test_BLAS_ddot_s_s(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");

			    printf("a\n");
			    ssy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    sge_print_matrix(b, m, n, ldb, order_type, "B");
			    dge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    dge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_zsymm_z_c
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_z_c";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_zsymm_z_c_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_zsymm_z_c(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
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

			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  zsy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      cge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      cge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin[0] = c_gen[cij];
			    rin[1] = c_gen[cij + 1];
			    rout[0] = c[cij];
			    rout[1] = c[cij + 1];
			    head_r_true_elem[0] = head_r_true[cij];
			    head_r_true_elem[1] = head_r_true[cij + 1];
			    tail_r_true_elem[0] = tail_r_true[cij];
			    tail_r_true_elem[1] = tail_r_true[cij + 1];

			    test_BLAS_zdot_z_c(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

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
			    zsy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    cge_print_matrix(b, m, n, ldb, order_type, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_zsymm_c_z
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_c_z";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_zsymm_c_z_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_zsymm_c_z(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
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

			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  csy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      zge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      zge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin[0] = c_gen[cij];
			    rin[1] = c_gen[cij + 1];
			    rout[0] = c[cij];
			    rout[1] = c[cij + 1];
			    head_r_true_elem[0] = head_r_true[cij];
			    head_r_true_elem[1] = head_r_true[cij + 1];
			    tail_r_true_elem[0] = tail_r_true[cij];
			    tail_r_true_elem[1] = tail_r_true[cij + 1];

			    test_BLAS_zdot_c_z(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

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
			    csy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    zge_print_matrix(b, m, n, ldb, order_type, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_zsymm_c_c
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_c_c";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_zsymm_c_c_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_zsymm_c_c(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
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

			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  csy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      cge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      cge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin[0] = c_gen[cij];
			    rin[1] = c_gen[cij + 1];
			    rout[0] = c[cij];
			    rout[1] = c[cij + 1];
			    head_r_true_elem[0] = head_r_true[cij];
			    head_r_true_elem[1] = head_r_true[cij + 1];
			    tail_r_true_elem[0] = tail_r_true[cij];
			    tail_r_true_elem[1] = tail_r_true[cij + 1];

			    test_BLAS_zdot_c_c(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

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
			    csy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    cge_print_matrix(b, m, n, ldb, order_type, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_csymm_c_s
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_csymm_c_s";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_csymm_c_s_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			cge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_csymm_c_s(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
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

			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  csy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      sge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      sge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin[0] = c_gen[cij];
			    rin[1] = c_gen[cij + 1];
			    rout[0] = c[cij];
			    rout[1] = c[cij + 1];
			    head_r_true_elem[0] = head_r_true[cij];
			    head_r_true_elem[1] = head_r_true[cij + 1];
			    tail_r_true_elem[0] = tail_r_true[cij];
			    tail_r_true_elem[1] = tail_r_true[cij + 1];

			    test_BLAS_cdot_c_s(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

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
			    csy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    sge_print_matrix(b, m, n, ldb, order_type, "B");
			    cge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    cge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_csymm_s_c
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_csymm_s_c";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_csymm_s_c_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			cge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_csymm_s_c(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
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

			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  ssy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      cge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      cge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin[0] = c_gen[cij];
			    rin[1] = c_gen[cij + 1];
			    rout[0] = c[cij];
			    rout[1] = c[cij + 1];
			    head_r_true_elem[0] = head_r_true[cij];
			    head_r_true_elem[1] = head_r_true[cij + 1];
			    tail_r_true_elem[0] = tail_r_true[cij];
			    tail_r_true_elem[1] = tail_r_true[cij + 1];

			    test_BLAS_cdot_s_c(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

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
			    ssy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    cge_print_matrix(b, m, n, ldb, order_type, "B");
			    cge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    cge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_csymm_s_s
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_csymm_s_s";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_csymm_s_s_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			cge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_csymm_s_s(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
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

			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  ssy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      sge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      sge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin[0] = c_gen[cij];
			    rin[1] = c_gen[cij + 1];
			    rout[0] = c[cij];
			    rout[1] = c[cij + 1];
			    head_r_true_elem[0] = head_r_true[cij];
			    head_r_true_elem[1] = head_r_true[cij + 1];
			    tail_r_true_elem[0] = tail_r_true[cij];
			    tail_r_true_elem[1] = tail_r_true[cij + 1];

			    test_BLAS_cdot_s_s(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

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
			    ssy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    sge_print_matrix(b, m, n, ldb, order_type, "B");
			    cge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    cge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_zsymm_z_d
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_z_d";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_zsymm_z_d_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_zsymm_z_d(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
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

			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  zsy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      dge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      dge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin[0] = c_gen[cij];
			    rin[1] = c_gen[cij + 1];
			    rout[0] = c[cij];
			    rout[1] = c[cij + 1];
			    head_r_true_elem[0] = head_r_true[cij];
			    head_r_true_elem[1] = head_r_true[cij + 1];
			    tail_r_true_elem[0] = tail_r_true[cij];
			    tail_r_true_elem[1] = tail_r_true[cij + 1];

			    test_BLAS_zdot_z_d(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

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
			    zsy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    dge_print_matrix(b, m, n, ldb, order_type, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_zsymm_d_z
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_d_z";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_zsymm_d_z_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_zsymm_d_z(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
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

			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  dsy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      zge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      zge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin[0] = c_gen[cij];
			    rin[1] = c_gen[cij + 1];
			    rout[0] = c[cij];
			    rout[1] = c[cij + 1];
			    head_r_true_elem[0] = head_r_true[cij];
			    head_r_true_elem[1] = head_r_true[cij + 1];
			    tail_r_true_elem[0] = tail_r_true[cij];
			    tail_r_true_elem[1] = tail_r_true[cij + 1];

			    test_BLAS_zdot_d_z(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

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
			    dsy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    zge_print_matrix(b, m, n, ldb, order_type, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_zsymm_d_d
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_d_d";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && b_vec == NULL) {
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

	    order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	    /* vary upper / lower variation */
	    for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

	      uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

	      /* vary left / right multiplication */
	      for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		side_type = (side_val == 0) ?
		  blas_left_side : blas_right_side;

		if (side_type == blas_left_side) {
		  m_i = m;
		  n_i = n;
		} else {
		  m_i = n;
		  n_i = m;
		}

		/* vary lda = m_i, m_i+1, 2*m_i */
		for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : 2 * m_i;

		  /* vary ldb = n_i, n_i+1, 2*n_i */
		  for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		    if (order_type == blas_colmajor)
		      ldb = (ldb_val == 0) ? m :
			(ldb_val == 1) ? m + 1 : 2 * m;
		    else
		      ldb = (ldb_val == 0) ? n :
			(ldb_val == 1) ? n + 1 : 2 * n;

		    /* vary ldc = k, k+1, 2*k */
		    for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

		      if (order_type == blas_colmajor)
			ldc = (ldc_val == 0) ? m :
			  (ldc_val == 1) ? m + 1 : 2 * m;
		      else
			ldc = (ldc_val == 0) ? n :
			  (ldc_val == 1) ? n + 1 : 2 * n;

		      for (randomize_val = RANDOMIZE_START;
			   randomize_val <= RANDOMIZE_END; randomize_val++) {

			/* For the sake of speed, we throw out this case at random */
			if (xrand(seed) >= test_prob)
			  continue;

			saved_seed = *seed;

			/* finally we are here to generate the test case */
			BLAS_zsymm_d_d_testgen(norm, order_type,
					       uplo_type, side_type, m, n,
					       randomize_val, &alpha,
					       alpha_flag, &beta, beta_flag,
					       a, lda, b, ldb, c, ldc, seed,
					       head_r_true, tail_r_true);
			test_count++;

			/* copy generated C matrix since this will be
			   over written */
			zge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

			/* call symm routines to be tested */
			FPU_FIX_STOP;
			BLAS_zsymm_d_d(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc);
			FPU_FIX_START;

			/* now compute the ratio using test_c_xdot */
			/* copy a row from A, a column from B, run 
			   dot test */

			if ((order_type == blas_colmajor &&
			     side_type == blas_left_side) ||
			    (order_type == blas_rowmajor &&
			     side_type == blas_right_side)) {
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

			for (i = 0, ci = 0, ri = 0;
			     i < m_i; i++, ci += incci, ri += incri) {
			  dsy_copy_row(order_type, uplo_type, m_i, a, lda,
				       a_vec, i);
			  for (j = 0, cij = ci, rij = ri; j < n_i;
			       j++, cij += inccij, rij += incrij) {
			    /* copy i-th row of A and j-th col of B */
			    if (side_type == blas_left_side)
			      dge_copy_col(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    else
			      dge_copy_row(order_type, blas_no_trans,
					   m, n, b, ldb, b_vec, j);
			    rin[0] = c_gen[cij];
			    rin[1] = c_gen[cij + 1];
			    rout[0] = c[cij];
			    rout[1] = c[cij + 1];
			    head_r_true_elem[0] = head_r_true[cij];
			    head_r_true_elem[1] = head_r_true[cij + 1];
			    tail_r_true_elem[0] = tail_r_true[cij];
			    tail_r_true_elem[1] = tail_r_true[cij + 1];

			    test_BLAS_zdot_d_d(m_i,
					       blas_no_conj,
					       alpha, beta, rin, rout,
					       head_r_true_elem,
					       tail_r_true_elem, a_vec, 1,
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
			    printf("Seed = %d\n", saved_seed);
			    printf("m %d   n %d\n", m, n);
			    printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    if (uplo_type == blas_upper)
			      printf("upper ");
			    else
			      printf("lower");

			    if (side_type == blas_left_side)
			      printf(" left\n");
			    else
			      printf(" right\n");

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
			    dsy_print_matrix(a, m_i, lda, order_type,
					     uplo_type);
			    dge_print_matrix(b, m, n, ldb, order_type, "B");
			    zge_print_matrix(c_gen, m, n, ldc, order_type,
					     "C_gen");
			    zge_print_matrix(c, m, n, ldc, order_type, "C");

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

		      }		/* end of randmize loop */

		    }		/* end of ldc loop */

		  }		/* end of ldb loop */

		}		/* end of lda loop */

	      }			/* end of side loop */

	    }			/* end of uplo loop */

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
void do_test_ssymm_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_ssymm_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin;
  float rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_ssymm_testgen(norm, order_type,
					     uplo_type, side_type, m, n,
					     randomize_val, &alpha,
					     alpha_flag, &beta, beta_flag, a,
					     lda, b, ldb, c, ldc, seed,
					     head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  sge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_ssymm_x(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;



			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    ssy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				sge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				sge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin = c_gen[cij];
			      rout = c[cij];
			      head_r_true_elem = head_r_true[cij];
			      tail_r_true_elem = tail_r_true[cij];

			      test_BLAS_sdot(m_i,
					     blas_no_conj,
					     alpha, beta, rin, rout,
					     head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("%16.8e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%16.8e", beta);;
			      printf("\n");

			      printf("a\n");
			      ssy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      sge_print_matrix(b, m, n, ldb, order_type, "B");
			      sge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      sge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_dsymm_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dsymm_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_dsymm_testgen(norm, order_type,
					     uplo_type, side_type, m, n,
					     randomize_val, &alpha,
					     alpha_flag, &beta, beta_flag, a,
					     lda, b, ldb, c, ldc, seed,
					     head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  dge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_dsymm_x(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;



			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    dsy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				dge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				dge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin = c_gen[cij];
			      rout = c[cij];
			      head_r_true_elem = head_r_true[cij];
			      tail_r_true_elem = tail_r_true[cij];

			      test_BLAS_ddot(m_i,
					     blas_no_conj,
					     alpha, beta, rin, rout,
					     head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("%24.16e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%24.16e", beta);;
			      printf("\n");

			      printf("a\n");
			      dsy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      dge_print_matrix(b, m, n, ldb, order_type, "B");
			      dge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      dge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_csymm_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_csymm_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_csymm_testgen(norm, order_type,
					     uplo_type, side_type, m, n,
					     randomize_val, &alpha,
					     alpha_flag, &beta, beta_flag, a,
					     lda, b, ldb, c, ldc, seed,
					     head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  cge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_csymm_x(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    csy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				cge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				cge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_cdot(m_i,
					     blas_no_conj,
					     alpha, beta, rin, rout,
					     head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

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
			      csy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      cge_print_matrix(b, m, n, ldb, order_type, "B");
			      cge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      cge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_zsymm_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_zsymm_testgen(norm, order_type,
					     uplo_type, side_type, m, n,
					     randomize_val, &alpha,
					     alpha_flag, &beta, beta_flag, a,
					     lda, b, ldb, c, ldc, seed,
					     head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zsymm_x(order_type, side_type,
				       uplo_type, m, n, alpha, a, lda, b, ldb,
				       beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    zsy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				zge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				zge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_zdot(m_i,
					     blas_no_conj,
					     alpha, beta, rin, rout,
					     head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      printf("a\n");
			      zsy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      zge_print_matrix(b, m, n, ldb, order_type, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_dsymm_d_s_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dsymm_d_s_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_dsymm_d_s_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  dge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_dsymm_d_s_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;



			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    dsy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				sge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				sge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin = c_gen[cij];
			      rout = c[cij];
			      head_r_true_elem = head_r_true[cij];
			      tail_r_true_elem = tail_r_true[cij];

			      test_BLAS_ddot_d_s(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("%24.16e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%24.16e", beta);;
			      printf("\n");

			      printf("a\n");
			      dsy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      sge_print_matrix(b, m, n, ldb, order_type, "B");
			      dge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      dge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_dsymm_s_d_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dsymm_s_d_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_dsymm_s_d_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  dge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_dsymm_s_d_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;



			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    ssy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				dge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				dge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin = c_gen[cij];
			      rout = c[cij];
			      head_r_true_elem = head_r_true[cij];
			      tail_r_true_elem = tail_r_true[cij];

			      test_BLAS_ddot_s_d(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("%24.16e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%24.16e", beta);;
			      printf("\n");

			      printf("a\n");
			      ssy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      dge_print_matrix(b, m, n, ldb, order_type, "B");
			      dge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      dge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_dsymm_s_s_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dsymm_s_s_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_dsymm_s_s_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  dge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_dsymm_s_s_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
			    incci = 1;
			    inccij = ldc;
			  } else {
			    incci = ldc;
			    inccij = 1;
			  }

			  incri = incci;
			  incrij = inccij;



			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    ssy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				sge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				sge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin = c_gen[cij];
			      rout = c[cij];
			      head_r_true_elem = head_r_true[cij];
			      tail_r_true_elem = tail_r_true[cij];

			      test_BLAS_ddot_s_s(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("%24.16e", alpha);;
			      printf("   ");
			      printf("beta = ");
			      printf("%24.16e", beta);;
			      printf("\n");

			      printf("a\n");
			      ssy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      sge_print_matrix(b, m, n, ldb, order_type, "B");
			      dge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      dge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_zsymm_z_c_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_z_c_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_zsymm_z_c_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zsymm_z_c_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    zsy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				cge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				cge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_zdot_z_c(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      printf("a\n");
			      zsy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      cge_print_matrix(b, m, n, ldb, order_type, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_zsymm_c_z_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_c_z_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_zsymm_c_z_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zsymm_c_z_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    csy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				zge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				zge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_zdot_c_z(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      printf("a\n");
			      csy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      zge_print_matrix(b, m, n, ldb, order_type, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_zsymm_c_c_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_c_c_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_zsymm_c_c_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zsymm_c_c_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    csy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				cge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				cge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_zdot_c_c(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      printf("a\n");
			      csy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      cge_print_matrix(b, m, n, ldb, order_type, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_csymm_c_s_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_csymm_c_s_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_csymm_c_s_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  cge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_csymm_c_s_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    csy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				sge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				sge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_cdot_c_s(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

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
			      csy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      sge_print_matrix(b, m, n, ldb, order_type, "B");
			      cge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      cge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_csymm_s_c_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_csymm_s_c_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_csymm_s_c_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  cge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_csymm_s_c_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    ssy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				cge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				cge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_cdot_s_c(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

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
			      ssy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      cge_print_matrix(b, m, n, ldb, order_type, "B");
			      cge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      cge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_csymm_s_s_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_csymm_s_s_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (float *) blas_malloc(2 * max_mn * max_mn * sizeof(float));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (float *) blas_malloc(2 * m * n * sizeof(float));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_csymm_s_s_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  cge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_csymm_s_s_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    ssy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				sge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				sge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_cdot_s_s(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

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
			      ssy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      sge_print_matrix(b, m, n, ldb, order_type, "B");
			      cge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      cge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_zsymm_z_d_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_z_d_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double) * 2);
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_zsymm_z_d_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zsymm_z_d_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    zsy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				dge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				dge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_zdot_z_d(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      printf("a\n");
			      zsy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      dge_print_matrix(b, m, n, ldb, order_type, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_zsymm_d_z_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_d_z_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double) * 2);
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_zsymm_d_z_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zsymm_d_z_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    dsy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				zge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				zge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_zdot_d_z(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      printf("a\n");
			      dsy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      zge_print_matrix(b, m, n, ldb, order_type, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
void do_test_zsymm_d_d_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zsymm_d_d_x";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;

  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

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

  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;

  max_mn = (m > n) ? m : n;

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
  a = (double *) blas_malloc(2 * max_mn * max_mn * sizeof(double));
  if (2 * max_mn * max_mn > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b = (double *) blas_malloc(2 * m * n * sizeof(double));
  if (2 * m * n > 0 && b == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  b_vec = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && b_vec == NULL) {
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

	      order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

	      /* vary upper / lower variation */
	      for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

		uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;

		/* vary left / right multiplication */
		for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

		  side_type = (side_val == 0) ?
		    blas_left_side : blas_right_side;

		  if (side_type == blas_left_side) {
		    m_i = m;
		    n_i = n;
		  } else {
		    m_i = n;
		    n_i = m;
		  }

		  /* vary lda = m_i, m_i+1, 2*m_i */
		  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		    lda = (lda_val == 0) ? m_i :
		      (lda_val == 1) ? m_i + 1 : 2 * m_i;

		    /* vary ldb = n_i, n_i+1, 2*n_i */
		    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {

		      if (order_type == blas_colmajor)
			ldb = (ldb_val == 0) ? m :
			  (ldb_val == 1) ? m + 1 : 2 * m;
		      else
			ldb = (ldb_val == 0) ? n :
			  (ldb_val == 1) ? n + 1 : 2 * n;

		      /* vary ldc = k, k+1, 2*k */
		      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {

			if (order_type == blas_colmajor)
			  ldc = (ldc_val == 0) ? m :
			    (ldc_val == 1) ? m + 1 : 2 * m;
			else
			  ldc = (ldc_val == 0) ? n :
			    (ldc_val == 1) ? n + 1 : 2 * n;

			for (randomize_val = RANDOMIZE_START;
			     randomize_val <= RANDOMIZE_END;
			     randomize_val++) {

			  /* For the sake of speed, we throw out this case at random */
			  if (xrand(seed) >= test_prob)
			    continue;

			  saved_seed = *seed;

			  /* finally we are here to generate the test case */
			  BLAS_zsymm_d_d_testgen(norm, order_type,
						 uplo_type, side_type, m, n,
						 randomize_val, &alpha,
						 alpha_flag, &beta, beta_flag,
						 a, lda, b, ldb, c, ldc, seed,
						 head_r_true, tail_r_true);
			  test_count++;

			  /* copy generated C matrix since this will be
			     over written */
			  zge_copy_matrix(order_type, m, n, c_gen, ldc, c,
					  ldc);

			  /* call symm routines to be tested */
			  FPU_FIX_STOP;
			  BLAS_zsymm_d_d_x(order_type, side_type,
					   uplo_type, m, n, alpha, a, lda, b,
					   ldb, beta, c, ldc, prec);
			  FPU_FIX_START;

			  /* now compute the ratio using test_c_xdot */
			  /* copy a row from A, a column from B, run 
			     dot test */

			  if ((order_type == blas_colmajor &&
			       side_type == blas_left_side) ||
			      (order_type == blas_rowmajor &&
			       side_type == blas_right_side)) {
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

			  for (i = 0, ci = 0, ri = 0;
			       i < m_i; i++, ci += incci, ri += incri) {
			    dsy_copy_row(order_type, uplo_type, m_i, a, lda,
					 a_vec, i);
			    for (j = 0, cij = ci, rij = ri; j < n_i;
				 j++, cij += inccij, rij += incrij) {
			      /* copy i-th row of A and j-th col of B */
			      if (side_type == blas_left_side)
				dge_copy_col(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      else
				dge_copy_row(order_type, blas_no_trans,
					     m, n, b, ldb, b_vec, j);
			      rin[0] = c_gen[cij];
			      rin[1] = c_gen[cij + 1];
			      rout[0] = c[cij];
			      rout[1] = c[cij + 1];
			      head_r_true_elem[0] = head_r_true[cij];
			      head_r_true_elem[1] = head_r_true[cij + 1];
			      tail_r_true_elem[0] = tail_r_true[cij];
			      tail_r_true_elem[1] = tail_r_true[cij + 1];

			      test_BLAS_zdot_d_d(m_i,
						 blas_no_conj,
						 alpha, beta, rin, rout,
						 head_r_true_elem,
						 tail_r_true_elem, a_vec, 1,
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
			      printf("Seed = %d\n", saved_seed);
			      printf("m %d   n %d\n", m, n);
			      printf("LDA %d  LDB %d  LDC %d\n", lda, ldb,
				     ldc);

			      if (order_type == blas_rowmajor)
				printf("row ");
			      else
				printf("col ");

			      if (uplo_type == blas_upper)
				printf("upper ");
			      else
				printf("lower");

			      if (side_type == blas_left_side)
				printf(" left\n");
			      else
				printf(" right\n");

			      printf("NORM %d, ALPHA %d, BETA %d\n",
				     norm, alpha_val, beta_val);

			      /* print out info */
			      printf("alpha = ");
			      printf("(%24.16e, %24.16e)", alpha[0],
				     alpha[1]);;
			      printf("   ");
			      printf("beta = ");
			      printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			      printf("\n");

			      printf("a\n");
			      dsy_print_matrix(a, m_i, lda, order_type,
					       uplo_type);
			      dge_print_matrix(b, m, n, ldb, order_type, "B");
			      zge_print_matrix(c_gen, m, n, ldc, order_type,
					       "C_gen");
			      zge_print_matrix(c, m, n, ldc, order_type, "C");

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

			}	/* end of randmize loop */

		      }		/* end of ldc loop */

		    }		/* end of ldb loop */

		  }		/* end of lda loop */

		}		/* end of side loop */

	      }			/* end of uplo loop */

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
  const char *base_routine = "symm";
  char *fname;
  int n;

  int m, i;
  int mn_data[NUM_DATA][2] =
    { {4, 4}, {2, 3}, {4, 2}, {8, 8}, {10, 1}, {1, 1}, {1, 7} };

  if (argc != 6) {
    printf("Usage:\n");
    printf("do_test_symm <nsizes> <ntests> <thresh> <debug> <test_prob>\n");
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
    BLAS_error("Testing symm", 0, 0, NULL);

  printf("Testing %s...\n", base_routine);
  printf("INPUT: nsizes = %d, ntests = %d, thresh = %4.2f, debug = %d\n\n",
	 nsizes, ntests, thresh, debug);





  fname = "BLAS_dsymm_d_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_dsymm_d_s(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_dsymm_s_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_dsymm_s_d(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_dsymm_s_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_dsymm_s_s(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_z_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_z_c(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_c_z";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_c_z(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_c_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_c_c(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_csymm_c_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_csymm_c_s(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_csymm_s_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_csymm_s_c(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_csymm_s_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_csymm_s_s(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_z_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_z_d(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_d_z";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_d_z(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_d_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_d_d(m, n, ntests, &seed, thresh, debug,
		      test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		      &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_ssymm_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_ssymm_x(m, n, ntests, &seed, thresh, debug,
		    test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		    &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_dsymm_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_dsymm_x(m, n, ntests, &seed, thresh, debug,
		    test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		    &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_csymm_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_csymm_x(m, n, ntests, &seed, thresh, debug,
		    test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		    &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_x(m, n, ntests, &seed, thresh, debug,
		    test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
		    &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_dsymm_d_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_dsymm_d_s_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_dsymm_s_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_dsymm_s_d_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_dsymm_s_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_dsymm_s_s_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_z_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_z_c_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_c_z_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_c_z_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_c_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_c_c_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_csymm_c_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_csymm_c_s_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_csymm_s_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_csymm_s_c_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_csymm_s_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_csymm_s_s_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_z_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_z_d_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_d_z_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_d_z_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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

  fname = "BLAS_zsymm_d_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    do_test_zsymm_d_d_x(m, n, ntests, &seed, thresh, debug,
			test_prob, &min_ratio, &max_ratio, &num_bad_ratio,
			&num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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
