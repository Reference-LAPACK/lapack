#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

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










void do_test_dge_sum_mv_d_s
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dge_sum_mv_d_s";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha;
  double beta;
  double beta_zero_fake;
  double alpha_use;
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  float *x;
  double *y;
  double *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  beta_zero_fake = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;




  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double));
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double));
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_dge_sum_mv_d_s_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    scopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_dge_sum_mv_d_s(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;



		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      scopy_vector(x, n_i, incx, x_vec, 1);
		      scopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin = 0.0;
			rout = y[yi];
			head_r_true_elem = head_r_true[ri];
			tail_r_true_elem = tail_r_true[ri];

			test_BLAS_ddot_d_s(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : d, a type : d, x type : s\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("%24.16e", alpha);;
			  printf("   ");
			  printf("beta = ");
			  printf("%24.16e", beta);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("%24.16e", alpha_use);;
			  printf("\n");

			  dge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  dge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  sprint_vector(x, n_i, incx, "x");

			  dprint_vector(y, m_i, incy, "y");

			  dprint_vector(head_r_true, m_i, 1, "head_r_true");

			  dge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  dge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dge_sum_mv_s_d
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dge_sum_mv_s_d";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha;
  double beta;
  double beta_zero_fake;
  double alpha_use;
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  double *x;
  double *y;
  float *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  beta_zero_fake = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;




  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double));
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double));
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_dge_sum_mv_s_d_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    dcopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_dge_sum_mv_s_d(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;



		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      dcopy_vector(x, n_i, incx, x_vec, 1);
		      dcopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin = 0.0;
			rout = y[yi];
			head_r_true_elem = head_r_true[ri];
			tail_r_true_elem = tail_r_true[ri];

			test_BLAS_ddot_s_d(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : d, a type : s, x type : d\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("%24.16e", alpha);;
			  printf("   ");
			  printf("beta = ");
			  printf("%24.16e", beta);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("%24.16e", alpha_use);;
			  printf("\n");

			  sge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  sge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  dprint_vector(x, n_i, incx, "x");

			  dprint_vector(y, m_i, incy, "y");

			  dprint_vector(head_r_true, m_i, 1, "head_r_true");

			  sge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  sge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dge_sum_mv_s_s
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dge_sum_mv_s_s";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha;
  double beta;
  double beta_zero_fake;
  double alpha_use;
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  double *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  beta_zero_fake = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;




  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double));
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double));
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_dge_sum_mv_s_s_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    scopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_dge_sum_mv_s_s(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;



		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      scopy_vector(x, n_i, incx, x_vec, 1);
		      scopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin = 0.0;
			rout = y[yi];
			head_r_true_elem = head_r_true[ri];
			tail_r_true_elem = tail_r_true[ri];

			test_BLAS_ddot_s_s(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : d, a type : s, x type : s\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("%24.16e", alpha);;
			  printf("   ");
			  printf("beta = ");
			  printf("%24.16e", beta);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("%24.16e", alpha_use);;
			  printf("\n");

			  sge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  sge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  sprint_vector(x, n_i, incx, "x");

			  dprint_vector(y, m_i, incy, "y");

			  dprint_vector(head_r_true, m_i, 1, "head_r_true");

			  sge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  sge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_z_c
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_z_c";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  float *x;
  double *y;
  double *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;
  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_zge_sum_mv_z_c_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    ccopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zge_sum_mv_z_c(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;
		      incx_veci *= 2;
		      incyi *= 2;
		      incri *= 2;
		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      ccopy_vector(x, n_i, incx, x_vec, 1);
		      ccopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin[0] = rin[1] = 0.0;
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot_z_c(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : z, x type : c\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("(%24.16e, %24.16e)", alpha_use[0],
				 alpha_use[1]);;
			  printf("\n");

			  zge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  zge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  cprint_vector(x, n_i, incx, "x");

			  zprint_vector(y, m_i, incy, "y");

			  zprint_vector(head_r_true, m_i, 1, "head_r_true");

			  zge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  zge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_c_z
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_c_z";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  double *x;
  double *y;
  float *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;
  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_zge_sum_mv_c_z_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    zcopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zge_sum_mv_c_z(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;
		      incx_veci *= 2;
		      incyi *= 2;
		      incri *= 2;
		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      zcopy_vector(x, n_i, incx, x_vec, 1);
		      zcopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin[0] = rin[1] = 0.0;
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot_c_z(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : c, x type : z\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("(%24.16e, %24.16e)", alpha_use[0],
				 alpha_use[1]);;
			  printf("\n");

			  cge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  cge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  zprint_vector(x, n_i, incx, "x");

			  zprint_vector(y, m_i, incy, "y");

			  zprint_vector(head_r_true, m_i, 1, "head_r_true");

			  cge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  cge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_c_c
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_c_c";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  double *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;
  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_zge_sum_mv_c_c_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    ccopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zge_sum_mv_c_c(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;
		      incx_veci *= 2;
		      incyi *= 2;
		      incri *= 2;
		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      ccopy_vector(x, n_i, incx, x_vec, 1);
		      ccopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin[0] = rin[1] = 0.0;
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot_c_c(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : c, x type : c\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("(%24.16e, %24.16e)", alpha_use[0],
				 alpha_use[1]);;
			  printf("\n");

			  cge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  cge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  cprint_vector(x, n_i, incx, "x");

			  zprint_vector(y, m_i, incy, "y");

			  zprint_vector(head_r_true, m_i, 1, "head_r_true");

			  cge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  cge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cge_sum_mv_c_s
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_cge_sum_mv_c_s";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  float alpha[2];
  float beta[2];
  float beta_zero_fake[2];
  float alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  float *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;

  incy *= 2;

  /* allocate memory for arrays */
  y = (float *) blas_malloc(4 * m_i * sizeof(float) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_cge_sum_mv_c_s_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    scopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_cge_sum_mv_c_s(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;

		      incyi *= 2;
		      incri *= 2;
		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      scopy_vector(x, n_i, incx, x_vec, 1);
		      scopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin[0] = rin[1] = 0.0;
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_cdot_c_s(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : c, a type : c, x type : s\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("(%16.8e, %16.8e)", alpha_use[0],
				 alpha_use[1]);;
			  printf("\n");

			  cge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  cge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  sprint_vector(x, n_i, incx, "x");

			  cprint_vector(y, m_i, incy, "y");

			  zprint_vector(head_r_true, m_i, 1, "head_r_true");

			  cge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  cge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cge_sum_mv_s_c
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_cge_sum_mv_s_c";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  float alpha[2];
  float beta[2];
  float beta_zero_fake[2];
  float alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  float *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;

  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (float *) blas_malloc(4 * m_i * sizeof(float) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_cge_sum_mv_s_c_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    ccopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_cge_sum_mv_s_c(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;
		      incx_veci *= 2;
		      incyi *= 2;
		      incri *= 2;
		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      ccopy_vector(x, n_i, incx, x_vec, 1);
		      ccopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin[0] = rin[1] = 0.0;
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_cdot_s_c(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : c, a type : s, x type : c\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("(%16.8e, %16.8e)", alpha_use[0],
				 alpha_use[1]);;
			  printf("\n");

			  sge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  sge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  cprint_vector(x, n_i, incx, "x");

			  cprint_vector(y, m_i, incy, "y");

			  zprint_vector(head_r_true, m_i, 1, "head_r_true");

			  sge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  sge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cge_sum_mv_s_s
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_cge_sum_mv_s_s";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  float alpha[2];
  float beta[2];
  float beta_zero_fake[2];
  float alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  float *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;


  incy *= 2;

  /* allocate memory for arrays */
  y = (float *) blas_malloc(4 * m_i * sizeof(float) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_cge_sum_mv_s_s_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    scopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_cge_sum_mv_s_s(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;

		      incyi *= 2;
		      incri *= 2;
		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      scopy_vector(x, n_i, incx, x_vec, 1);
		      scopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin[0] = rin[1] = 0.0;
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_cdot_s_s(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : c, a type : s, x type : s\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("(%16.8e, %16.8e)", alpha_use[0],
				 alpha_use[1]);;
			  printf("\n");

			  sge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  sge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  sprint_vector(x, n_i, incx, "x");

			  cprint_vector(y, m_i, incy, "y");

			  zprint_vector(head_r_true, m_i, 1, "head_r_true");

			  sge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  sge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_z_d
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_z_d";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  double *x;
  double *y;
  double *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;

  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_zge_sum_mv_z_d_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    dcopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zge_sum_mv_z_d(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;

		      incyi *= 2;
		      incri *= 2;
		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      dcopy_vector(x, n_i, incx, x_vec, 1);
		      dcopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin[0] = rin[1] = 0.0;
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot_z_d(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : z, x type : d\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("(%24.16e, %24.16e)", alpha_use[0],
				 alpha_use[1]);;
			  printf("\n");

			  zge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  zge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  dprint_vector(x, n_i, incx, "x");

			  zprint_vector(y, m_i, incy, "y");

			  zprint_vector(head_r_true, m_i, 1, "head_r_true");

			  zge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  zge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_d_z
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_d_z";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  double *x;
  double *y;
  double *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;

  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_zge_sum_mv_d_z_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    zcopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zge_sum_mv_d_z(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;
		      incx_veci *= 2;
		      incyi *= 2;
		      incri *= 2;
		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      zcopy_vector(x, n_i, incx, x_vec, 1);
		      zcopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin[0] = rin[1] = 0.0;
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot_d_z(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : d, x type : z\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("(%24.16e, %24.16e)", alpha_use[0],
				 alpha_use[1]);;
			  printf("\n");

			  dge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  dge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  zprint_vector(x, n_i, incx, "x");

			  zprint_vector(y, m_i, incy, "y");

			  zprint_vector(head_r_true, m_i, 1, "head_r_true");

			  dge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  dge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_d_d
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_d_d";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;



  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  double *x;
  double *y;
  double *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;


  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	    /* vary lda = n_i, n_i+1, 2*n_i */
	    for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

	      if (order_type == blas_rowmajor) {
		lda = (lda_val == 0) ? n_i :
		  (lda_val == 1) ? n_i + 1 : n_i * n_i;
	      } else {
		lda = (lda_val == 0) ? m_i :
		  (lda_val == 1) ? m_i + 1 : m_i * m_i;
	      }

	      /* vary ldb = n_i, n_i+1, 2*n_i */
	      for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		if (order_type == blas_rowmajor) {
		  ldb = (ldb_val == 0) ? n_i :
		    (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  ldb = (ldb_val == 0) ? m_i :
		    (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		}

		for (randomize_val = RANDOMIZE_START;
		     randomize_val <= RANDOMIZE_END; randomize_val++) {

		  /* For the sake of speed, we throw out this case at random */
		  if (xrand(seed) >= test_prob)
		    continue;

		  /* finally we are here to generate the test case */
		  /* alpha_use, a_use, B_use are the generated alpha, a, B
		   *  before any scaling.  
		   *  That is, in the generator, alpha == beta == alpha_use 
		   *  before scaling. */

		  saved_seed = *seed;
		  BLAS_zge_sum_mv_d_d_testgen(norm, order_type,
					      m, n, randomize_val, &alpha,
					      alpha_flag, &beta, beta_flag, a,
					      lda, B, ldb, x_vec, 1,
					      &alpha_use, a_use, B_use, seed,
					      head_r_true, tail_r_true);

		  /* vary incx = 1, 2 */
		  for (incx_val = INCX_START; incx_val <= INCX_END;
		       incx_val++) {

		    incx = incx_val;
		    if (0 == incx)
		      continue;

		    dcopy_vector(x_vec, n_i, 1, x, incx);

		    /* vary incy = 1, 2 */
		    for (incy_val = INCY_START; incy_val <= INCY_END;
			 incy_val++) {

		      incy = incy_val;
		      if (0 == incy)
			continue;

		      test_count++;

		      /* call ge_sum_mv routines to be tested */
		      FPU_FIX_STOP;
		      BLAS_zge_sum_mv_d_d(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy);
		      FPU_FIX_START;

		      /* now compute the ratio using test_BLAS_xdot */
		      /* copy a row from A, use x, run 
		         dot test */

		      incyi = incy;

		      incri = 1;
		      incx_veci = 1;

		      incyi *= 2;
		      incri *= 2;
		      if (incy < 0) {
			y_starti = (-m_i + 1) * incyi;
		      } else {
			y_starti = 0;
		      }
		      /* make two copies of x into x_vec. redundant */
		      dcopy_vector(x, n_i, incx, x_vec, 1);
		      dcopy_vector(x, n_i, incx, (x_vec + (n_i * incx_veci)),
				   1);
		      for (i = 0, yi = y_starti, ri = 0; i < m_i;
			   i++, yi += incyi, ri += incri) {
			dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     a_use, lda, a_vec, i);
			dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				     B_use, ldb, (a_vec + inca_veci * n_i),
				     i);

			rin[0] = rin[1] = 0.0;
			rout[0] = y[yi];
			rout[1] = y[yi + 1];
			head_r_true_elem[0] = head_r_true[ri];
			head_r_true_elem[1] = head_r_true[ri + 1];
			tail_r_true_elem[0] = tail_r_true[ri];
			tail_r_true_elem[1] = tail_r_true[ri + 1];

			test_BLAS_zdot_d_d(2 * n_i,
					   blas_no_conj,
					   alpha_use, beta_zero_fake, rin,
					   rout, head_r_true_elem,
					   tail_r_true_elem, a_vec, 1, x_vec,
					   1, eps_int, un_int, &ratios[i]);

			/* take the max ratio */
			if (i == 0) {
			  ratio = ratios[0];
			  /* The !<= below causes NaN errors
			   *  to be included.
			   * Note that (NaN > 0) is false */
			} else if (!(ratios[i] <= ratio)) {
			  ratio = ratios[i];
			}
		      }		/* end of dot-test loop */

		      /* The !<= below causes NaN errors
		       *  to be included.
		       * Note that (NaN > 0) is false */
		      if (!(ratio <= thresh)) {

			if (debug == 3) {
			  printf("\n\t\tTest # %d\n", test_count);
			  printf("y type : z, a type : d, x type : d\n");
			  printf("Seed = %d\t", saved_seed);
			  printf("n %d, m %d\n", n, m);
			  printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				 ldb, incx, incx);

			  if (order_type == blas_rowmajor)
			    printf("row ");
			  else
			    printf("col ");

			  printf("NORM %d, ALPHA %d, BETA %d\n",
				 norm, alpha_val, beta_val);
			  printf("randomize %d\n", randomize_val);

			  /* print out info */
			  printf("alpha = ");
			  printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			  printf("   ");
			  printf("beta = ");
			  printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			  printf("\n");
			  printf("alpha_use = ");
			  printf("(%24.16e, %24.16e)", alpha_use[0],
				 alpha_use[1]);;
			  printf("\n");

			  dge_print_matrix(a, m_i, n_i, lda, order_type, "A");
			  dge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
			  dprint_vector(x, n_i, incx, "x");

			  zprint_vector(y, m_i, incy, "y");

			  zprint_vector(head_r_true, m_i, 1, "head_r_true");

			  dge_print_matrix(a_use, m_i, n_i, lda, order_type,
					   "A_use");
			  dge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					   "B_use");

			  dprint_vector(ratios, m_i, 1, "ratios");
			  printf("ratio = %g\n", ratio);
			  fflush(stdout);
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

		}		/* end of randmize loop */

	      }			/* end of ldb loop */

	    }			/* end of lda loop */

	  }			/* end of order loop */

	}			/* end of nr test loop */

      }				/* end of norm loop */



    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_sge_sum_mv_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_sge_sum_mv_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin;
  float rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  float alpha;
  float beta;
  float beta_zero_fake;
  float alpha_use;
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  float *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  beta_zero_fake = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;




  /* allocate memory for arrays */
  y = (float *) blas_malloc(4 * m_i * sizeof(float));
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double));
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_sge_sum_mv_testgen(norm, order_type,
					    m, n, randomize_val, &alpha,
					    alpha_flag, &beta, beta_flag, a,
					    lda, B, ldb, x_vec, 1, &alpha_use,
					    a_use, B_use, seed, head_r_true,
					    tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      scopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_sge_sum_mv_x(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;



			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			scopy_vector(x, n_i, incx, x_vec, 1);
			scopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin = 0.0;
			  rout = y[yi];
			  head_r_true_elem = head_r_true[ri];
			  tail_r_true_elem = tail_r_true[ri];

			  test_BLAS_sdot(2 * n_i,
					 blas_no_conj,
					 alpha_use, beta_zero_fake, rin, rout,
					 head_r_true_elem, tail_r_true_elem,
					 a_vec, 1, x_vec, 1, eps_int, un_int,
					 &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : s, a type : s, x type : s\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("%16.8e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%16.8e", beta);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("%16.8e", alpha_use);;
			    printf("\n");

			    sge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    sge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    sprint_vector(x, n_i, incx, "x");

			    sprint_vector(y, m_i, incy, "y");

			    dprint_vector(head_r_true, m_i, 1, "head_r_true");

			    sge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    sge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dge_sum_mv_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dge_sum_mv_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha;
  double beta;
  double beta_zero_fake;
  double alpha_use;
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  double *x;
  double *y;
  double *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  beta_zero_fake = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;




  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double));
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double));
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_dge_sum_mv_testgen(norm, order_type,
					    m, n, randomize_val, &alpha,
					    alpha_flag, &beta, beta_flag, a,
					    lda, B, ldb, x_vec, 1, &alpha_use,
					    a_use, B_use, seed, head_r_true,
					    tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      dcopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_dge_sum_mv_x(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;



			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			dcopy_vector(x, n_i, incx, x_vec, 1);
			dcopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin = 0.0;
			  rout = y[yi];
			  head_r_true_elem = head_r_true[ri];
			  tail_r_true_elem = tail_r_true[ri];

			  test_BLAS_ddot(2 * n_i,
					 blas_no_conj,
					 alpha_use, beta_zero_fake, rin, rout,
					 head_r_true_elem, tail_r_true_elem,
					 a_vec, 1, x_vec, 1, eps_int, un_int,
					 &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : d, a type : d, x type : d\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("%24.16e", alpha_use);;
			    printf("\n");

			    dge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    dge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    dprint_vector(x, n_i, incx, "x");

			    dprint_vector(y, m_i, incy, "y");

			    dprint_vector(head_r_true, m_i, 1, "head_r_true");

			    dge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    dge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cge_sum_mv_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_cge_sum_mv_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  float alpha[2];
  float beta[2];
  float beta_zero_fake[2];
  float alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  float *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;
  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (float *) blas_malloc(4 * m_i * sizeof(float) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_cge_sum_mv_testgen(norm, order_type,
					    m, n, randomize_val, &alpha,
					    alpha_flag, &beta, beta_flag, a,
					    lda, B, ldb, x_vec, 1, &alpha_use,
					    a_use, B_use, seed, head_r_true,
					    tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      ccopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_cge_sum_mv_x(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;
			incx_veci *= 2;
			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			ccopy_vector(x, n_i, incx, x_vec, 1);
			ccopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_cdot(2 * n_i,
					 blas_no_conj,
					 alpha_use, beta_zero_fake, rin, rout,
					 head_r_true_elem, tail_r_true_elem,
					 a_vec, 1, x_vec, 1, eps_int, un_int,
					 &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : c, a type : c, x type : c\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%16.8e, %16.8e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    cge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    cge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    cprint_vector(x, n_i, incx, "x");

			    cprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    cge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    cge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  double *x;
  double *y;
  double *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;
  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_zge_sum_mv_testgen(norm, order_type,
					    m, n, randomize_val, &alpha,
					    alpha_flag, &beta, beta_flag, a,
					    lda, B, ldb, x_vec, 1, &alpha_use,
					    a_use, B_use, seed, head_r_true,
					    tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      zcopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_zge_sum_mv_x(order_type,
					  m, n, alpha, a, lda, x, incx, beta,
					  B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;
			incx_veci *= 2;
			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			zcopy_vector(x, n_i, incx, x_vec, 1);
			zcopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_zdot(2 * n_i,
					 blas_no_conj,
					 alpha_use, beta_zero_fake, rin, rout,
					 head_r_true_elem, tail_r_true_elem,
					 a_vec, 1, x_vec, 1, eps_int, un_int,
					 &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : z, a type : z, x type : z\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%24.16e, %24.16e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    zge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    zge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    zprint_vector(x, n_i, incx, "x");

			    zprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    zge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    zge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dge_sum_mv_d_s_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dge_sum_mv_d_s_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha;
  double beta;
  double beta_zero_fake;
  double alpha_use;
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  float *x;
  double *y;
  double *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  beta_zero_fake = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;




  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double));
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double));
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_dge_sum_mv_d_s_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      scopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_dge_sum_mv_d_s_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;



			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			scopy_vector(x, n_i, incx, x_vec, 1);
			scopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin = 0.0;
			  rout = y[yi];
			  head_r_true_elem = head_r_true[ri];
			  tail_r_true_elem = tail_r_true[ri];

			  test_BLAS_ddot_d_s(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : d, a type : d, x type : s\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("%24.16e", alpha_use);;
			    printf("\n");

			    dge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    dge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    sprint_vector(x, n_i, incx, "x");

			    dprint_vector(y, m_i, incy, "y");

			    dprint_vector(head_r_true, m_i, 1, "head_r_true");

			    dge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    dge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dge_sum_mv_s_d_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dge_sum_mv_s_d_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha;
  double beta;
  double beta_zero_fake;
  double alpha_use;
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  double *x;
  double *y;
  float *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  beta_zero_fake = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;




  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double));
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double));
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_dge_sum_mv_s_d_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      dcopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_dge_sum_mv_s_d_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;



			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			dcopy_vector(x, n_i, incx, x_vec, 1);
			dcopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin = 0.0;
			  rout = y[yi];
			  head_r_true_elem = head_r_true[ri];
			  tail_r_true_elem = tail_r_true[ri];

			  test_BLAS_ddot_s_d(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : d, a type : s, x type : d\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("%24.16e", alpha_use);;
			    printf("\n");

			    sge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    sge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    dprint_vector(x, n_i, incx, "x");

			    dprint_vector(y, m_i, incy, "y");

			    dprint_vector(head_r_true, m_i, 1, "head_r_true");

			    sge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    sge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_dge_sum_mv_s_s_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_dge_sum_mv_s_s_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin;
  double rout;
  double head_r_true_elem, tail_r_true_elem;

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha;
  double beta;
  double beta_zero_fake;
  double alpha_use;
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  double *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;

  FPU_FIX_DECL;

  beta_zero_fake = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;




  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double));
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double));
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_dge_sum_mv_s_s_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      scopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_dge_sum_mv_s_s_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;



			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			scopy_vector(x, n_i, incx, x_vec, 1);
			scopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin = 0.0;
			  rout = y[yi];
			  head_r_true_elem = head_r_true[ri];
			  tail_r_true_elem = tail_r_true[ri];

			  test_BLAS_ddot_s_s(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : d, a type : s, x type : s\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("%24.16e", alpha);;
			    printf("   ");
			    printf("beta = ");
			    printf("%24.16e", beta);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("%24.16e", alpha_use);;
			    printf("\n");

			    sge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    sge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    sprint_vector(x, n_i, incx, "x");

			    dprint_vector(y, m_i, incy, "y");

			    dprint_vector(head_r_true, m_i, 1, "head_r_true");

			    sge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    sge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_z_c_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_z_c_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  float *x;
  double *y;
  double *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;
  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_zge_sum_mv_z_c_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      ccopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_zge_sum_mv_z_c_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;
			incx_veci *= 2;
			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			ccopy_vector(x, n_i, incx, x_vec, 1);
			ccopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_zdot_z_c(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : z, a type : z, x type : c\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%24.16e, %24.16e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    zge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    zge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    cprint_vector(x, n_i, incx, "x");

			    zprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    zge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    zge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_c_z_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_c_z_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  double *x;
  double *y;
  float *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;
  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_zge_sum_mv_c_z_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      zcopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_zge_sum_mv_c_z_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;
			incx_veci *= 2;
			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			zcopy_vector(x, n_i, incx, x_vec, 1);
			zcopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_zdot_c_z(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : z, a type : c, x type : z\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%24.16e, %24.16e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    cge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    cge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    zprint_vector(x, n_i, incx, "x");

			    zprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    cge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    cge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_c_c_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_c_c_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  double *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;
  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_zge_sum_mv_c_c_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      ccopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_zge_sum_mv_c_c_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;
			incx_veci *= 2;
			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			ccopy_vector(x, n_i, incx, x_vec, 1);
			ccopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_zdot_c_c(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : z, a type : c, x type : c\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%24.16e, %24.16e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    cge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    cge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    cprint_vector(x, n_i, incx, "x");

			    zprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    cge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    cge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cge_sum_mv_c_s_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_cge_sum_mv_c_s_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  float alpha[2];
  float beta[2];
  float beta_zero_fake[2];
  float alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  float *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;

  incy *= 2;

  /* allocate memory for arrays */
  y = (float *) blas_malloc(4 * m_i * sizeof(float) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_cge_sum_mv_c_s_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      scopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_cge_sum_mv_c_s_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;

			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			scopy_vector(x, n_i, incx, x_vec, 1);
			scopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  cge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_cdot_c_s(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : c, a type : c, x type : s\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%16.8e, %16.8e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    cge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    cge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    sprint_vector(x, n_i, incx, "x");

			    cprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    cge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    cge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cge_sum_mv_s_c_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_cge_sum_mv_s_c_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  float alpha[2];
  float beta[2];
  float beta_zero_fake[2];
  float alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  float *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;

  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (float *) blas_malloc(4 * m_i * sizeof(float) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_cge_sum_mv_s_c_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      ccopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_cge_sum_mv_s_c_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;
			incx_veci *= 2;
			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			ccopy_vector(x, n_i, incx, x_vec, 1);
			ccopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_cdot_s_c(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : c, a type : s, x type : c\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%16.8e, %16.8e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    sge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    sge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    cprint_vector(x, n_i, incx, "x");

			    cprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    sge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    sge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_cge_sum_mv_s_s_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_cge_sum_mv_s_s_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  float rin[2];
  float rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  float alpha[2];
  float beta[2];
  float beta_zero_fake[2];
  float alpha_use[2];
  float *a;
  float *a_use;
  float *B;
  float *B_use;
  float *x;
  float *y;
  float *a_vec;
  float *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;


  incy *= 2;

  /* allocate memory for arrays */
  y = (float *) blas_malloc(4 * m_i * sizeof(float) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (float *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(float));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (float *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(float));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (float *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(float));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (float *) blas_malloc(4 * n_i * sizeof(float));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (float *) blas_malloc(2 * n_i * sizeof(float));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_cge_sum_mv_s_s_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      scopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_cge_sum_mv_s_s_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;

			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			scopy_vector(x, n_i, incx, x_vec, 1);
			scopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  sge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_cdot_s_s(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : c, a type : s, x type : s\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%16.8e, %16.8e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%16.8e, %16.8e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%16.8e, %16.8e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    sge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    sge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    sprint_vector(x, n_i, incx, "x");

			    cprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    sge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    sge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_z_d_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_z_d_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  double *x;
  double *y;
  double *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  inca *= 2;

  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double) * 2);
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use =
    (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double) * 2);
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use =
    (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double) * 2);
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;
  inca_veci *= 2;
  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_zge_sum_mv_z_d_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      dcopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_zge_sum_mv_z_d_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;

			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			dcopy_vector(x, n_i, incx, x_vec, 1);
			dcopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  zge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_zdot_z_d(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : z, a type : z, x type : d\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%24.16e, %24.16e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    zge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    zge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    dprint_vector(x, n_i, incx, "x");

			    zprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    zge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    zge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_d_z_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_d_z_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  double *x;
  double *y;
  double *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;

  incx *= 2;
  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double) * 2);
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double) * 2);
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_zge_sum_mv_d_z_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      zcopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_zge_sum_mv_d_z_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;
			incx_veci *= 2;
			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			zcopy_vector(x, n_i, incx, x_vec, 1);
			zcopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_zdot_d_z(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : z, a type : d, x type : z\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%24.16e, %24.16e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    dge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    dge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    zprint_vector(x, n_i, incx, "x");

			    zprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    dge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    dge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;

}
void do_test_zge_sum_mv_d_d_x
  (int m, int n,
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests) {

  /* Function name */
  const char fname[] = "BLAS_zge_sum_mv_d_d_x";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;		/* internal machine epsilon     */
  double un_int;		/* internal underflow threshold */

  double rin[2];
  double rout[2];
  double head_r_true_elem[2], tail_r_true_elem[2];

  enum blas_order_type order_type;
  enum blas_prec_type prec;

  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;

  int prec_val;

  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  double alpha[2];
  double beta[2];
  double beta_zero_fake[2];
  double alpha_use[2];
  double *a;
  double *a_use;
  double *B;
  double *B_use;
  double *x;
  double *y;
  double *a_vec;
  double *x_vec;


  double *ratios;

  /* true result calculated by testgen, in double-double */
  double *head_r_true, *tail_r_true;


  FPU_FIX_DECL;

  beta_zero_fake[0] = beta_zero_fake[1] = 0.0;

  if (n < 0 || ntests < 0)
    BLAS_error(fname, -3, n, NULL);

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

  FPU_FIX_START;

  n_i = n;
  m_i = m;

  inca = incx = incy = 1;


  incy *= 2;

  /* allocate memory for arrays */
  y = (double *) blas_malloc(4 * m_i * sizeof(double) * 2);
  if (4 * m_i > 0 && y == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a = (double *) blas_malloc(2 * n_i * m_i * m_i * n_i * sizeof(double));
  if (2 * n_i * m_i * m_i * n_i > 0 && a == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  a_use = (double *) blas_malloc(2 * m_i * n_i * m_i * n_i * sizeof(double));
  if (2 * m_i * n_i * m_i * n_i > 0 && a_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  B_use = (double *) blas_malloc(2 * n_i * n_i * m_i * m_i * sizeof(double));
  if (2 * n_i * n_i * m_i * m_i > 0 && B_use == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x = (double *) blas_malloc(4 * n_i * sizeof(double));
  if (4 * n_i > 0 && x == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  inca_veci = 1;

  a_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && a_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_vec = (double *) blas_malloc(2 * n_i * sizeof(double));
  if (2 * n_i > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  tail_r_true = (double *) blas_malloc(m_i * sizeof(double) * 2);
  if (m_i > 0 && (head_r_true == NULL || tail_r_true == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  ratios = (double *) blas_malloc(m_i * sizeof(double));
  if (m_i > 0 && ratios == NULL) {
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

	      /* vary lda = n_i, n_i+1, 2*n_i */
	      for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {

		if (order_type == blas_rowmajor) {
		  lda = (lda_val == 0) ? n_i :
		    (lda_val == 1) ? n_i + 1 : n_i * n_i;
		} else {
		  lda = (lda_val == 0) ? m_i :
		    (lda_val == 1) ? m_i + 1 : m_i * m_i;
		}

		/* vary ldb = n_i, n_i+1, 2*n_i */
		for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {

		  if (order_type == blas_rowmajor) {
		    ldb = (ldb_val == 0) ? n_i :
		      (ldb_val == 1) ? n_i + 1 : n_i * n_i;
		  } else {
		    ldb = (ldb_val == 0) ? m_i :
		      (ldb_val == 1) ? m_i + 1 : m_i * m_i;
		  }

		  for (randomize_val = RANDOMIZE_START;
		       randomize_val <= RANDOMIZE_END; randomize_val++) {

		    /* For the sake of speed, we throw out this case at random */
		    if (xrand(seed) >= test_prob)
		      continue;

		    /* finally we are here to generate the test case */
		    /* alpha_use, a_use, B_use are the generated alpha, a, B
		     *  before any scaling.  
		     *  That is, in the generator, alpha == beta == alpha_use 
		     *  before scaling. */

		    saved_seed = *seed;
		    BLAS_zge_sum_mv_d_d_testgen(norm, order_type,
						m, n, randomize_val, &alpha,
						alpha_flag, &beta, beta_flag,
						a, lda, B, ldb, x_vec, 1,
						&alpha_use, a_use, B_use,
						seed, head_r_true,
						tail_r_true);

		    /* vary incx = 1, 2 */
		    for (incx_val = INCX_START; incx_val <= INCX_END;
			 incx_val++) {

		      incx = incx_val;
		      if (0 == incx)
			continue;

		      dcopy_vector(x_vec, n_i, 1, x, incx);

		      /* vary incy = 1, 2 */
		      for (incy_val = INCY_START; incy_val <= INCY_END;
			   incy_val++) {

			incy = incy_val;
			if (0 == incy)
			  continue;

			test_count++;

			/* call ge_sum_mv routines to be tested */
			FPU_FIX_STOP;
			BLAS_zge_sum_mv_d_d_x(order_type,
					      m, n, alpha, a, lda, x, incx,
					      beta, B, ldb, y, incy, prec);
			FPU_FIX_START;

			/* now compute the ratio using test_BLAS_xdot */
			/* copy a row from A, use x, run 
			   dot test */

			incyi = incy;

			incri = 1;
			incx_veci = 1;

			incyi *= 2;
			incri *= 2;
			if (incy < 0) {
			  y_starti = (-m_i + 1) * incyi;
			} else {
			  y_starti = 0;
			}
			/* make two copies of x into x_vec. redundant */
			dcopy_vector(x, n_i, incx, x_vec, 1);
			dcopy_vector(x, n_i, incx,
				     (x_vec + (n_i * incx_veci)), 1);
			for (i = 0, yi = y_starti, ri = 0; i < m_i;
			     i++, yi += incyi, ri += incri) {
			  dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       a_use, lda, a_vec, i);
			  dge_copy_row(order_type, blas_no_trans, m_i, n_i,
				       B_use, ldb, (a_vec + inca_veci * n_i),
				       i);

			  rin[0] = rin[1] = 0.0;
			  rout[0] = y[yi];
			  rout[1] = y[yi + 1];
			  head_r_true_elem[0] = head_r_true[ri];
			  head_r_true_elem[1] = head_r_true[ri + 1];
			  tail_r_true_elem[0] = tail_r_true[ri];
			  tail_r_true_elem[1] = tail_r_true[ri + 1];

			  test_BLAS_zdot_d_d(2 * n_i,
					     blas_no_conj,
					     alpha_use, beta_zero_fake, rin,
					     rout, head_r_true_elem,
					     tail_r_true_elem, a_vec, 1,
					     x_vec, 1, eps_int, un_int,
					     &ratios[i]);

			  /* take the max ratio */
			  if (i == 0) {
			    ratio = ratios[0];
			    /* The !<= below causes NaN errors
			     *  to be included.
			     * Note that (NaN > 0) is false */
			  } else if (!(ratios[i] <= ratio)) {
			    ratio = ratios[i];
			  }
			}	/* end of dot-test loop */

			/* The !<= below causes NaN errors
			 *  to be included.
			 * Note that (NaN > 0) is false */
			if (!(ratio <= thresh)) {

			  if (debug == 3) {
			    printf("\n\t\tTest # %d\n", test_count);
			    printf("y type : z, a type : d, x type : d\n");
			    printf("Seed = %d\t", saved_seed);
			    printf("n %d, m %d\n", n, m);
			    printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda,
				   ldb, incx, incx);

			    if (order_type == blas_rowmajor)
			      printf("row ");
			    else
			      printf("col ");

			    printf("NORM %d, ALPHA %d, BETA %d\n",
				   norm, alpha_val, beta_val);
			    printf("randomize %d\n", randomize_val);

			    /* print out info */
			    printf("alpha = ");
			    printf("(%24.16e, %24.16e)", alpha[0], alpha[1]);;
			    printf("   ");
			    printf("beta = ");
			    printf("(%24.16e, %24.16e)", beta[0], beta[1]);;
			    printf("\n");
			    printf("alpha_use = ");
			    printf("(%24.16e, %24.16e)", alpha_use[0],
				   alpha_use[1]);;
			    printf("\n");

			    dge_print_matrix(a, m_i, n_i, lda, order_type,
					     "A");
			    dge_print_matrix(B, m_i, n_i, ldb, order_type,
					     "B");
			    dprint_vector(x, n_i, incx, "x");

			    zprint_vector(y, m_i, incy, "y");

			    zprint_vector(head_r_true, m_i, 1, "head_r_true");

			    dge_print_matrix(a_use, m_i, n_i, lda, order_type,
					     "A_use");
			    dge_print_matrix(B_use, m_i, n_i, ldb, order_type,
					     "B_use");

			    dprint_vector(ratios, m_i, 1, "ratios");
			    printf("ratio = %g\n", ratio);
			    fflush(stdout);
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

		  }		/* end of randmize loop */

		}		/* end of ldb loop */

	      }			/* end of lda loop */

	    }			/* end of order loop */

	  }			/* end of nr test loop */

	}			/* end of norm loop */


      }				/* end of prec loop */

    }				/* end of beta loop */

  }				/* end of alpha loop */

  FPU_FIX_STOP;

end:
  blas_free(y);
  blas_free(a);
  blas_free(a_use);
  blas_free(B);
  blas_free(B_use);
  blas_free(x);
  blas_free(head_r_true);
  blas_free(tail_r_true);
  blas_free(ratios);
  blas_free(a_vec);
  blas_free(x_vec);

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
  const char *base_routine = "ge_sum_mv";
  char *fname;
  int n;

  int m, i;
  int n_data[NUM_DATA][2] =
    { {1, 1}, {1, 2}, {3, 2}, {8, 6}, {9, 10}, {4, 4}, {7, 7} };

  if (argc != 6) {
    printf("Usage:\n");
    printf
      ("do_test_ge_sum_mv <nsizes> <ntests> <thresh> <debug> <test_prob>\n");
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
    BLAS_error("Testing ge_sum_mv", 0, 0, NULL);

  printf("Testing %s...\n", base_routine);
  printf("INPUT: nsizes = %d, ntests = %d, thresh = %4.2f, debug = %d\n\n",
	 nsizes, ntests, thresh, debug);





  fname = "BLAS_dge_sum_mv_d_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_dge_sum_mv_d_s(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_dge_sum_mv_s_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_dge_sum_mv_s_d(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_dge_sum_mv_s_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_dge_sum_mv_s_s(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_zge_sum_mv_z_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_z_c(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_zge_sum_mv_c_z";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_c_z(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_zge_sum_mv_c_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_c_c(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_cge_sum_mv_c_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_cge_sum_mv_c_s(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_cge_sum_mv_s_c";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_cge_sum_mv_s_c(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_cge_sum_mv_s_s";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_cge_sum_mv_s_s(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_zge_sum_mv_z_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_z_d(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_zge_sum_mv_d_z";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_d_z(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_zge_sum_mv_d_d";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_d_d(m, n,
			   ntests, &seed, thresh, debug,
			   test_prob,
			   &min_ratio, &max_ratio, &num_bad_ratio,
			   &num_tests);

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

  fname = "BLAS_sge_sum_mv_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_sge_sum_mv_x(m, n,
			 ntests, &seed, thresh, debug,
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

  fname = "BLAS_dge_sum_mv_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_dge_sum_mv_x(m, n,
			 ntests, &seed, thresh, debug,
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

  fname = "BLAS_cge_sum_mv_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_cge_sum_mv_x(m, n,
			 ntests, &seed, thresh, debug,
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

  fname = "BLAS_zge_sum_mv_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_x(m, n,
			 ntests, &seed, thresh, debug,
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

  fname = "BLAS_dge_sum_mv_d_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_dge_sum_mv_d_s_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_dge_sum_mv_s_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_dge_sum_mv_s_d_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_dge_sum_mv_s_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_dge_sum_mv_s_s_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_zge_sum_mv_z_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_z_c_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_zge_sum_mv_c_z_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_c_z_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_zge_sum_mv_c_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_c_c_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_cge_sum_mv_c_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_cge_sum_mv_c_s_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_cge_sum_mv_s_c_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_cge_sum_mv_s_c_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_cge_sum_mv_s_s_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_cge_sum_mv_s_s_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_zge_sum_mv_z_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_z_d_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_zge_sum_mv_d_z_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_d_z_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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

  fname = "BLAS_zge_sum_mv_d_d_x";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    do_test_zge_sum_mv_d_d_x(m, n,
			     ntests, &seed, thresh, debug,
			     test_prob,
			     &min_ratio, &max_ratio, &num_bad_ratio,
			     &num_tests);

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
