#include "blas_extended.h"
#include "blas_extended_test.h"


void BLAS_strmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, float *alpha,
			int alpha_flag, float *T, int ldt, float *x,
			int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, T and x, where T is a triangular matrix; and 
 * computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) float*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) float*
 *
 * x            (input/output) float*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  float *x_i = x;
  float *T_i = T;
  float *alpha_i = alpha;
  float *x_vec;
  float *t_vec;
  float beta;
  float r;
  double head_r_true_elem, tail_r_true_elem;
  float x_elem;
  float t_elem;

  int inc_tvec = 1, inc_xvec = 1;
  int xvec_i, tvec_j;
  int xi;
  int ti, tij;
  int inc_ti, inc_tij;
  int inc_xi;
  int i, j;

  r = 0.0;
  beta = 0.0;




  t_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && t_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  x_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -ldt;
	inc_tij = -1;
      } else {
	inc_ti = -1;
	inc_tij = -ldt;
      }
    } else {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = ldt;
	inc_tij = 1;
      } else {
	inc_ti = 1;
	inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = 1;
	inc_tij = ldt;
      } else {
	inc_ti = ldt;
	inc_tij = 1;
      }
    } else {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -1;
	inc_tij = -ldt;
      } else {
	inc_ti = -ldt;
	inc_tij = -1;
      }
    }
  }






  /* Call dot_testgen n times.  Each call will generate
   * one row of T and one element of x.                   */
  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
  xi = (inc_xi > 0 ? 0 : -(n - 1) * inc_xi);
  xvec_i = 0;
  for (i = 0; i < n; i++) {

    /* Generate the i-th element of x_vec and all of t_vec. */
    if (diag == blas_unit_diag) {
      /* Since we need alpha = beta, we fix alpha if alpha_flag = 0. */
      if (i == 0 && alpha_flag == 0) {
	*alpha_i = xrand(seed);
      }
      BLAS_sdot_testgen(i, 0, i, norm, blas_no_conj, alpha_i,
			1, alpha_i, 1, x_vec, t_vec,
			seed, &r, &head_r_true_elem, &tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j < i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Set the diagonal element to 1. */
      t_elem = 1.0;
      T_i[tij] = t_elem;

      /* Set x[i] to be r. */
      x_i[xi] = r;
      x_vec[xvec_i] = r;

    } else {
      BLAS_sdot_testgen(i + 1, 0, i, norm, blas_no_conj, alpha,
			(i == 0 ? alpha_flag : 1), &beta, 1, x_vec, t_vec,
			seed, &r, &head_r_true_elem, &tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j <= i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Copy generated x_vec[i] to appropriate position in x. */
      x_elem = x_vec[xvec_i];
      x_i[xi] = x_elem;
    }

    /* Copy r_true */
    head_r_true[xi] = head_r_true_elem;
    tail_r_true[xi] = tail_r_true_elem;

    xvec_i += inc_xvec;
    xi += inc_xi;

    ti += inc_ti;
  }

  blas_free(x_vec);
  blas_free(t_vec);
}

  /* end of BLAS_strmv_testgen */


void BLAS_dtrmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, double *alpha,
			int alpha_flag, double *T, int ldt, double *x,
			int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, T and x, where T is a triangular matrix; and 
 * computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) double*
 *
 * x            (input/output) double*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *x_i = x;
  double *T_i = T;
  double *alpha_i = alpha;
  double *x_vec;
  double *t_vec;
  double beta;
  double r;
  double head_r_true_elem, tail_r_true_elem;
  double x_elem;
  double t_elem;

  int inc_tvec = 1, inc_xvec = 1;
  int xvec_i, tvec_j;
  int xi;
  int ti, tij;
  int inc_ti, inc_tij;
  int inc_xi;
  int i, j;

  r = 0.0;
  beta = 0.0;




  t_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && t_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  x_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -ldt;
	inc_tij = -1;
      } else {
	inc_ti = -1;
	inc_tij = -ldt;
      }
    } else {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = ldt;
	inc_tij = 1;
      } else {
	inc_ti = 1;
	inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = 1;
	inc_tij = ldt;
      } else {
	inc_ti = ldt;
	inc_tij = 1;
      }
    } else {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -1;
	inc_tij = -ldt;
      } else {
	inc_ti = -ldt;
	inc_tij = -1;
      }
    }
  }






  /* Call dot_testgen n times.  Each call will generate
   * one row of T and one element of x.                   */
  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
  xi = (inc_xi > 0 ? 0 : -(n - 1) * inc_xi);
  xvec_i = 0;
  for (i = 0; i < n; i++) {

    /* Generate the i-th element of x_vec and all of t_vec. */
    if (diag == blas_unit_diag) {
      /* Since we need alpha = beta, we fix alpha if alpha_flag = 0. */
      if (i == 0 && alpha_flag == 0) {
	*alpha_i = xrand(seed);
      }
      BLAS_ddot_testgen(i, 0, i, norm, blas_no_conj, alpha_i,
			1, alpha_i, 1, x_vec, t_vec,
			seed, &r, &head_r_true_elem, &tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j < i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Set the diagonal element to 1. */
      t_elem = 1.0;
      T_i[tij] = t_elem;

      /* Set x[i] to be r. */
      x_i[xi] = r;
      x_vec[xvec_i] = r;

    } else {
      BLAS_ddot_testgen(i + 1, 0, i, norm, blas_no_conj, alpha,
			(i == 0 ? alpha_flag : 1), &beta, 1, x_vec, t_vec,
			seed, &r, &head_r_true_elem, &tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j <= i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Copy generated x_vec[i] to appropriate position in x. */
      x_elem = x_vec[xvec_i];
      x_i[xi] = x_elem;
    }

    /* Copy r_true */
    head_r_true[xi] = head_r_true_elem;
    tail_r_true[xi] = tail_r_true_elem;

    xvec_i += inc_xvec;
    xi += inc_xi;

    ti += inc_ti;
  }

  blas_free(x_vec);
  blas_free(t_vec);
}

  /* end of BLAS_dtrmv_testgen */


void BLAS_dtrmv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, double *alpha,
			  int alpha_flag, float *T, int ldt, double *x,
			  int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, T and x, where T is a triangular matrix; and 
 * computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) float*
 *
 * x            (input/output) double*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *x_i = x;
  float *T_i = T;
  double *alpha_i = alpha;
  double *x_vec;
  float *t_vec;
  double beta;
  double r;
  double head_r_true_elem, tail_r_true_elem;
  double x_elem;
  float t_elem;

  int inc_tvec = 1, inc_xvec = 1;
  int xvec_i, tvec_j;
  int xi;
  int ti, tij;
  int inc_ti, inc_tij;
  int inc_xi;
  int i, j;

  r = 0.0;
  beta = 0.0;




  t_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && t_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  x_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -ldt;
	inc_tij = -1;
      } else {
	inc_ti = -1;
	inc_tij = -ldt;
      }
    } else {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = ldt;
	inc_tij = 1;
      } else {
	inc_ti = 1;
	inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = 1;
	inc_tij = ldt;
      } else {
	inc_ti = ldt;
	inc_tij = 1;
      }
    } else {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -1;
	inc_tij = -ldt;
      } else {
	inc_ti = -ldt;
	inc_tij = -1;
      }
    }
  }






  /* Call dot_testgen n times.  Each call will generate
   * one row of T and one element of x.                   */
  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
  xi = (inc_xi > 0 ? 0 : -(n - 1) * inc_xi);
  xvec_i = 0;
  for (i = 0; i < n; i++) {

    /* Generate the i-th element of x_vec and all of t_vec. */
    if (diag == blas_unit_diag) {
      /* Since we need alpha = beta, we fix alpha if alpha_flag = 0. */
      if (i == 0 && alpha_flag == 0) {
	*alpha_i = xrand(seed);
      }
      BLAS_ddot_d_s_testgen(i, 0, i, norm, blas_no_conj, alpha_i,
			    1, alpha_i, 1, x_vec, t_vec,
			    seed, &r, &head_r_true_elem, &tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j < i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Set the diagonal element to 1. */
      t_elem = 1.0;
      T_i[tij] = t_elem;

      /* Set x[i] to be r. */
      x_i[xi] = r;
      x_vec[xvec_i] = r;

    } else {
      BLAS_ddot_d_s_testgen(i + 1, 0, i, norm, blas_no_conj, alpha,
			    (i == 0 ? alpha_flag : 1), &beta, 1, x_vec, t_vec,
			    seed, &r, &head_r_true_elem, &tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j <= i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Copy generated x_vec[i] to appropriate position in x. */
      x_elem = x_vec[xvec_i];
      x_i[xi] = x_elem;
    }

    /* Copy r_true */
    head_r_true[xi] = head_r_true_elem;
    tail_r_true[xi] = tail_r_true_elem;

    xvec_i += inc_xvec;
    xi += inc_xi;

    ti += inc_ti;
  }

  blas_free(x_vec);
  blas_free(t_vec);
}

  /* end of BLAS_dtrmv_s_testgen */


void BLAS_ctrmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *T, int ldt, void *x, int *seed,
			double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, T and x, where T is a triangular matrix; and 
 * computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) void*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  float *x_i = (float *) x;
  float *T_i = (float *) T;
  float *alpha_i = (float *) alpha;
  float *x_vec;
  float *t_vec;
  float beta[2];
  float r[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  float x_elem[2];
  float t_elem[2];

  int inc_tvec = 1, inc_xvec = 1;
  int xvec_i, tvec_j;
  int xi;
  int ti, tij;
  int inc_ti, inc_tij;
  int inc_xi;
  int i, j;

  r[0] = r[1] = 0.0;
  beta[0] = beta[1] = 0.0;

  inc_tvec *= 2;
  inc_xvec *= 2;

  t_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && t_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  x_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -ldt;
	inc_tij = -1;
      } else {
	inc_ti = -1;
	inc_tij = -ldt;
      }
    } else {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = ldt;
	inc_tij = 1;
      } else {
	inc_ti = 1;
	inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = 1;
	inc_tij = ldt;
      } else {
	inc_ti = ldt;
	inc_tij = 1;
      }
    } else {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -1;
	inc_tij = -ldt;
      } else {
	inc_ti = -ldt;
	inc_tij = -1;
      }
    }
  }

  inc_xi *= 2;

  inc_ti *= 2;
  inc_tij *= 2;

  /* Call dot_testgen n times.  Each call will generate
   * one row of T and one element of x.                   */
  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
  xi = (inc_xi > 0 ? 0 : -(n - 1) * inc_xi);
  xvec_i = 0;
  for (i = 0; i < n; i++) {

    /* Generate the i-th element of x_vec and all of t_vec. */
    if (diag == blas_unit_diag) {
      /* Since we need alpha = beta, we fix alpha if alpha_flag = 0. */
      if (i == 0 && alpha_flag == 0) {
	alpha_i[0] = xrand(seed);
	alpha_i[1] = xrand(seed);
      }
      BLAS_cdot_testgen(i, 0, i, norm, blas_no_conj, alpha_i,
			1, alpha_i, 1, x_vec, t_vec,
			seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j < i; j++) {
	t_elem[0] = t_vec[tvec_j];
	t_elem[1] = t_vec[tvec_j + 1];

	if (trans == blas_conj_trans) {
	  t_elem[1] = -t_elem[1];
	}
	T_i[tij] = t_elem[0];
	T_i[tij + 1] = t_elem[1];
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Set the diagonal element to 1. */
      t_elem[0] = 1.0;
      t_elem[1] = 0.0;
      T_i[tij] = t_elem[0];
      T_i[tij + 1] = t_elem[1];

      /* Set x[i] to be r. */
      x_i[xi] = r[0];
      x_i[xi + 1] = r[1];
      x_vec[xvec_i] = r[0];
      x_vec[xvec_i + 1] = r[1];

    } else {
      BLAS_cdot_testgen(i + 1, 0, i, norm, blas_no_conj, alpha,
			(i == 0 ? alpha_flag : 1), beta, 1, x_vec, t_vec,
			seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j <= i; j++) {
	t_elem[0] = t_vec[tvec_j];
	t_elem[1] = t_vec[tvec_j + 1];

	if (trans == blas_conj_trans) {
	  t_elem[1] = -t_elem[1];
	}
	T_i[tij] = t_elem[0];
	T_i[tij + 1] = t_elem[1];
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Copy generated x_vec[i] to appropriate position in x. */
      x_elem[0] = x_vec[xvec_i];
      x_elem[1] = x_vec[xvec_i + 1];
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* Copy r_true */
    head_r_true[xi] = head_r_true_elem[0];
    head_r_true[xi + 1] = head_r_true_elem[1];
    tail_r_true[xi] = tail_r_true_elem[0];
    tail_r_true[xi + 1] = tail_r_true_elem[1];

    xvec_i += inc_xvec;
    xi += inc_xi;

    ti += inc_ti;
  }

  blas_free(x_vec);
  blas_free(t_vec);
}

  /* end of BLAS_ctrmv_testgen */


void BLAS_ztrmv_c_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, void *T, int ldt, void *x,
			  int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, T and x, where T is a triangular matrix; and 
 * computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) void*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *x_i = (double *) x;
  float *T_i = (float *) T;
  double *alpha_i = (double *) alpha;
  double *x_vec;
  float *t_vec;
  double beta[2];
  double r[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  double x_elem[2];
  float t_elem[2];

  int inc_tvec = 1, inc_xvec = 1;
  int xvec_i, tvec_j;
  int xi;
  int ti, tij;
  int inc_ti, inc_tij;
  int inc_xi;
  int i, j;

  r[0] = r[1] = 0.0;
  beta[0] = beta[1] = 0.0;

  inc_tvec *= 2;
  inc_xvec *= 2;

  t_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && t_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  x_vec = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -ldt;
	inc_tij = -1;
      } else {
	inc_ti = -1;
	inc_tij = -ldt;
      }
    } else {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = ldt;
	inc_tij = 1;
      } else {
	inc_ti = 1;
	inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = 1;
	inc_tij = ldt;
      } else {
	inc_ti = ldt;
	inc_tij = 1;
      }
    } else {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -1;
	inc_tij = -ldt;
      } else {
	inc_ti = -ldt;
	inc_tij = -1;
      }
    }
  }

  inc_xi *= 2;

  inc_ti *= 2;
  inc_tij *= 2;

  /* Call dot_testgen n times.  Each call will generate
   * one row of T and one element of x.                   */
  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
  xi = (inc_xi > 0 ? 0 : -(n - 1) * inc_xi);
  xvec_i = 0;
  for (i = 0; i < n; i++) {

    /* Generate the i-th element of x_vec and all of t_vec. */
    if (diag == blas_unit_diag) {
      /* Since we need alpha = beta, we fix alpha if alpha_flag = 0. */
      if (i == 0 && alpha_flag == 0) {
	alpha_i[0] = xrand(seed);
	alpha_i[1] = xrand(seed);
      }
      BLAS_zdot_z_c_testgen(i, 0, i, norm, blas_no_conj, alpha_i,
			    1, alpha_i, 1, x_vec, t_vec,
			    seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j < i; j++) {
	t_elem[0] = t_vec[tvec_j];
	t_elem[1] = t_vec[tvec_j + 1];

	if (trans == blas_conj_trans) {
	  t_elem[1] = -t_elem[1];
	}
	T_i[tij] = t_elem[0];
	T_i[tij + 1] = t_elem[1];
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Set the diagonal element to 1. */
      t_elem[0] = 1.0;
      t_elem[1] = 0.0;
      T_i[tij] = t_elem[0];
      T_i[tij + 1] = t_elem[1];

      /* Set x[i] to be r. */
      x_i[xi] = r[0];
      x_i[xi + 1] = r[1];
      x_vec[xvec_i] = r[0];
      x_vec[xvec_i + 1] = r[1];

    } else {
      BLAS_zdot_z_c_testgen(i + 1, 0, i, norm, blas_no_conj, alpha,
			    (i == 0 ? alpha_flag : 1), beta, 1, x_vec, t_vec,
			    seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j <= i; j++) {
	t_elem[0] = t_vec[tvec_j];
	t_elem[1] = t_vec[tvec_j + 1];

	if (trans == blas_conj_trans) {
	  t_elem[1] = -t_elem[1];
	}
	T_i[tij] = t_elem[0];
	T_i[tij + 1] = t_elem[1];
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Copy generated x_vec[i] to appropriate position in x. */
      x_elem[0] = x_vec[xvec_i];
      x_elem[1] = x_vec[xvec_i + 1];
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* Copy r_true */
    head_r_true[xi] = head_r_true_elem[0];
    head_r_true[xi + 1] = head_r_true_elem[1];
    tail_r_true[xi] = tail_r_true_elem[0];
    tail_r_true[xi + 1] = tail_r_true_elem[1];

    xvec_i += inc_xvec;
    xi += inc_xi;

    ti += inc_ti;
  }

  blas_free(x_vec);
  blas_free(t_vec);
}

  /* end of BLAS_ztrmv_c_testgen */


void BLAS_ztrmv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *T, int ldt, void *x, int *seed,
			double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, T and x, where T is a triangular matrix; and 
 * computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) void*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *x_i = (double *) x;
  double *T_i = (double *) T;
  double *alpha_i = (double *) alpha;
  double *x_vec;
  double *t_vec;
  double beta[2];
  double r[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  double x_elem[2];
  double t_elem[2];

  int inc_tvec = 1, inc_xvec = 1;
  int xvec_i, tvec_j;
  int xi;
  int ti, tij;
  int inc_ti, inc_tij;
  int inc_xi;
  int i, j;

  r[0] = r[1] = 0.0;
  beta[0] = beta[1] = 0.0;

  inc_tvec *= 2;
  inc_xvec *= 2;

  t_vec = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && t_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  x_vec = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -ldt;
	inc_tij = -1;
      } else {
	inc_ti = -1;
	inc_tij = -ldt;
      }
    } else {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = ldt;
	inc_tij = 1;
      } else {
	inc_ti = 1;
	inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = 1;
	inc_tij = ldt;
      } else {
	inc_ti = ldt;
	inc_tij = 1;
      }
    } else {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -1;
	inc_tij = -ldt;
      } else {
	inc_ti = -ldt;
	inc_tij = -1;
      }
    }
  }

  inc_xi *= 2;

  inc_ti *= 2;
  inc_tij *= 2;

  /* Call dot_testgen n times.  Each call will generate
   * one row of T and one element of x.                   */
  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
  xi = (inc_xi > 0 ? 0 : -(n - 1) * inc_xi);
  xvec_i = 0;
  for (i = 0; i < n; i++) {

    /* Generate the i-th element of x_vec and all of t_vec. */
    if (diag == blas_unit_diag) {
      /* Since we need alpha = beta, we fix alpha if alpha_flag = 0. */
      if (i == 0 && alpha_flag == 0) {
	alpha_i[0] = xrand(seed);
	alpha_i[1] = xrand(seed);
      }
      BLAS_zdot_testgen(i, 0, i, norm, blas_no_conj, alpha_i,
			1, alpha_i, 1, x_vec, t_vec,
			seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j < i; j++) {
	t_elem[0] = t_vec[tvec_j];
	t_elem[1] = t_vec[tvec_j + 1];

	if (trans == blas_conj_trans) {
	  t_elem[1] = -t_elem[1];
	}
	T_i[tij] = t_elem[0];
	T_i[tij + 1] = t_elem[1];
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Set the diagonal element to 1. */
      t_elem[0] = 1.0;
      t_elem[1] = 0.0;
      T_i[tij] = t_elem[0];
      T_i[tij + 1] = t_elem[1];

      /* Set x[i] to be r. */
      x_i[xi] = r[0];
      x_i[xi + 1] = r[1];
      x_vec[xvec_i] = r[0];
      x_vec[xvec_i + 1] = r[1];

    } else {
      BLAS_zdot_testgen(i + 1, 0, i, norm, blas_no_conj, alpha,
			(i == 0 ? alpha_flag : 1), beta, 1, x_vec, t_vec,
			seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j <= i; j++) {
	t_elem[0] = t_vec[tvec_j];
	t_elem[1] = t_vec[tvec_j + 1];

	if (trans == blas_conj_trans) {
	  t_elem[1] = -t_elem[1];
	}
	T_i[tij] = t_elem[0];
	T_i[tij + 1] = t_elem[1];
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Copy generated x_vec[i] to appropriate position in x. */
      x_elem[0] = x_vec[xvec_i];
      x_elem[1] = x_vec[xvec_i + 1];
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* Copy r_true */
    head_r_true[xi] = head_r_true_elem[0];
    head_r_true[xi + 1] = head_r_true_elem[1];
    tail_r_true[xi] = tail_r_true_elem[0];
    tail_r_true[xi + 1] = tail_r_true_elem[1];

    xvec_i += inc_xvec;
    xi += inc_xi;

    ti += inc_ti;
  }

  blas_free(x_vec);
  blas_free(t_vec);
}

  /* end of BLAS_ztrmv_testgen */


void BLAS_ctrmv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, float *T, int ldt, void *x,
			  int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, T and x, where T is a triangular matrix; and 
 * computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) float*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  float *x_i = (float *) x;
  float *T_i = T;
  float *alpha_i = (float *) alpha;
  float *x_vec;
  float *t_vec;
  float beta[2];
  float r[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  float x_elem[2];
  float t_elem;

  int inc_tvec = 1, inc_xvec = 1;
  int xvec_i, tvec_j;
  int xi;
  int ti, tij;
  int inc_ti, inc_tij;
  int inc_xi;
  int i, j;

  r[0] = r[1] = 0.0;
  beta[0] = beta[1] = 0.0;


  inc_xvec *= 2;

  t_vec = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && t_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  x_vec = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -ldt;
	inc_tij = -1;
      } else {
	inc_ti = -1;
	inc_tij = -ldt;
      }
    } else {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = ldt;
	inc_tij = 1;
      } else {
	inc_ti = 1;
	inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = 1;
	inc_tij = ldt;
      } else {
	inc_ti = ldt;
	inc_tij = 1;
      }
    } else {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -1;
	inc_tij = -ldt;
      } else {
	inc_ti = -ldt;
	inc_tij = -1;
      }
    }
  }

  inc_xi *= 2;




  /* Call dot_testgen n times.  Each call will generate
   * one row of T and one element of x.                   */
  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
  xi = (inc_xi > 0 ? 0 : -(n - 1) * inc_xi);
  xvec_i = 0;
  for (i = 0; i < n; i++) {

    /* Generate the i-th element of x_vec and all of t_vec. */
    if (diag == blas_unit_diag) {
      /* Since we need alpha = beta, we fix alpha if alpha_flag = 0. */
      if (i == 0 && alpha_flag == 0) {
	alpha_i[0] = xrand(seed);
	alpha_i[1] = xrand(seed);
      }
      BLAS_cdot_c_s_testgen(i, 0, i, norm, blas_no_conj, alpha_i,
			    1, alpha_i, 1, x_vec, t_vec,
			    seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j < i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Set the diagonal element to 1. */
      t_elem = 1.0;
      T_i[tij] = t_elem;

      /* Set x[i] to be r. */
      x_i[xi] = r[0];
      x_i[xi + 1] = r[1];
      x_vec[xvec_i] = r[0];
      x_vec[xvec_i + 1] = r[1];

    } else {
      BLAS_cdot_c_s_testgen(i + 1, 0, i, norm, blas_no_conj, alpha,
			    (i == 0 ? alpha_flag : 1), beta, 1, x_vec, t_vec,
			    seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j <= i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Copy generated x_vec[i] to appropriate position in x. */
      x_elem[0] = x_vec[xvec_i];
      x_elem[1] = x_vec[xvec_i + 1];
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* Copy r_true */
    head_r_true[xi] = head_r_true_elem[0];
    head_r_true[xi + 1] = head_r_true_elem[1];
    tail_r_true[xi] = tail_r_true_elem[0];
    tail_r_true[xi + 1] = tail_r_true_elem[1];

    xvec_i += inc_xvec;
    xi += inc_xi;

    ti += inc_ti;
  }

  blas_free(x_vec);
  blas_free(t_vec);
}

  /* end of BLAS_ctrmv_s_testgen */


void BLAS_ztrmv_d_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, double *T, int ldt, void *x,
			  int *seed, double *head_r_true, double *tail_r_true)

/*
 * Purpose
 * =======
 *
 * Generates alpha, T and x, where T is a triangular matrix; and 
 * computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) double*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * head_r_true     (output) double*
 *              The leading part of the truth in double-double.
 *
 * tail_r_true     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  double *x_i = (double *) x;
  double *T_i = T;
  double *alpha_i = (double *) alpha;
  double *x_vec;
  double *t_vec;
  double beta[2];
  double r[2];
  double head_r_true_elem[2], tail_r_true_elem[2];
  double x_elem[2];
  double t_elem;

  int inc_tvec = 1, inc_xvec = 1;
  int xvec_i, tvec_j;
  int xi;
  int ti, tij;
  int inc_ti, inc_tij;
  int inc_xi;
  int i, j;

  r[0] = r[1] = 0.0;
  beta[0] = beta[1] = 0.0;


  inc_xvec *= 2;

  t_vec = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && t_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };
  x_vec = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && x_vec == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  };

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -ldt;
	inc_tij = -1;
      } else {
	inc_ti = -1;
	inc_tij = -ldt;
      }
    } else {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = ldt;
	inc_tij = 1;
      } else {
	inc_ti = 1;
	inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_xi = 1;
      if (order == blas_rowmajor) {
	inc_ti = 1;
	inc_tij = ldt;
      } else {
	inc_ti = ldt;
	inc_tij = 1;
      }
    } else {
      inc_xi = -1;
      if (order == blas_rowmajor) {
	inc_ti = -1;
	inc_tij = -ldt;
      } else {
	inc_ti = -ldt;
	inc_tij = -1;
      }
    }
  }

  inc_xi *= 2;




  /* Call dot_testgen n times.  Each call will generate
   * one row of T and one element of x.                   */
  ti = (inc_ti > 0 ? 0 : -(n - 1) * inc_ti);
  xi = (inc_xi > 0 ? 0 : -(n - 1) * inc_xi);
  xvec_i = 0;
  for (i = 0; i < n; i++) {

    /* Generate the i-th element of x_vec and all of t_vec. */
    if (diag == blas_unit_diag) {
      /* Since we need alpha = beta, we fix alpha if alpha_flag = 0. */
      if (i == 0 && alpha_flag == 0) {
	alpha_i[0] = xrand(seed);
	alpha_i[1] = xrand(seed);
      }
      BLAS_zdot_z_d_testgen(i, 0, i, norm, blas_no_conj, alpha_i,
			    1, alpha_i, 1, x_vec, t_vec,
			    seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j < i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Set the diagonal element to 1. */
      t_elem = 1.0;
      T_i[tij] = t_elem;

      /* Set x[i] to be r. */
      x_i[xi] = r[0];
      x_i[xi + 1] = r[1];
      x_vec[xvec_i] = r[0];
      x_vec[xvec_i + 1] = r[1];

    } else {
      BLAS_zdot_z_d_testgen(i + 1, 0, i, norm, blas_no_conj, alpha,
			    (i == 0 ? alpha_flag : 1), beta, 1, x_vec, t_vec,
			    seed, r, head_r_true_elem, tail_r_true_elem);

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n - 1) * inc_tij);
      for (j = 0; j <= i; j++) {
	t_elem = t_vec[tvec_j];

	T_i[tij] = t_elem;
	tvec_j += inc_tvec;
	tij += inc_tij;
      }

      /* Copy generated x_vec[i] to appropriate position in x. */
      x_elem[0] = x_vec[xvec_i];
      x_elem[1] = x_vec[xvec_i + 1];
      x_i[xi] = x_elem[0];
      x_i[xi + 1] = x_elem[1];
    }

    /* Copy r_true */
    head_r_true[xi] = head_r_true_elem[0];
    head_r_true[xi + 1] = head_r_true_elem[1];
    tail_r_true[xi] = tail_r_true_elem[0];
    tail_r_true[xi + 1] = tail_r_true_elem[1];

    xvec_i += inc_xvec;
    xi += inc_xi;

    ti += inc_ti;
  }

  blas_free(x_vec);
  blas_free(t_vec);
}

  /* end of BLAS_ztrmv_d_testgen */
