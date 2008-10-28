#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

void BLAS_stbsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, int k, int randomize,
			float *alpha, int alpha_flag, float *T, int ldt,
			float *x, int *seed, double *head_r_true,
			double *tail_r_true, int row,
			enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular Banded matrix; 
 *      and computes r_true.
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
 *              Dimension of T and the length of vector x
 *
 * k            (input) int
 *              Number of super(sub)diagonals of T
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
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  int start;
  int length;
  int i, j;
  int tempi, inc_tempi;
  float alpha_i;
  float minus_one;
  float Tii;
  float *temp;
  float *xtemp2;

  temp = (float *) blas_malloc(n * 2 * sizeof(float));
  if (n * 2 > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* always allocate, not always needed */
  xtemp2 = (float *) blas_malloc(n * 2 * sizeof(float));
  if (n * 2 > 0 && xtemp2 == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  minus_one = -1.0;

  /* if alpha_flag=0, gives a random value to alpha */
  if (alpha_flag == 0) {
    alpha_i = xrand(seed);
    *alpha = alpha_i;
    alpha_flag = 1;
  }

  for (i = 0; i < 4 * n * ldt; i++) {
    T[i] = 0.0;
  }

  inc_tempi = 1;


  for (i = 0; i < n; i++) {

    if (i != row) {
      if (diag == blas_non_unit_diag) {
	Tii = xrand(seed);
      } else {
	Tii = 1.0;
      }
      for (j = 0, tempi = 0; j < n; j++, tempi += inc_tempi) {
	if (j != i) {
	  temp[tempi] = 0.0;
	} else {
	  temp[tempi] = Tii;
	}
      }
      stbsv_commit(order, uplo, trans, n, k, T, ldt, temp, i);


      x[i] = xrand(seed);

      switch (prec) {
      case blas_prec_single:
	{
	  float multemp;
	  float divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  head_r_true[i] = divtemp;
	  tail_r_true[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_indigenous:
      case blas_prec_double:
	{
	  double multemp;
	  double divtemp;

	  multemp = (double) x[i] * *alpha;
	  divtemp = (double) multemp / Tii;
	  head_r_true[i] = divtemp;
	  tail_r_true[i] = 0.0;
	  break;
	}
      case blas_prec_extra:
	{
	  double head_multemp, tail_multemp;
	  double head_divtemp, tail_divtemp;

	  head_multemp = (double) x[i] * *alpha;
	  tail_multemp = 0.0;
	  {
	    double dt = (double) Tii;
	    {
	      /* Compute double-double = double-double / double,
	         using a Newton iteration scheme. */
	      double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

	      /* Compute a DP approximation to the quotient. */
	      t1 = head_multemp / dt;

	      /* Split t1 and b into two parts with at most 26 bits each,
	         using the Dekker-Veltkamp method. */
	      con = t1 * split;
	      t11 = con - (con - t1);
	      t21 = t1 - t11;
	      con = dt * split;
	      b1 = con - (con - dt);
	      b2 = dt - b1;

	      /* Compute t1 * b using Dekker method. */
	      t12 = t1 * dt;
	      t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

	      /* Compute dda - (t12, t22) using Knuth trick. */
	      t11 = head_multemp - t12;
	      e = t11 - head_multemp;
	      t21 =
		((-t12 - e) + (head_multemp - (t11 - e))) + tail_multemp -
		t22;

	      /* Compute high-order word of (t11, t21) and divide by b. */
	      t2 = (t11 + t21) / dt;

	      /* The result is t1 + t2, after normalization. */
	      head_divtemp = t1 + t2;
	      tail_divtemp = t2 - (head_divtemp - t1);
	    }
	  }
	  head_r_true[i] = head_divtemp;
	  tail_r_true[i] = tail_divtemp;
	  break;
	}
      }				/* case */
    }				/* if */
  }				/* for */

  for (j = 0; j < n; j++) {
    temp[j] = 0.0;
  }

  /* now set the important row - */
  Tii = 1.0;
  for (j = 0, tempi = 0; j < n; j++, tempi += inc_tempi) {
    if (j != row) {
      temp[tempi] = 0.0;
    } else {
      temp[tempi] = Tii;
    }
  }
  /* this is extra, it will really be committed later */
  stbsv_commit(order, uplo, trans, n, k, T, ldt, temp, row);


  if ((uplo == blas_lower && (trans == blas_no_trans || trans == blas_conj))
      || (uplo == blas_upper
	  && (trans == blas_trans || trans == blas_conj_trans))) {
    length = MIN(row, k);
    start = MAX(0, row - k);
  } else {
    length = MIN(n - row - 1, k);
    start = row + 1;
  }

  if (length != 0) {

    switch (prec) {
    case blas_prec_single:
      BLAS_sdot_testgen(length, 0, length, norm,
			blas_no_conj, &minus_one, 1, alpha, 1,
			&xtemp2[start], &temp[start], seed, &x[row],
			&head_r_true[row], &tail_r_true[row]);
      break;
    case blas_prec_indigenous:
    case blas_prec_double:
    case blas_prec_extra:
      BLAS_sdot_x_testgen(length, 0, length, norm,
			  blas_no_conj, &minus_one, 1, alpha, 1,
			  &head_r_true[start], &tail_r_true[start],
			  &temp[start],
			  seed, &x[row], &head_r_true[row],
			  &tail_r_true[row]);
      break;
    }
    stbsv_commit(order, uplo, trans, n, k, T, ldt, temp, row);
  } else {
    /* probably low k- must use tricks to test */
    /* tricks are not yet implemented */
    x[row] = xrand(seed);

    switch (prec) {
    case blas_prec_single:
      {
	float multemp;

	multemp = x[row] * *alpha;
	head_r_true[row] = multemp;
	tail_r_true[row] = 0.0;
	break;
      }
    case blas_prec_indigenous:
    case blas_prec_double:
      {
	double multemp;

	multemp = (double) x[row] * *alpha;
	head_r_true[row] = multemp;
	tail_r_true[row] = 0.0;
	break;
      }
    case blas_prec_extra:
      {
	double head_multemp, tail_multemp;

	head_multemp = (double) x[row] * *alpha;
	tail_multemp = 0.0;
	head_r_true[row] = head_multemp;
	tail_r_true[row] = tail_multemp;
	break;
      }
    }
  }

  blas_free(temp);
  blas_free(xtemp2);
}

void BLAS_dtbsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, int k, int randomize,
			double *alpha, int alpha_flag, double *T, int ldt,
			double *x, int *seed, double *head_r_true,
			double *tail_r_true, int row,
			enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular Banded matrix; 
 *      and computes r_true.
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
 *              Dimension of T and the length of vector x
 *
 * k            (input) int
 *              Number of super(sub)diagonals of T
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
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  int start;
  int length;
  int i, j;
  int tempi, inc_tempi;
  float alpha_i;
  double minus_one;
  double Tii;
  double *temp;
  double *xtemp2;

  temp = (double *) blas_malloc(n * 2 * sizeof(double));
  if (n * 2 > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* always allocate, not always needed */
  xtemp2 = (double *) blas_malloc(n * 2 * sizeof(double));
  if (n * 2 > 0 && xtemp2 == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  minus_one = -1.0;

  /* if alpha_flag=0, gives a random value to alpha */
  if (alpha_flag == 0) {
    alpha_i = xrand(seed);
    *alpha = alpha_i;
    alpha_flag = 1;
  }

  for (i = 0; i < 4 * n * ldt; i++) {
    T[i] = 0.0;
  }

  inc_tempi = 1;


  for (i = 0; i < n; i++) {

    if (i != row) {
      if (diag == blas_non_unit_diag) {
	Tii = xrand(seed);
      } else {
	Tii = 1.0;
      }
      for (j = 0, tempi = 0; j < n; j++, tempi += inc_tempi) {
	if (j != i) {
	  temp[tempi] = 0.0;
	} else {
	  temp[tempi] = Tii;
	}
      }
      dtbsv_commit(order, uplo, trans, n, k, T, ldt, temp, i);


      x[i] = xrand(seed);

      switch (prec) {
      case blas_prec_single:
	{
	  double multemp;
	  double divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  head_r_true[i] = divtemp;
	  tail_r_true[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_indigenous:
      case blas_prec_double:
	{
	  double multemp;
	  double divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  head_r_true[i] = divtemp;
	  tail_r_true[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_extra:
	{
	  double head_multemp, tail_multemp;
	  double head_divtemp, tail_divtemp;

	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = x[i] * split;
	    a1 = con - x[i];
	    a1 = con - a1;
	    a2 = x[i] - a1;
	    con = *alpha * split;
	    b1 = con - *alpha;
	    b1 = con - b1;
	    b2 = *alpha - b1;

	    head_multemp = x[i] * *alpha;
	    tail_multemp =
	      (((a1 * b1 - head_multemp) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double-double = double-double / double,
	       using a Newton iteration scheme. */
	    double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

	    /* Compute a DP approximation to the quotient. */
	    t1 = head_multemp / Tii;

	    /* Split t1 and b into two parts with at most 26 bits each,
	       using the Dekker-Veltkamp method. */
	    con = t1 * split;
	    t11 = con - (con - t1);
	    t21 = t1 - t11;
	    con = Tii * split;
	    b1 = con - (con - Tii);
	    b2 = Tii - b1;

	    /* Compute t1 * b using Dekker method. */
	    t12 = t1 * Tii;
	    t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

	    /* Compute dda - (t12, t22) using Knuth trick. */
	    t11 = head_multemp - t12;
	    e = t11 - head_multemp;
	    t21 =
	      ((-t12 - e) + (head_multemp - (t11 - e))) + tail_multemp - t22;

	    /* Compute high-order word of (t11, t21) and divide by b. */
	    t2 = (t11 + t21) / Tii;

	    /* The result is t1 + t2, after normalization. */
	    head_divtemp = t1 + t2;
	    tail_divtemp = t2 - (head_divtemp - t1);
	  }
	  head_r_true[i] = head_divtemp;
	  tail_r_true[i] = tail_divtemp;
	  break;
	}
      }				/* case */
    }				/* if */
  }				/* for */

  for (j = 0; j < n; j++) {
    temp[j] = 0.0;
  }

  /* now set the important row - */
  Tii = 1.0;
  for (j = 0, tempi = 0; j < n; j++, tempi += inc_tempi) {
    if (j != row) {
      temp[tempi] = 0.0;
    } else {
      temp[tempi] = Tii;
    }
  }
  /* this is extra, it will really be committed later */
  dtbsv_commit(order, uplo, trans, n, k, T, ldt, temp, row);


  if ((uplo == blas_lower && (trans == blas_no_trans || trans == blas_conj))
      || (uplo == blas_upper
	  && (trans == blas_trans || trans == blas_conj_trans))) {
    length = MIN(row, k);
    start = MAX(0, row - k);
  } else {
    length = MIN(n - row - 1, k);
    start = row + 1;
  }

  if (length != 0) {

    switch (prec) {
    case blas_prec_single:
      BLAS_ddot_testgen(length, 0, length, norm,
			blas_no_conj, &minus_one, 1, alpha, 1,
			&xtemp2[start], &temp[start], seed, &x[row],
			&head_r_true[row], &tail_r_true[row]);
      break;
    case blas_prec_indigenous:
    case blas_prec_double:
    case blas_prec_extra:
      BLAS_ddot_x_testgen(length, 0, length, norm,
			  blas_no_conj, &minus_one, 1, alpha, 1,
			  &head_r_true[start], &tail_r_true[start],
			  &temp[start],
			  seed, &x[row], &head_r_true[row],
			  &tail_r_true[row]);
      break;
    }
    dtbsv_commit(order, uplo, trans, n, k, T, ldt, temp, row);
  } else {
    /* probably low k- must use tricks to test */
    /* tricks are not yet implemented */
    x[row] = xrand(seed);

    switch (prec) {
    case blas_prec_single:
      {
	double multemp;

	multemp = x[row] * *alpha;
	head_r_true[row] = multemp;
	tail_r_true[row] = 0.0;
	break;
      }
    case blas_prec_indigenous:
    case blas_prec_double:
      {
	double multemp;

	multemp = x[row] * *alpha;
	head_r_true[row] = multemp;
	tail_r_true[row] = 0.0;
	break;
      }
    case blas_prec_extra:
      {
	double head_multemp, tail_multemp;

	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x[row] * split;
	  a1 = con - x[row];
	  a1 = con - a1;
	  a2 = x[row] - a1;
	  con = *alpha * split;
	  b1 = con - *alpha;
	  b1 = con - b1;
	  b2 = *alpha - b1;

	  head_multemp = x[row] * *alpha;
	  tail_multemp =
	    (((a1 * b1 - head_multemp) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	head_r_true[row] = head_multemp;
	tail_r_true[row] = tail_multemp;
	break;
      }
    }
  }

  blas_free(temp);
  blas_free(xtemp2);
}

void BLAS_dtbsv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, int k,
			  int randomize, double *alpha, int alpha_flag,
			  float *T, int ldt, double *x, int *seed,
			  double *head_r_true, double *tail_r_true, int row,
			  enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular Banded matrix; 
 *      and computes r_true.
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
 *              Dimension of T and the length of vector x
 *
 * k            (input) int
 *              Number of super(sub)diagonals of T
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
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  int start;
  int length;
  int i, j;
  int tempi, inc_tempi;
  float alpha_i;
  double minus_one;
  float Tii;
  float *temp;
  double *xtemp2;

  temp = (float *) blas_malloc(n * 2 * sizeof(float));
  if (n * 2 > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  /* always allocate, not always needed */
  xtemp2 = (double *) blas_malloc(n * 2 * sizeof(double));
  if (n * 2 > 0 && xtemp2 == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  minus_one = -1.0;

  /* if alpha_flag=0, gives a random value to alpha */
  if (alpha_flag == 0) {
    alpha_i = xrand(seed);
    *alpha = alpha_i;
    alpha_flag = 1;
  }

  for (i = 0; i < 4 * n * ldt; i++) {
    T[i] = 0.0;
  }

  inc_tempi = 1;


  for (i = 0; i < n; i++) {

    if (i != row) {
      if (diag == blas_non_unit_diag) {
	Tii = xrand(seed);
      } else {
	Tii = 1.0;
      }
      for (j = 0, tempi = 0; j < n; j++, tempi += inc_tempi) {
	if (j != i) {
	  temp[tempi] = 0.0;
	} else {
	  temp[tempi] = Tii;
	}
      }
      stbsv_commit(order, uplo, trans, n, k, T, ldt, temp, i);


      x[i] = xrand(seed);

      switch (prec) {
      case blas_prec_single:
	{
	  double multemp;
	  double divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  head_r_true[i] = divtemp;
	  tail_r_true[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_indigenous:
      case blas_prec_double:
	{
	  double multemp;
	  double divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  head_r_true[i] = divtemp;
	  tail_r_true[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_extra:
	{
	  double head_multemp, tail_multemp;
	  double head_divtemp, tail_divtemp;

	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = x[i] * split;
	    a1 = con - x[i];
	    a1 = con - a1;
	    a2 = x[i] - a1;
	    con = *alpha * split;
	    b1 = con - *alpha;
	    b1 = con - b1;
	    b2 = *alpha - b1;

	    head_multemp = x[i] * *alpha;
	    tail_multemp =
	      (((a1 * b1 - head_multemp) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    double dt = (double) Tii;
	    {
	      /* Compute double-double = double-double / double,
	         using a Newton iteration scheme. */
	      double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

	      /* Compute a DP approximation to the quotient. */
	      t1 = head_multemp / dt;

	      /* Split t1 and b into two parts with at most 26 bits each,
	         using the Dekker-Veltkamp method. */
	      con = t1 * split;
	      t11 = con - (con - t1);
	      t21 = t1 - t11;
	      con = dt * split;
	      b1 = con - (con - dt);
	      b2 = dt - b1;

	      /* Compute t1 * b using Dekker method. */
	      t12 = t1 * dt;
	      t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

	      /* Compute dda - (t12, t22) using Knuth trick. */
	      t11 = head_multemp - t12;
	      e = t11 - head_multemp;
	      t21 =
		((-t12 - e) + (head_multemp - (t11 - e))) + tail_multemp -
		t22;

	      /* Compute high-order word of (t11, t21) and divide by b. */
	      t2 = (t11 + t21) / dt;

	      /* The result is t1 + t2, after normalization. */
	      head_divtemp = t1 + t2;
	      tail_divtemp = t2 - (head_divtemp - t1);
	    }
	  }
	  head_r_true[i] = head_divtemp;
	  tail_r_true[i] = tail_divtemp;
	  break;
	}
      }				/* case */
    }				/* if */
  }				/* for */

  for (j = 0; j < n; j++) {
    temp[j] = 0.0;
  }

  /* now set the important row - */
  Tii = 1.0;
  for (j = 0, tempi = 0; j < n; j++, tempi += inc_tempi) {
    if (j != row) {
      temp[tempi] = 0.0;
    } else {
      temp[tempi] = Tii;
    }
  }
  /* this is extra, it will really be committed later */
  stbsv_commit(order, uplo, trans, n, k, T, ldt, temp, row);


  if ((uplo == blas_lower && (trans == blas_no_trans || trans == blas_conj))
      || (uplo == blas_upper
	  && (trans == blas_trans || trans == blas_conj_trans))) {
    length = MIN(row, k);
    start = MAX(0, row - k);
  } else {
    length = MIN(n - row - 1, k);
    start = row + 1;
  }

  if (length != 0) {

    switch (prec) {
    case blas_prec_single:
    case blas_prec_indigenous:
    case blas_prec_double:	/*BLAS_ddot_s_x_testgen(length, 0, length, norm, 
				   blas_no_conj, &minus_one, 1, alpha, 1, 
				   &head_r_true[start], &tail_r_true[start], 
				   &temp[start], 
				   seed, &x[row], &head_r_true[row], &tail_r_true[row]);
				   break; */
    case blas_prec_extra:
      BLAS_ddot_s_x_testgen(length, 0, length, norm,
			    blas_no_conj, &minus_one, 1, alpha, 1,
			    &head_r_true[start], &tail_r_true[start],
			    &temp[start],
			    seed, &x[row], &head_r_true[row],
			    &tail_r_true[row]);
      break;
    }
    stbsv_commit(order, uplo, trans, n, k, T, ldt, temp, row);
  } else {
    /* probably low k- must use tricks to test */
    /* tricks are not yet implemented */
    x[row] = xrand(seed);

    switch (prec) {
    case blas_prec_single:
      {
	double multemp;

	multemp = x[row] * *alpha;
	head_r_true[row] = multemp;
	tail_r_true[row] = 0.0;
	break;
      }
    case blas_prec_indigenous:
    case blas_prec_double:
      {
	double multemp;

	multemp = x[row] * *alpha;
	head_r_true[row] = multemp;
	tail_r_true[row] = 0.0;
	break;
      }
    case blas_prec_extra:
      {
	double head_multemp, tail_multemp;

	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x[row] * split;
	  a1 = con - x[row];
	  a1 = con - a1;
	  a2 = x[row] - a1;
	  con = *alpha * split;
	  b1 = con - *alpha;
	  b1 = con - b1;
	  b2 = *alpha - b1;

	  head_multemp = x[row] * *alpha;
	  tail_multemp =
	    (((a1 * b1 - head_multemp) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	head_r_true[row] = head_multemp;
	tail_r_true[row] = tail_multemp;
	break;
      }
    }
  }

  blas_free(temp);
  blas_free(xtemp2);
}

void BLAS_ctbsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, int k, int randomize,
			void *alpha, int alpha_flag, void *T, int ldt,
			void *x, int *seed, double *head_r_true,
			double *tail_r_true, int row,
			enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular Banded matrix; 
 *      and computes r_true.
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
 *              Dimension of T and the length of vector x
 *
 * k            (input) int
 *              Number of super(sub)diagonals of T
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
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  float *x_i = (float *) x;
  float *alpha_i = (float *) alpha;
  float *T_i = (float *) T;
  float alpha_r;
  float *T_r;
  float *x_r;
  float *T_temp_r;
  float *T_temp_c;
  double *head_r_true_r, *tail_r_true_r;
  int i, inc = 2, length, j;

  T_r = (float *) blas_malloc(8 * n * ldt * sizeof(float));
  if (8 * n * ldt > 0 && T_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_r = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true_r = (double *) blas_malloc(n * sizeof(double));
  tail_r_true_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && (head_r_true_r == NULL || tail_r_true_r == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  T_temp_c = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && T_temp_c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  T_temp_r = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && T_temp_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  if ((uplo == blas_lower && trans == blas_no_trans) ||
      (uplo == blas_upper && trans != blas_no_trans)) {
    length = row;
  } else {
    length = n - row - 1;
  }

  BLAS_stbsv_testgen(norm, order, uplo, trans, diag, n, k, randomize,
		     &alpha_r,
		     alpha_flag, T_r, ldt, x_r, seed, head_r_true_r,
		     tail_r_true_r, row, prec);

  /* multiply alpha_r by 1+i */
  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  if (diag == blas_non_unit_diag) {
    for (i = 0; i < n; i++) {
      /* multiply x_r by i */
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      /* multiply non major rows truth by i */
      if (i != row) {
	head_r_true[i * inc] = 0.0;
	head_r_true[i * inc + 1] = head_r_true_r[i];
	tail_r_true[i * inc] = 0.0;
	tail_r_true[i * inc + 1] = tail_r_true_r[i];
      } else {
	/* multiply major rows truth by (-1+i) */
	head_r_true[i * inc] = -head_r_true_r[i];
	head_r_true[i * inc + 1] = head_r_true_r[i];
	tail_r_true[i * inc] = -tail_r_true_r[i];
	tail_r_true[i * inc + 1] = tail_r_true_r[i];
      }
    }

    /* copy T real to complex (and multiply by (1+i)) */
    for (i = 0; i < n; i++) {
      stbsv_copy(order, uplo, trans, n, k, T_r, ldt, T_temp_r, i);
      for (j = 0; j < n; j++) {
	T_temp_c[j * inc] = T_temp_r[j];
	/* conjugation is handled by the commit later - */
	T_temp_c[j * inc + 1] = T_temp_r[j];
      }
      if (i == row) {
	/* zero out imaginary part of diagonal on important row. */
	T_temp_c[row * inc + 1] = 0.0;
      }
      ctbsv_commit(order, uplo, trans, n, k, T_i, ldt, T_temp_c, i);
    }
  } else {
    for (i = 0; i < n; i++) {
      /* multiply x by i */
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      /* multiply non-major rows by -1+i */
      if (i != row || length == 0) {
	head_r_true[i * inc] = -head_r_true_r[i];
	head_r_true[i * inc + 1] = head_r_true_r[i];
	tail_r_true[i * inc] = -tail_r_true_r[i];
	tail_r_true[i * inc + 1] = tail_r_true_r[i];
      } else {
	/* multiply x by 1+i (overwrite multiply by i earlier) */
	x_i[i * inc] = x_r[i];
	x_i[i * inc + 1] = x_r[i];

	head_r_true[i * inc] = 0.0;
	head_r_true[i * inc + 1] = 2 * head_r_true_r[i];
	tail_r_true[i * inc] = 0.0;
	tail_r_true[i * inc + 1] = 2 * tail_r_true_r[i];
      }
    }

    /* copy T real to complex (and multiply by (1-i)) */
    for (i = 0; i < n; i++) {
      stbsv_copy(order, uplo, trans, n, k, T_r, ldt, T_temp_r, i);
      for (j = 0; j < n; j++) {
	T_temp_c[j * inc] = T_temp_r[j];
	/* conjugation is handled by the commit later - */
	T_temp_c[j * inc + 1] = -T_temp_r[j];
      }
      /* zero out imaginary part of diagonal on every row. */
      T_temp_c[i * inc + 1] = 0.0;
      ctbsv_commit(order, uplo, trans, n, k, T_i, ldt, T_temp_c, i);
    }
  }

  blas_free(T_temp_c);
  blas_free(T_temp_r);
  blas_free(T_r);
  blas_free(x_r);
  blas_free(head_r_true_r);
  blas_free(tail_r_true_r);
}

void BLAS_ztbsv_c_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, int k,
			  int randomize, void *alpha, int alpha_flag, void *T,
			  int ldt, void *x, int *seed, double *head_r_true,
			  double *tail_r_true, int row,
			  enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular Banded matrix; 
 *      and computes r_true.
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
 *              Dimension of T and the length of vector x
 *
 * k            (input) int
 *              Number of super(sub)diagonals of T
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
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  double *x_i = (double *) x;
  double *alpha_i = (double *) alpha;
  float *T_i = (float *) T;
  double alpha_r;
  float *T_r;
  double *x_r;
  float *T_temp_r;
  float *T_temp_c;
  double *head_r_true_r, *tail_r_true_r;
  int i, inc = 2, length, j;

  T_r = (float *) blas_malloc(8 * n * ldt * sizeof(float));
  if (8 * n * ldt > 0 && T_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true_r = (double *) blas_malloc(n * sizeof(double));
  tail_r_true_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && (head_r_true_r == NULL || tail_r_true_r == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  T_temp_c = (float *) blas_malloc(n * sizeof(float) * 2);
  if (n > 0 && T_temp_c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  T_temp_r = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && T_temp_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  if ((uplo == blas_lower && trans == blas_no_trans) ||
      (uplo == blas_upper && trans != blas_no_trans)) {
    length = row;
  } else {
    length = n - row - 1;
  }

  BLAS_dtbsv_s_testgen(norm, order, uplo, trans, diag, n, k, randomize,
		       &alpha_r,
		       alpha_flag, T_r, ldt, x_r, seed, head_r_true_r,
		       tail_r_true_r, row, prec);

  /* multiply alpha_r by 1+i */
  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  if (diag == blas_non_unit_diag) {
    for (i = 0; i < n; i++) {
      /* multiply x_r by i */
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      /* multiply non major rows truth by i */
      if (i != row) {
	head_r_true[i * inc] = 0.0;
	head_r_true[i * inc + 1] = head_r_true_r[i];
	tail_r_true[i * inc] = 0.0;
	tail_r_true[i * inc + 1] = tail_r_true_r[i];
      } else {
	/* multiply major rows truth by (-1+i) */
	head_r_true[i * inc] = -head_r_true_r[i];
	head_r_true[i * inc + 1] = head_r_true_r[i];
	tail_r_true[i * inc] = -tail_r_true_r[i];
	tail_r_true[i * inc + 1] = tail_r_true_r[i];
      }
    }

    /* copy T real to complex (and multiply by (1+i)) */
    for (i = 0; i < n; i++) {
      stbsv_copy(order, uplo, trans, n, k, T_r, ldt, T_temp_r, i);
      for (j = 0; j < n; j++) {
	T_temp_c[j * inc] = T_temp_r[j];
	/* conjugation is handled by the commit later - */
	T_temp_c[j * inc + 1] = T_temp_r[j];
      }
      if (i == row) {
	/* zero out imaginary part of diagonal on important row. */
	T_temp_c[row * inc + 1] = 0.0;
      }
      ctbsv_commit(order, uplo, trans, n, k, T_i, ldt, T_temp_c, i);
    }
  } else {
    for (i = 0; i < n; i++) {
      /* multiply x by i */
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      /* multiply non-major rows by -1+i */
      if (i != row || length == 0) {
	head_r_true[i * inc] = -head_r_true_r[i];
	head_r_true[i * inc + 1] = head_r_true_r[i];
	tail_r_true[i * inc] = -tail_r_true_r[i];
	tail_r_true[i * inc + 1] = tail_r_true_r[i];
      } else {
	/* multiply x by 1+i (overwrite multiply by i earlier) */
	x_i[i * inc] = x_r[i];
	x_i[i * inc + 1] = x_r[i];

	head_r_true[i * inc] = 0.0;
	head_r_true[i * inc + 1] = 2 * head_r_true_r[i];
	tail_r_true[i * inc] = 0.0;
	tail_r_true[i * inc + 1] = 2 * tail_r_true_r[i];
      }
    }

    /* copy T real to complex (and multiply by (1-i)) */
    for (i = 0; i < n; i++) {
      stbsv_copy(order, uplo, trans, n, k, T_r, ldt, T_temp_r, i);
      for (j = 0; j < n; j++) {
	T_temp_c[j * inc] = T_temp_r[j];
	/* conjugation is handled by the commit later - */
	T_temp_c[j * inc + 1] = -T_temp_r[j];
      }
      /* zero out imaginary part of diagonal on every row. */
      T_temp_c[i * inc + 1] = 0.0;
      ctbsv_commit(order, uplo, trans, n, k, T_i, ldt, T_temp_c, i);
    }
  }

  blas_free(T_temp_c);
  blas_free(T_temp_r);
  blas_free(T_r);
  blas_free(x_r);
  blas_free(head_r_true_r);
  blas_free(tail_r_true_r);
}

void BLAS_ztbsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, int k, int randomize,
			void *alpha, int alpha_flag, void *T, int ldt,
			void *x, int *seed, double *head_r_true,
			double *tail_r_true, int row,
			enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular Banded matrix; 
 *      and computes r_true.
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
 *              Dimension of T and the length of vector x
 *
 * k            (input) int
 *              Number of super(sub)diagonals of T
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
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  double *x_i = (double *) x;
  double *alpha_i = (double *) alpha;
  double *T_i = (double *) T;
  double alpha_r;
  double *T_r;
  double *x_r;
  double *T_temp_r;
  double *T_temp_c;
  double *head_r_true_r, *tail_r_true_r;
  int i, inc = 2, length, j;

  T_r = (double *) blas_malloc(8 * n * ldt * sizeof(double));
  if (8 * n * ldt > 0 && T_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  x_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true_r = (double *) blas_malloc(n * sizeof(double));
  tail_r_true_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && (head_r_true_r == NULL || tail_r_true_r == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  T_temp_c = (double *) blas_malloc(n * sizeof(double) * 2);
  if (n > 0 && T_temp_c == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  T_temp_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && T_temp_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  if ((uplo == blas_lower && trans == blas_no_trans) ||
      (uplo == blas_upper && trans != blas_no_trans)) {
    length = row;
  } else {
    length = n - row - 1;
  }

  BLAS_dtbsv_testgen(norm, order, uplo, trans, diag, n, k, randomize,
		     &alpha_r,
		     alpha_flag, T_r, ldt, x_r, seed, head_r_true_r,
		     tail_r_true_r, row, prec);

  /* multiply alpha_r by 1+i */
  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  if (diag == blas_non_unit_diag) {
    for (i = 0; i < n; i++) {
      /* multiply x_r by i */
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      /* multiply non major rows truth by i */
      if (i != row) {
	head_r_true[i * inc] = 0.0;
	head_r_true[i * inc + 1] = head_r_true_r[i];
	tail_r_true[i * inc] = 0.0;
	tail_r_true[i * inc + 1] = tail_r_true_r[i];
      } else {
	/* multiply major rows truth by (-1+i) */
	head_r_true[i * inc] = -head_r_true_r[i];
	head_r_true[i * inc + 1] = head_r_true_r[i];
	tail_r_true[i * inc] = -tail_r_true_r[i];
	tail_r_true[i * inc + 1] = tail_r_true_r[i];
      }
    }

    /* copy T real to complex (and multiply by (1+i)) */
    for (i = 0; i < n; i++) {
      dtbsv_copy(order, uplo, trans, n, k, T_r, ldt, T_temp_r, i);
      for (j = 0; j < n; j++) {
	T_temp_c[j * inc] = T_temp_r[j];
	/* conjugation is handled by the commit later - */
	T_temp_c[j * inc + 1] = T_temp_r[j];
      }
      if (i == row) {
	/* zero out imaginary part of diagonal on important row. */
	T_temp_c[row * inc + 1] = 0.0;
      }
      ztbsv_commit(order, uplo, trans, n, k, T_i, ldt, T_temp_c, i);
    }
  } else {
    for (i = 0; i < n; i++) {
      /* multiply x by i */
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      /* multiply non-major rows by -1+i */
      if (i != row || length == 0) {
	head_r_true[i * inc] = -head_r_true_r[i];
	head_r_true[i * inc + 1] = head_r_true_r[i];
	tail_r_true[i * inc] = -tail_r_true_r[i];
	tail_r_true[i * inc + 1] = tail_r_true_r[i];
      } else {
	/* multiply x by 1+i (overwrite multiply by i earlier) */
	x_i[i * inc] = x_r[i];
	x_i[i * inc + 1] = x_r[i];

	head_r_true[i * inc] = 0.0;
	head_r_true[i * inc + 1] = 2 * head_r_true_r[i];
	tail_r_true[i * inc] = 0.0;
	tail_r_true[i * inc + 1] = 2 * tail_r_true_r[i];
      }
    }

    /* copy T real to complex (and multiply by (1-i)) */
    for (i = 0; i < n; i++) {
      dtbsv_copy(order, uplo, trans, n, k, T_r, ldt, T_temp_r, i);
      for (j = 0; j < n; j++) {
	T_temp_c[j * inc] = T_temp_r[j];
	/* conjugation is handled by the commit later - */
	T_temp_c[j * inc + 1] = -T_temp_r[j];
      }
      /* zero out imaginary part of diagonal on every row. */
      T_temp_c[i * inc + 1] = 0.0;
      ztbsv_commit(order, uplo, trans, n, k, T_i, ldt, T_temp_c, i);
    }
  }

  blas_free(T_temp_c);
  blas_free(T_temp_r);
  blas_free(T_r);
  blas_free(x_r);
  blas_free(head_r_true_r);
  blas_free(tail_r_true_r);
}

void BLAS_ctbsv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, int k,
			  int randomize, void *alpha, int alpha_flag,
			  float *T, int ldt, void *x, int *seed,
			  double *head_r_true, double *tail_r_true, int row,
			  enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular Banded matrix; 
 *      and computes r_true.
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
 *              Dimension of T and the length of vector x
 *
 * k            (input) int
 *              Number of super(sub)diagonals of T
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
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  float *x_i = (float *) x;
  float *alpha_i = (float *) alpha;
  float *T_i = T;
  float alpha_r;
  float *x_r;
  double *head_r_true_r, *tail_r_true_r;
  int i, inc = 2;

  x_r = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true_r = (double *) blas_malloc(n * sizeof(double));
  tail_r_true_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && (head_r_true_r == NULL || tail_r_true_r == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  BLAS_stbsv_testgen(norm, order, uplo, trans, diag, n, k,
		     randomize, &alpha_r,
		     alpha_flag, T_i, ldt, x_r, seed, head_r_true_r,
		     tail_r_true_r, row, prec);

  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  for (i = 0; i < n; i++) {
    x_i[i * inc] = 0.0;
    x_i[i * inc + 1] = x_r[i];

    head_r_true[i * inc] = -head_r_true_r[i];
    head_r_true[i * inc + 1] = head_r_true_r[i];
    tail_r_true[i * inc] = -tail_r_true_r[i];
    tail_r_true[i * inc + 1] = tail_r_true_r[i];
  }

  blas_free(x_r);
  blas_free(head_r_true_r);
  blas_free(tail_r_true_r);
}

void BLAS_ztbsv_d_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, int k,
			  int randomize, void *alpha, int alpha_flag,
			  double *T, int ldt, void *x, int *seed,
			  double *head_r_true, double *tail_r_true, int row,
			  enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular Banded matrix; 
 *      and computes r_true.
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
 *              Dimension of T and the length of vector x
 *
 * k            (input) int
 *              Number of super(sub)diagonals of T
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
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  double *x_i = (double *) x;
  double *alpha_i = (double *) alpha;
  double *T_i = T;
  double alpha_r;
  double *x_r;
  double *head_r_true_r, *tail_r_true_r;
  int i, inc = 2;

  x_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_r == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  head_r_true_r = (double *) blas_malloc(n * sizeof(double));
  tail_r_true_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && (head_r_true_r == NULL || tail_r_true_r == NULL)) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  BLAS_dtbsv_testgen(norm, order, uplo, trans, diag, n, k,
		     randomize, &alpha_r,
		     alpha_flag, T_i, ldt, x_r, seed, head_r_true_r,
		     tail_r_true_r, row, prec);

  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  for (i = 0; i < n; i++) {
    x_i[i * inc] = 0.0;
    x_i[i * inc + 1] = x_r[i];

    head_r_true[i * inc] = -head_r_true_r[i];
    head_r_true[i * inc + 1] = head_r_true_r[i];
    tail_r_true[i * inc] = -tail_r_true_r[i];
    tail_r_true[i * inc + 1] = tail_r_true_r[i];
  }

  blas_free(x_r);
  blas_free(head_r_true_r);
  blas_free(tail_r_true_r);
}
