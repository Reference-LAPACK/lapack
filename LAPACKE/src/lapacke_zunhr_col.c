#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_zunhr_col)( int matrix_layout, lapack_int m, lapack_int n,
                              lapack_int nb, lapack_complex_double* a,
                              lapack_int lda, lapack_complex_double* t,
                              lapack_int ldt, lapack_complex_double* d)
{
  lapack_int info = 0;
  if( matrix_layout != LAPACK_COL_MAJOR && matrix_layout != LAPACK_ROW_MAJOR ) {
    API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_zunhr_col", -1 );
    return -1;
  }
#ifndef LAPACK_DISABLE_NAN_CHECK
  if( LAPACKE_get_nancheck() ) {
    /* Optionally check input matrices for NaNs */
    if( API_SUFFIX(LAPACKE_zge_nancheck)( matrix_layout, m, n, a, lda ) ) {
      return -5;
    }
  }
#endif
  /* Call middle-level interface */
  info = API_SUFFIX(LAPACKE_zunhr_col_work)( matrix_layout, m, n, nb, a, lda, t, ldt, d );
  return info;
}
