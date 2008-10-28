      SUBROUTINE CLAPLL( N, X, INCX, Y, INCY, SSMIN )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      REAL               SSMIN
*     ..
*     .. Array Arguments ..
      COMPLEX            X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  Given two column vectors X and Y, let
*
*                       A = ( X Y ).
*
*  The subroutine first computes the QR factorization of A = Q*R,
*  and then computes the SVD of the 2-by-2 upper triangular matrix R.
*  The smaller singular value of R is returned in SSMIN, which is used
*  as the measurement of the linear dependency of the vectors X and Y.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vectors X and Y.
*
*  X       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)
*          On entry, X contains the N-vector X.
*          On exit, X is overwritten.
*
*  INCX    (input) INTEGER
*          The increment between successive elements of X. INCX > 0.
*
*  Y       (input/output) COMPLEX array, dimension (1+(N-1)*INCY)
*          On entry, Y contains the N-vector Y.
*          On exit, Y is overwritten.
*
*  INCY    (input) INTEGER
*          The increment between successive elements of Y. INCY > 0.
*
*  SSMIN   (output) REAL
*          The smallest singular value of the N-by-2 matrix A = ( X Y ).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      REAL               SSMAX
      COMPLEX            A11, A12, A22, C, TAU
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG
*     ..
*     .. External Functions ..
      COMPLEX            CDOTC
      EXTERNAL           CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CLARFG, SLAS2
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 ) THEN
         SSMIN = ZERO
         RETURN
      END IF
*
*     Compute the QR factorization of the N-by-2 matrix ( X Y )
*
      CALL CLARFG( N, X( 1 ), X( 1+INCX ), INCX, TAU )
      A11 = X( 1 )
      X( 1 ) = CONE
*
      C = -CONJG( TAU )*CDOTC( N, X, INCX, Y, INCY )
      CALL CAXPY( N, C, X, INCX, Y, INCY )
*
      CALL CLARFG( N-1, Y( 1+INCY ), Y( 1+2*INCY ), INCY, TAU )
*
      A12 = Y( 1 )
      A22 = Y( 1+INCY )
*
*     Compute the SVD of 2-by-2 Upper triangular matrix.
*
      CALL SLAS2( ABS( A11 ), ABS( A12 ), ABS( A22 ), SSMIN, SSMAX )
*
      RETURN
*
*     End of CLAPLL
*
      END
