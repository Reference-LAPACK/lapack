      SUBROUTINE CLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INCV, LDC, N
      COMPLEX            TAU
*     ..
*     .. Array Arguments ..
      COMPLEX            C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CLARFY applies an elementary reflector, or Householder matrix, H,
*  to an n x n Hermitian matrix C, from both the left and the right.
*
*  H is represented in the form
*
*     H = I - tau * v * v'
*
*  where  tau  is a scalar and  v  is a vector.
*
*  If  tau  is  zero, then  H  is taken to be the unit matrix.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix C is stored.
*          = 'U':  Upper triangle
*          = 'L':  Lower triangle
*
*  N       (input) INTEGER
*          The number of rows and columns of the matrix C.  N >= 0.
*
*  V       (input) COMPLEX array, dimension
*                  (1 + (N-1)*abs(INCV))
*          The vector v as described above.
*
*  INCV    (input) INTEGER
*          The increment between successive elements of v.  INCV must
*          not be zero.
*
*  TAU     (input) COMPLEX
*          The value tau as described above.
*
*  C       (input/output) COMPLEX array, dimension (LDC, N)
*          On entry, the matrix C.
*          On exit, C is overwritten by H * C * H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C.  LDC >= max( 1, N ).
*
*  WORK    (workspace) COMPLEX array, dimension (N)
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ),
     $                   HALF = ( 0.5E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX            ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CHEMV, CHER2
*     ..
*     .. External Functions ..
      COMPLEX            CDOTC
      EXTERNAL           CDOTC
*     ..
*     .. Executable Statements ..
*
      IF( TAU.EQ.ZERO )
     $   RETURN
*
*     Form  w:= C * v
*
      CALL CHEMV( UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )
*
      ALPHA = -HALF*TAU*CDOTC( N, WORK, 1, V, INCV )
      CALL CAXPY( N, ALPHA, V, INCV, WORK, 1 )
*
*     C := C - v * w' - w * v'
*
      CALL CHER2( UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC )
*
      RETURN
*
*     End of CLARFY
*
      END
