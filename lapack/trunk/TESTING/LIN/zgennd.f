      LOGICAL FUNCTION ZGENND (M, N, A, LDA)
      IMPLICIT NONE
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     February 2008
*
*     .. Scalar Arguments ..
      INTEGER M, N, LDA
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*     ZGENND tests that its argument has a real, non-negative diagonal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows in A.
*
*  N       (input) INTEGER
*          The number of columns in A.
*
*  A       (input) COMPLEX*16 array, dimension (LDA, N)
*          The matrix.
*
*  LDA     (input) INTEGER
*          Leading dimension of A.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL OUT
      INTEGER I, K
      COMPLEX*16 AII
*     ..
*     .. Intrinsics ..
      INTRINSIC MIN, DBLE, DIMAG
*     ..
*     .. Executable Statements ..
      K = MIN( M, N )
      DO I = 1, K
         AII = A( I, I )
         IF( DBLE( AII ).LT.ZERO.OR.DIMAG( AII ).NE.ZERO ) THEN
            ZGENND = .FALSE.
            RETURN
         END IF
      END DO
      ZGENND = .TRUE.
      RETURN
      END
