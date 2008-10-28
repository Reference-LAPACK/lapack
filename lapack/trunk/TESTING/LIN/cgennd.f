      LOGICAL FUNCTION CGENND (M, N, A, LDA)
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
      COMPLEX A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*     CGENND tests that its argument has a real, non-negative diagonal.
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
*  A       (input) COMPLEX array, dimension (LDA, N)
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
      COMPLEX AII
*     ..
*     .. Intrinsics ..
      INTRINSIC MIN, REAL, AIMAG
*     ..
*     .. Executable Statements ..
      K = MIN( M, N )
      DO I = 1, K
         AII = A( I, I )
         IF( REAL( AII ).LT.ZERO.OR.AIMAG( AII ).NE.ZERO ) THEN
            CGENND = .FALSE.
            RETURN
         END IF
      END DO
      CGENND = .TRUE.
      RETURN
      END
