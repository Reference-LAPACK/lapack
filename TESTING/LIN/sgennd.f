      LOGICAL FUNCTION SGENND (M, N, A, LDA)
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
      REAL A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*     SGENND tests that its argument has a non-negative diagonal.
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
*  A       (input) REAL array, dimension (LDA, N)
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
      INTEGER I, K
*     ..
*     .. Intrinsics ..
      INTRINSIC MIN
*     ..
*     .. Executable Statements ..
      K = MIN( M, N )
      DO I = 1, K
         IF( A( I, I ).LT.ZERO ) THEN
            SGENND = .FALSE.
            RETURN
         END IF
      END DO
      SGENND = .TRUE.
      RETURN
      END
