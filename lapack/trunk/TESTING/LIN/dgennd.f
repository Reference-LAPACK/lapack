      LOGICAL FUNCTION DGENND (M, N, A, LDA)
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
      DOUBLE PRECISION A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*     DGENND tests that its argument has a non-negative diagonal.
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
*  A       (input) DOUBLE PRECISION array, dimension (LDA, N)
*          The matrix.
*
*  LDA     (input) INTEGER
*          Leading dimension of A.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
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
            DGENND = .FALSE.
            RETURN
         END IF
      END DO
      DGENND = .TRUE.
      RETURN
      END
