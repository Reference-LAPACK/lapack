      INTEGER FUNCTION ILAZLC(M, N, A, LDA)
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2.1)                        --
*
*  -- April 2009                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ILAZLC scans A for its last non-zero column.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16       ZERO
      PARAMETER ( ZERO = (0.0D+0, 0.0D+0) )
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILAZLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILAZLC = N
      ELSE
*     Now scan each column from the end, returning with the first non-zero.
         DO ILAZLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILAZLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END FUNCTION
