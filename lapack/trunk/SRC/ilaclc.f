      INTEGER FUNCTION ILACLC(M, N, A, LDA)
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     December 2007
!
!     .. Scalar Arguments ..
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ILACLC scans A for its last non-zero column.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input) COMPLEX array, dimension (LDA,N)
!          The m by n matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX          ZERO
      PARAMETER ( ZERO = (0.0E+0, 0.0E+0) )
!     ..
!     .. Local Scalars ..
      INTEGER I, J
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 .OR. A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILACLC = N
      ELSE
!     Now scan each column from the end, returning with the first non-zero.
         DO ILACLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILACLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END FUNCTION
