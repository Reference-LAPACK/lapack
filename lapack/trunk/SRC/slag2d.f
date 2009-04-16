      SUBROUTINE SLAG2D( M, N, SA, LDSA, A, LDA, INFO )
*
*  -- LAPACK PROTOTYPE auxiliary routine (version 3.1.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     August 2007
*
*     ..
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDSA, M, N
*     ..
*     .. Array Arguments ..
      REAL               SA( LDSA, * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  SLAG2D converts a SINGLE PRECISION matrix, SA, to a DOUBLE
*  PRECISION matrix, A.
*
*  Note that while it is possible to overflow while converting
*  from double to single, it is not possible to overflow when
*  converting from single to double.
*
*  This is an auxiliary routine so there is no argument checking.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of lines of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  SA      (input) REAL array, dimension (LDSA,N)
*          On entry, the M-by-N coefficient matrix SA.
*
*  LDSA    (input) INTEGER
*          The leading dimension of the array SA.  LDSA >= max(1,M).
*
*  A       (output) DOUBLE PRECISION array, dimension (LDA,N)
*          On exit, the M-by-N coefficient matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*  =========
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      DO 20 J = 1, N
         DO 10 I = 1, M
            A( I, J ) = SA( I, J )
   10    CONTINUE
   20 CONTINUE
      RETURN
*
*     End of SLAG2D
*
      END
