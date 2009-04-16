      SUBROUTINE DLAG2S( M, N, A, LDA, SA, LDSA, INFO )
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
*  DLAG2S converts a DOUBLE PRECISION matrix, SA, to a SINGLE
*  PRECISION matrix, A.
*
*  RMAX is the overflow for the SINGLE PRECISION arithmetic
*  DLAG2S checks that all the entries of A are between -RMAX and
*  RMAX. If not the convertion is aborted and a flag is raised.
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
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N coefficient matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  SA      (output) REAL array, dimension (LDSA,N)
*          On exit, if INFO=0, the M-by-N coefficient matrix SA; if
*          INFO>0, the content of SA is unspecified.
*
*  LDSA    (input) INTEGER
*          The leading dimension of the array SA.  LDSA >= max(1,M).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          = 1:  an entry of the matrix A is greater than the SINGLE
*                PRECISION overflow threshold, in this case, the content
*                of SA in exit is unspecified.
*
*  =========
*
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   RMAX
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Executable Statements ..
*
      RMAX = SLAMCH( 'O' )
      DO 20 J = 1, N
         DO 10 I = 1, M
            IF( ( A( I, J ).LT.-RMAX ) .OR. ( A( I, J ).GT.RMAX ) ) THEN
               INFO = 1
               GO TO 30
            END IF
            SA( I, J ) = A( I, J )
   10    CONTINUE
   20 CONTINUE
      INFO = 0
   30 CONTINUE
      RETURN
*
*     End of DLAG2S
*
      END
