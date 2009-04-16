      SUBROUTINE ZLAG2C( M, N, A, LDA, SA, LDSA, INFO )
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
      COMPLEX            SA( LDSA, * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLAG2C converts a COMPLEX*16 matrix, SA, to a COMPLEX matrix, A.
*
*  RMAX is the overflow for the SINGLE PRECISION arithmetic
*  ZLAG2C checks that all the entries of A are between -RMAX and
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
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N coefficient matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  SA      (output) COMPLEX array, dimension (LDSA,N)
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
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DIMAG
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
            IF( ( DBLE( A( I, J ) ).LT.-RMAX ) .OR.
     +          ( DBLE( A( I, J ) ).GT.RMAX ) .OR.
     +          ( DIMAG( A( I, J ) ).LT.-RMAX ) .OR.
     +          ( DIMAG( A( I, J ) ).GT.RMAX ) ) THEN
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
*     End of ZLAG2C
*
      END
