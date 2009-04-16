      SUBROUTINE DLAT2S( UPLO, N, A, LDA, SA, LDSA, INFO )
*
*  -- LAPACK PROTOTYPE auxiliary routine (version 3.1.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     May 2007
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDSA, N
*     ..
*     .. Array Arguments ..
      REAL               SA( LDSA, * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLAT2S converts a DOUBLE PRECISION triangular matrix, SA, to a SINGLE
*  PRECISION triangular matrix, A.
*
*  RMAX is the overflow for the SINGLE PRECISION arithmetic
*  DLAS2S checks that all the entries of A are between -RMAX and
*  RMAX. If not the convertion is aborted and a flag is raised.
*
*  This is an auxiliary routine so there is no argument checking.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  N       (input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N triangular coefficient matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  SA      (output) REAL array, dimension (LDSA,N)
*          Only the UPLO part of SA is referenced.  On exit, if INFO=0,
*          the N-by-N coefficient matrix SA; if INFO>0, the content of
*          the UPLO part of SA is unspecified.
*
*  LDSA    (input) INTEGER
*          The leading dimension of the array SA.  LDSA >= max(1,M).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          = 1:  an entry of the matrix A is greater than the SINGLE
*                PRECISION overflow threshold, in this case, the content
*                of the UPLO part of SA in exit is unspecified.
*
*  =========
*
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   RMAX
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      LOGICAL            LSAME
      EXTERNAL           SLAMCH, LSAME
*     ..
*     .. Executable Statements ..
*
      RMAX = SLAMCH( 'O' )
      UPPER = LSAME( UPLO, 'U' )
      IF( UPPER ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, J
               IF( ( A( I, J ).LT.-RMAX ) .OR. ( A( I, J ).GT.RMAX ) )
     +             THEN
                  INFO = 1
                  GO TO 50
               END IF
               SA( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE
         DO 40 J = 1, N
            DO 30 I = J, N
               IF( ( A( I, J ).LT.-RMAX ) .OR. ( A( I, J ).GT.RMAX ) )
     +             THEN
                  INFO = 1
                  GO TO 50
               END IF
               SA( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      END IF
   50 CONTINUE
*
      RETURN
*
*     End of DLAT2S
*
      END
