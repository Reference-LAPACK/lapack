      SUBROUTINE CQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, M, N, SCALE
      REAL               NORMA
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  CQRT13 generates a full-rank matrix that may be scaled to have large
*  or small norm.
*
*  Arguments
*  =========
*
*  SCALE   (input) INTEGER
*          SCALE = 1: normally scaled matrix
*          SCALE = 2: matrix scaled up
*          SCALE = 3: matrix scaled down
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  N       (input) INTEGER
*          The number of columns of A.
*
*  A       (output) COMPLEX array, dimension (LDA,N)
*          The M-by-N matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  NORMA   (output) REAL
*          The one-norm of A.
*
*  ISEED   (input/output) integer array, dimension (4)
*          Seed for random number generator
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, J
      REAL               BIGNUM, SMLNUM
*     ..
*     .. External Functions ..
      REAL               CLANGE, SCASUM, SLAMCH
      EXTERNAL           CLANGE, SCASUM, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARNV, CLASCL, SLABAD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, REAL, SIGN
*     ..
*     .. Local Arrays ..
      REAL               DUMMY( 1 )
*     ..
*     .. Executable Statements ..
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*     benign matrix
*
      DO 10 J = 1, N
         CALL CLARNV( 2, ISEED, M, A( 1, J ) )
         IF( J.LE.M ) THEN
            A( J, J ) = A( J, J ) + CMPLX( SIGN( SCASUM( M, A( 1, J ),
     $                  1 ), REAL( A( J, J ) ) ) )
         END IF
   10 CONTINUE
*
*     scaled versions
*
      IF( SCALE.NE.1 ) THEN
         NORMA = CLANGE( 'Max', M, N, A, LDA, DUMMY )
         SMLNUM = SLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
         CALL SLABAD( SMLNUM, BIGNUM )
         SMLNUM = SMLNUM / SLAMCH( 'Epsilon' )
         BIGNUM = ONE / SMLNUM
*
         IF( SCALE.EQ.2 ) THEN
*
*           matrix scaled up
*
            CALL CLASCL( 'General', 0, 0, NORMA, BIGNUM, M, N, A, LDA,
     $                   INFO )
         ELSE IF( SCALE.EQ.3 ) THEN
*
*           matrix scaled down
*
            CALL CLASCL( 'General', 0, 0, NORMA, SMLNUM, M, N, A, LDA,
     $                   INFO )
         END IF
      END IF
*
      NORMA = CLANGE( 'One-norm', M, N, A, LDA, DUMMY )
      RETURN
*
*     End of CQRT13
*
      END
