      SUBROUTINE ZPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK,
     $                   RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            KD, LDA, LDAFAC, N
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * )
*     ..
*
*  Purpose
*  =======
*
*  ZPBT01 reconstructs a Hermitian positive definite band matrix A from
*  its L*L' or U'*U factorization and computes the residual
*     norm( L*L' - A ) / ( N * norm(A) * EPS ) or
*     norm( U'*U - A ) / ( N * norm(A) * EPS ),
*  where EPS is the machine epsilon, L' is the conjugate transpose of
*  L, and U' is the conjugate transpose of U.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of super-diagonals of the matrix A if UPLO = 'U',
*          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The original Hermitian band matrix A.  If UPLO = 'U', the
*          upper triangular part of A is stored as a band matrix; if
*          UPLO = 'L', the lower triangular part of A is stored.  The
*          columns of the appropriate triangle are stored in the columns
*          of A and the diagonals of the triangle are stored in the rows
*          of A.  See ZPBTRF for further details.
*
*  LDA     (input) INTEGER.
*          The leading dimension of the array A.  LDA >= max(1,KD+1).
*
*  AFAC    (input) COMPLEX*16 array, dimension (LDAFAC,N)
*          The factored form of the matrix A.  AFAC contains the factor
*          L or U from the L*L' or U'*U factorization in band storage
*          format, as computed by ZPBTRF.
*
*  LDAFAC  (input) INTEGER
*          The leading dimension of the array AFAC.
*          LDAFAC >= max(1,KD+1).
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RESID   (output) DOUBLE PRECISION
*          If UPLO = 'L', norm(L*L' - A) / ( N * norm(A) * EPS )
*          If UPLO = 'U', norm(U'*U - A) / ( N * norm(A) * EPS )
*
*  =====================================================================
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K, KC, KLEN, ML, MU
      DOUBLE PRECISION   AKK, ANORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHB
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, DLAMCH, ZLANHB, ZDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZDSCAL, ZHER, ZTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DIMAG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHB( '1', UPLO, N, KD, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Check the imaginary parts of the diagonal elements and return with
*     an error code if any are nonzero.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 10 J = 1, N
            IF( DIMAG( AFAC( KD+1, J ) ).NE.ZERO ) THEN
               RESID = ONE / EPS
               RETURN
            END IF
   10    CONTINUE
      ELSE
         DO 20 J = 1, N
            IF( DIMAG( AFAC( 1, J ) ).NE.ZERO ) THEN
               RESID = ONE / EPS
               RETURN
            END IF
   20    CONTINUE
      END IF
*
*     Compute the product U'*U, overwriting U.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 30 K = N, 1, -1
            KC = MAX( 1, KD+2-K )
            KLEN = KD + 1 - KC
*
*           Compute the (K,K) element of the result.
*
            AKK = ZDOTC( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 )
            AFAC( KD+1, K ) = AKK
*
*           Compute the rest of column K.
*
            IF( KLEN.GT.0 )
     $         CALL ZTRMV( 'Upper', 'Conjugate', 'Non-unit', KLEN,
     $                     AFAC( KD+1, K-KLEN ), LDAFAC-1,
     $                     AFAC( KC, K ), 1 )
*
   30    CONTINUE
*
*     UPLO = 'L':  Compute the product L*L', overwriting L.
*
      ELSE
         DO 40 K = N, 1, -1
            KLEN = MIN( KD, N-K )
*
*           Add a multiple of column K of the factor L to each of
*           columns K+1 through N.
*
            IF( KLEN.GT.0 )
     $         CALL ZHER( 'Lower', KLEN, ONE, AFAC( 2, K ), 1,
     $                    AFAC( 1, K+1 ), LDAFAC-1 )
*
*           Scale column K by the diagonal element.
*
            AKK = AFAC( 1, K )
            CALL ZDSCAL( KLEN+1, AKK, AFAC( 1, K ), 1 )
*
   40    CONTINUE
      END IF
*
*     Compute the difference  L*L' - A  or  U'*U - A.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 60 J = 1, N
            MU = MAX( 1, KD+2-J )
            DO 50 I = MU, KD + 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   50       CONTINUE
   60    CONTINUE
      ELSE
         DO 80 J = 1, N
            ML = MIN( KD+1, N-J+1 )
            DO 70 I = 1, ML
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   70       CONTINUE
   80    CONTINUE
      END IF
*
*     Compute norm( L*L' - A ) / ( N * norm(A) * EPS )
*
      RESID = ZLANHB( '1', UPLO, N, KD, AFAC, LDAFAC, RWORK )
*
      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
*
      RETURN
*
*     End of ZPBT01
*
      END
