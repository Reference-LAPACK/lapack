      SUBROUTINE ZPOT01( UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDAFAC, N
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
*  ZPOT01 reconstructs a Hermitian positive definite matrix  A  from
*  its L*L' or U'*U factorization and computes the residual
*     norm( L*L' - A ) / ( N * norm(A) * EPS ) or
*     norm( U'*U - A ) / ( N * norm(A) * EPS ),
*  where EPS is the machine epsilon, L' is the conjugate transpose of L,
*  and U' is the conjugate transpose of U.
*
*  Arguments
*  ==========
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
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The original Hermitian matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N)
*
*  AFAC    (input/output) COMPLEX*16 array, dimension (LDAFAC,N)
*          On entry, the factor L or U from the L*L' or U'*U
*          factorization of A.
*          Overwritten with the reconstructed matrix, and then with the
*          difference L*L' - A (or U'*U - A).
*
*  LDAFAC  (input) INTEGER
*          The leading dimension of the array AFAC.  LDAFAC >= max(1,N).
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RESID   (output) DOUBLE PRECISION
*          If UPLO = 'L', norm(L*L' - A) / ( N * norm(A) * EPS )
*          If UPLO = 'U', norm(U'*U - A) / ( N * norm(A) * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K
      DOUBLE PRECISION   ANORM, EPS, TR
      COMPLEX*16         TC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHE
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, DLAMCH, ZLANHE, ZDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZHER, ZSCAL, ZTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DIMAG
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
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Check the imaginary parts of the diagonal elements and return with
*     an error code if any are nonzero.
*
      DO 10 J = 1, N
         IF( DIMAG( AFAC( J, J ) ).NE.ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
   10 CONTINUE
*
*     Compute the product U'*U, overwriting U.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 K = N, 1, -1
*
*           Compute the (K,K) element of the result.
*
            TR = ZDOTC( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 )
            AFAC( K, K ) = TR
*
*           Compute the rest of column K.
*
            CALL ZTRMV( 'Upper', 'Conjugate', 'Non-unit', K-1, AFAC,
     $                  LDAFAC, AFAC( 1, K ), 1 )
*
   20    CONTINUE
*
*     Compute the product L*L', overwriting L.
*
      ELSE
         DO 30 K = N, 1, -1
*
*           Add a multiple of column K of the factor L to each of
*           columns K+1 through N.
*
            IF( K+1.LE.N )
     $         CALL ZHER( 'Lower', N-K, ONE, AFAC( K+1, K ), 1,
     $                    AFAC( K+1, K+1 ), LDAFAC )
*
*           Scale column K by the diagonal element.
*
            TC = AFAC( K, K )
            CALL ZSCAL( N-K+1, TC, AFAC( K, K ), 1 )
*
   30    CONTINUE
      END IF
*
*     Compute the difference  L*L' - A (or U'*U - A).
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 50 J = 1, N
            DO 40 I = 1, J - 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   40       CONTINUE
            AFAC( J, J ) = AFAC( J, J ) - DBLE( A( J, J ) )
   50    CONTINUE
      ELSE
         DO 70 J = 1, N
            AFAC( J, J ) = AFAC( J, J ) - DBLE( A( J, J ) )
            DO 60 I = J + 1, N
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   60       CONTINUE
   70    CONTINUE
      END IF
*
*     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
*
      RESID = ZLANHE( '1', UPLO, N, AFAC, LDAFAC, RWORK )
*
      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
*
      RETURN
*
*     End of ZPOT01
*
      END
