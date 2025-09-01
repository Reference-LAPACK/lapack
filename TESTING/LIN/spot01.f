*> \brief \b SPOT01
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SPOT01( UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDA, LDAFAC, N
*       REAL               RESID
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SPOT01 reconstructs a symmetric positive definite matrix  A  from
*> its L*L' or U'*U factorization and computes the residual
*>    norm( L*L' - A ) / ( N * norm(A) * EPS ) or
*>    norm( U'*U - A ) / ( N * norm(A) * EPS ),
*> where EPS is the machine epsilon.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows and columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The original symmetric matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N)
*> \endverbatim
*>
*> \param[in,out] AFAC
*> \verbatim
*>          AFAC is REAL array, dimension (LDAFAC,N)
*>          On entry, the factor L or U from the L * L**T or U**T * U
*>          factorization of A.
*>          Overwritten with the reconstructed matrix, and then with
*>          the difference L * L**T - A (or U**T * U - A).
*> \endverbatim
*>
*> \param[in] LDAFAC
*> \verbatim
*>          LDAFAC is INTEGER
*>          The leading dimension of the array AFAC.  LDAFAC >= max(1,N).
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (N)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is REAL
*>          If UPLO = 'L', norm(L * L**T - A) / ( N * norm(A) * EPS )
*>          If UPLO = 'U', norm(U**T * U - A) / ( N * norm(A) * EPS )
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup single_lin
*
*  =====================================================================
      SUBROUTINE SPOT01( UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDAFAC, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K
      REAL               ANORM, EPS, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT, SLAMCH, SLANSY
      EXTERNAL           LSAME, SDOT, SLAMCH, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SSYR, STRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
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
      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute the product U**T * U, overwriting U.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 10 K = N, 1, -1
*
*           Compute the (K,K) element of the result.
*
            T = SDOT( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 )
            AFAC( K, K ) = T
*
*           Compute the rest of column K.
*
            CALL STRMV( 'Upper', 'Transpose', 'Non-unit', K-1, AFAC,
     $                  LDAFAC, AFAC( 1, K ), 1 )
*
   10    CONTINUE
*
*     Compute the product L * L**T, overwriting L.
*
      ELSE
         DO 20 K = N, 1, -1
*
*           Add a multiple of column K of the factor L to each of
*           columns K+1 through N.
*
            IF( K+1.LE.N )
     $         CALL SSYR( 'Lower', N-K, ONE, AFAC( K+1, K ), 1,
     $                    AFAC( K+1, K+1 ), LDAFAC )
*
*           Scale column K by the diagonal element.
*
            T = AFAC( K, K )
            CALL SSCAL( N-K+1, T, AFAC( K, K ), 1 )
*
   20    CONTINUE
      END IF
*
*     Compute the difference L * L**T - A (or U**T * U - A).
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = 1, J
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = J, N
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
*     Compute norm(L*U - A) / ( N * norm(A) * EPS )
*
      RESID = SLANSY( '1', UPLO, N, AFAC, LDAFAC, RWORK )
*
      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
*
      RETURN
*
*     End of SPOT01
*
      END
