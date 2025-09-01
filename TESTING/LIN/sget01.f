*> \brief \b SGET01
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGET01( M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK,
*                          RESID )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDAFAC, M, N
*       REAL               RESID
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL               A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGET01 reconstructs a matrix A from its L*U factorization and
*> computes the residual
*>    norm(L*U - A) / ( N * norm(A) * EPS ),
*> where EPS is the machine epsilon.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The original M x N matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] AFAC
*> \verbatim
*>          AFAC is REAL array, dimension (LDAFAC,N)
*>          The factored form of the matrix A.  AFAC contains the factors
*>          L and U from the L*U factorization as computed by SGETRF.
*>          Overwritten with the reconstructed matrix, and then with the
*>          difference L*U - A.
*> \endverbatim
*>
*> \param[in] LDAFAC
*> \verbatim
*>          LDAFAC is INTEGER
*>          The leading dimension of the array AFAC.  LDAFAC >= max(1,M).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices from SGETRF.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (M)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is REAL
*>          norm(L*U - A) / ( N * norm(A) * EPS )
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
      SUBROUTINE SGET01( M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK,
     $                   RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDAFAC, M, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
*     ..
*
*  =====================================================================
*
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
      REAL               SDOT, SLAMCH, SLANGE
      EXTERNAL           SDOT, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV, SLASWP, SSCAL, STRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick exit if M = 0 or N = 0.
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Determine EPS and the norm of A.
*
      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANGE( '1', M, N, A, LDA, RWORK )
*
*     Compute the product L*U and overwrite AFAC with the result.
*     A column at a time of the product is obtained, starting with
*     column N.
*
      DO 10 K = N, 1, -1
         IF( K.GT.M ) THEN
            CALL STRMV( 'Lower', 'No transpose', 'Unit', M, AFAC,
     $                  LDAFAC, AFAC( 1, K ), 1 )
         ELSE
*
*           Compute elements (K+1:M,K)
*
            T = AFAC( K, K )
            IF( K+1.LE.M ) THEN
               CALL SSCAL( M-K, T, AFAC( K+1, K ), 1 )
               CALL SGEMV( 'No transpose', M-K, K-1, ONE,
     $                     AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, ONE,
     $                     AFAC( K+1, K ), 1 )
            END IF
*
*           Compute the (K,K) element
*
            AFAC( K, K ) = T + SDOT( K-1, AFAC( K, 1 ), LDAFAC,
     $                     AFAC( 1, K ), 1 )
*
*           Compute elements (1:K-1,K)
*
            CALL STRMV( 'Lower', 'No transpose', 'Unit', K-1, AFAC,
     $                  LDAFAC, AFAC( 1, K ), 1 )
         END IF
   10 CONTINUE
      CALL SLASWP( N, AFAC, LDAFAC, 1, MIN( M, N ), IPIV, -1 )
*
*     Compute the difference  L*U - A  and store in AFAC.
*
      DO 30 J = 1, N
         DO 20 I = 1, M
            AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   20    CONTINUE
   30 CONTINUE
*
*     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
*
      RESID = SLANGE( '1', M, N, AFAC, LDAFAC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      END IF
*
      RETURN
*
*     End of SGET01
*
      END
