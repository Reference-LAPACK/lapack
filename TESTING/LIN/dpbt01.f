*> \brief \b DPBT01
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK,
*                          RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            KD, LDA, LDAFAC, N
*       DOUBLE PRECISION   RESID
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DPBT01 reconstructs a symmetric positive definite band matrix A from
*> its L*L' or U'*U factorization and computes the residual
*>    norm( L*L' - A ) / ( N * norm(A) * EPS ) or
*>    norm( U'*U - A ) / ( N * norm(A) * EPS ),
*> where EPS is the machine epsilon, L' is the conjugate transpose of
*> L, and U' is the conjugate transpose of U.
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
*> \param[in] KD
*> \verbatim
*>          KD is INTEGER
*>          The number of super-diagonals of the matrix A if UPLO = 'U',
*>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The original symmetric band matrix A.  If UPLO = 'U', the
*>          upper triangular part of A is stored as a band matrix; if
*>          UPLO = 'L', the lower triangular part of A is stored.  The
*>          columns of the appropriate triangle are stored in the columns
*>          of A and the diagonals of the triangle are stored in the rows
*>          of A.  See DPBTRF for further details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER.
*>          The leading dimension of the array A.  LDA >= max(1,KD+1).
*> \endverbatim
*>
*> \param[in] AFAC
*> \verbatim
*>          AFAC is DOUBLE PRECISION array, dimension (LDAFAC,N)
*>          The factored form of the matrix A.  AFAC contains the factor
*>          L or U from the L*L' or U'*U factorization in band storage
*>          format, as computed by DPBTRF.
*> \endverbatim
*>
*> \param[in] LDAFAC
*> \verbatim
*>          LDAFAC is INTEGER
*>          The leading dimension of the array AFAC.
*>          LDAFAC >= max(1,KD+1).
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is DOUBLE PRECISION
*>          If UPLO = 'L', norm(L*L' - A) / ( N * norm(A) * EPS )
*>          If UPLO = 'U', norm(U'*U - A) / ( N * norm(A) * EPS )
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
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE DPBT01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK,
     $                   RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            KD, LDA, LDAFAC, N
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
*     ..
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
      DOUBLE PRECISION   ANORM, EPS, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLANSB
      EXTERNAL           LSAME, DDOT, DLAMCH, DLANSB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSYR, DTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
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
      ANORM = DLANSB( '1', UPLO, N, KD, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute the product U'*U, overwriting U.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 10 K = N, 1, -1
            KC = MAX( 1, KD+2-K )
            KLEN = KD + 1 - KC
*
*           Compute the (K,K) element of the result.
*
            T = DDOT( KLEN+1, AFAC( KC, K ), 1, AFAC( KC, K ), 1 )
            AFAC( KD+1, K ) = T
*
*           Compute the rest of column K.
*
            IF( KLEN.GT.0 )
     $         CALL DTRMV( 'Upper', 'Transpose', 'Non-unit', KLEN,
     $                     AFAC( KD+1, K-KLEN ), LDAFAC-1,
     $                     AFAC( KC, K ), 1 )
*
   10    CONTINUE
*
*     UPLO = 'L':  Compute the product L*L', overwriting L.
*
      ELSE
         DO 20 K = N, 1, -1
            KLEN = MIN( KD, N-K )
*
*           Add a multiple of column K of the factor L to each of
*           columns K+1 through N.
*
            IF( KLEN.GT.0 )
     $         CALL DSYR( 'Lower', KLEN, ONE, AFAC( 2, K ), 1,
     $                    AFAC( 1, K+1 ), LDAFAC-1 )
*
*           Scale column K by the diagonal element.
*
            T = AFAC( 1, K )
            CALL DSCAL( KLEN+1, T, AFAC( 1, K ), 1 )
*
   20    CONTINUE
      END IF
*
*     Compute the difference  L*L' - A  or  U'*U - A.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 40 J = 1, N
            MU = MAX( 1, KD+2-J )
            DO 30 I = MU, KD + 1
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            ML = MIN( KD+1, N-J+1 )
            DO 50 I = 1, ML
               AFAC( I, J ) = AFAC( I, J ) - A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
*     Compute norm( L*L' - A ) / ( N * norm(A) * EPS )
*
      RESID = DLANSB( 'I', UPLO, N, KD, AFAC, LDAFAC, RWORK )
*
      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
*
      RETURN
*
*     End of DPBT01
*
      END
