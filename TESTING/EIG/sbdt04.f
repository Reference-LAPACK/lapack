*> \brief \b SBDT04
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SBDT04( UPLO, N, D, E, S, NS, U, LDU, VT, LDVT,
*                          WORK, RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDU, LDVT, N, NS
*       REAL               RESID
*       ..
*       .. Array Arguments ..
*       REAL               D( * ), E( * ), S( * ), U( LDU, * ),
*      $                   VT( LDVT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SBDT04 reconstructs a bidiagonal matrix B from its (partial) SVD:
*>    S = U' * B * V
*> where U and V are orthogonal matrices and S is diagonal.
*>
*> The test ratio to test the singular value decomposition is
*>    RESID = norm( S - U' * B * V ) / ( n * norm(B) * EPS )
*> where VT = V' and EPS is the machine precision.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the matrix B is upper or lower bidiagonal.
*>          = 'U':  Upper bidiagonal
*>          = 'L':  Lower bidiagonal
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix B.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is REAL array, dimension (N)
*>          The n diagonal elements of the bidiagonal matrix B.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is REAL array, dimension (N-1)
*>          The (n-1) superdiagonal elements of the bidiagonal matrix B
*>          if UPLO = 'U', or the (n-1) subdiagonal elements of B if
*>          UPLO = 'L'.
*> \endverbatim
*>
*> \param[in] S
*> \verbatim
*>          S is REAL array, dimension (NS)
*>          The singular values from the (partial) SVD of B, sorted in
*>          decreasing order.
*> \endverbatim
*>
*> \param[in] NS
*> \verbatim
*>          NS is INTEGER
*>          The number of singular values/vectors from the (partial)
*>          SVD of B.
*> \endverbatim
*>
*> \param[in] U
*> \verbatim
*>          U is REAL array, dimension (LDU,NS)
*>          The n by ns orthogonal matrix U in S = U' * B * V.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of the array U.  LDU >= max(1,N)
*> \endverbatim
*>
*> \param[in] VT
*> \verbatim
*>          VT is REAL array, dimension (LDVT,N)
*>          The n by ns orthogonal matrix V in S = U' * B * V.
*> \endverbatim
*>
*> \param[in] LDVT
*> \verbatim
*>          LDVT is INTEGER
*>          The leading dimension of the array VT.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is REAL
*>          The test ratio:  norm(S - U' * B * V) / ( n * norm(B) * EPS )
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
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE SBDT04( UPLO, N, D, E, S, NS, U, LDU, VT, LDVT, WORK,
     $                   RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDU, LDVT, N, NS
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), S( * ), U( LDU, * ),
     $                   VT( LDVT, * ), WORK( * )
*     ..
*
* ======================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K
      REAL               BNORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ISAMAX
      REAL               SASUM, SLAMCH
      EXTERNAL           LSAME, ISAMAX, SASUM, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      RESID = ZERO
      IF( N.LE.0 .OR. NS.LE.0 )
     $   RETURN
*
      EPS = SLAMCH( 'Precision' )
*
*     Compute S - U' * B * V.
*
      BNORM = ZERO
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        B is upper bidiagonal.
*
         K = 0
         DO 20 I = 1, NS
            DO 10 J = 1, N-1
               K = K + 1
               WORK( K ) = D( J )*VT( I, J ) + E( J )*VT( I, J+1 )
   10       CONTINUE
            K = K + 1
            WORK( K ) = D( N )*VT( I, N )
   20    CONTINUE
         BNORM = ABS( D( 1 ) )
         DO 30 I = 2, N
            BNORM = MAX( BNORM, ABS( D( I ) )+ABS( E( I-1 ) ) )
   30    CONTINUE
      ELSE
*
*        B is lower bidiagonal.
*
         K = 0
         DO 50 I = 1, NS
            K = K + 1
            WORK( K ) = D( 1 )*VT( I, 1 )
            DO 40 J = 1, N-1
               K = K + 1
               WORK( K ) = E( J )*VT( I, J ) + D( J+1 )*VT( I, J+1 )
   40       CONTINUE
   50    CONTINUE
         BNORM = ABS( D( N ) )
         DO 60 I = 1, N-1
            BNORM = MAX( BNORM, ABS( D( I ) )+ABS( E( I ) ) )
   60    CONTINUE
      END IF
*
      CALL SGEMM( 'T', 'N', NS, NS, N, -ONE, U, LDU, WORK( 1 ),
     $            N, ZERO, WORK( 1+N*NS ), NS )
*
*     norm(S - U' * B * V)
*
      K = N*NS
      DO 70 I = 1, NS
         WORK( K+I ) =  WORK( K+I ) + S( I )
         RESID = MAX( RESID, SASUM( NS, WORK( K+1 ), 1 ) )
         K = K + NS
   70 CONTINUE
*
      IF( BNORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         IF( BNORM.GE.RESID ) THEN
            RESID = ( RESID / BNORM ) / ( REAL( N )*EPS )
         ELSE
            IF( BNORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REAL( N )*BNORM ) / BNORM ) /
     $                 ( REAL( N )*EPS )
            ELSE
               RESID = MIN( RESID / BNORM, REAL( N ) ) /
     $                 ( REAL( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of SBDT04
*
      END
