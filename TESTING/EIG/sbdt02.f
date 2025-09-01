*> \brief \b SBDT02
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RESID )
*
*       .. Scalar Arguments ..
*       INTEGER            LDB, LDC, LDU, M, N
*       REAL               RESID
*       ..
*       .. Array Arguments ..
*       REAL               B( LDB, * ), C( LDC, * ), U( LDU, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SBDT02 tests the change of basis C = U**H * B by computing the
*> residual
*>
*>    RESID = norm(B - U * C) / ( max(m,n) * norm(B) * EPS ),
*>
*> where B and C are M by N matrices, U is an M by M orthogonal matrix,
*> and EPS is the machine precision.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrices B and C and the order of
*>          the matrix Q.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices B and C.
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is REAL array, dimension (LDB,N)
*>          The m by n matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,M).
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is REAL array, dimension (LDC,N)
*>          The m by n matrix C, assumed to contain U**H * B.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C.  LDC >= max(1,M).
*> \endverbatim
*>
*> \param[in] U
*> \verbatim
*>          U is REAL array, dimension (LDU,M)
*>          The m by m orthogonal matrix U.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of the array U.  LDU >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (M)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is REAL
*>          RESID = norm(B - U * C) / ( max(m,n) * norm(B) * EPS ),
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
*> \ingroup single_eig
*
*  =====================================================================
      SUBROUTINE SBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDB, LDC, LDU, M, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               B( LDB, * ), C( LDC, * ), U( LDU, * ),
     $                   WORK( * )
*     ..
*
* ======================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      REAL               BNORM, EPS, REALMN
*     ..
*     .. External Functions ..
      REAL               SASUM, SLAMCH, SLANGE
      EXTERNAL           SASUM, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      RESID = ZERO
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
      REALMN = REAL( MAX( M, N ) )
      EPS = SLAMCH( 'Precision' )
*
*     Compute norm(B - U * C)
*
      DO 10 J = 1, N
         CALL SCOPY( M, B( 1, J ), 1, WORK, 1 )
         CALL SGEMV( 'No transpose', M, M, -ONE, U, LDU, C( 1, J ), 1,
     $               ONE, WORK, 1 )
         RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
   10 CONTINUE
*
*     Compute norm of B.
*
      BNORM = SLANGE( '1', M, N, B, LDB, WORK )
*
      IF( BNORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         IF( BNORM.GE.RESID ) THEN
            RESID = ( RESID / BNORM ) / ( REALMN*EPS )
         ELSE
            IF( BNORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REALMN*BNORM ) / BNORM ) /
     $                 ( REALMN*EPS )
            ELSE
               RESID = MIN( RESID / BNORM, REALMN ) / ( REALMN*EPS )
            END IF
         END IF
      END IF
      RETURN
*
*     End of SBDT02
*
      END
