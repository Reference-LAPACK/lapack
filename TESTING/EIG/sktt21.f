*> \brief \b SKTT21
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKTT21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK,
*                          RESULT )
*
*       .. Scalar Arguments ..
*       INTEGER            KBAND, LDU, N
*       ..
*       .. Array Arguments ..
*       REAL               AD( * ), AE( * ), RESULT( 2 ), SD( * ),
*      $                   SE( * ), U( LDU, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKTT21 checks a decomposition of the form
*>
*>    A = U S U'
*>
*> where ' means transpose, A is skew-symmetric tridiagonal, U is orthogonal,
*> and S is diagonal (if KBAND=0) or skew-symmetric tridiagonal (if KBAND=1).
*> Two tests are performed:
*>
*>    RESULT(1) = | A - U S U' | / ( |A| n ulp )
*>
*>    RESULT(2) = | I - UU' | / ( n ulp )
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The size of the matrix.  If it is zero, SKTT21 does nothing.
*>          It must be at least zero.
*> \endverbatim
*>
*> \param[in] KBAND
*> \verbatim
*>          KBAND is INTEGER
*>          The bandwidth of the matrix S.  It may only be zero or one.
*>          If zero, then S is diagonal, and SE is not referenced.  If
*>          one, then S is skew-symmetric tri-diagonal.
*> \endverbatim
*>
*> \param[in] AD
*> \verbatim
*>          AD is REAL array, dimension (N)
*>          AD is not referenced.
*> \endverbatim
*>
*> \param[in] AE
*> \verbatim
*>          AE is REAL array, dimension (N-1)
*>          The off-diagonal of the original (unfactored) matrix A.  A
*>          is assumed to be skew-symmetric tridiagonal.  AE(1) is the (1,2)
*>          and (2,1) element, AE(2) is the (2,3) and (3,2) element, etc.
*> \endverbatim
*>
*> \param[in] SD
*> \verbatim
*>          SD is REAL array, dimension (N)
*>          SD is not referenced.
*> \endverbatim
*>
*> \param[in] SE
*> \verbatim
*>          SE is REAL array, dimension (N-1)
*>          The off-diagonal of the (skew-symmetric tri-) diagonal matrix S.
*>          Not referenced if KBSND=0.  If KBAND=1, then AE(1) is the
*>          (1,2) and (2,1) element, SE(2) is the (2,3) and (3,2)
*>          element, etc.
*> \endverbatim
*>
*> \param[in] U
*> \verbatim
*>          U is REAL array, dimension (LDU, N)
*>          The orthogonal matrix in the decomposition.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of U.  LDU must be at least N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (N*(N+1))
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is REAL array, dimension (2)
*>          The values computed by the two tests described above.  The
*>          values are currently limited to 1/ulp, to avoid overflow.
*>          RESULT(1) is always modified.
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
      SUBROUTINE SKTT21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK,
     $                   RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            KBAND, LDU, N
*     ..
*     .. Array Arguments ..
      REAL               AD( * ), AE( * ), RESULT( 2 ), SD( * ),
     $                   SE( * ), U( LDU, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      REAL               ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANKY
      EXTERNAL           SLAMCH, SLANGE, SLANKY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SLASET, SKYR2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     1)      Constants
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )
*
*     Do Test 1
*
*     Copy A & Compute its 1-Norm:
*
      CALL SLASET( 'Full', N, N, ZERO, ZERO, WORK, N )
*
      ANORM = ZERO
      TEMP1 = ZERO
*
      DO 10 J = 1, N - 1
         WORK( ( N+1 )*( J-1 )+1 ) = ZERO
         WORK( ( N+1 )*( J-1 )+2 ) = AE( J )
         TEMP2 = ABS( AE( J ) )
         ANORM = MAX( ANORM, ABS( ZERO )+TEMP1+TEMP2 )
         TEMP1 = TEMP2
   10 CONTINUE
*
      WORK( N**2 ) = ZERO
      ANORM = MAX( ANORM, ABS( ZERO )+TEMP1, UNFL )
*
*     Norm of A - USU'
*
      IF( N.GT.1 .AND. KBAND.EQ.1 ) THEN
         DO 30 J = 1, N - 1
            CALL SKYR2( 'L', N, -SE( J ), U( 1, J ), 1, U( 1, J+1 ), 1,
     $                  WORK, N )
   30    CONTINUE
      END IF
*
      WNORM = SLANKY( '1', 'L', N, WORK, N, WORK( N**2+1 ) )
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         ELSE
            RESULT( 1 ) = MIN( WNORM / ANORM, REAL( N ) ) / ( N*ULP )
         END IF
      END IF
*
*     Do Test 2
*
*     Compute  UU' - I
*
      CALL SGEMM( 'N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK,
     $            N )
*
      DO 40 J = 1, N
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
   40 CONTINUE
*
      RESULT( 2 ) = MIN( REAL( N ), SLANGE( '1', N, N, WORK, N,
     $              WORK( N**2+1 ) ) ) / ( N*ULP )
*
      RETURN
*
*     End of SKTT21
*
      END
