*> \brief \b ZGLMTS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U,
*                          WORK, LWORK, RWORK, RESULT )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDB, LWORK, M, N, P
*       DOUBLE PRECISION   RESULT
*       ..
*       .. Array Arguments ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGLMTS tests ZGGGLM - a subroutine for solving the generalized
*> linear model problem.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of columns of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,M)
*>          The N-by-M matrix A.
*> \endverbatim
*>
*> \param[out] AF
*> \verbatim
*>          AF is COMPLEX*16 array, dimension (LDA,M)
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the arrays A, AF. LDA >= max(M,N).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,P)
*>          The N-by-P matrix A.
*> \endverbatim
*>
*> \param[out] BF
*> \verbatim
*>          BF is COMPLEX*16 array, dimension (LDB,P)
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the arrays B, BF. LDB >= max(P,N).
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is COMPLEX*16 array, dimension( N )
*>          On input, the left hand side of the GLM.
*> \endverbatim
*>
*> \param[out] DF
*> \verbatim
*>          DF is COMPLEX*16 array, dimension( N )
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension( M )
*>          solution vector X in the GLM problem.
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is COMPLEX*16 array, dimension( P )
*>          solution vector U in the GLM problem.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (M)
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is DOUBLE PRECISION
*>          The test ratio:
*>                           norm( d - A*x - B*u )
*>            RESULT = -----------------------------------------
*>                     (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
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
*> \ingroup complex16_eig
*
*  =====================================================================
      SUBROUTINE ZGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U,
     $                   WORK, LWORK, RWORK, RESULT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, N, P
      DOUBLE PRECISION   RESULT
*     ..
*     .. Array Arguments ..
*
*  ====================================================================
*
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ),
     $                   BF( LDB, * ), D( * ), DF( * ), U( * ),
     $                   WORK( LWORK ), X( * )
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      INTEGER            MAXEXP
      PARAMETER          ( MAXEXP = MAXEXPONENT( ZERO ) - 2 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, J
      DOUBLE PRECISION   ANORM, BNORM, DNORM, EPS, UNFL, XNORM, YNORM
      DOUBLE PRECISION    SCL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANGE
      LOGICAL            DISNAN
      EXTERNAL           DISNAN, DLAMCH, DZASUM, ZLANGE
*     ..
*     .. External Subroutines ..
*
      EXTERNAL           ZCOPY, ZGEMV, ZGGGLM, ZLACPY, ZDSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          EXPONENT, MAX, MAXEXPONENT, SCALE
*     ..
*     .. Executable Statements ..
*
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
      ANORM = MAX( ZLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( ZLANGE( '1', N, P, B, LDB, RWORK ), UNFL )
*
*     Copy the matrices A and B to the arrays AF and BF,
*     and the vector D the array DF.
*
      CALL ZLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL ZLACPY( 'Full', N, P, B, LDB, BF, LDB )
      CALL ZCOPY( N, D, 1, DF, 1 )
*
*     Solve GLM problem
*
      CALL ZGGGLM( N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK,
     $             INFO )
*
*     Test the residual for the solution of LSE
*
*                       norm( d - A*x - B*u )
*       RESULT = -----------------------------------------
*                (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
*
      CALL ZCOPY( N, D, 1, DF, 1 )
      CALL ZGEMV( 'No transpose', N, M, -CONE, A, LDA, X, 1, CONE, DF,
     $            1 )
*
      CALL ZGEMV( 'No transpose', N, P, -CONE, B, LDB, U, 1, CONE, DF,
     $            1 )
*
      DNORM = DZASUM( N, DF, 1 )
      XNORM = DZASUM( M, X, 1 ) + DZASUM( P, U, 1 )
      YNORM = ANORM + BNORM
*
      IF( XNORM.LE.ZERO ) THEN
         RESULT = ZERO
      ELSE
         RESULT = ( ( DNORM / YNORM ) / XNORM ) / EPS
      END IF
*
*     The problem is exactly invariant under scaling A, B and d by one
*     power of two, so solving it again with the largest entry near the
*     overflow threshold has to give the same residual.
*
      SCL = MAX( ZLANGE( 'M', N, M, A, LDA, RWORK ),
     $           ZLANGE( 'M', N, P, B, LDB, RWORK ),
     $           ZLANGE( 'M', N, 1, D, N, RWORK ) )
      IF( SCL.GT.ZERO .AND. SCL.LE.DLAMCH( 'Overflow' ) ) THEN
         SCL = SCALE( ONE, MAXEXP-EXPONENT( SCL ) )
         CALL ZLACPY( 'Full', N, M, A, LDA, AF, LDA )
         CALL ZLACPY( 'Full', N, P, B, LDB, BF, LDB )
         CALL ZCOPY( N, D, 1, DF, 1 )
         DO 10 J = 1, M
            CALL ZDSCAL( N, SCL, AF( 1, J ), 1 )
   10    CONTINUE
         DO 20 J = 1, P
            CALL ZDSCAL( N, SCL, BF( 1, J ), 1 )
   20    CONTINUE
         CALL ZDSCAL( N, SCL, DF, 1 )
*
         CALL ZGGGLM( N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK,
     $                INFO )
*
         CALL ZCOPY( N, D, 1, DF, 1 )
         CALL ZGEMV( 'No transpose', N, M, -CONE, A, LDA, X, 1, CONE,
     $               DF, 1 )
         CALL ZGEMV( 'No transpose', N, P, -CONE, B, LDB, U, 1, CONE,
     $               DF, 1 )
         DNORM = DZASUM( N, DF, 1 )
         XNORM = DZASUM( M, X, 1 ) + DZASUM( P, U, 1 )
         IF( DISNAN( DNORM ) .OR. DISNAN( XNORM ) ) THEN
            RESULT = ONE / EPS
         ELSE IF( XNORM.GT.ZERO ) THEN
            RESULT = MAX( RESULT, ( ( DNORM / YNORM ) / XNORM ) / EPS )
         END IF
      END IF
*
      RETURN
*
*     End of ZGLMTS
*
      END
