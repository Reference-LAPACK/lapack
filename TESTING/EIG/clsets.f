*> \brief \b CLSETS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLSETS( M, P, N, A, AF, LDA, B, BF, LDB, C, CF,
*                          D, DF, X, WORK, LWORK, RWORK, RESULT )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDB, LWORK, M, P, N
*       ..
*       .. Array Arguments ..
*       REAL               RESULT( 2 ), RWORK( * )
*       COMPLEX            A( LDA, * ), AF( LDA, * ), B( LDB, * ),
*      $                   BF( LDB, * ), C( * ), D( * ), CF( * ),
*      $                   DF( * ), WORK( LWORK ), X( * )
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLSETS tests CGGLSE - a subroutine for solving linear equality
*> constrained least square problem (LSE).
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
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          The M-by-N matrix A.
*> \endverbatim
*>
*> \param[out] AF
*> \verbatim
*>          AF is COMPLEX array, dimension (LDA,N)
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the arrays A, AF, Q and R.
*>          LDA >= max(M,N).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,N)
*>          The P-by-N matrix A.
*> \endverbatim
*>
*> \param[out] BF
*> \verbatim
*>          BF is COMPLEX array, dimension (LDB,N)
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the arrays B, BF, V and S.
*>          LDB >= max(P,N).
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is COMPLEX array, dimension( M )
*>          the vector C in the LSE problem.
*> \endverbatim
*>
*> \param[out] CF
*> \verbatim
*>          CF is COMPLEX array, dimension( M )
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is COMPLEX array, dimension( P )
*>          the vector D in the LSE problem.
*> \endverbatim
*>
*> \param[out] DF
*> \verbatim
*>          DF is COMPLEX array, dimension( P )
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX array, dimension( N )
*>          solution vector X in the LSE problem.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (LWORK)
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
*>          RWORK is REAL array, dimension (M)
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is REAL array, dimension (2)
*>          The test ratios:
*>            RESULT(1) = norm( A*x - c )/ norm(A)*norm(X)*EPS
*>            RESULT(2) = norm( B*x - d )/ norm(B)*norm(X)*EPS
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
*> \ingroup complex_eig
*
*  =====================================================================
      SUBROUTINE CLSETS( M, P, N, A, AF, LDA, B, BF, LDB, C, CF,
     $                   D, DF, X, WORK, LWORK, RWORK, RESULT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, P, N
*     ..
*     .. Array Arguments ..
      REAL               RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), B( LDB, * ),
     $                   BF( LDB, * ), C( * ), D( * ), CF( * ),
     $                   DF( * ), WORK( LWORK ), X( * )
*
*  ====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, J
      REAL              RESID, SCL
*     ..
*     .. Parameters ..
      REAL              ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      INTEGER            MAXEXP
      PARAMETER          ( MAXEXP = MAXEXPONENT( ZERO ) - 2 )
*     ..
*     .. External Functions ..
      LOGICAL            SISNAN
      REAL              CLANGE, SLAMCH
      EXTERNAL           SISNAN, CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CGGLSE, CLACPY, CGET02, CSSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          EXPONENT, MAX, MAXEXPONENT, SCALE
*     ..
*     .. Executable Statements ..
*
*     Copy the matrices A and B to the arrays AF and BF,
*     and the vectors C and D to the arrays CF and DF,
*
      CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL CLACPY( 'Full', P, N, B, LDB, BF, LDB )
      CALL CCOPY( M, C, 1, CF, 1 )
      CALL CCOPY( P, D, 1, DF, 1 )
*
*     Solve LSE problem
*
      CALL CGGLSE( M, N, P, AF, LDA, BF, LDB, CF, DF, X,
     $             WORK, LWORK, INFO )
*
*     Test the residual for the solution of LSE
*
*     Compute RESULT(1) = norm( A*x - c ) / norm(A)*norm(X)*EPS
*
      CALL CCOPY( M, C, 1, CF, 1 )
      CALL CCOPY( P, D, 1, DF, 1 )
      CALL CGET02( 'No transpose', M, N, 1, A, LDA, X, N, CF, M,
     $             RWORK, RESULT( 1 ) )
*
*     Compute result(2) = norm( B*x - d ) / norm(B)*norm(X)*EPS
*
      CALL CGET02( 'No transpose', P, N, 1, B, LDB, X, N, DF, P,
     $             RWORK, RESULT( 2 ) )
*
*     The problem is exactly invariant under scaling A, B, c and d by
*     one power of two, so solving it again with the largest entry near
*     the overflow threshold has to give the same residuals.
*
      SCL = MAX( CLANGE( 'M', M, N, A, LDA, RWORK ),
     $           CLANGE( 'M', P, N, B, LDB, RWORK ),
     $           CLANGE( 'M', M, 1, C, M, RWORK ),
     $           CLANGE( 'M', P, 1, D, P, RWORK ) )
      IF( SCL.GT.ZERO .AND. SCL.LE.SLAMCH( 'Overflow' ) ) THEN
         SCL = SCALE( ONE, MAXEXP-EXPONENT( SCL ) )
         CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )
         CALL CLACPY( 'Full', P, N, B, LDB, BF, LDB )
         CALL CCOPY( M, C, 1, CF, 1 )
         CALL CCOPY( P, D, 1, DF, 1 )
         DO 10 J = 1, N
            CALL CSSCAL( M, SCL, AF( 1, J ), 1 )
            CALL CSSCAL( P, SCL, BF( 1, J ), 1 )
   10    CONTINUE
         CALL CSSCAL( M, SCL, CF, 1 )
         CALL CSSCAL( P, SCL, DF, 1 )
*
         CALL CGGLSE( M, N, P, AF, LDA, BF, LDB, CF, DF, X, WORK,
     $                LWORK, INFO )
*
         CALL CCOPY( M, C, 1, CF, 1 )
         CALL CCOPY( P, D, 1, DF, 1 )
         CALL CGET02( 'No transpose', M, N, 1, A, LDA, X, N, CF, M,
     $                RWORK, RESID )
         IF( SISNAN( RESID ) ) THEN
            RESULT( 1 ) = ONE / SLAMCH( 'Epsilon' )
         ELSE
            RESULT( 1 ) = MAX( RESULT( 1 ), RESID )
         END IF
         CALL CGET02( 'No transpose', P, N, 1, B, LDB, X, N, DF, P,
     $                RWORK, RESID )
         IF( SISNAN( RESID ) ) THEN
            RESULT( 2 ) = ONE / SLAMCH( 'Epsilon' )
         ELSE
            RESULT( 2 ) = MAX( RESULT( 2 ), RESID )
         END IF
      END IF
*
      RETURN
*
*     End of CLSETS
*
      END
