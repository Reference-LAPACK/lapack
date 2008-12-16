      SUBROUTINE CLSETS( M, P, N, A, AF, LDA, B, BF, LDB, C, CF,
     $                   D, DF, X, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
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
*  Purpose
*  =======
*
*  CLSETS tests CGGLSE - a subroutine for solving linear equality
*  constrained least square problem (LSE).
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The M-by-N matrix A.
*
*  AF      (workspace) COMPLEX array, dimension (LDA,N)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and R.
*          LDA >= max(M,N).
*
*  B       (input) COMPLEX array, dimension (LDB,N)
*          The P-by-N matrix A.
*
*  BF      (workspace) COMPLEX array, dimension (LDB,N)
*
*  LDB     (input) INTEGER
*          The leading dimension of the arrays B, BF, V and S.
*          LDB >= max(P,N).
*
*  C       (input) COMPLEX array, dimension( M )
*          the vector C in the LSE problem.
*
*  CF      (workspace) COMPLEX array, dimension( M )
*
*  D       (input) COMPLEX array, dimension( P )
*          the vector D in the LSE problem.
*
*  DF      (workspace) COMPLEX array, dimension( P )
*
*  X       (output) COMPLEX array, dimension( N )
*          solution vector X in the LSE problem.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) REAL array, dimension (M)
*
*  RESULT  (output) REAL array, dimension (2)
*          The test ratios:
*            RESULT(1) = norm( A*x - c )/ norm(A)*norm(X)*EPS
*            RESULT(2) = norm( B*x - d )/ norm(B)*norm(X)*EPS
*
*  ====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGGLSE, CLACPY, CGET02
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
      RETURN
*
*     End of CLSETS
*
      END
