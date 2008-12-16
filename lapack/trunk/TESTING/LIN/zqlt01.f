      SUBROUTINE ZQLT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK,
     $                   RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RESULT( * ), RWORK( * )
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), L( LDA, * ),
     $                   Q( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  ZQLT01 tests ZGEQLF, which computes the QL factorization of an m-by-n
*  matrix A, and partially tests ZUNGQL which forms the m-by-m
*  orthogonal matrix Q.
*
*  ZQLT01 compares L with Q'*A, and checks that Q is orthogonal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The m-by-n matrix A.
*
*  AF      (output) COMPLEX*16 array, dimension (LDA,N)
*          Details of the QL factorization of A, as returned by ZGEQLF.
*          See ZGEQLF for further details.
*
*  Q       (output) COMPLEX*16 array, dimension (LDA,M)
*          The m-by-m orthogonal matrix Q.
*
*  L       (workspace) COMPLEX*16 array, dimension (LDA,max(M,N))
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and R.
*          LDA >= max(M,N).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors, as returned
*          by ZGEQLF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (M)
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (2)
*          The test ratios:
*          RESULT(1) = norm( L - Q'*A ) / ( M * norm(A) * EPS )
*          RESULT(2) = norm( I - Q'*Q ) / ( M * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         ROGUE
      PARAMETER          ( ROGUE = ( -1.0D+10, -1.0D+10 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, MINMN
      DOUBLE PRECISION   ANORM, EPS, RESID
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANSY
      EXTERNAL           DLAMCH, ZLANGE, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMM, ZGEQLF, ZHERK, ZLACPY, ZLASET, ZUNGQL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      MINMN = MIN( M, N )
      EPS = DLAMCH( 'Epsilon' )
*
*     Copy the matrix A to the array AF.
*
      CALL ZLACPY( 'Full', M, N, A, LDA, AF, LDA )
*
*     Factorize the matrix A in the array AF.
*
      SRNAMT = 'ZGEQLF'
      CALL ZGEQLF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy details of Q
*
      CALL ZLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      IF( M.GE.N ) THEN
         IF( N.LT.M .AND. N.GT.0 )
     $      CALL ZLACPY( 'Full', M-N, N, AF, LDA, Q( 1, M-N+1 ), LDA )
         IF( N.GT.1 )
     $      CALL ZLACPY( 'Upper', N-1, N-1, AF( M-N+1, 2 ), LDA,
     $                   Q( M-N+1, M-N+2 ), LDA )
      ELSE
         IF( M.GT.1 )
     $      CALL ZLACPY( 'Upper', M-1, M-1, AF( 1, N-M+2 ), LDA,
     $                   Q( 1, 2 ), LDA )
      END IF
*
*     Generate the m-by-m matrix Q
*
      SRNAMT = 'ZUNGQL'
      CALL ZUNGQL( M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy L
*
      CALL ZLASET( 'Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), L,
     $             LDA )
      IF( M.GE.N ) THEN
         IF( N.GT.0 )
     $      CALL ZLACPY( 'Lower', N, N, AF( M-N+1, 1 ), LDA,
     $                   L( M-N+1, 1 ), LDA )
      ELSE
         IF( N.GT.M .AND. M.GT.0 )
     $      CALL ZLACPY( 'Full', M, N-M, AF, LDA, L, LDA )
         IF( M.GT.0 )
     $      CALL ZLACPY( 'Lower', M, M, AF( 1, N-M+1 ), LDA,
     $                   L( 1, N-M+1 ), LDA )
      END IF
*
*     Compute L - Q'*A
*
      CALL ZGEMM( 'Conjugate transpose', 'No transpose', M, N, M,
     $            DCMPLX( -ONE ), Q, LDA, A, LDA, DCMPLX( ONE ), L,
     $            LDA )
*
*     Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      RESID = ZLANGE( '1', M, N, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL ZLASET( 'Full', M, M, DCMPLX( ZERO ), DCMPLX( ONE ), L, LDA )
      CALL ZHERK( 'Upper', 'Conjugate transpose', M, M, -ONE, Q, LDA,
     $            ONE, L, LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = ZLANSY( '1', 'Upper', M, L, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of ZQLT01
*
      END
