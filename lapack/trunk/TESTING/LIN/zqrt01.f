      SUBROUTINE ZQRT01( M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK,
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
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
     $                   R( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  ZQRT01 tests ZGEQRF, which computes the QR factorization of an m-by-n
*  matrix A, and partially tests ZUNGQR which forms the m-by-m
*  orthogonal matrix Q.
*
*  ZQRT01 compares R with Q'*A, and checks that Q is orthogonal.
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
*          Details of the QR factorization of A, as returned by ZGEQRF.
*          See ZGEQRF for further details.
*
*  Q       (output) COMPLEX*16 array, dimension (LDA,M)
*          The m-by-m orthogonal matrix Q.
*
*  R       (workspace) COMPLEX*16 array, dimension (LDA,max(M,N))
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF, Q and R.
*          LDA >= max(M,N).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors, as returned
*          by ZGEQRF.
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
*          RESULT(1) = norm( R - Q'*A ) / ( M * norm(A) * EPS )
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
      EXTERNAL           ZGEMM, ZGEQRF, ZHERK, ZLACPY, ZLASET, ZUNGQR
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
      SRNAMT = 'ZGEQRF'
      CALL ZGEQRF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy details of Q
*
      CALL ZLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      CALL ZLACPY( 'Lower', M-1, N, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )
*
*     Generate the m-by-m matrix Q
*
      SRNAMT = 'ZUNGQR'
      CALL ZUNGQR( M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL ZLASET( 'Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), R,
     $             LDA )
      CALL ZLACPY( 'Upper', M, N, AF, LDA, R, LDA )
*
*     Compute R - Q'*A
*
      CALL ZGEMM( 'Conjugate transpose', 'No transpose', M, N, M,
     $            DCMPLX( -ONE ), Q, LDA, A, LDA, DCMPLX( ONE ), R,
     $            LDA )
*
*     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
*
      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      RESID = ZLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL ZLASET( 'Full', M, M, DCMPLX( ZERO ), DCMPLX( ONE ), R, LDA )
      CALL ZHERK( 'Upper', 'Conjugate transpose', M, M, -ONE, Q, LDA,
     $            ONE, R, LDA )
*
*     Compute norm( I - Q'*Q ) / ( M * EPS ) .
*
      RESID = ZLANSY( '1', 'Upper', M, R, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / EPS
*
      RETURN
*
*     End of ZQRT01
*
      END
