      SUBROUTINE CHST01( N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK,
     $                   LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, LDA, LDH, LDQ, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), H( LDH, * ), Q( LDQ, * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  CHST01 tests the reduction of a general matrix A to upper Hessenberg
*  form:  A = Q*H*Q'.  Two test ratios are computed;
*
*  RESULT(1) = norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
*  RESULT(2) = norm( I - Q'*Q ) / ( N * EPS )
*
*  The matrix Q is assumed to be given explicitly as it would be
*  following CGEHRD + CUNGHR.
*
*  In this version, ILO and IHI are not used, but they could be used
*  to save some work if this is desired.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          A is assumed to be upper triangular in rows and columns
*          1:ILO-1 and IHI+1:N, so Q differs from the identity only in
*          rows and columns ILO+1:IHI.
*
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The original n by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  H       (input) COMPLEX array, dimension (LDH,N)
*          The upper Hessenberg matrix H from the reduction A = Q*H*Q'
*          as computed by CGEHRD.  H is assumed to be zero below the
*          first subdiagonal.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H.  LDH >= max(1,N).
*
*  Q       (input) COMPLEX array, dimension (LDQ,N)
*          The orthogonal matrix Q from the reduction A = Q*H*Q' as
*          computed by CGEHRD + CUNGHR.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= 2*N*N.
*
*  RWORK   (workspace) REAL array, dimension (N)
*
*  RESULT  (output) REAL array, dimension (2)
*          RESULT(1) = norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
*          RESULT(2) = norm( I - Q'*Q ) / ( N * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            LDWORK
      REAL               ANORM, EPS, OVFL, SMLNUM, UNFL, WNORM
*     ..
*     .. External Functions ..
      REAL               CLANGE, SLAMCH
      EXTERNAL           CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CLACPY, CUNT01, SLABAD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         RESULT( 1 ) = ZERO
         RESULT( 2 ) = ZERO
         RETURN
      END IF
*
      UNFL = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      OVFL = ONE / UNFL
      CALL SLABAD( UNFL, OVFL )
      SMLNUM = UNFL*N / EPS
*
*     Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
*
*     Copy A to WORK
*
      LDWORK = MAX( 1, N )
      CALL CLACPY( ' ', N, N, A, LDA, WORK, LDWORK )
*
*     Compute Q*H
*
      CALL CGEMM( 'No transpose', 'No transpose', N, N, N, CMPLX( ONE ),
     $            Q, LDQ, H, LDH, CMPLX( ZERO ), WORK( LDWORK*N+1 ),
     $            LDWORK )
*
*     Compute A - Q*H*Q'
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', N, N, N,
     $            CMPLX( -ONE ), WORK( LDWORK*N+1 ), LDWORK, Q, LDQ,
     $            CMPLX( ONE ), WORK, LDWORK )
*
      ANORM = MAX( CLANGE( '1', N, N, A, LDA, RWORK ), UNFL )
      WNORM = CLANGE( '1', N, N, WORK, LDWORK, RWORK )
*
*     Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)
*
      RESULT( 1 ) = MIN( WNORM, ANORM ) / MAX( SMLNUM, ANORM*EPS ) / N
*
*     Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )
*
      CALL CUNT01( 'Columns', N, N, Q, LDQ, WORK, LWORK, RWORK,
     $             RESULT( 2 ) )
*
      RETURN
*
*     End of CHST01
*
      END
