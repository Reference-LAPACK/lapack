      SUBROUTINE DHST01( N, ILO, IHI, A, LDA, H, LDH, Q, LDQ, WORK,
     $                   LWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, LDA, LDH, LDQ, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), H( LDH, * ), Q( LDQ, * ),
     $                   RESULT( 2 ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DHST01 tests the reduction of a general matrix A to upper Hessenberg
*  form:  A = Q*H*Q'.  Two test ratios are computed;
*
*  RESULT(1) = norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
*  RESULT(2) = norm( I - Q'*Q ) / ( N * EPS )
*
*  The matrix Q is assumed to be given explicitly as it would be
*  following DGEHRD + DORGHR.
*
*  In this version, ILO and IHI are not used and are assumed to be 1 and
*  N, respectively.
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
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The original n by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  H       (input) DOUBLE PRECISION array, dimension (LDH,N)
*          The upper Hessenberg matrix H from the reduction A = Q*H*Q'
*          as computed by DGEHRD.  H is assumed to be zero below the
*          first subdiagonal.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H.  LDH >= max(1,N).
*
*  Q       (input) DOUBLE PRECISION array, dimension (LDQ,N)
*          The orthogonal matrix Q from the reduction A = Q*H*Q' as
*          computed by DGEHRD + DORGHR.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= 2*N*N.
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (2)
*          RESULT(1) = norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
*          RESULT(2) = norm( I - Q'*Q ) / ( N * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            LDWORK
      DOUBLE PRECISION   ANORM, EPS, OVFL, SMLNUM, UNFL, WNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLABAD, DLACPY, DORT01
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
      UNFL = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      SMLNUM = UNFL*N / EPS
*
*     Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
*
*     Copy A to WORK
*
      LDWORK = MAX( 1, N )
      CALL DLACPY( ' ', N, N, A, LDA, WORK, LDWORK )
*
*     Compute Q*H
*
      CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE, Q, LDQ,
     $            H, LDH, ZERO, WORK( LDWORK*N+1 ), LDWORK )
*
*     Compute A - Q*H*Q'
*
      CALL DGEMM( 'No transpose', 'Transpose', N, N, N, -ONE,
     $            WORK( LDWORK*N+1 ), LDWORK, Q, LDQ, ONE, WORK,
     $            LDWORK )
*
      ANORM = MAX( DLANGE( '1', N, N, A, LDA, WORK( LDWORK*N+1 ) ),
     $        UNFL )
      WNORM = DLANGE( '1', N, N, WORK, LDWORK, WORK( LDWORK*N+1 ) )
*
*     Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)
*
      RESULT( 1 ) = MIN( WNORM, ANORM ) / MAX( SMLNUM, ANORM*EPS ) / N
*
*     Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )
*
      CALL DORT01( 'Columns', N, N, Q, LDQ, WORK, LWORK, RESULT( 2 ) )
*
      RETURN
*
*     End of DHST01
*
      END
