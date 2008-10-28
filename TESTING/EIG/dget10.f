      SUBROUTINE DGET10( M, N, A, LDA, B, LDB, WORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, M, N
      DOUBLE PRECISION   RESULT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGET10 compares two matrices A and B and computes the ratio
*  RESULT = norm( A - B ) / ( norm(A) * M * EPS )
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrices A and B.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
*          The m by n matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
*
*  RESULT  (output) DOUBLE PRECISION
*          RESULT = norm( A - B ) / ( norm(A) * M * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   ANORM, EPS, UNFL, WNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DASUM, DLAMCH, DLANGE
      EXTERNAL           DASUM, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RESULT = ZERO
         RETURN
      END IF
*
      UNFL = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
*
      WNORM = ZERO
      DO 10 J = 1, N
         CALL DCOPY( M, A( 1, J ), 1, WORK, 1 )
         CALL DAXPY( M, -ONE, B( 1, J ), 1, WORK, 1 )
         WNORM = MAX( WNORM, DASUM( N, WORK, 1 ) )
   10 CONTINUE
*
      ANORM = MAX( DLANGE( '1', M, N, A, LDA, WORK ), UNFL )
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT = ( WNORM / ANORM ) / ( M*EPS )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*EPS )
         ELSE
            RESULT = MIN( WNORM / ANORM, DBLE( M ) ) / ( M*EPS )
         END IF
      END IF
*
      RETURN
*
*     End of DGET10
*
      END
