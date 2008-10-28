      SUBROUTINE ZGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK,
     $                   RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            ITYPE, LDA, LDB, LDU, LDV, N
      DOUBLE PRECISION   RESULT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), U( LDU, * ),
     $                   V( LDV, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*       ZGET51  generally checks a decomposition of the form
*
*               A = U B V*
*
*       where * means conjugate transpose and U and V are unitary.
*
*       Specifically, if ITYPE=1
*
*               RESULT = | A - U B V* | / ( |A| n ulp )
*
*       If ITYPE=2, then:
*
*               RESULT = | A - B | / ( |A| n ulp )
*
*       If ITYPE=3, then:
*
*               RESULT = | I - UU* | / ( n ulp )
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          Specifies the type of tests to be performed.
*          =1: RESULT = | A - U B V* | / ( |A| n ulp )
*          =2: RESULT = | A - B | / ( |A| n ulp )
*          =3: RESULT = | I - UU* | / ( n ulp )
*
*  N       (input) INTEGER
*          The size of the matrix.  If it is zero, ZGET51 does nothing.
*          It must be at least zero.
*
*  A       (input) COMPLEX*16 array, dimension (LDA, N)
*          The original (unfactored) matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  It must be at least 1
*          and at least N.
*
*  B       (input) COMPLEX*16 array, dimension (LDB, N)
*          The factored matrix.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  It must be at least 1
*          and at least N.
*
*  U       (input) COMPLEX*16 array, dimension (LDU, N)
*          The unitary matrix on the left-hand side in the
*          decomposition.
*          Not referenced if ITYPE=2
*
*  LDU     (input) INTEGER
*          The leading dimension of U.  LDU must be at least N and
*          at least 1.
*
*  V       (input) COMPLEX*16 array, dimension (LDV, N)
*          The unitary matrix on the left-hand side in the
*          decomposition.
*          Not referenced if ITYPE=2
*
*  LDV     (input) INTEGER
*          The leading dimension of V.  LDV must be at least N and
*          at least 1.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (2*N**2)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RESULT  (output) DOUBLE PRECISION
*          The values computed by the test specified by ITYPE.  The
*          value is currently limited to 1/ulp, to avoid overflow.
*          Errors are flagged by RESULT=10/ulp.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TEN
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 10.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            JCOL, JDIAG, JROW
      DOUBLE PRECISION   ANORM, ULP, UNFL, WNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMM, ZLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      RESULT = ZERO
      IF( N.LE.0 )
     $   RETURN
*
*     Constants
*
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
*
*     Some Error Checks
*
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         RESULT = TEN / ULP
         RETURN
      END IF
*
      IF( ITYPE.LE.2 ) THEN
*
*        Tests scaled by the norm(A)
*
         ANORM = MAX( ZLANGE( '1', N, N, A, LDA, RWORK ), UNFL )
*
         IF( ITYPE.EQ.1 ) THEN
*
*           ITYPE=1: Compute W = A - UBV'
*
            CALL ZLACPY( ' ', N, N, A, LDA, WORK, N )
            CALL ZGEMM( 'N', 'N', N, N, N, CONE, U, LDU, B, LDB, CZERO,
     $                  WORK( N**2+1 ), N )
*
            CALL ZGEMM( 'N', 'C', N, N, N, -CONE, WORK( N**2+1 ), N, V,
     $                  LDV, CONE, WORK, N )
*
         ELSE
*
*           ITYPE=2: Compute W = A - B
*
            CALL ZLACPY( ' ', N, N, B, LDB, WORK, N )
*
            DO 20 JCOL = 1, N
               DO 10 JROW = 1, N
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) )
     $                - A( JROW, JCOL )
   10          CONTINUE
   20       CONTINUE
         END IF
*
*        Compute norm(W)/ ( ulp*norm(A) )
*
         WNORM = ZLANGE( '1', N, N, WORK, N, RWORK )
*
         IF( ANORM.GT.WNORM ) THEN
            RESULT = ( WNORM / ANORM ) / ( N*ULP )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESULT = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
            ELSE
               RESULT = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
            END IF
         END IF
*
      ELSE
*
*        Tests not scaled by norm(A)
*
*        ITYPE=3: Compute  UU' - I
*
         CALL ZGEMM( 'N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO,
     $               WORK, N )
*
         DO 30 JDIAG = 1, N
            WORK( ( N+1 )*( JDIAG-1 )+1 ) = WORK( ( N+1 )*( JDIAG-1 )+
     $         1 ) - CONE
   30    CONTINUE
*
         RESULT = MIN( ZLANGE( '1', N, N, WORK, N, RWORK ),
     $            DBLE( N ) ) / ( N*ULP )
      END IF
*
      RETURN
*
*     End of ZGET51
*
      END
