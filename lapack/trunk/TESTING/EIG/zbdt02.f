      SUBROUTINE ZBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RWORK,
     $                   RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDB, LDC, LDU, M, N
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         B( LDB, * ), C( LDC, * ), U( LDU, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZBDT02 tests the change of basis C = U' * B by computing the residual
*
*     RESID = norm( B - U * C ) / ( max(m,n) * norm(B) * EPS ),
*
*  where B and C are M by N matrices, U is an M by M orthogonal matrix,
*  and EPS is the machine precision.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrices B and C and the order of
*          the matrix Q.
*
*  N       (input) INTEGER
*          The number of columns of the matrices B and C.
*
*  B       (input) COMPLEX*16 array, dimension (LDB,N)
*          The m by n matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  C       (input) COMPLEX*16 array, dimension (LDC,N)
*          The m by n matrix C, assumed to contain U' * B.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C.  LDC >= max(1,M).
*
*  U       (input) COMPLEX*16 array, dimension (LDU,M)
*          The m by m orthogonal matrix U.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= max(1,M).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (M)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (M)
*
*  RESID   (output) DOUBLE PRECISION
*          RESID = norm( B - U * C ) / ( max(m,n) * norm(B) * EPS ),
*
* ======================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   BNORM, EPS, REALMN
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANGE
      EXTERNAL           DLAMCH, DZASUM, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      RESID = ZERO
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
      REALMN = DBLE( MAX( M, N ) )
      EPS = DLAMCH( 'Precision' )
*
*     Compute norm( B - U * C )
*
      DO 10 J = 1, N
         CALL ZCOPY( M, B( 1, J ), 1, WORK, 1 )
         CALL ZGEMV( 'No transpose', M, M, -DCMPLX( ONE ), U, LDU,
     $               C( 1, J ), 1, DCMPLX( ONE ), WORK, 1 )
         RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   10 CONTINUE
*
*     Compute norm of B.
*
      BNORM = ZLANGE( '1', M, N, B, LDB, RWORK )
*
      IF( BNORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         IF( BNORM.GE.RESID ) THEN
            RESID = ( RESID / BNORM ) / ( REALMN*EPS )
         ELSE
            IF( BNORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REALMN*BNORM ) / BNORM ) /
     $                 ( REALMN*EPS )
            ELSE
               RESID = MIN( RESID / BNORM, REALMN ) / ( REALMN*EPS )
            END IF
         END IF
      END IF
      RETURN
*
*     End of ZBDT02
*
      END
