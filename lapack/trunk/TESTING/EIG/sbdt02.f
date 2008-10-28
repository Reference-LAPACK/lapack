      SUBROUTINE SBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDB, LDC, LDU, M, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               B( LDB, * ), C( LDC, * ), U( LDU, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SBDT02 tests the change of basis C = U' * B by computing the residual
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
*  B       (input) REAL array, dimension (LDB,N)
*          The m by n matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  C       (input) REAL array, dimension (LDC,N)
*          The m by n matrix C, assumed to contain U' * B.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C.  LDC >= max(1,M).
*
*  U       (input) REAL array, dimension (LDU,M)
*          The m by m orthogonal matrix U.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= max(1,M).
*
*  WORK    (workspace) REAL array, dimension (M)
*
*  RESID   (output) REAL
*          RESID = norm( B - U * C ) / ( max(m,n) * norm(B) * EPS ),
*
* ======================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      REAL               BNORM, EPS, REALMN
*     ..
*     .. External Functions ..
      REAL               SASUM, SLAMCH, SLANGE
      EXTERNAL           SASUM, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      RESID = ZERO
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
      REALMN = REAL( MAX( M, N ) )
      EPS = SLAMCH( 'Precision' )
*
*     Compute norm( B - U * C )
*
      DO 10 J = 1, N
         CALL SCOPY( M, B( 1, J ), 1, WORK, 1 )
         CALL SGEMV( 'No transpose', M, M, -ONE, U, LDU, C( 1, J ), 1,
     $               ONE, WORK, 1 )
         RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
   10 CONTINUE
*
*     Compute norm of B.
*
      BNORM = SLANGE( '1', M, N, B, LDB, WORK )
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
*     End of SBDT02
*
      END
