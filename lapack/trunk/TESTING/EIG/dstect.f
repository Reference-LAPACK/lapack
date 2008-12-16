      SUBROUTINE DSTECT( N, A, B, SHIFT, NUM )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            N, NUM
      DOUBLE PRECISION   SHIFT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*     DSTECT counts the number NUM of eigenvalues of a tridiagonal
*     matrix T which are less than or equal to SHIFT. T has
*     diagonal entries A(1), ... , A(N), and offdiagonal entries
*     B(1), ..., B(N-1).
*     See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
*     Matrix", Report CS41, Computer Science Dept., Stanford
*     University, July 21, 1966
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The dimension of the tridiagonal matrix T.
*
*  A       (input) DOUBLE PRECISION array, dimension (N)
*          The diagonal entries of the tridiagonal matrix T.
*
*  B       (input) DOUBLE PRECISION array, dimension (N-1)
*          The offdiagonal entries of the tridiagonal matrix T.
*
*  SHIFT   (input) DOUBLE PRECISION
*          The shift, used as described under Purpose.
*
*  NUM     (output) INTEGER
*          The number of eigenvalues of T less than or equal
*          to SHIFT.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, THREE = 3.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   M1, M2, MX, OVFL, SOV, SSHIFT, SSUN, SUN, TMP,
     $                   TOM, U, UNFL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Get machine constants
*
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
*
*     Find largest entry
*
      MX = ABS( A( 1 ) )
      DO 10 I = 1, N - 1
         MX = MAX( MX, ABS( A( I+1 ) ), ABS( B( I ) ) )
   10 CONTINUE
*
*     Handle easy cases, including zero matrix
*
      IF( SHIFT.GE.THREE*MX ) THEN
         NUM = N
         RETURN
      END IF
      IF( SHIFT.LT.-THREE*MX ) THEN
         NUM = 0
         RETURN
      END IF
*
*     Compute scale factors as in Kahan's report
*     At this point, MX .NE. 0 so we can divide by it
*
      SUN = SQRT( UNFL )
      SSUN = SQRT( SUN )
      SOV = SQRT( OVFL )
      TOM = SSUN*SOV
      IF( MX.LE.ONE ) THEN
         M1 = ONE / MX
         M2 = TOM
      ELSE
         M1 = ONE
         M2 = TOM / MX
      END IF
*
*     Begin counting
*
      NUM = 0
      SSHIFT = ( SHIFT*M1 )*M2
      U = ( A( 1 )*M1 )*M2 - SSHIFT
      IF( U.LE.SUN ) THEN
         IF( U.LE.ZERO ) THEN
            NUM = NUM + 1
            IF( U.GT.-SUN )
     $         U = -SUN
         ELSE
            U = SUN
         END IF
      END IF
      DO 20 I = 2, N
         TMP = ( B( I-1 )*M1 )*M2
         U = ( ( A( I )*M1 )*M2-TMP*( TMP / U ) ) - SSHIFT
         IF( U.LE.SUN ) THEN
            IF( U.LE.ZERO ) THEN
               NUM = NUM + 1
               IF( U.GT.-SUN )
     $            U = -SUN
            ELSE
               U = SUN
            END IF
         END IF
   20 CONTINUE
      RETURN
*
*     End of DSTECT
*
      END
