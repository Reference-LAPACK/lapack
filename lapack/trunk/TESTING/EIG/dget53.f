      SUBROUTINE DGET53( A, LDA, B, LDB, SCALE, WR, WI, RESULT, INFO )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB
      DOUBLE PRECISION   RESULT, SCALE, WI, WR
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGET53  checks the generalized eigenvalues computed by DLAG2.
*
*  The basic test for an eigenvalue is:
*
*                               | det( s A - w B ) |
*      RESULT =  ---------------------------------------------------
*                ulp max( s norm(A), |w| norm(B) )*norm( s A - w B )
*
*  Two "safety checks" are performed:
*
*  (1)  ulp*max( s*norm(A), |w|*norm(B) )  must be at least
*       safe_minimum.  This insures that the test performed is
*       not essentially  det(0*A + 0*B)=0.
*
*  (2)  s*norm(A) + |w|*norm(B) must be less than 1/safe_minimum.
*       This insures that  s*A - w*B  will not overflow.
*
*  If these tests are not passed, then  s  and  w  are scaled and
*  tested anyway, if this is possible.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA, 2)
*          The 2x2 matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  It must be at least 2.
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB, N)
*          The 2x2 upper-triangular matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  It must be at least 2.
*
*  SCALE   (input) DOUBLE PRECISION
*          The "scale factor" s in the formula  s A - w B .  It is
*          assumed to be non-negative.
*
*  WR      (input) DOUBLE PRECISION
*          The real part of the eigenvalue  w  in the formula
*          s A - w B .
*
*  WI      (input) DOUBLE PRECISION
*          The imaginary part of the eigenvalue  w  in the formula
*          s A - w B .
*
*  RESULT  (output) DOUBLE PRECISION
*          If INFO is 2 or less, the value computed by the test
*             described above.
*          If INFO=3, this will just be 1/ulp.
*
*  INFO    (output) INTEGER
*          =0:  The input data pass the "safety checks".
*          =1:  s*norm(A) + |w|*norm(B) > 1/safe_minimum.
*          =2:  ulp*max( s*norm(A), |w|*norm(B) ) < safe_minimum
*          =3:  same as INFO=2, but  s  and  w  could not be scaled so
*               as to compute the test.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ABSW, ANORM, BNORM, CI11, CI12, CI22, CNORM,
     $                   CR11, CR12, CR21, CR22, CSCALE, DETI, DETR, S1,
     $                   SAFMIN, SCALES, SIGMIN, TEMP, ULP, WIS, WRS
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
*     Initialize
*
      INFO = 0
      RESULT = ZERO
      SCALES = SCALE
      WRS = WR
      WIS = WI
*
*     Machine constants and norms
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ABSW = ABS( WRS ) + ABS( WIS )
      ANORM = MAX( ABS( A( 1, 1 ) )+ABS( A( 2, 1 ) ),
     $        ABS( A( 1, 2 ) )+ABS( A( 2, 2 ) ), SAFMIN )
      BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 1, 2 ) )+ABS( B( 2, 2 ) ),
     $        SAFMIN )
*
*     Check for possible overflow.
*
      TEMP = ( SAFMIN*BNORM )*ABSW + ( SAFMIN*ANORM )*SCALES
      IF( TEMP.GE.ONE ) THEN
*
*        Scale down to avoid overflow
*
         INFO = 1
         TEMP = ONE / TEMP
         SCALES = SCALES*TEMP
         WRS = WRS*TEMP
         WIS = WIS*TEMP
         ABSW = ABS( WRS ) + ABS( WIS )
      END IF
      S1 = MAX( ULP*MAX( SCALES*ANORM, ABSW*BNORM ),
     $     SAFMIN*MAX( SCALES, ABSW ) )
*
*     Check for W and SCALE essentially zero.
*
      IF( S1.LT.SAFMIN ) THEN
         INFO = 2
         IF( SCALES.LT.SAFMIN .AND. ABSW.LT.SAFMIN ) THEN
            INFO = 3
            RESULT = ONE / ULP
            RETURN
         END IF
*
*        Scale up to avoid underflow
*
         TEMP = ONE / MAX( SCALES*ANORM+ABSW*BNORM, SAFMIN )
         SCALES = SCALES*TEMP
         WRS = WRS*TEMP
         WIS = WIS*TEMP
         ABSW = ABS( WRS ) + ABS( WIS )
         S1 = MAX( ULP*MAX( SCALES*ANORM, ABSW*BNORM ),
     $        SAFMIN*MAX( SCALES, ABSW ) )
         IF( S1.LT.SAFMIN ) THEN
            INFO = 3
            RESULT = ONE / ULP
            RETURN
         END IF
      END IF
*
*     Compute C = s A - w B
*
      CR11 = SCALES*A( 1, 1 ) - WRS*B( 1, 1 )
      CI11 = -WIS*B( 1, 1 )
      CR21 = SCALES*A( 2, 1 )
      CR12 = SCALES*A( 1, 2 ) - WRS*B( 1, 2 )
      CI12 = -WIS*B( 1, 2 )
      CR22 = SCALES*A( 2, 2 ) - WRS*B( 2, 2 )
      CI22 = -WIS*B( 2, 2 )
*
*     Compute the smallest singular value of s A - w B:
*
*                 |det( s A - w B )|
*     sigma_min = ------------------
*                 norm( s A - w B )
*
      CNORM = MAX( ABS( CR11 )+ABS( CI11 )+ABS( CR21 ),
     $        ABS( CR12 )+ABS( CI12 )+ABS( CR22 )+ABS( CI22 ), SAFMIN )
      CSCALE = ONE / SQRT( CNORM )
      DETR = ( CSCALE*CR11 )*( CSCALE*CR22 ) -
     $       ( CSCALE*CI11 )*( CSCALE*CI22 ) -
     $       ( CSCALE*CR12 )*( CSCALE*CR21 )
      DETI = ( CSCALE*CR11 )*( CSCALE*CI22 ) +
     $       ( CSCALE*CI11 )*( CSCALE*CR22 ) -
     $       ( CSCALE*CI12 )*( CSCALE*CR21 )
      SIGMIN = ABS( DETR ) + ABS( DETI )
      RESULT = SIGMIN / S1
      RETURN
*
*     End of DGET53
*
      END
