      SUBROUTINE SLARFGP( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               ALPHA, TAU
*     ..
*     .. Array Arguments ..
      REAL               X( * )
*     ..
*
*  Purpose
*  =======
*
*  SLARFGP generates a real elementary reflector H of order n, such
*  that
*
*        H * ( alpha ) = ( beta ),   H**T * H = I.
*            (   x   )   (   0  )
*
*  where alpha and beta are scalars, beta is non-negative, and x is
*  an (n-1)-element real vector.  H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v**T ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (n-1)-element
*  vector.
*
*  If the elements of x are all zero, then tau = 0 and H is taken to be
*  the unit matrix.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) REAL
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) REAL array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) REAL
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               TWO, ONE, ZERO
      PARAMETER          ( TWO = 2.0E+0, ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      REAL               BETA, BIGNUM, SAVEALPHA, SMLNUM, XNORM
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLAPY2, SNRM2
      EXTERNAL           SLAMCH, SLAPY2, SNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = SNRM2( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0.
*
         IF( ALPHA.GE.ZERO ) THEN
*           When TAU.eq.ZERO, the vector is special-cased to be
*           all zeros in the application routines.  We do not need
*           to clear it.
            TAU = ZERO
         ELSE
*           However, the application routines rely on explicit
*           zero checks when TAU.ne.ZERO, and we must clear X.
            TAU = TWO
            DO J = 1, N-1
               X( 1 + (J-1)*INCX ) = 0
            END DO
            ALPHA = -ALPHA
         END IF
      ELSE
*
*        general case
*
         BETA = SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
         SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SMLNUM ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            BIGNUM = ONE / SMLNUM
   10       CONTINUE
            KNT = KNT + 1
            CALL SSCAL( N-1, BIGNUM, X, INCX )
            BETA = BETA*BIGNUM
            ALPHA = ALPHA*BIGNUM
            IF( ABS( BETA ).LT.SMLNUM )
     $         GO TO 10
*
*           New BETA is at most 1, at least SMLNUM
*
            XNORM = SNRM2( N-1, X, INCX )
            BETA = SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         SAVEALPHA = ALPHA
         ALPHA = ALPHA + BETA
         IF( BETA.LT.ZERO ) THEN
            BETA = -BETA
            TAU = -ALPHA / BETA
         ELSE
            ALPHA = XNORM * (XNORM/ALPHA)
            TAU = ALPHA / BETA
            ALPHA = -ALPHA
         END IF
*
         IF ( ABS(TAU).LE.SMLNUM ) THEN
*
*           In the case where the computed TAU ends up being a denormalized number,
*           it loses relative accuracy. This is a BIG problem. Solution: flush TAU 
*           to ZERO. This explains the next IF statement.
*
*           (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
*           (Thanks Pat. Thanks MathWorks.)
*
            IF( SAVEALPHA.GE.ZERO ) THEN
               TAU = ZERO
            ELSE
               TAU = TWO
               DO J = 1, N-1
                  X( 1 + (J-1)*INCX ) = 0
               END DO
               BETA = -SAVEALPHA
            END IF
*
         ELSE 
*
*           This is the general case.
*
            CALL SSCAL( N-1, ONE / ALPHA, X, INCX )
*
         END IF
*
*        If BETA is subnormal, it may lose relative accuracy
*
         DO 20 J = 1, KNT
            BETA = BETA*SMLNUM
 20      CONTINUE
         ALPHA = BETA
      END IF
*
      RETURN
*
*     End of SLARFGP
*
      END
