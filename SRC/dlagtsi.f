*> \brief \b DLAGTSI solves a tridiagonal inverse-iteration system.
*>
*> \par Purpose:
*> =============
*>
*> \verbatim
*> Given the pivoted factors returned by DLAGTF, DLAGTSI overwrites Y
*> with the solution of (T - shift I) X = SCALE Y. It may scale the
*> right hand side to preserve the direction of X without overflowing.
*> Only exact zero pivots are perturbed. A nonpositive TOL selects a
*> pivot perturbation relative to the upper triangular factor; on exit
*> TOL holds the perturbation used and may be retained across solves.
*> The product of the scaling factors can underflow when the exact
*> solution is outside the exponent range by more than one such scale.
*> \endverbatim
*>
*  Arguments:
*  ==========
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the tridiagonal matrix.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (N)
*>          Diagonal of U, as returned by DLAGTF.
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (N-1)
*>          First superdiagonal of U.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (N-1)
*>          Elimination multipliers from DLAGTF.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N-2)
*>          Second superdiagonal of U.
*> \endverbatim
*>
*> \param[in] IN
*> \verbatim
*>          IN is INTEGER array, dimension (N)
*>          Pivot information from DLAGTF.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array, dimension (N)
*>          On entry, the right hand side b. On exit, the scaled
*>          solution x of (T - shift I) x = SCALE b, up to any
*>          perturbations to exact zero pivots.
*> \endverbatim
*>
*> \param[in,out] TOL
*> \verbatim
*>          TOL is DOUBLE PRECISION
*>          A positive perturbation is used as supplied. Otherwise
*>          choose one relative to the largest entry of U and return it.
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>          The cumulative right hand side scale, at most one.
*> \endverbatim
*>
*> \ingroup lagts
*
      SUBROUTINE DLAGTSI( N, A, B, C, D, IN, Y, TOL, SCALE )
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N
      DOUBLE PRECISION   TOL, SCALE
*     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), Y( * )
*     .. Local Scalars ..
      INTEGER            K
      DOUBLE PRECISION   TEMP, PIV, BOUND, REC, PMAX
      DOUBLE PRECISION   ONE, ZERO, DLAMCH
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
      INTEGER            IDAMAX
      EXTERNAL           DLAMCH, IDAMAX, DSCAL
*
*     Keep headroom in the updates and in the subsequent dot products.
*
      SCALE = ONE
      IF( N.LE.0 ) RETURN
      BOUND = DLAMCH( 'Precision' ) / DLAMCH( 'Safe minimum' )
      IF( TOL.LE.ZERO ) THEN
         TOL = ABS( A( 1 ) )
         DO K = 2, N
            TOL = MAX( TOL, ABS( A( K ) ), ABS( B( K-1 ) ) )
         END DO
         DO K = 3, N
            TOL = MAX( TOL, ABS( D( K-2 ) ) )
         END DO
         TOL = TOL*DLAMCH( 'Epsilon' )
         IF( TOL.EQ.ZERO ) TOL = DLAMCH( 'Epsilon' )
      END IF
      PMAX = ABS( Y( IDAMAX( N, Y, 1 ) ) )
      IF( PMAX.GT.BOUND ) THEN
         REC = BOUND / PMAX
         CALL DSCAL( N, REC, Y, 1 )
         SCALE = REC
      END IF
*
*     Apply row interchanges and forward substitution.
*
      DO K = 1, N-1
         IF( IN( K ).NE.0 ) THEN
            TEMP = Y( K+1 )
            Y( K+1 ) = Y( K )
            Y( K ) = TEMP
         END IF
         IF( ABS( C( K ) ).GT.ONE ) THEN
            PMAX = BOUND / ABS( C( K ) )
            IF( ABS( Y( K ) ).GT.PMAX ) THEN
               REC = PMAX / ABS( Y( K ) )
               CALL DSCAL( N, REC, Y, 1 )
               SCALE = SCALE*REC
            END IF
         END IF
         Y( K+1 ) = Y( K+1 ) - C( K )*Y( K )
         IF( ABS( Y( K+1 ) ).GT.BOUND ) THEN
            REC = BOUND / ABS( Y( K+1 ) )
            CALL DSCAL( N, REC, Y, 1 )
            SCALE = SCALE*REC
         END IF
      END DO
*
*     Back substitute, rescaling the iterate instead of perturbing a
*     nonzero diagonal of U when division would overflow.
*
      DO K = N, 1, -1
         REC = ONE
         IF( K.LT.N ) THEN
            IF( ABS( B( K ) ).GT.ONE ) THEN
               PMAX = BOUND / ABS( B( K ) )
               IF( ABS( Y( K+1 ) ).GT.PMAX )
     $            REC = PMAX / ABS( Y( K+1 ) )
            END IF
         END IF
         IF( K.LT.N-1 ) THEN
            IF( ABS( D( K ) ).GT.ONE ) THEN
               PMAX = BOUND / ABS( D( K ) )
               IF( ABS( Y( K+2 ) ).GT.PMAX )
     $            REC = MIN( REC, PMAX / ABS( Y( K+2 ) ) )
            END IF
         END IF
         IF( REC.LT.ONE ) THEN
            CALL DSCAL( N, REC, Y, 1 )
            SCALE = SCALE*REC
         END IF
         TEMP = Y( K )
         IF( K.LT.N ) TEMP = TEMP - B( K )*Y( K+1 )
         IF( K.LT.N-1 ) TEMP = TEMP - D( K )*Y( K+2 )
         PIV = A( K )
         IF( PIV.EQ.ZERO ) PIV = SIGN( TOL, PIV )
         IF( ABS( PIV ).LT.ONE ) THEN
            PMAX = BOUND*ABS( PIV )
            IF( ABS( TEMP ).GT.PMAX ) THEN
*              Do not form a ratio that could be smaller than SAFMIN.
               IF( ABS( TEMP ).GT.ONE ) THEN
                  REC = ONE / ABS( TEMP )
                  CALL DSCAL( N, REC, Y, 1 )
                  SCALE = SCALE*REC
                  TEMP = TEMP*REC
               END IF
               IF( ABS( TEMP ).GT.PMAX ) THEN
                  REC = PMAX / ABS( TEMP )
                  CALL DSCAL( N, REC, Y, 1 )
                  SCALE = SCALE*REC
                  TEMP = TEMP*REC
               END IF
            END IF
         END IF
         Y( K ) = TEMP / PIV
         IF( ABS( Y( K ) ).GT.BOUND ) THEN
            REC = BOUND / ABS( Y( K ) )
            CALL DSCAL( N, REC, Y, 1 )
            SCALE = SCALE*REC
         END IF
      END DO
      RETURN
      END
