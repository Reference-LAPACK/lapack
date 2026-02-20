*> \brief \b DLAGTSK solves the system of equations (T-λI)x = y
*> or (T-λI)^Tx = y, where T is a general tridiagonal matrix
*> with purely imaginary diagonal, and λ a purely imaginary scalar,
*> using the LU factorization computed by slagtfk.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DLAGTSK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlagtsk.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlagtsk.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlagtsk.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAGTSK( JOB, N, A, B, C, D, IN, YR, YI, TOL, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, JOB, N
*       DOUBLE PRECISION   TOL
*       ..
*       .. Array Arguments ..
*       INTEGER            IN( * )
*       DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), YR( * ), YI( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAGTSK may be used to solve one of the systems of equations
*>
*>    (T - lambda*I)*x = y   or   (T - lambda*I)**T*x = y,
*>
*> where T is an n by n tridiagonal matrix, for x, following the
*> factorization of (T - lambda*I) as
*>
*>    (T - lambda*I) = P*L*U ,
*>
*> by routine SLAGTF. The choice of equation to be solved is
*> controlled by the argument JOB, and in each case there is an option
*> to perturb zero or very small diagonal elements of U, this option
*> being intended for use in applications such as inverse iteration.
*>
*> Purely imaginary values in this subroutine are stored with imaginary
*> part, and every value in this subroutine are either real, or
*> purely imaginary, which can be determined by the permutation matrix
*> stored in IN with the rule in slagtfk.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is INTEGER
*>          Specifies the job to be performed by DLAGTSK as follows:
*>          =  1: The equations  (T - lambda*I)x = y  are to be solved,
*>                but diagonal elements of U are not to be perturbed.
*>          = -1: The equations  (T - lambda*I)x = y  are to be solved
*>                and, if overflow would otherwise occur, the diagonal
*>                elements of U are to be perturbed. See argument TOL
*>                below.
*>          =  2: The equations  (T - lambda*I)**Tx = y  are to be solved,
*>                but diagonal elements of U are not to be perturbed.
*>                **T here is conjugate transpose.
*>          = -2: The equations  (T - lambda*I)**Tx = y  are to be solved
*>                and, if overflow would otherwise occur, the diagonal
*>                elements of U are to be perturbed. See argument TOL
*>                below. **T here is conjugate transpose.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix T.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (N)
*>          On entry, A must contain the diagonal elements of U as
*>          returned from SLAGTF.
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (N-1)
*>          On entry, B must contain the first super-diagonal elements of
*>          U as returned from SLAGTF.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (N-1)
*>          On entry, C must contain the sub-diagonal elements of L as
*>          returned from SLAGTF.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N-2)
*>          On entry, D must contain the second super-diagonal elements
*>          of U as returned from SLAGTF.
*> \endverbatim
*>
*> \param[in] IN
*> \verbatim
*>          IN is INTEGER array, dimension (N+1)
*>          On entry, IN must contain details of the matrix P as returned
*>          from SLAGTF.
*> \endverbatim
*>
*> \param[in,out] YR
*> \verbatim
*>          YR is DOUBLE PRECISION array, dimension (N)
*>          On entry, the real part of right hand side vector y.
*>          On exit, YR is overwritten by the real part of solution vector x.
*> \endverbatim
*>
*> \param[in,out] YI
*> \verbatim
*>          YI is DOUBLE PRECISION array, dimension (N)
*>          On entry, the imaginary part of right hand side vector y.
*>          On exit, YI is overwritten by the imaginary part of solution vector x.
*> \endverbatim
*>
*> \param[in,out] TOL
*> \verbatim
*>          TOL is DOUBLE PRECISION
*>          On entry, with  JOB < 0, TOL should be the minimum
*>          perturbation to be made to very small diagonal elements of U.
*>          TOL should normally be chosen as about eps*norm(U), where eps
*>          is the relative machine precision, but if TOL is supplied as
*>          non-positive, then it is reset to eps*max( abs( u(i,j) ) ).
*>          If  JOB > 0  then TOL is not referenced.
*>
*>          On exit, TOL is changed as described above, only if TOL is
*>          non-positive on entry. Otherwise TOL is unchanged.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: overflow would occur when computing the INFO(th)
*>               element of the solution vector x. This can only occur
*>               when JOB is supplied as positive and either means
*>               that a diagonal element of U is very small, or that
*>               the elements of the right-hand side vector y are very
*>               large.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup lagts
*
*  =====================================================================
      SUBROUTINE DLAGTSK( JOB, N, A, B, C, D, IN, YR, YI, TOL,
     $                    INFO )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, JOB, N
      DOUBLE PRECISION   TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), YR( * ),
     $                   YI( * ) 
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            K
      DOUBLE PRECISION   ABSAK, AK, BIGNUM, EPS, PERT, SFMIN,
     $                   TEMPR, TEMPI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( ( ABS( JOB ).GT.2 ) .OR. ( JOB.EQ.0 ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAGTSK', -INFO )
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
*
      EPS = DLAMCH( 'Epsilon' )
      SFMIN = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SFMIN
*
      IF( JOB.LT.0 ) THEN
         IF( TOL.LE.ZERO ) THEN
            TOL = ABS( A( 1 ) )
            IF( N.GT.1 )
     $         TOL = MAX( TOL, ABS( A( 2 ) ), ABS( B( 1 ) ) )
            DO 10 K = 3, N
               TOL = MAX( TOL, ABS( A( K ) ), ABS( B( K-1 ) ),
     $               ABS( D( K-2 ) ) )
   10       CONTINUE
            TOL = TOL*EPS
            IF( TOL.EQ.ZERO )
     $         TOL = EPS
         END IF
      END IF
*
      IF( ABS( JOB ).EQ.1 ) THEN
         DO 20 K = 2, N
            IF( MOD( IN( K-1 ), 2 ).EQ.0 ) THEN
               IF ( IN( K-1 ) .LT. 2 ) THEN
                  YR( K ) = YR( K ) + C( K-1 )*YI( K-1 )
                  YI( K ) = YI( K ) - C( K-1 )*YR( K-1 )
               ELSE
                  YR( K ) = YR( K ) - C( K-1 )*YR( K-1 )
                  YI( K ) = YI( K ) - C( K-1 )*YI( K-1 )
               END IF
            ELSE
               TEMPR = YR( K-1 )
               TEMPI = YI( K-1 )
               YR( K-1 ) = YR( K )
               YI( K-1 ) = YI( K )
               IF ( IN( K-1 ) .LT. 2 ) THEN
                  YR( K ) = TEMPR + C( K-1 )*YI( K-1 )
                  YI( K ) = TEMPI - C( K-1 )*YR( K-1 )
               ELSE
                  YR( K ) = TEMPR - C( K-1 )*YR( K-1 )
                  YI( K ) = TEMPI - C( K-1 )*YI( K-1 )
               END IF
            END IF
   20    CONTINUE
         IF( JOB.EQ.1 ) THEN
            DO 30 K = N, 1, -1
               IF( K.LE.N-2 ) THEN
                  IF ( IN( K ) .EQ. 0 ) THEN
                     TEMPR = YR( K ) - B( K )*YR( K+1 ) -
     $                       D( K )*YR( K+2 )
                     TEMPI = YI( K ) - B( K )*YI( K+1 ) -
     $                       D( K )*YI( K+2 )
                  ELSE
                     TEMPR = YR( K ) + B( K )*YI( K+1 ) -
     $                       D( K )*YR( K+2 )
                     TEMPI = YI( K ) - B( K )*YR( K+1 ) -
     $                       D( K )*YI( K+2 )
                  END IF
               ELSE IF( K.EQ.N-1 ) THEN
                  IF ( IN( K ) .EQ. 0 ) THEN
                     TEMPR = YR( K ) - B( K )*YR( K+1 )
                     TEMPI = YI( K ) - B( K )*YI( K+1 )
                  ELSE
                     TEMPR = YR( K ) + B( K )*YI( K+1 )
                     TEMPI = YI( K ) - B( K )*YR( K+1 )
                  END IF
               ELSE
                  TEMPR = YR( K )
                  TEMPI = YI( K )
               END IF
               AK = A( K )
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR.
     $               ( ABS( TEMPR )+ABS( TEMPI ) )*SFMIN.GT.ABSAK )
     $                    THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMPR = TEMPR*BIGNUM
                        TEMPI = TEMPI*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ( ABS( TEMPR )+ABS( TEMPI ) ).GT.
     $                     ABSAK*BIGNUM ) THEN
                     INFO = K
                     RETURN
                  END IF
               END IF
               IF ( IN( K ) .EQ. 0 ) THEN
                  YR( K ) = TEMPI / AK
                  YI( K ) = -TEMPR / AK
               ELSE
                  YR( K ) = TEMPR / AK
                  YI( K ) = TEMPI / AK
               END IF
   30       CONTINUE
         ELSE
            DO 50 K = N, 1, -1
               IF( K.LE.N-2 ) THEN
                  IF ( IN( K ) .EQ. 0 ) THEN
                     TEMPR = YR( K ) - B( K )*YR( K+1 ) -
     $                       D( K )*YR( K+2 )
                     TEMPI = YI( K ) - B( K )*YI( K+1 ) -
     $                       D( K )*YI( K+2 )
                  ELSE
                     TEMPR = YR( K ) + B( K )*YI( K+1 ) -
     $                       D( K )*YR( K+2 )
                     TEMPI = YI( K ) - B( K )*YR( K+1 ) -
     $                       D( K )*YI( K+2 )
                  END IF
               ELSE IF( K.EQ.N-1 ) THEN
                  IF ( IN( K ) .EQ. 0 ) THEN
                     TEMPR = YR( K ) - B( K )*YR( K+1 )
                     TEMPI = YI( K ) - B( K )*YI( K+1 )
                  ELSE
                     TEMPR = YR( K ) + B( K )*YI( K+1 )
                     TEMPI = YI( K ) - B( K )*YR( K+1 )
                  END IF
               ELSE
                  TEMPR = YR( K )
                  TEMPI = YI( K )
               END IF
               AK = A( K )
               PERT = SIGN( TOL, AK )
   40          CONTINUE
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR.
     $               ( ABS( TEMPR )+ABS( TEMPI ) )*SFMIN.GT.ABSAK )
     $                    THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 40
                     ELSE
                        TEMPR = TEMPR*BIGNUM
                        TEMPI = TEMPI*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ( ABS( TEMPR )+ABS( TEMPI ) ).GT.
     $                     ABSAK*BIGNUM ) THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 40
                  END IF
               END IF
               IF ( IN( K ) .EQ. 0 ) THEN
                  YR( K ) = TEMPI / AK
                  YI( K ) = -TEMPR / AK
               ELSE
                  YR( K ) = TEMPR / AK
                  YI( K ) = TEMPI / AK
               END IF
   50       CONTINUE
         END IF
      ELSE
*
*        Come to here if  JOB = 2 or -2
*
         IF( JOB.EQ.2 ) THEN
            DO 60 K = 1, N
               IF( K.GE.3 ) THEN
                  IF ( IN( K-1 ) .EQ. 0 ) THEN
                     TEMPR = YR( K ) - B( K-1 )*YR( K-1 ) -
     $                       D( K-2 )*YR( K-2 )
                     TEMPI = YI( K ) - B( K-1 )*YI( K-1 ) -
     $                     D( K-2 )*YI( K-2 )
                  ELSE
                     TEMPR = YR( K ) - B( K-1 )*YI( K-1 ) -
     $                       D( K-2 )*YR( K-2 )
                     TEMPI = YI( K ) + B( K-1 )*YR( K-1 ) -
     $                       D( K-2 )*YI( K-2 )
                  END IF
               ELSE IF( K.EQ.2 ) THEN
                  IF ( IN( K-1 ) .EQ. 0 ) THEN
                     TEMPR = YR( K ) - B( K-1 )*YR( K-1 )
                     TEMPI = YI( K ) - B( K-1 )*YI( K-1 )
                  ELSE
                     TEMPR = YR( K ) - B( K-1 )*YI( K-1 )
                     TEMPI = YI( K ) + B( K-1 )*YR( K-1 )
                  END IF
               ELSE
                  TEMPR = YR( K )
                  TEMPI = YI( K )
               END IF
               AK = A( K )
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR.
     $               ( ABS( TEMPR )+ABS( TEMPI ) )*SFMIN.GT.ABSAK )
     $                    THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMPR = TEMPR*BIGNUM
                        TEMPI = TEMPI*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ( ABS( TEMPR )+ABS( TEMPI ) ).GT.
     $                     ABSAK*BIGNUM ) THEN
                     INFO = K
                     RETURN
                  END IF
               END IF
               IF ( IN( K ) .EQ. 0 ) THEN
                  YR( K ) = -TEMPI / AK
                  YI( K ) = TEMPR / AK
               ELSE
                  YR( K ) = TEMPR / AK
                  YI( K ) = TEMPI / AK
               END IF
   60       CONTINUE
         ELSE
            DO 80 K = 1, N
               IF( K.GE.3 ) THEN
                  IF ( IN( K-1 ) .EQ. 0 ) THEN
                     TEMPR = YR( K ) - B( K-1 )*YR( K-1 ) -
     $                     D( K-2 )*YR( K-2 )
                     TEMPI = YI( K ) - B( K-1 )*YI( K-1 ) -
     $                       D( K-2 )*YI( K-2 )
                  ELSE
                     TEMPR = YR( K ) - B( K-1 )*YI( K-1 ) -
     $                       D( K-2 )*YR( K-2 )
                     TEMPI = YI( K ) + B( K-1 )*YR( K-1 ) -
     $                       D( K-2 )*YI( K-2 )
                  END IF
               ELSE IF( K.EQ.2 ) THEN
                  IF ( IN( K-1 ) .EQ. 0 ) THEN
                     TEMPR = YR( K ) - B( K-1 )*YR( K-1 )
                     TEMPI = YI( K ) - B( K-1 )*YI( K-1 )
                  ELSE
                     TEMPR = YR( K ) - B( K-1 )*YI( K-1 )
                     TEMPI = YI( K ) + B( K-1 )*YR( K-1 )
                  END IF
               ELSE
                  TEMPR = YR( K )
                  TEMPI = YI( K )
               END IF
               AK = A( K )
               PERT = SIGN( TOL, AK )
   70          CONTINUE
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR.
     $               ( ABS( TEMPR )+ABS( TEMPI ) )*SFMIN.GT.ABSAK )
     $                    THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 70
                     ELSE
                        TEMPR = TEMPR*BIGNUM
                        TEMPI = TEMPI*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ( ABS( TEMPR )+ABS( TEMPI ) ).GT.
     $                     ABSAK*BIGNUM ) THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 70
                  END IF
               END IF
               IF ( IN( K ) .EQ. 0 ) THEN
                  YR( K ) = -TEMPI / AK
                  YI( K ) = TEMPR / AK
               ELSE
                  YR( K ) = TEMPR / AK
                  YI( K ) = TEMPI / AK
               END IF
   80       CONTINUE
         END IF
*
         DO 90 K = N, 2, -1
            IF( MOD( IN( K-1 ), 2 ).EQ.0 ) THEN
               IF ( IN( K-1 ) .LT. 2 ) THEN
                  YR( K-1 ) = YR( K-1 ) - C( K-1 )*YI( K )
                  YI( K-1 ) = YI( K-1 ) + C( K-1 )*YR( K )
               ELSE
                  YR( K-1 ) = YR( K-1 ) - C( K-1 )*YR( K )
                  YI( K-1 ) = YI( K-1 ) - C( K-1 )*YI( K )
               END IF
            ELSE
               TEMPR = YR( K-1 )
               TEMPI = YI( K-1 )
               YR( K-1 ) = YR( K )
               YI( K-1 ) = YI( K )
               IF ( IN( K-1 ) .LT. 2 ) THEN
                  YR( K ) = TEMPR - C( K-1 )*YI( K-1 )
                  YI( K ) = TEMPI + C( K-1 )*YR( K-1 )
               ELSE
                  YR( K ) = TEMPR - C( K-1 )*YR( K-1 )
                  YI( K ) = TEMPI - C( K-1 )*YI( K-1 )
               END IF
            END IF
   90    CONTINUE
      END IF
*
*     End of DLAGTSK
*
      END
