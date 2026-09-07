*> \brief \b CCHKSY_2X2 checks extreme-scale 2-by-2 pivots.
*>
*> \par Purpose:
*> =============
*>
*> \verbatim
*> Test factorization and solve with subnormal and near-overflow pivots,
*> using both triangles and the classic, rook, RK and packed routines.
*> Orders 3 and 130 exercise unblocked code and panels with NB = 64.
*> Each right-hand side is a column of A; the solution is a unit vector.
*> Subnormal cases are skipped if gradual underflow is unavailable.
*> Failures and test counts are added to NFAIL and NRUN, respectively.
*> \endverbatim
*>
*> \param[in] NOUT
*>          Output unit for failure diagnostics.
*> \param[in,out] NFAIL
*>          Number of failed tests.
*> \param[in,out] NRUN
*>          Number of tests run.
*>
      SUBROUTINE CCHKSY_2X2( NOUT, NFAIL, NRUN )
      IMPLICIT NONE
      INTEGER            NOUT, NFAIL, NRUN
*
*     .. Parameters ..
      INTEGER            NMAX, NB, LWORK
      PARAMETER          ( NMAX = 130, NB = 64, LWORK = NMAX*NB )
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
*     .. Local Scalars ..
      CHARACTER          UPLO
      CHARACTER*2        KIND
      INTEGER            I, J, II, JJ, IC, IFAM, IMETH, ISIZE,
     $                   IUPLO, K, N, NBOLD, ICOL, INFO
      LOGICAL            BAD
      REAL               BIG, SMALL, TOL, ERR, MAG
*     .. Local Arrays ..
      INTEGER            IPIV( NMAX )
      COMPLEX            G( 3, 3 ), A( NMAX, NMAX ), B( NMAX ),
     $                   AP( NMAX*( NMAX+1 ) / 2 ), E( NMAX ),
     $                   WORK( LWORK )
      COMPLEX            PHASE( 3 )
*     .. External Functions ..
      REAL               SLAMCH
      INTEGER            ILAENV
      EXTERNAL           SLAMCH, ILAENV, XLAENV
      EXTERNAL           CSYTRF, CSYTRF_ROOK, CSYTRF_RK
      EXTERNAL           CSYTRS, CSYTRS_ROOK, CSYTRS_3
      EXTERNAL           CHETRF, CHETRF_ROOK, CHETRF_RK
      EXTERNAL           CHETRS, CHETRS_ROOK, CHETRS_3
      EXTERNAL           CSPTRF, CSPTRS, CHPTRF
      EXTERNAL           CHPTRS
*
      BIG = HUGE( ONE )
      SMALL = SCALE( TINY( ONE ), -10 )
      TOL = 128*SLAMCH( 'Epsilon' )
      NBOLD = ILAENV( 1, 'CSYTRF', 'L', NMAX, -1, -1, -1 )
      CALL XLAENV( 1, NB )
*
      DO IFAM = 1, 2
         KIND = 'SY'
         IF( IFAM.EQ.2 ) KIND = 'HE'
         PHASE( 1 ) = ONE
         PHASE( 2 ) = CMPLX( ZERO, ONE )
         PHASE( 3 ) = -ONE
         DO IC = 1, 7
            IF( IC.EQ.1 .AND. ( SMALL.LE.ZERO .OR.
     $          SMALL.GE.TINY( ONE ) ) ) CYCLE
            IF( IFAM.EQ.2 .AND. IC.GE.3 ) CYCLE
            G = ZERO
            IF( IC.EQ.1 ) THEN
*              1 / (4*SMALL) overflows, but the multipliers are bounded.
               G( 1, 1 ) = SMALL
               G( 2, 2 ) = SMALL
               G( 1, 2 ) = 4*SMALL
               G( 1, 3 ) = SMALL
               G( 2, 3 ) = SMALL
               G( 3, 3 ) = ONE
            ELSE IF( IC.EQ.2 ) THEN
*              T*W overflows before division by the off-diagonal pivot.
               G( 1, 1 ) = 0.5E0*BIG
               G( 2, 2 ) = 0.5E0*BIG
               G( 1, 2 ) = BIG
               G( 2, 3 ) = 0.9E0*BIG
            ELSE
*              The intrinsic complex division can overflow internally:
*              (.53125 + .53125*i) / (-.5 - .5*i), both scaled by BIG.
*              Also cover ordinary values and both reciprocal guards.
               MAG = BIG
               IF( IC.EQ.4 ) MAG = ONE
               IF( IC.EQ.5 ) MAG = SQRT( SLAMCH( 'S' ) ) / 16
               IF( IC.EQ.6 ) MAG = 16 / SQRT( SLAMCH( 'S' ) )
               G( 1, 2 ) = CMPLX( MAG / 8, -MAG / 8 )
               G( 1, 3 ) = CMPLX( -MAG / 2, -MAG / 2 )
               G( 2, 3 ) = G( 1, 3 )
               G( 3, 3 ) = G( 1, 2 )
               IF( IC.EQ.7 ) THEN
*                 Moderate pivot, large numerator: a pivot-only guard
*                 would form (-1.2-.4*i)*(.9+.34*i)*BIG and overflow.
                  G( 1, 2 ) = CMPLX( .75E0, -.25E0 )
                  G( 1, 3 ) = CMPLX( .4E0, .2E0 )
                  G( 2, 2 ) = CMPLX( .3E0, -.2E0 )*BIG
                  G( 2, 3 ) = CMPLX( -.7E0, -.3E0 )*BIG
                  G( 3, 3 ) = CMPLX( -.2E0, -.3E0 )*BIG
               END IF
            END IF
            DO J = 1, 3
               DO I = 1, J - 1
                  G( J, I ) = G( I, J )
               END DO
            END DO
            IF( IFAM.EQ.2 ) THEN
*              Unitary diagonal scaling gives Hermitian imaginary pivots.
               DO J = 1, 3
                  DO I = 1, 3
                     G( I, J ) = PHASE( I )*G( I, J )*
     $                           CONJG( PHASE( J ) )
                  END DO
               END DO
            END IF
            DO ISIZE = 1, 2
               N = 3
               IF( ISIZE.EQ.2 ) N = NMAX
               DO IUPLO = 1, 2
                  UPLO = 'L'
                  ICOL = 3
                  IF( IUPLO.EQ.2 ) THEN
                     UPLO = 'U'
                     ICOL = N - 2
                  END IF
                  DO IMETH = 1, 4
                     IF( IC.EQ.7 .AND.
     $                   ( IMETH.EQ.2 .OR. IMETH.EQ.3 ) ) CYCLE
                     A = ZERO
                     E = ZERO
                     DO I = 1, N
                        A( I, I ) = ONE
                     END DO
                     DO J = 1, 3
                        DO I = 1, 3
                           II = I
                           JJ = J
                           IF( IUPLO.EQ.2 ) THEN
                              II = N + 1 - I
                              JJ = N + 1 - J
                           END IF
                           A( II, JJ ) = G( I, J )
                        END DO
                     END DO
                     B( 1:N ) = A( 1:N, ICOL )
                     K = 0
                     DO J = 1, N
                        DO I = 1, N
                           IF( ( UPLO.EQ.'U' .AND. I.LE.J ) .OR.
     $                         ( UPLO.EQ.'L' .AND. I.GE.J ) ) THEN
                              K = K + 1
                              AP( K ) = A( I, J )
                           END IF
                        END DO
                     END DO
                     IF( IFAM.EQ.1 ) THEN
                        IF( IMETH.EQ.1 ) THEN
                           CALL CSYTRF( UPLO, N, A, NMAX, IPIV, WORK,
     $                        LWORK, INFO )
                        ELSE IF( IMETH.EQ.2 ) THEN
                           CALL CSYTRF_ROOK( UPLO, N, A, NMAX, IPIV,
     $                        WORK, LWORK, INFO )
                        ELSE IF( IMETH.EQ.3 ) THEN
                           CALL CSYTRF_RK( UPLO, N, A, NMAX, E, IPIV,
     $                        WORK, LWORK, INFO )
                        ELSE
                           CALL CSPTRF( UPLO, N, AP, IPIV, INFO )
                        END IF
                     ELSE
                        IF( IMETH.EQ.1 ) THEN
                           CALL CHETRF( UPLO, N, A, NMAX, IPIV, WORK,
     $                        LWORK, INFO )
                        ELSE IF( IMETH.EQ.2 ) THEN
                           CALL CHETRF_ROOK( UPLO, N, A, NMAX, IPIV,
     $                        WORK, LWORK, INFO )
                        ELSE IF( IMETH.EQ.3 ) THEN
                           CALL CHETRF_RK( UPLO, N, A, NMAX, E, IPIV,
     $                        WORK, LWORK, INFO )
                        ELSE
                           CALL CHPTRF( UPLO, N, AP, IPIV, INFO )
                        END IF
                     END IF
                     BAD = INFO.NE.0
*                    Confirm a 2-by-2 pivot and reject nonfinite factors
*                    before passing them to the solve routine.
                     IF( UPLO.EQ.'L' ) THEN
                        BAD = BAD .OR. IPIV( 1 ).GE.0
                     ELSE
                        BAD = BAD .OR. IPIV( N ).GE.0
                     END IF
                     DO J = 1, N
                        DO I = 1, N
                           IF( .NOT.( ABS( REAL( A( I, J ) ) ).LE.
     $                          BIG .AND. ABS( AIMAG( A( I, J ) ) ).LE.
     $                          BIG ) ) BAD = .TRUE.
                        END DO
                     END DO
                     DO I = 1, K
                        IF( .NOT.( ABS( REAL( AP( I ) ) ).LE.
     $                          BIG .AND. ABS( AIMAG( AP( I ) ) ).LE.
     $                          BIG ) ) BAD = .TRUE.
                     END DO
                     DO I = 1, N
                        IF( .NOT.( ABS( REAL( E( I ) ) ).LE.
     $                          BIG .AND. ABS( AIMAG( E( I ) ) ).LE.
     $                          BIG ) ) BAD = .TRUE.
                     END DO
                     IF( IC.EQ.7 .AND. .NOT.BAD ) THEN
*                       This matrix is ill-conditioned; check the known
*                       multiplier instead of a forward solve error.
                        II = 3
                        JJ = 1
                        K = 3
                        IF( UPLO.EQ.'U' ) THEN
                           II = N - 2
                           JJ = N
                           K = N*( N-1 ) / 2 + N - 2
                        END IF
                        WORK( 1 ) = A( II, JJ ) / BIG
                        IF( IMETH.EQ.4 ) WORK( 1 ) = AP( K ) / BIG
                        ERR = ABS( WORK( 1 )-
     $                        CMPLX( -.944E0, -.768E0 ) )
                        BAD = .NOT.( ERR.LE.TOL )
                     END IF
                     IF( .NOT.BAD .AND. IC.NE.7 ) THEN
                        IF( IFAM.EQ.1 ) THEN
                           IF( IMETH.EQ.1 ) THEN
                              CALL CSYTRS( UPLO, N, 1, A, NMAX, IPIV, B,
     $                        NMAX, INFO )
                           ELSE IF( IMETH.EQ.2 ) THEN
                              CALL CSYTRS_ROOK( UPLO, N, 1, A, NMAX,
     $                        IPIV, B, NMAX, INFO )
                           ELSE IF( IMETH.EQ.3 ) THEN
                              CALL CSYTRS_3( UPLO, N, 1, A, NMAX, E,
     $                        IPIV, B, NMAX, INFO )
                           ELSE
                              CALL CSPTRS( UPLO, N, 1, AP, IPIV, B,
     $                        NMAX, INFO )
                           END IF
                        ELSE
                           IF( IMETH.EQ.1 ) THEN
                              CALL CHETRS( UPLO, N, 1, A, NMAX, IPIV, B,
     $                        NMAX, INFO )
                           ELSE IF( IMETH.EQ.2 ) THEN
                              CALL CHETRS_ROOK( UPLO, N, 1, A, NMAX,
     $                        IPIV, B, NMAX, INFO )
                           ELSE IF( IMETH.EQ.3 ) THEN
                              CALL CHETRS_3( UPLO, N, 1, A, NMAX, E,
     $                        IPIV, B, NMAX, INFO )
                           ELSE
                              CALL CHPTRS( UPLO, N, 1, AP, IPIV, B,
     $                        NMAX, INFO )
                           END IF
                        END IF
                        BAD = INFO.NE.0
                        B( ICOL ) = B( ICOL ) - ONE
                        ERR = ZERO
                        DO I = 1, N
                           IF( .NOT.( ABS( REAL( B( I ) ) ).LE.
     $                          BIG .AND. ABS( AIMAG( B( I ) ) ).LE.
     $                          BIG ) ) BAD = .TRUE.
                           ERR = MAX( ERR, ABS( B( I ) ) )
                        END DO
                        BAD = BAD .OR. .NOT.( ERR.LE.TOL )
                     END IF
                     NRUN = NRUN + 1
                     IF( BAD ) THEN
                        NFAIL = NFAIL + 1
                        WRITE( NOUT, 9999 ) KIND, UPLO, N, IC,
     $                                     IMETH, INFO
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL XLAENV( 1, NBOLD )
      RETURN
 9999 FORMAT( ' CCHKSY_2X2: ', A2, ' UPLO=', A1, ' N=', I3,
     $        ' case=', I1, ' method=', I1, ' INFO=', I3 )
      END
