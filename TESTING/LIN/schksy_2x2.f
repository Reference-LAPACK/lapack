*> \brief \b SCHKSY_2X2 checks extreme-scale 2-by-2 pivots.
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
      SUBROUTINE SCHKSY_2X2( NOUT, NFAIL, NRUN )
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
      REAL               BIG, SMALL, TOL, ERR
*     .. Local Arrays ..
      INTEGER            IPIV( NMAX )
      REAL               G( 3, 3 ), A( NMAX, NMAX ), B( NMAX ),
     $                   AP( NMAX*( NMAX+1 ) / 2 ), E( NMAX ),
     $                   WORK( LWORK )
*     .. External Functions ..
      REAL               SLAMCH
      INTEGER            ILAENV
      EXTERNAL           SLAMCH, ILAENV, XLAENV
      EXTERNAL           SSYTRF, SSYTRF_ROOK, SSYTRF_RK
      EXTERNAL           SSYTRS, SSYTRS_ROOK, SSYTRS_3
      EXTERNAL           SSPTRF, SSPTRS
*
      BIG = HUGE( ONE )
      SMALL = SCALE( TINY( ONE ), -10 )
      TOL = 128*SLAMCH( 'Epsilon' )
      NBOLD = ILAENV( 1, 'SSYTRF', 'L', NMAX, -1, -1, -1 )
      CALL XLAENV( 1, NB )
*
      DO IFAM = 1, 1
         KIND = 'SY'
         DO IC = 1, 2
            IF( IC.EQ.1 .AND. ( SMALL.LE.ZERO .OR.
     $          SMALL.GE.TINY( ONE ) ) ) CYCLE
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
            END IF
            DO J = 1, 3
               DO I = 1, J - 1
                  G( J, I ) = G( I, J )
               END DO
            END DO
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
                     IF( IMETH.EQ.1 ) THEN
                        CALL SSYTRF( UPLO, N, A, NMAX, IPIV, WORK,
     $                        LWORK, INFO )
                     ELSE IF( IMETH.EQ.2 ) THEN
                        CALL SSYTRF_ROOK( UPLO, N, A, NMAX, IPIV, WORK,
     $                        LWORK, INFO )
                     ELSE IF( IMETH.EQ.3 ) THEN
                        CALL SSYTRF_RK( UPLO, N, A, NMAX, E, IPIV, WORK,
     $                        LWORK, INFO )
                     ELSE
                        CALL SSPTRF( UPLO, N, AP, IPIV, INFO )
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
                           IF( .NOT.( ABS( A( I, J ) ).LE.BIG ) )
     $                          BAD = .TRUE.
                        END DO
                     END DO
                     DO I = 1, K
                        IF( .NOT.( ABS( AP( I ) ).LE.BIG ) )
     $                          BAD = .TRUE.
                     END DO
                     DO I = 1, N
                        IF( .NOT.( ABS( E( I ) ).LE.BIG ) )
     $                          BAD = .TRUE.
                     END DO
                     IF( .NOT.BAD ) THEN
                        IF( IMETH.EQ.1 ) THEN
                           CALL SSYTRS( UPLO, N, 1, A, NMAX, IPIV, B,
     $                        NMAX, INFO )
                        ELSE IF( IMETH.EQ.2 ) THEN
                           CALL SSYTRS_ROOK( UPLO, N, 1, A, NMAX, IPIV,
     $                        B, NMAX, INFO )
                        ELSE IF( IMETH.EQ.3 ) THEN
                           CALL SSYTRS_3( UPLO, N, 1, A, NMAX, E, IPIV,
     $                        B, NMAX, INFO )
                        ELSE
                           CALL SSPTRS( UPLO, N, 1, AP, IPIV, B, NMAX,
     $                        INFO )
                        END IF
                        BAD = INFO.NE.0
                        B( ICOL ) = B( ICOL ) - ONE
                        ERR = ZERO
                        DO I = 1, N
                           IF( .NOT.( ABS( B( I ) ).LE.BIG ) )
     $                          BAD = .TRUE.
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
 9999 FORMAT( ' SCHKSY_2X2: ', A2, ' UPLO=', A1, ' N=', I3,
     $        ' case=', I1, ' method=', I1, ' INFO=', I3 )
      END
