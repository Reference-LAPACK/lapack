*> \brief \b CDRVSY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CDRVSY( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
*                          A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK,
*                          NOUT )
*
*       .. Scalar Arguments ..
*       LOGICAL            TSTERR
*       INTEGER            NMAX, NN, NOUT, NRHS
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            IWORK( * ), NVAL( * )
*       REAL               RWORK( * )
*       COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ),
*      $                   WORK( * ), X( * ), XACT( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CDRVSY tests the driver routines CSYSV and -SVX.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DOTYPE
*> \verbatim
*>          DOTYPE is LOGICAL array, dimension (NTYPES)
*>          The matrix types to be used for testing.  Matrices of type j
*>          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
*>          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
*> \endverbatim
*>
*> \param[in] NN
*> \verbatim
*>          NN is INTEGER
*>          The number of values of N contained in the vector NVAL.
*> \endverbatim
*>
*> \param[in] NVAL
*> \verbatim
*>          NVAL is INTEGER array, dimension (NN)
*>          The values of the matrix dimension N.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand side vectors to be generated for
*>          each linear system.
*> \endverbatim
*>
*> \param[in] THRESH
*> \verbatim
*>          THRESH is REAL
*>          The threshold value for the test ratios.  A result is
*>          included in the output file if RESULT >= THRESH.  To have
*>          every test ratio printed, use THRESH = 0.
*> \endverbatim
*>
*> \param[in] TSTERR
*> \verbatim
*>          TSTERR is LOGICAL
*>          Flag that indicates whether error exits are to be tested.
*> \endverbatim
*>
*> \param[in] NMAX
*> \verbatim
*>          NMAX is INTEGER
*>          The maximum value permitted for N, used in dimensioning the
*>          work arrays.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AFAC
*> \verbatim
*>          AFAC is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AINV
*> \verbatim
*>          AINV is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is COMPLEX array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] XACT
*> \verbatim
*>          XACT is COMPLEX array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (NMAX*max(2,NRHS))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (NMAX+2*NRHS)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (NMAX)
*> \endverbatim
*>
*> \param[in] NOUT
*> \verbatim
*>          NOUT is INTEGER
*>          The unit number for output.
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
*> \ingroup complex_lin
*
*  =====================================================================
      SUBROUTINE CDRVSY( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
     $                   A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK,
     $                   NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NMAX, NN, NOUT, NRHS
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            IWORK( * ), NVAL( * )
      REAL               RWORK( * )
      COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ),
     $                   WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NTYPES, NTESTS
      PARAMETER          ( NTYPES = 11, NTESTS = 6 )
      INTEGER            NFACT
      PARAMETER          ( NFACT = 2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ZEROT
      CHARACTER          DIST, FACT, TYPE, UPLO, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO,
     $                   IZERO, J, K, K1, KL, KU, LDA, LWORK, MODE, N,
     $                   NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT
      REAL               AINVNM, ANORM, CNDNUM, RCOND, RCONDC
*     ..
*     .. Local Arrays ..
      CHARACTER          FACTS( NFACT ), UPLOS( 2 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Functions ..
      REAL               CLANSY, SGET06
      EXTERNAL           CLANSY, SGET06
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALADHD, ALAERH, ALASVM, CERRVX, CGET04, CLACPY,
     $                   CLARHS, CLASET, CLATB4, CLATMS, CLATSY, CPOT05,
     $                   CSYSV, CSYSVX, CSYT01, CSYT02, CSYTRF, CSYTRI2,
     $                   XLAENV
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'SY'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      LWORK = MAX( 2*NMAX, NMAX*NRHS )
*
*     Test the error exits
*
      IF( TSTERR )
     $   CALL CERRVX( PATH, NOUT )
      INFOT = 0
*
*     Set the block size and minimum block size for testing.
*
      NB = 1
      NBMIN = 2
      CALL XLAENV( 1, NB )
      CALL XLAENV( 2, NBMIN )
*
*     Do for each value of N in NVAL
*
      DO 180 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         DO 170 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 170
*
*           Skip types 3, 4, 5, or 6 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            IF( ZEROT .AND. N.LT.IMAT-2 )
     $         GO TO 170
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 160 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
*
               IF( IMAT.NE.NTYPES ) THEN
*
*                 Set up parameters with CLATB4 and generate a test
*                 matrix with CLATMS.
*
                  CALL CLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM,
     $                         MODE, CNDNUM, DIST )
*
                  SRNAMT = 'CLATMS'
                  CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                         CNDNUM, ANORM, KL, KU, UPLO, A, LDA,
     $                         WORK, INFO )
*
*                 Check error code from CLATMS.
*
                  IF( INFO.NE.0 ) THEN
                     CALL ALAERH( PATH, 'CLATMS', INFO, 0, UPLO, N, N,
     $                            -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 160
                  END IF
*
*                 For types 3-6, zero one or more rows and columns of
*                 the matrix to test that INFO is returned correctly.
*
                  IF( ZEROT ) THEN
                     IF( IMAT.EQ.3 ) THEN
                        IZERO = 1
                     ELSE IF( IMAT.EQ.4 ) THEN
                        IZERO = N
                     ELSE
                        IZERO = N / 2 + 1
                     END IF
*
                     IF( IMAT.LT.6 ) THEN
*
*                       Set row and column IZERO to zero.
*
                        IF( IUPLO.EQ.1 ) THEN
                           IOFF = ( IZERO-1 )*LDA
                           DO 20 I = 1, IZERO - 1
                              A( IOFF+I ) = ZERO
   20                      CONTINUE
                           IOFF = IOFF + IZERO
                           DO 30 I = IZERO, N
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   30                      CONTINUE
                        ELSE
                           IOFF = IZERO
                           DO 40 I = 1, IZERO - 1
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   40                      CONTINUE
                           IOFF = IOFF - IZERO
                           DO 50 I = IZERO, N
                              A( IOFF+I ) = ZERO
   50                      CONTINUE
                        END IF
                     ELSE
                        IF( IUPLO.EQ.1 ) THEN
*
*                          Set the first IZERO rows to zero.
*
                           IOFF = 0
                           DO 70 J = 1, N
                              I2 = MIN( J, IZERO )
                              DO 60 I = 1, I2
                                 A( IOFF+I ) = ZERO
   60                         CONTINUE
                              IOFF = IOFF + LDA
   70                      CONTINUE
                        ELSE
*
*                          Set the last IZERO rows to zero.
*
                           IOFF = 0
                           DO 90 J = 1, N
                              I1 = MAX( J, IZERO )
                              DO 80 I = I1, N
                                 A( IOFF+I ) = ZERO
   80                         CONTINUE
                              IOFF = IOFF + LDA
   90                      CONTINUE
                        END IF
                     END IF
                  ELSE
                     IZERO = 0
                  END IF
               ELSE
*
*                 IMAT = NTYPES:  Use a special block diagonal matrix to
*                 test alternate code for the 2-by-2 blocks.
*
                  CALL CLATSY( UPLO, N, A, LDA, ISEED )
               END IF
*
               DO 150 IFACT = 1, NFACT
*
*                 Do first for FACT = 'F', then for other values.
*
                  FACT = FACTS( IFACT )
*
*                 Compute the condition number for comparison with
*                 the value returned by CSYSVX.
*
                  IF( ZEROT ) THEN
                     IF( IFACT.EQ.1 )
     $                  GO TO 150
                     RCONDC = ZERO
*
                  ELSE IF( IFACT.EQ.1 ) THEN
*
*                    Compute the 1-norm of A.
*
                     ANORM = CLANSY( '1', UPLO, N, A, LDA, RWORK )
*
*                    Factor the matrix A.
*
                     CALL CLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL CSYTRF( UPLO, N, AFAC, LDA, IWORK, WORK,
     $                            LWORK, INFO )
*
*                    Compute inv(A) and take its norm.
*
                     CALL CLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
                     LWORK = (N+NB+1)*(NB+3)
                     CALL CSYTRI2( UPLO, N, AINV, LDA, IWORK, WORK,
     $                            LWORK, INFO )
                     AINVNM = CLANSY( '1', UPLO, N, AINV, LDA, RWORK )
*
*                    Compute the 1-norm condition number of A.
*
                     IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                        RCONDC = ONE
                     ELSE
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     END IF
                  END IF
*
*                 Form an exact solution and set the right hand side.
*
                  SRNAMT = 'CLARHS'
                  CALL CLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU,
     $                         NRHS, A, LDA, XACT, LDA, B, LDA, ISEED,
     $                         INFO )
                  XTYPE = 'C'
*
*                 --- Test CSYSV  ---
*
                  IF( IFACT.EQ.2 ) THEN
                     CALL CLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
*                    Factor the matrix and solve the system using CSYSV.
*
                     SRNAMT = 'CSYSV '
                     CALL CSYSV( UPLO, N, NRHS, AFAC, LDA, IWORK, X,
     $                           LDA, WORK, LWORK, INFO )
*
*                    Adjust the expected value of INFO to account for
*                    pivoting.
*
                     K = IZERO
                     IF( K.GT.0 ) THEN
  100                   CONTINUE
                        IF( IWORK( K ).LT.0 ) THEN
                           IF( IWORK( K ).NE.-K ) THEN
                              K = -IWORK( K )
                              GO TO 100
                           END IF
                        ELSE IF( IWORK( K ).NE.K ) THEN
                           K = IWORK( K )
                           GO TO 100
                        END IF
                     END IF
*
*                    Check error code from CSYSV .
*
                     IF( INFO.NE.K ) THEN
                        CALL ALAERH( PATH, 'CSYSV ', INFO, K, UPLO, N,
     $                               N, -1, -1, NRHS, IMAT, NFAIL,
     $                               NERRS, NOUT )
                        GO TO 120
                     ELSE IF( INFO.NE.0 ) THEN
                        GO TO 120
                     END IF
*
*                    Reconstruct matrix from factors and compute
*                    residual.
*
                     CALL CSYT01( UPLO, N, A, LDA, AFAC, LDA, IWORK,
     $                            AINV, LDA, RWORK, RESULT( 1 ) )
*
*                    Compute residual of the computed solution.
*
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL CSYT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK,
     $                            LDA, RWORK, RESULT( 2 ) )
*
*                    Check solution from generated exact solution.
*
                     CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                            RESULT( 3 ) )
                     NT = 3
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     DO 110 K = 1, NT
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALADHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9999 )'CSYSV ', UPLO, N,
     $                        IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
  110                CONTINUE
                     NRUN = NRUN + NT
  120                CONTINUE
                  END IF
*
*                 --- Test CSYSVX ---
*
                  IF( IFACT.EQ.2 )
     $               CALL CLASET( UPLO, N, N, CMPLX( ZERO ),
     $                            CMPLX( ZERO ), AFAC, LDA )
                  CALL CLASET( 'Full', N, NRHS, CMPLX( ZERO ),
     $                         CMPLX( ZERO ), X, LDA )
*
*                 Solve the system and compute the condition number and
*                 error bounds using CSYSVX.
*
                  SRNAMT = 'CSYSVX'
                  CALL CSYSVX( FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA,
     $                         IWORK, B, LDA, X, LDA, RCOND, RWORK,
     $                         RWORK( NRHS+1 ), WORK, LWORK,
     $                         RWORK( 2*NRHS+1 ), INFO )
*
*                 Adjust the expected value of INFO to account for
*                 pivoting.
*
                  K = IZERO
                  IF( K.GT.0 ) THEN
  130                CONTINUE
                     IF( IWORK( K ).LT.0 ) THEN
                        IF( IWORK( K ).NE.-K ) THEN
                           K = -IWORK( K )
                           GO TO 130
                        END IF
                     ELSE IF( IWORK( K ).NE.K ) THEN
                        K = IWORK( K )
                        GO TO 130
                     END IF
                  END IF
*
*                 Check the error code from CSYSVX.
*
                  IF( INFO.NE.K ) THEN
                     CALL ALAERH( PATH, 'CSYSVX', INFO, K, FACT // UPLO,
     $                            N, N, -1, -1, NRHS, IMAT, NFAIL,
     $                            NERRS, NOUT )
                     GO TO 150
                  END IF
*
                  IF( INFO.EQ.0 ) THEN
                     IF( IFACT.GE.2 ) THEN
*
*                       Reconstruct matrix from factors and compute
*                       residual.
*
                        CALL CSYT01( UPLO, N, A, LDA, AFAC, LDA, IWORK,
     $                               AINV, LDA, RWORK( 2*NRHS+1 ),
     $                               RESULT( 1 ) )
                        K1 = 1
                     ELSE
                        K1 = 2
                     END IF
*
*                    Compute residual of the computed solution.
*
                     CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL CSYT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK,
     $                            LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )
*
*                    Check solution from generated exact solution.
*
                     CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                            RESULT( 3 ) )
*
*                    Check the error bounds from iterative refinement.
*
                     CALL CPOT05( UPLO, N, NRHS, A, LDA, B, LDA, X, LDA,
     $                            XACT, LDA, RWORK, RWORK( NRHS+1 ),
     $                            RESULT( 4 ) )
                  ELSE
                     K1 = 6
                  END IF
*
*                 Compare RCOND from CSYSVX with the computed value
*                 in RCONDC.
*
                  RESULT( 6 ) = SGET06( RCOND, RCONDC )
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 140 K = K1, 6
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALADHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9998 )'CSYSVX', FACT, UPLO,
     $                     N, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  140             CONTINUE
                  NRUN = NRUN + 7 - K1
*
  150          CONTINUE
*
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2,
     $      ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N =', I5,
     $      ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
      RETURN
*
*     End of CDRVSY
*
      END
