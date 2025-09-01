*> \brief \b CCHKPB
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CCHKPB( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
*                          THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
*                          XACT, WORK, RWORK, NOUT )
*
*       .. Scalar Arguments ..
*       LOGICAL            TSTERR
*       INTEGER            NMAX, NN, NNB, NNS, NOUT
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            NBVAL( * ), NSVAL( * ), NVAL( * )
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
*> CCHKPB tests CPBTRF, -TRS, -RFS, and -CON.
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
*> \param[in] NNB
*> \verbatim
*>          NNB is INTEGER
*>          The number of values of NB contained in the vector NBVAL.
*> \endverbatim
*>
*> \param[in] NBVAL
*> \verbatim
*>          NBVAL is INTEGER array, dimension (NNB)
*>          The values of the blocksize NB.
*> \endverbatim
*>
*> \param[in] NNS
*> \verbatim
*>          NNS is INTEGER
*>          The number of values of NRHS contained in the vector NSVAL.
*> \endverbatim
*>
*> \param[in] NSVAL
*> \verbatim
*>          NSVAL is INTEGER array, dimension (NNS)
*>          The values of the number of right hand sides NRHS.
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
*>          A is REAL array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AFAC
*> \verbatim
*>          AFAC is REAL array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AINV
*> \verbatim
*>          AINV is REAL array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is REAL array, dimension (NMAX*NSMAX)
*>          where NSMAX is the largest entry in NSVAL.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is REAL array, dimension (NMAX*NSMAX)
*> \endverbatim
*>
*> \param[out] XACT
*> \verbatim
*>          XACT is REAL array, dimension (NMAX*NSMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension
*>                      (NMAX*max(3,NSMAX))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension
*>                      (max(NMAX,2*NSMAX))
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
      SUBROUTINE CCHKPB( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
     $                   THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
     $                   XACT, WORK, RWORK, NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NMAX, NN, NNB, NNS, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            NBVAL( * ), NSVAL( * ), NVAL( * )
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
      PARAMETER          ( NTYPES = 8, NTESTS = 7 )
      INTEGER            NBW
      PARAMETER          ( NBW = 4 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ZEROT
      CHARACTER          DIST, PACKIT, TYPE, UPLO, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, I1, I2, IKD, IMAT, IN, INB, INFO, IOFF,
     $                   IRHS, IUPLO, IW, IZERO, K, KD, KL, KOFF, KU,
     $                   LDA, LDAB, MODE, N, NB, NERRS, NFAIL, NIMAT,
     $                   NKD, NRHS, NRUN
      REAL               AINVNM, ANORM, CNDNUM, RCOND, RCONDC
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 ), ISEEDY( 4 ), KDVAL( NBW )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Functions ..
      REAL               CLANGE, CLANHB, SGET06
      EXTERNAL           CLANGE, CLANHB, SGET06
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, CCOPY, CERRPO, CGET04,
     $                   CLACPY, CLAIPD, CLARHS, CLASET, CLATB4, CLATMS,
     $                   CPBCON, CPBRFS, CPBT01, CPBT02, CPBT05, CPBTRF,
     $                   CPBTRS, CSWAP, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN
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
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'PB'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
*     Test the error exits
*
      IF( TSTERR )
     $   CALL CERRPO( PATH, NOUT )
      INFOT = 0
      KDVAL( 1 ) = 0
*
*     Do for each value of N in NVAL
*
      DO 90 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
*
*        Set limits on the number of loop iterations.
*
         NKD = MAX( 1, MIN( N, 4 ) )
         NIMAT = NTYPES
         IF( N.EQ.0 )
     $      NIMAT = 1
*
         KDVAL( 2 ) = N + ( N+1 ) / 4
         KDVAL( 3 ) = ( 3*N-1 ) / 4
         KDVAL( 4 ) = ( N+1 ) / 4
*
         DO 80 IKD = 1, NKD
*
*           Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
*           makes it easier to skip redundant values for small values
*           of N.
*
            KD = KDVAL( IKD )
            LDAB = KD + 1
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 70 IUPLO = 1, 2
               KOFF = 1
               IF( IUPLO.EQ.1 ) THEN
                  UPLO = 'U'
                  KOFF = MAX( 1, KD+2-N )
                  PACKIT = 'Q'
               ELSE
                  UPLO = 'L'
                  PACKIT = 'B'
               END IF
*
               DO 60 IMAT = 1, NIMAT
*
*                 Do the tests only if DOTYPE( IMAT ) is true.
*
                  IF( .NOT.DOTYPE( IMAT ) )
     $               GO TO 60
*
*                 Skip types 2, 3, or 4 if the matrix size is too small.
*
                  ZEROT = IMAT.GE.2 .AND. IMAT.LE.4
                  IF( ZEROT .AND. N.LT.IMAT-1 )
     $               GO TO 60
*
                  IF( .NOT.ZEROT .OR. .NOT.DOTYPE( 1 ) ) THEN
*
*                    Set up parameters with CLATB4 and generate a test
*                    matrix with CLATMS.
*
                     CALL CLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM,
     $                            MODE, CNDNUM, DIST )
*
                     SRNAMT = 'CLATMS'
                     CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                            CNDNUM, ANORM, KD, KD, PACKIT,
     $                            A( KOFF ), LDAB, WORK, INFO )
*
*                    Check error code from CLATMS.
*
                     IF( INFO.NE.0 ) THEN
                        CALL ALAERH( PATH, 'CLATMS', INFO, 0, UPLO, N,
     $                               N, KD, KD, -1, IMAT, NFAIL, NERRS,
     $                               NOUT )
                        GO TO 60
                     END IF
                  ELSE IF( IZERO.GT.0 ) THEN
*
*                    Use the same matrix for types 3 and 4 as for type
*                    2 by copying back the zeroed out column,
*
                     IW = 2*LDA + 1
                     IF( IUPLO.EQ.1 ) THEN
                        IOFF = ( IZERO-1 )*LDAB + KD + 1
                        CALL CCOPY( IZERO-I1, WORK( IW ), 1,
     $                              A( IOFF-IZERO+I1 ), 1 )
                        IW = IW + IZERO - I1
                        CALL CCOPY( I2-IZERO+1, WORK( IW ), 1,
     $                              A( IOFF ), MAX( LDAB-1, 1 ) )
                     ELSE
                        IOFF = ( I1-1 )*LDAB + 1
                        CALL CCOPY( IZERO-I1, WORK( IW ), 1,
     $                              A( IOFF+IZERO-I1 ),
     $                              MAX( LDAB-1, 1 ) )
                        IOFF = ( IZERO-1 )*LDAB + 1
                        IW = IW + IZERO - I1
                        CALL CCOPY( I2-IZERO+1, WORK( IW ), 1,
     $                              A( IOFF ), 1 )
                     END IF
                  END IF
*
*                 For types 2-4, zero one row and column of the matrix
*                 to test that INFO is returned correctly.
*
                  IZERO = 0
                  IF( ZEROT ) THEN
                     IF( IMAT.EQ.2 ) THEN
                        IZERO = 1
                     ELSE IF( IMAT.EQ.3 ) THEN
                        IZERO = N
                     ELSE
                        IZERO = N / 2 + 1
                     END IF
*
*                    Save the zeroed out row and column in WORK(*,3)
*
                     IW = 2*LDA
                     DO 20 I = 1, MIN( 2*KD+1, N )
                        WORK( IW+I ) = ZERO
   20                CONTINUE
                     IW = IW + 1
                     I1 = MAX( IZERO-KD, 1 )
                     I2 = MIN( IZERO+KD, N )
*
                     IF( IUPLO.EQ.1 ) THEN
                        IOFF = ( IZERO-1 )*LDAB + KD + 1
                        CALL CSWAP( IZERO-I1, A( IOFF-IZERO+I1 ), 1,
     $                              WORK( IW ), 1 )
                        IW = IW + IZERO - I1
                        CALL CSWAP( I2-IZERO+1, A( IOFF ),
     $                              MAX( LDAB-1, 1 ), WORK( IW ), 1 )
                     ELSE
                        IOFF = ( I1-1 )*LDAB + 1
                        CALL CSWAP( IZERO-I1, A( IOFF+IZERO-I1 ),
     $                              MAX( LDAB-1, 1 ), WORK( IW ), 1 )
                        IOFF = ( IZERO-1 )*LDAB + 1
                        IW = IW + IZERO - I1
                        CALL CSWAP( I2-IZERO+1, A( IOFF ), 1,
     $                              WORK( IW ), 1 )
                     END IF
                  END IF
*
*                 Set the imaginary part of the diagonals.
*
                  IF( IUPLO.EQ.1 ) THEN
                     CALL CLAIPD( N, A( KD+1 ), LDAB, 0 )
                  ELSE
                     CALL CLAIPD( N, A( 1 ), LDAB, 0 )
                  END IF
*
*                 Do for each value of NB in NBVAL
*
                  DO 50 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
*
*                    Compute the L*L' or U'*U factorization of the band
*                    matrix.
*
                     CALL CLACPY( 'Full', KD+1, N, A, LDAB, AFAC, LDAB )
                     SRNAMT = 'CPBTRF'
                     CALL CPBTRF( UPLO, N, KD, AFAC, LDAB, INFO )
*
*                    Check error code from CPBTRF.
*
                     IF( INFO.NE.IZERO ) THEN
                        CALL ALAERH( PATH, 'CPBTRF', INFO, IZERO, UPLO,
     $                               N, N, KD, KD, NB, IMAT, NFAIL,
     $                               NERRS, NOUT )
                        GO TO 50
                     END IF
*
*                    Skip the tests if INFO is not 0.
*
                     IF( INFO.NE.0 )
     $                  GO TO 50
*
*+    TEST 1
*                    Reconstruct matrix from factors and compute
*                    residual.
*
                     CALL CLACPY( 'Full', KD+1, N, AFAC, LDAB, AINV,
     $                            LDAB )
                     CALL CPBT01( UPLO, N, KD, A, LDAB, AINV, LDAB,
     $                            RWORK, RESULT( 1 ) )
*
*                    Print the test ratio if it is .GE. THRESH.
*
                     IF( RESULT( 1 ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 )UPLO, N, KD, NB, IMAT,
     $                     1, RESULT( 1 )
                        NFAIL = NFAIL + 1
                     END IF
                     NRUN = NRUN + 1
*
*                    Only do other tests if this is the first blocksize.
*
                     IF( INB.GT.1 )
     $                  GO TO 50
*
*                    Form the inverse of A so we can get a good estimate
*                    of RCONDC = 1/(norm(A) * norm(inv(A))).
*
                     CALL CLASET( 'Full', N, N, CMPLX( ZERO ),
     $                            CMPLX( ONE ), AINV, LDA )
                     SRNAMT = 'CPBTRS'
                     CALL CPBTRS( UPLO, N, KD, N, AFAC, LDAB, AINV, LDA,
     $                            INFO )
*
*                    Compute RCONDC = 1/(norm(A) * norm(inv(A))).
*
                     ANORM = CLANHB( '1', UPLO, N, KD, A, LDAB, RWORK )
                     AINVNM = CLANGE( '1', N, N, AINV, LDA, RWORK )
                     IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                        RCONDC = ONE
                     ELSE
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     END IF
*
                     DO 40 IRHS = 1, NNS
                        NRHS = NSVAL( IRHS )
*
*+    TEST 2
*                    Solve and compute residual for A * X = B.
*
                        SRNAMT = 'CLARHS'
                        CALL CLARHS( PATH, XTYPE, UPLO, ' ', N, N, KD,
     $                               KD, NRHS, A, LDAB, XACT, LDA, B,
     $                               LDA, ISEED, INFO )
                        CALL CLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                        SRNAMT = 'CPBTRS'
                        CALL CPBTRS( UPLO, N, KD, NRHS, AFAC, LDAB, X,
     $                               LDA, INFO )
*
*                    Check error code from CPBTRS.
*
                        IF( INFO.NE.0 )
     $                     CALL ALAERH( PATH, 'CPBTRS', INFO, 0, UPLO,
     $                                  N, N, KD, KD, NRHS, IMAT, NFAIL,
     $                                  NERRS, NOUT )
*
                        CALL CLACPY( 'Full', N, NRHS, B, LDA, WORK,
     $                               LDA )
                        CALL CPBT02( UPLO, N, KD, NRHS, A, LDAB, X, LDA,
     $                               WORK, LDA, RWORK, RESULT( 2 ) )
*
*+    TEST 3
*                    Check solution from generated exact solution.
*
                        CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                               RESULT( 3 ) )
*
*+    TESTS 4, 5, and 6
*                    Use iterative refinement to improve the solution.
*
                        SRNAMT = 'CPBRFS'
                        CALL CPBRFS( UPLO, N, KD, NRHS, A, LDAB, AFAC,
     $                               LDAB, B, LDA, X, LDA, RWORK,
     $                               RWORK( NRHS+1 ), WORK,
     $                               RWORK( 2*NRHS+1 ), INFO )
*
*                    Check error code from CPBRFS.
*
                        IF( INFO.NE.0 )
     $                     CALL ALAERH( PATH, 'CPBRFS', INFO, 0, UPLO,
     $                                  N, N, KD, KD, NRHS, IMAT, NFAIL,
     $                                  NERRS, NOUT )
*
                        CALL CGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                               RESULT( 4 ) )
                        CALL CPBT05( UPLO, N, KD, NRHS, A, LDAB, B, LDA,
     $                               X, LDA, XACT, LDA, RWORK,
     $                               RWORK( NRHS+1 ), RESULT( 5 ) )
*
*                       Print information about the tests that did not
*                       pass the threshold.
*
                        DO 30 K = 2, 6
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                           CALL ALAHD( NOUT, PATH )
                              WRITE( NOUT, FMT = 9998 )UPLO, N, KD,
     $                           NRHS, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           END IF
   30                   CONTINUE
                        NRUN = NRUN + 5
   40                CONTINUE
*
*+    TEST 7
*                    Get an estimate of RCOND = 1/CNDNUM.
*
                     SRNAMT = 'CPBCON'
                     CALL CPBCON( UPLO, N, KD, AFAC, LDAB, ANORM, RCOND,
     $                            WORK, RWORK, INFO )
*
*                    Check error code from CPBCON.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'CPBCON', INFO, 0, UPLO, N,
     $                               N, KD, KD, -1, IMAT, NFAIL, NERRS,
     $                               NOUT )
*
                     RESULT( 7 ) = SGET06( RCOND, RCONDC )
*
*                    Print the test ratio if it is .GE. THRESH.
*
                     IF( RESULT( 7 ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9997 )UPLO, N, KD, IMAT, 7,
     $                     RESULT( 7 )
                        NFAIL = NFAIL + 1
                     END IF
                     NRUN = NRUN + 1
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ', NB=', I4,
     $      ', type ', I2, ', test ', I2, ', ratio= ', G12.5 )
 9998 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ', NRHS=', I3,
     $      ', type ', I2, ', test(', I2, ') = ', G12.5 )
 9997 FORMAT( ' UPLO=''', A1, ''', N=', I5, ', KD=', I5, ',', 10X,
     $      ' type ', I2, ', test(', I2, ') = ', G12.5 )
      RETURN
*
*     End of CCHKPB
*
      END
