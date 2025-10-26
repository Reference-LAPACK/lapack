*> \brief \b SCHKKY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SCHKKY( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
*                          THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
*                          XACT, WORK, RWORK, IWORK, NOUT )
*
*       .. Scalar Arguments ..
*       LOGICAL            TSTERR
*       INTEGER            NMAX, NN, NNB, NNS, NOUT
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
*       REAL               A( * ), AFAC( * ), AINV( * ), B( * ),
*      $                   RWORK( * ), WORK( * ), X( * ), XACT( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SCHKKY tests SKYTRF, -TRI2, -TRS, -TRS2, -RFS, and -CON.
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
*>          WORK is REAL array, dimension (NMAX*max(3,NSMAX))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (max(NMAX,2*NSMAX))
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (2*NMAX)
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
*> \ingroup single_lin
*
*  =====================================================================
      SUBROUTINE SCHKKY( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
     $                   THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
     $                   XACT, WORK, RWORK, IWORK, NOUT )
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
      INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
      REAL               A( * ), AFAC( * ), AINV( * ), B( * ),
     $                   RWORK( * ), WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 10 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TRFCON, ZEROT, LSAME
      CHARACTER          DIST, TYPE, UPLO, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS,
     $                   IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE,
     $                   N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT
      REAL               ANORM, CNDNUM, RCOND, RCONDC
*     ..
*     .. Local Arrays ..
      CHARACTER          UPLOS( 2 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Functions ..
      REAL               SGET06, SLANKY
      EXTERNAL           SGET06, SLANKY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, SERRSY, SGET04, SLACPY,
     $                   SLARHS, SLATB4, SLATMS, SPOT08, SPOT07, SPOT09,
     $                   SKYCON, SKYCONV, SKYRFS, SKYT01, SKYTRF,
     $                   SKYTRI2, SKYTRS, SKYTRS2, LSAME, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
      DATA               UPLOS / 'U', 'L' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'KY'
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
     $   CALL SERRKY( PATH, NOUT )
      INFOT = 0
*
*     Set the minimum block size for which the block routine should
*     be used, which will be later returned by ILAENV
*
      CALL XLAENV( 2, 2 )
*
*     Do for each value of N in NVAL
*
      DO 180 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 2
*
         IZERO = 0
*
*        Do for each value of matrix type IMAT, except IMAT.EQ.1
*
         DO 170 IMAT = 2, NIMAT
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
            ZEROT = MOD(N,2).NE.0
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 160 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
*
*              Begin generate the test matrix A.
*
*              Set up parameters with SLATB4 for the matrix generator
*              based on the type of matrix to be generated.
*
               CALL SLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE,
     $                      CNDNUM, DIST )
*
*              Generate a matrix with SLATMS.
*
               SRNAMT = 'SLATMS'
               CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                      CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK,
     $                      INFO )
*
*              Check error code from SLATMS and handle error.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'SLATMS', INFO, 0, UPLO, N, N, -1,
     $                         -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
*                 Skip all tests for this generated matrix
*
                  GO TO 160
               END IF
*
*              For matrix types 3-6, zero one or more rows and
*              columns of the matrix to test that INFO is returned
*              correctly.
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
*                    Set row and column IZERO to zero.
*
                     IF( IUPLO.EQ.1 ) THEN
                        IOFF = ( IZERO-1 )*LDA
                        DO 20 I = 1, IZERO - 1
                           A( IOFF+I ) = ZERO
   20                   CONTINUE
                        IOFF = IOFF + IZERO
                        DO 30 I = IZERO, N
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
   30                   CONTINUE
                     ELSE
                        IOFF = IZERO
                        DO 40 I = 1, IZERO - 1
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
   40                   CONTINUE
                        IOFF = IOFF - IZERO
                        DO 50 I = IZERO, N
                           A( IOFF+I ) = ZERO
   50                   CONTINUE
                     END IF
                  ELSE
                     IF( IUPLO.EQ.1 ) THEN
*
*                       Set the first IZERO rows and columns to zero.
*
                        IOFF = 0
                        DO 70 J = 1, N
                           I2 = MIN( J, IZERO )
                           DO 60 I = 1, I2
                              A( IOFF+I ) = ZERO
   60                      CONTINUE
                           IOFF = IOFF + LDA
   70                   CONTINUE
                     ELSE
*
*                       Set the last IZERO rows and columns to zero.
*
                        IOFF = 0
                        DO 90 J = 1, N
                           I1 = MAX( J, IZERO )
                           DO 80 I = I1, N
                              A( IOFF+I ) = ZERO
   80                      CONTINUE
                           IOFF = IOFF + LDA
   90                   CONTINUE
                     END IF
                  END IF
               ELSE
                  IZERO = 0
               END IF
*
*              End generate the test matrix A.
*
*
*              Do for each value of NB in NBVAL
*
               DO 150 INB = 1, NNB
*
*                 Set the optimal blocksize, which will be later
*                 returned by ILAENV.
*
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
*
*                 Copy the test matrix A into matrix AFAC which
*                 will be factorized in place. This is needed to
*                 preserve the test matrix A for subsequent tests.
*
                  CALL SLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
*
*                 Compute the L*D*L**T or U*D*U**T factorization of the
*                 matrix. IWORK stores details of the interchanges and
*                 the block structure of D. AINV is a work array for
*                 block factorization, LWORK is the length of AINV.
*
                  LWORK = MAX( 2, NB )*LDA
                  SRNAMT = 'SKYTRF'
                  CALL SKYTRF( UPLO, N, AFAC, LDA, IWORK, AINV, LWORK,
     $                         INFO )
*
*                 Adjust the expected value of INFO to account for
*                 pivoting.
*
                  K = IZERO
                  IF (N.EQ.3 .AND. IMAT.GE.7 .AND. IMAT.LE.10) THEN
                     IF (LSAME( UPLO, 'U' )) THEN
                        K = N
                     ELSEIF (LSAME( UPLO, 'L' )) THEN
                        K = 1
                     END IF
                  ELSEIF (N.EQ.5 .AND. IMAT.GE.6 .AND. IMAT.LE.10)
     $            THEN
                     K = (N + 1) / 2
                  ELSEIF (MOD(N,2).NE.0 .AND. LSAME( UPLO, 'U' ))
     $            THEN
                     K = 1
                  ELSEIF (MOD(N,2).NE.0 .AND. LSAME( UPLO, 'L' ))
     $            THEN
                     K = N
                  ELSEIF( K.GT.0 ) THEN
  100                CONTINUE
                     IF(LSAME( UPLO, 'U' )) THEN
                        IF(MOD(N-K+1,2).NE.0 .AND. IWORK(K).LT.0)
     $                  THEN
                           K = -IWORK( K )
                           GO TO 100
                        ELSEIF(MOD(N-K+1,2).EQ.0 .AND.
     $                  IWORK(K+1).GT.0) THEN
                           K = IWORK( K+1 )
                           GO TO 100
                        ELSEIF(MOD(N-K+1,2).EQ.0 .AND.
     $                  IWORK(K+1).EQ.0) THEN
                           K = K+1
                        END IF
                     ELSE IF(LSAME( UPLO, 'L' )) THEN
                        IF(MOD(K,2).NE.0 .AND. IWORK(K).LT.0)
     $                  THEN
                           K = -IWORK( K )
                           GO TO 100
                        ELSEIF(MOD(K,2).EQ.0 .AND. IWORK(K-1).GT.0)
     $                  THEN
                           K = IWORK( K-1 )
                           GO TO 100
                        ELSEIF(MOD(K,2).EQ.0 .AND. IWORK(K-1).EQ.0)
     $                  THEN
                           K = K-1  
                        END IF
                     END IF
                  END IF
*
*                 Check error code from SKYTRF and handle error.
*
                  IF( INFO.NE.K )
     $               CALL ALAERH( PATH, 'SKYTRF', INFO, K, UPLO, N, N,
     $                            -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )
*
*                 Set the condition estimate flag if the INFO is not 0.
*
                  IF( INFO.NE.0 ) THEN
                     TRFCON = .TRUE.
                  ELSE
                     TRFCON = .FALSE.
                  END IF
*
*+    TEST 1
*                 Reconstruct matrix from factors and compute residual.
*
                  CALL SKYT01( UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV,
     $                         LDA, RWORK, RESULT( 1 ) )
                  NT = 1
*
*+    TEST 2
*                 Form the inverse and compute the residual,
*                 if the factorization was competed without INFO > 0
*                 (i.e. there is no zero rows and columns).
*
                  IF( INB.EQ.1 .AND. .NOT.TRFCON ) THEN
                     CALL SLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
                     SRNAMT = 'SKYTRI2'
                     LWORK = (N+NB+1)*(NB+3)
                     CALL SKYTRI2( UPLO, N, AINV, LDA, IWORK, WORK,
     $                            LWORK, INFO )
*
*                    Check error code from SKYTRI2 and handle error.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'SKYTRI2', INFO, -1, UPLO, N,
     $                               N, -1, -1, -1, IMAT, NFAIL, NERRS,
     $                               NOUT )
*
*                    Compute the residual for a skew-symmetric matrix times
*                    its inverse.
*
                     CALL SPOT08( UPLO, N, A, LDA, AINV, LDA, WORK, LDA,
     $                            RWORK, RCONDC, RESULT( 2 ) )
                     NT = 2
                  END IF
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 110 K = 1, NT
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K,
     $                     RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  110             CONTINUE
                  NRUN = NRUN + NT
*
*                 Skip the other tests if this is not the first block
*                 size.
*
                  IF( INB.GT.1 )
     $               GO TO 150
*
*                 Do only the condition estimate if INFO is not 0.
*
                  IF( TRFCON ) THEN
                     RCONDC = ZERO
                     GO TO 140
                  END IF
*
*                 Do for each value of NRHS in NSVAL.
*
                  DO 130 IRHS = 1, NNS
                     NRHS = NSVAL( IRHS )
*
*+    TEST 3 (Using DSYTRS)
*                 Solve and compute residual for  A * X = B.
*
*                    Choose a set of NRHS random solution vectors
*                    stored in XACT and set up the right hand side B
*
                     SRNAMT = 'SLARHS'
                     CALL SLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU,
     $                            NRHS, A, LDA, XACT, LDA, B, LDA,
     $                            ISEED, INFO )
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                     SRNAMT = 'SKYTRS'
                     CALL SKYTRS( UPLO, N, NRHS, AFAC, LDA, IWORK, X,
     $                            LDA, INFO )
*
*                    Check error code from SKYTRS and handle error.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'SKYTRS', INFO, 0, UPLO, N,
     $                               N, -1, -1, NRHS, IMAT, NFAIL,
     $                               NERRS, NOUT )
*
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
*
*                    Compute the residual for the solution
*
                     CALL SPOT07( UPLO, N, NRHS, A, LDA, X, LDA, WORK,
     $                            LDA, RWORK, RESULT( 3 ) )
*
*+    TEST 4 (Using DSYTRS2)
*                 Solve and compute residual for  A * X = B.
*
*                    Choose a set of NRHS random solution vectors
*                    stored in XACT and set up the right hand side B
*
                     SRNAMT = 'SLARHS'
                     CALL SLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU,
     $                            NRHS, A, LDA, XACT, LDA, B, LDA,
     $                            ISEED, INFO )
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                     SRNAMT = 'DSYTRS2'
                     CALL SKYTRS2( UPLO, N, NRHS, AFAC, LDA, IWORK, X,
     $                            LDA, WORK, INFO )
*
*                    Check error code from SKYTRS2 and handle error.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'SKYTRS2', INFO, 0, UPLO, N,
     $                               N, -1, -1, NRHS, IMAT, NFAIL,
     $                               NERRS, NOUT )
*
                     CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
*
*                    Compute the residual for the solution
*
                     CALL SPOT07( UPLO, N, NRHS, A, LDA, X, LDA, WORK,
     $                            LDA, RWORK, RESULT( 4 ) )
*
*+    TEST 5
*                 Check solution from generated exact solution.
*
                     CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                            RESULT( 5 ) )
*
*+    TESTS 6, 7, and 8
*                 Use iterative refinement to improve the solution.
*
                     SRNAMT = 'SKYRFS'
                     CALL SKYRFS( UPLO, N, NRHS, A, LDA, AFAC, LDA,
     $                            IWORK, B, LDA, X, LDA, RWORK,
     $                            RWORK( NRHS+1 ), WORK, IWORK( N+1 ),
     $                            INFO )
*
*                    Check error code from SKYRFS and handle error.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'SKYRFS', INFO, 0, UPLO, N,
     $                               N, -1, -1, NRHS, IMAT, NFAIL,
     $                               NERRS, NOUT )
*
                     CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                            RESULT( 6 ) )
                     CALL SPOT09( UPLO, N, NRHS, A, LDA, B, LDA, X, LDA,
     $                            XACT, LDA, RWORK, RWORK( NRHS+1 ),
     $                            RESULT( 7 ) )
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     DO 120 K = 3, 8
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALAHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS,
     $                        IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
  120                CONTINUE
                     NRUN = NRUN + 6
*
*                 End do for each value of NRHS in NSVAL.
*
  130             CONTINUE
*
*+    TEST 9
*                 Get an estimate of RCOND = 1/CNDNUM.
*
  140             CONTINUE
                  ANORM = SLANKY( '1', UPLO, N, A, LDA, RWORK )
                  SRNAMT = 'SKYCON'
                  CALL SKYCON( UPLO, N, AFAC, LDA, IWORK, ANORM, RCOND,
     $                         WORK, IWORK( N+1 ), INFO )
*
*                 Check error code from SKYCON and handle error.
*
                  IF( INFO.NE.0 )
     $               CALL ALAERH( PATH, 'SKYCON', INFO, 0, UPLO, N, N,
     $                            -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
*                 Compute the test ratio to compare to values of RCOND
*
                  RESULT( 9 ) = SGET06( RCOND, RCONDC )
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  IF( RESULT( 9 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                  CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9997 )UPLO, N, IMAT, 9,
     $                  RESULT( 9 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 1
  150          CONTINUE
*
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ',
     $      I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
 9997 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ',', 10X, ' type ', I2,
     $      ', test(', I2, ') =', G12.5 )
      RETURN
*
*     End of SCHKKY
*
      END
