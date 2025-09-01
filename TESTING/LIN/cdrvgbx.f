*> \brief \b CDRVGBX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, LA,
*                          AFB, LAFB, ASAV, B, BSAV, X, XACT, S, WORK,
*                          RWORK, IWORK, NOUT )
*
*       .. Scalar Arguments ..
*       LOGICAL            TSTERR
*       INTEGER            LA, LAFB, NN, NOUT, NRHS
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            IWORK( * ), NVAL( * )
*       REAL               RWORK( * ), S( * )
*       COMPLEX            A( * ), AFB( * ), ASAV( * ), B( * ), BSAV( * ),
*      $                   WORK( * ), X( * ), XACT( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CDRVGB tests the driver routines CGBSV, -SVX, and -SVXX.
*>
*> Note that this file is used only when the XBLAS are available,
*> otherwise cdrvgb.f defines this subroutine.
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
*>          The values of the matrix column dimension N.
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
*> \param[out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LA)
*> \endverbatim
*>
*> \param[in] LA
*> \verbatim
*>          LA is INTEGER
*>          The length of the array A.  LA >= (2*NMAX-1)*NMAX
*>          where NMAX is the largest entry in NVAL.
*> \endverbatim
*>
*> \param[out] AFB
*> \verbatim
*>          AFB is COMPLEX array, dimension (LAFB)
*> \endverbatim
*>
*> \param[in] LAFB
*> \verbatim
*>          LAFB is INTEGER
*>          The length of the array AFB.  LAFB >= (3*NMAX-2)*NMAX
*>          where NMAX is the largest entry in NVAL.
*> \endverbatim
*>
*> \param[out] ASAV
*> \verbatim
*>          ASAV is COMPLEX array, dimension (LA)
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is COMPLEX array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] BSAV
*> \verbatim
*>          BSAV is COMPLEX array, dimension (NMAX*NRHS)
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
*> \param[out] S
*> \verbatim
*>          S is REAL array, dimension (2*NMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension
*>                      (NMAX*max(3,NRHS,NMAX))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension
*>                      (max(2*NMAX,NMAX+2*NRHS))
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
      SUBROUTINE CDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, LA,
     $                   AFB, LAFB, ASAV, B, BSAV, X, XACT, S, WORK,
     $                   RWORK, IWORK, NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            LA, LAFB, NN, NOUT, NRHS
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            IWORK( * ), NVAL( * )
      REAL               RWORK( * ), S( * )
      COMPLEX            A( * ), AFB( * ), ASAV( * ), B( * ), BSAV( * ),
     $                   WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 8 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 7 )
      INTEGER            NTRAN
      PARAMETER          ( NTRAN = 3 )
*     ..
*     .. Local Scalars ..
      LOGICAL            EQUIL, NOFACT, PREFAC, TRFCON, ZEROT
      CHARACTER          DIST, EQUED, FACT, TRANS, TYPE, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, I1, I2, IEQUED, IFACT, IKL, IKU, IMAT, IN,
     $                   INFO, IOFF, ITRAN, IZERO, J, K, K1, KL, KU,
     $                   LDA, LDAFB, LDB, MODE, N, NB, NBMIN, NERRS,
     $                   NFACT, NFAIL, NIMAT, NKL, NKU, NRUN, NT,
     $                   N_ERR_BNDS
      REAL               AINVNM, AMAX, ANORM, ANORMI, ANORMO, ANRMPV,
     $                   CNDNUM, COLCND, RCOND, RCONDC, RCONDI, RCONDO,
     $                   ROLDC, ROLDI, ROLDO, ROWCND, RPVGRW,
     $                   RPVGRW_SVXX
*     ..
*     .. Local Arrays ..
      CHARACTER          EQUEDS( 4 ), FACTS( 3 ), TRANSS( NTRAN )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RDUM( 1 ), RESULT( NTESTS ), BERR( NRHS ),
     $                   ERRBNDS_N( NRHS,3 ), ERRBNDS_C( NRHS, 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANGB, CLANGE, CLANTB, SGET06, SLAMCH,
     $                   CLA_GBRPVGRW
      EXTERNAL           LSAME, CLANGB, CLANGE, CLANTB, SGET06, SLAMCH,
     $                   CLA_GBRPVGRW
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALADHD, ALAERH, ALASVM, CERRVX, CGBEQU, CGBSV,
     $                   CGBSVX, CGBT01, CGBT02, CGBT05, CGBTRF, CGBTRS,
     $                   CGET04, CLACPY, CLAQGB, CLARHS, CLASET, CLATB4,
     $                   CLATMS, XLAENV, CGBSVXX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, MAX, MIN
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
      DATA               TRANSS / 'N', 'T', 'C' /
      DATA               FACTS / 'F', 'N', 'E' /
      DATA               EQUEDS / 'N', 'R', 'C', 'B' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'GB'
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
      DO 150 IN = 1, NN
         N = NVAL( IN )
         LDB = MAX( N, 1 )
         XTYPE = 'N'
*
*        Set limits on the number of loop iterations.
*
         NKL = MAX( 1, MIN( N, 4 ) )
         IF( N.EQ.0 )
     $      NKL = 1
         NKU = NKL
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         DO 140 IKL = 1, NKL
*
*           Do for KL = 0, N-1, (3N-1)/4, and (N+1)/4. This order makes
*           it easier to skip redundant values for small values of N.
*
            IF( IKL.EQ.1 ) THEN
               KL = 0
            ELSE IF( IKL.EQ.2 ) THEN
               KL = MAX( N-1, 0 )
            ELSE IF( IKL.EQ.3 ) THEN
               KL = ( 3*N-1 ) / 4
            ELSE IF( IKL.EQ.4 ) THEN
               KL = ( N+1 ) / 4
            END IF
            DO 130 IKU = 1, NKU
*
*              Do for KU = 0, N-1, (3N-1)/4, and (N+1)/4. This order
*              makes it easier to skip redundant values for small
*              values of N.
*
               IF( IKU.EQ.1 ) THEN
                  KU = 0
               ELSE IF( IKU.EQ.2 ) THEN
                  KU = MAX( N-1, 0 )
               ELSE IF( IKU.EQ.3 ) THEN
                  KU = ( 3*N-1 ) / 4
               ELSE IF( IKU.EQ.4 ) THEN
                  KU = ( N+1 ) / 4
               END IF
*
*              Check that A and AFB are big enough to generate this
*              matrix.
*
               LDA = KL + KU + 1
               LDAFB = 2*KL + KU + 1
               IF( LDA*N.GT.LA .OR. LDAFB*N.GT.LAFB ) THEN
                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $               CALL ALADHD( NOUT, PATH )
                  IF( LDA*N.GT.LA ) THEN
                     WRITE( NOUT, FMT = 9999 )LA, N, KL, KU,
     $                  N*( KL+KU+1 )
                     NERRS = NERRS + 1
                  END IF
                  IF( LDAFB*N.GT.LAFB ) THEN
                     WRITE( NOUT, FMT = 9998 )LAFB, N, KL, KU,
     $                  N*( 2*KL+KU+1 )
                     NERRS = NERRS + 1
                  END IF
                  GO TO 130
               END IF
*
               DO 120 IMAT = 1, NIMAT
*
*                 Do the tests only if DOTYPE( IMAT ) is true.
*
                  IF( .NOT.DOTYPE( IMAT ) )
     $               GO TO 120
*
*                 Skip types 2, 3, or 4 if the matrix is too small.
*
                  ZEROT = IMAT.GE.2 .AND. IMAT.LE.4
                  IF( ZEROT .AND. N.LT.IMAT-1 )
     $               GO TO 120
*
*                 Set up parameters with CLATB4 and generate a
*                 test matrix with CLATMS.
*
                  CALL CLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM,
     $                         MODE, CNDNUM, DIST )
                  RCONDC = ONE / CNDNUM
*
                  SRNAMT = 'CLATMS'
                  CALL CLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                         CNDNUM, ANORM, KL, KU, 'Z', A, LDA, WORK,
     $                         INFO )
*
*                 Check the error code from CLATMS.
*
                  IF( INFO.NE.0 ) THEN
                     CALL ALAERH( PATH, 'CLATMS', INFO, 0, ' ', N, N,
     $                            KL, KU, -1, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 120
                  END IF
*
*                 For types 2, 3, and 4, zero one or more columns of
*                 the matrix to test that INFO is returned correctly.
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
                     IOFF = ( IZERO-1 )*LDA
                     IF( IMAT.LT.4 ) THEN
                        I1 = MAX( 1, KU+2-IZERO )
                        I2 = MIN( KL+KU+1, KU+1+( N-IZERO ) )
                        DO 20 I = I1, I2
                           A( IOFF+I ) = ZERO
   20                   CONTINUE
                     ELSE
                        DO 40 J = IZERO, N
                           DO 30 I = MAX( 1, KU+2-J ),
     $                             MIN( KL+KU+1, KU+1+( N-J ) )
                              A( IOFF+I ) = ZERO
   30                      CONTINUE
                           IOFF = IOFF + LDA
   40                   CONTINUE
                     END IF
                  END IF
*
*                 Save a copy of the matrix A in ASAV.
*
                  CALL CLACPY( 'Full', KL+KU+1, N, A, LDA, ASAV, LDA )
*
                  DO 110 IEQUED = 1, 4
                     EQUED = EQUEDS( IEQUED )
                     IF( IEQUED.EQ.1 ) THEN
                        NFACT = 3
                     ELSE
                        NFACT = 1
                     END IF
*
                     DO 100 IFACT = 1, NFACT
                        FACT = FACTS( IFACT )
                        PREFAC = LSAME( FACT, 'F' )
                        NOFACT = LSAME( FACT, 'N' )
                        EQUIL = LSAME( FACT, 'E' )
*
                        IF( ZEROT ) THEN
                           IF( PREFAC )
     $                        GO TO 100
                           RCONDO = ZERO
                           RCONDI = ZERO
*
                        ELSE IF( .NOT.NOFACT ) THEN
*
*                          Compute the condition number for comparison
*                          with the value returned by SGESVX (FACT =
*                          'N' reuses the condition number from the
*                          previous iteration with FACT = 'F').
*
                           CALL CLACPY( 'Full', KL+KU+1, N, ASAV, LDA,
     $                                  AFB( KL+1 ), LDAFB )
                           IF( EQUIL .OR. IEQUED.GT.1 ) THEN
*
*                             Compute row and column scale factors to
*                             equilibrate the matrix A.
*
                              CALL CGBEQU( N, N, KL, KU, AFB( KL+1 ),
     $                                     LDAFB, S, S( N+1 ), ROWCND,
     $                                     COLCND, AMAX, INFO )
                              IF( INFO.EQ.0 .AND. N.GT.0 ) THEN
                                 IF( LSAME( EQUED, 'R' ) ) THEN
                                    ROWCND = ZERO
                                    COLCND = ONE
                                 ELSE IF( LSAME( EQUED, 'C' ) ) THEN
                                    ROWCND = ONE
                                    COLCND = ZERO
                                 ELSE IF( LSAME( EQUED, 'B' ) ) THEN
                                    ROWCND = ZERO
                                    COLCND = ZERO
                                 END IF
*
*                                Equilibrate the matrix.
*
                                 CALL CLAQGB( N, N, KL, KU, AFB( KL+1 ),
     $                                        LDAFB, S, S( N+1 ),
     $                                        ROWCND, COLCND, AMAX,
     $                                        EQUED )
                              END IF
                           END IF
*
*                          Save the condition number of the
*                          non-equilibrated system for use in CGET04.
*
                           IF( EQUIL ) THEN
                              ROLDO = RCONDO
                              ROLDI = RCONDI
                           END IF
*
*                          Compute the 1-norm and infinity-norm of A.
*
                           ANORMO = CLANGB( '1', N, KL, KU, AFB( KL+1 ),
     $                              LDAFB, RWORK )
                           ANORMI = CLANGB( 'I', N, KL, KU, AFB( KL+1 ),
     $                              LDAFB, RWORK )
*
*                          Factor the matrix A.
*
                           CALL CGBTRF( N, N, KL, KU, AFB, LDAFB, IWORK,
     $                                  INFO )
*
*                          Form the inverse of A.
*
                           CALL CLASET( 'Full', N, N, CMPLX( ZERO ),
     $                                  CMPLX( ONE ), WORK, LDB )
                           SRNAMT = 'CGBTRS'
                           CALL CGBTRS( 'No transpose', N, KL, KU, N,
     $                                  AFB, LDAFB, IWORK, WORK, LDB,
     $                                  INFO )
*
*                          Compute the 1-norm condition number of A.
*
                           AINVNM = CLANGE( '1', N, N, WORK, LDB,
     $                              RWORK )
                           IF( ANORMO.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                              RCONDO = ONE
                           ELSE
                              RCONDO = ( ONE / ANORMO ) / AINVNM
                           END IF
*
*                          Compute the infinity-norm condition number
*                          of A.
*
                           AINVNM = CLANGE( 'I', N, N, WORK, LDB,
     $                              RWORK )
                           IF( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                              RCONDI = ONE
                           ELSE
                              RCONDI = ( ONE / ANORMI ) / AINVNM
                           END IF
                        END IF
*
                        DO 90 ITRAN = 1, NTRAN
*
*                          Do for each value of TRANS.
*
                           TRANS = TRANSS( ITRAN )
                           IF( ITRAN.EQ.1 ) THEN
                              RCONDC = RCONDO
                           ELSE
                              RCONDC = RCONDI
                           END IF
*
*                          Restore the matrix A.
*
                           CALL CLACPY( 'Full', KL+KU+1, N, ASAV, LDA,
     $                                  A, LDA )
*
*                          Form an exact solution and set the right hand
*                          side.
*
                           SRNAMT = 'CLARHS'
                           CALL CLARHS( PATH, XTYPE, 'Full', TRANS, N,
     $                                  N, KL, KU, NRHS, A, LDA, XACT,
     $                                  LDB, B, LDB, ISEED, INFO )
                           XTYPE = 'C'
                           CALL CLACPY( 'Full', N, NRHS, B, LDB, BSAV,
     $                                  LDB )
*
                           IF( NOFACT .AND. ITRAN.EQ.1 ) THEN
*
*                             --- Test CGBSV  ---
*
*                             Compute the LU factorization of the matrix
*                             and solve the system.
*
                              CALL CLACPY( 'Full', KL+KU+1, N, A, LDA,
     $                                     AFB( KL+1 ), LDAFB )
                              CALL CLACPY( 'Full', N, NRHS, B, LDB, X,
     $                                     LDB )
*
                              SRNAMT = 'CGBSV '
                              CALL CGBSV( N, KL, KU, NRHS, AFB, LDAFB,
     $                                    IWORK, X, LDB, INFO )
*
*                             Check error code from CGBSV .
*
                              IF( INFO.NE.IZERO )
     $                           CALL ALAERH( PATH, 'CGBSV ', INFO,
     $                                        IZERO, ' ', N, N, KL, KU,
     $                                        NRHS, IMAT, NFAIL, NERRS,
     $                                        NOUT )
*
*                             Reconstruct matrix from factors and
*                             compute residual.
*
                              CALL CGBT01( N, N, KL, KU, A, LDA, AFB,
     $                                     LDAFB, IWORK, WORK,
     $                                     RESULT( 1 ) )
                              NT = 1
                              IF( IZERO.EQ.0 ) THEN
*
*                                Compute residual of the computed
*                                solution.
*
                                 CALL CLACPY( 'Full', N, NRHS, B, LDB,
     $                                        WORK, LDB )
                                 CALL CGBT02( 'No transpose', N, N, KL,
     $                                        KU, NRHS, A, LDA, X, LDB,
     $                                        WORK, LDB, RWORK,
     $                                        RESULT( 2 ) )
*
*                                Check solution from generated exact
*                                solution.
*
                                 CALL CGET04( N, NRHS, X, LDB, XACT,
     $                                        LDB, RCONDC, RESULT( 3 ) )
                                 NT = 3
                              END IF
*
*                             Print information about the tests that did
*                             not pass the threshold.
*
                              DO 50 K = 1, NT
                                 IF( RESULT( K ).GE.THRESH ) THEN
                                    IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                                 CALL ALADHD( NOUT, PATH )
                                    WRITE( NOUT, FMT = 9997 )'CGBSV ',
     $                                 N, KL, KU, IMAT, K, RESULT( K )
                                    NFAIL = NFAIL + 1
                                 END IF
   50                         CONTINUE
                              NRUN = NRUN + NT
                           END IF
*
*                          --- Test CGBSVX ---
*
                           IF( .NOT.PREFAC )
     $                        CALL CLASET( 'Full', 2*KL+KU+1, N,
     $                                     CMPLX( ZERO ), CMPLX( ZERO ),
     $                                     AFB, LDAFB )
                           CALL CLASET( 'Full', N, NRHS, CMPLX( ZERO ),
     $                                  CMPLX( ZERO ), X, LDB )
                           IF( IEQUED.GT.1 .AND. N.GT.0 ) THEN
*
*                             Equilibrate the matrix if FACT = 'F' and
*                             EQUED = 'R', 'C', or 'B'.
*
                              CALL CLAQGB( N, N, KL, KU, A, LDA, S,
     $                                     S( N+1 ), ROWCND, COLCND,
     $                                     AMAX, EQUED )
                           END IF
*
*                          Solve the system and compute the condition
*                          number and error bounds using CGBSVX.
*
                           SRNAMT = 'CGBSVX'
                           CALL CGBSVX( FACT, TRANS, N, KL, KU, NRHS, A,
     $                                  LDA, AFB, LDAFB, IWORK, EQUED,
     $                                  S, S( LDB+1 ), B, LDB, X, LDB,
     $                                  RCOND, RWORK, RWORK( NRHS+1 ),
     $                                  WORK, RWORK( 2*NRHS+1 ), INFO )
*
*                          Check the error code from CGBSVX.
*
                           IF( INFO.NE.IZERO )
     $                        CALL ALAERH( PATH, 'CGBSVX', INFO, IZERO,
     $                                     FACT // TRANS, N, N, KL, KU,
     $                                     NRHS, IMAT, NFAIL, NERRS,
     $                                     NOUT )
*
*                          Compare RWORK(2*NRHS+1) from CGBSVX with the
*                          computed reciprocal pivot growth RPVGRW
*
                           IF( INFO.NE.0 ) THEN
                              ANRMPV = ZERO
                              DO 70 J = 1, INFO
                                 DO 60 I = MAX( KU+2-J, 1 ),
     $                                   MIN( N+KU+1-J, KL+KU+1 )
                                    ANRMPV = MAX( ANRMPV,
     $                                       ABS( A( I+( J-1 )*LDA ) ) )
   60                            CONTINUE
   70                         CONTINUE
                              RPVGRW = CLANTB( 'M', 'U', 'N', INFO,
     $                                 MIN( INFO-1, KL+KU ),
     $                                 AFB( MAX( 1, KL+KU+2-INFO ) ),
     $                                 LDAFB, RDUM )
                              IF( RPVGRW.EQ.ZERO ) THEN
                                 RPVGRW = ONE
                              ELSE
                                 RPVGRW = ANRMPV / RPVGRW
                              END IF
                           ELSE
                              RPVGRW = CLANTB( 'M', 'U', 'N', N, KL+KU,
     $                                 AFB, LDAFB, RDUM )
                              IF( RPVGRW.EQ.ZERO ) THEN
                                 RPVGRW = ONE
                              ELSE
                                 RPVGRW = CLANGB( 'M', N, KL, KU, A,
     $                                    LDA, RDUM ) / RPVGRW
                              END IF
                           END IF
                           RESULT( 7 ) = ABS( RPVGRW-RWORK( 2*NRHS+1 ) )
     $                                    / MAX( RWORK( 2*NRHS+1 ),
     $                                   RPVGRW ) / SLAMCH( 'E' )
*
                           IF( .NOT.PREFAC ) THEN
*
*                             Reconstruct matrix from factors and
*                             compute residual.
*
                              CALL CGBT01( N, N, KL, KU, A, LDA, AFB,
     $                                     LDAFB, IWORK, WORK,
     $                                     RESULT( 1 ) )
                              K1 = 1
                           ELSE
                              K1 = 2
                           END IF
*
                           IF( INFO.EQ.0 ) THEN
                              TRFCON = .FALSE.
*
*                             Compute residual of the computed solution.
*
                              CALL CLACPY( 'Full', N, NRHS, BSAV, LDB,
     $                                     WORK, LDB )
                              CALL CGBT02( TRANS, N, N, KL, KU, NRHS,
     $                                     ASAV, LDA, X, LDB, WORK, LDB,
     $                                     RWORK( 2*NRHS+1 ),
     $                                     RESULT( 2 ) )
*
*                             Check solution from generated exact
*                             solution.
*
                              IF( NOFACT .OR. ( PREFAC .AND.
     $                            LSAME( EQUED, 'N' ) ) ) THEN
                                 CALL CGET04( N, NRHS, X, LDB, XACT,
     $                                        LDB, RCONDC, RESULT( 3 ) )
                              ELSE
                                 IF( ITRAN.EQ.1 ) THEN
                                    ROLDC = ROLDO
                                 ELSE
                                    ROLDC = ROLDI
                                 END IF
                                 CALL CGET04( N, NRHS, X, LDB, XACT,
     $                                        LDB, ROLDC, RESULT( 3 ) )
                              END IF
*
*                             Check the error bounds from iterative
*                             refinement.
*
                              CALL CGBT05( TRANS, N, KL, KU, NRHS, ASAV,
     $                                     LDA, BSAV, LDB, X, LDB, XACT,
     $                                     LDB, RWORK, RWORK( NRHS+1 ),
     $                                     RESULT( 4 ) )
                           ELSE
                              TRFCON = .TRUE.
                           END IF
*
*                          Compare RCOND from CGBSVX with the computed
*                          value in RCONDC.
*
                           RESULT( 6 ) = SGET06( RCOND, RCONDC )
*
*                          Print information about the tests that did
*                          not pass the threshold.
*
                           IF( .NOT.TRFCON ) THEN
                              DO 80 K = K1, NTESTS
                                 IF( RESULT( K ).GE.THRESH ) THEN
                                    IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                                 CALL ALADHD( NOUT, PATH )
                                    IF( PREFAC ) THEN
                                       WRITE( NOUT, FMT = 9995 )
     $                                    'CGBSVX', FACT, TRANS, N, KL,
     $                                    KU, EQUED, IMAT, K,
     $                                    RESULT( K )
                                    ELSE
                                       WRITE( NOUT, FMT = 9996 )
     $                                    'CGBSVX', FACT, TRANS, N, KL,
     $                                    KU, IMAT, K, RESULT( K )
                                    END IF
                                    NFAIL = NFAIL + 1
                                 END IF
   80                         CONTINUE
                              NRUN = NRUN + 7 - K1
                           ELSE
                              IF( RESULT( 1 ).GE.THRESH .AND. .NOT.
     $                            PREFAC ) THEN
                                 IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                              CALL ALADHD( NOUT, PATH )
                                 IF( PREFAC ) THEN
                                    WRITE( NOUT, FMT = 9995 )'CGBSVX',
     $                                 FACT, TRANS, N, KL, KU, EQUED,
     $                                 IMAT, 1, RESULT( 1 )
                                 ELSE
                                    WRITE( NOUT, FMT = 9996 )'CGBSVX',
     $                                 FACT, TRANS, N, KL, KU, IMAT, 1,
     $                                 RESULT( 1 )
                                 END IF
                                 NFAIL = NFAIL + 1
                                 NRUN = NRUN + 1
                              END IF
                              IF( RESULT( 6 ).GE.THRESH ) THEN
                                 IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                              CALL ALADHD( NOUT, PATH )
                                 IF( PREFAC ) THEN
                                    WRITE( NOUT, FMT = 9995 )'CGBSVX',
     $                                 FACT, TRANS, N, KL, KU, EQUED,
     $                                 IMAT, 6, RESULT( 6 )
                                 ELSE
                                    WRITE( NOUT, FMT = 9996 )'CGBSVX',
     $                                 FACT, TRANS, N, KL, KU, IMAT, 6,
     $                                 RESULT( 6 )
                                 END IF
                                 NFAIL = NFAIL + 1
                                 NRUN = NRUN + 1
                              END IF
                              IF( RESULT( 7 ).GE.THRESH ) THEN
                                 IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                              CALL ALADHD( NOUT, PATH )
                                 IF( PREFAC ) THEN
                                    WRITE( NOUT, FMT = 9995 )'CGBSVX',
     $                                 FACT, TRANS, N, KL, KU, EQUED,
     $                                 IMAT, 7, RESULT( 7 )
                                 ELSE
                                    WRITE( NOUT, FMT = 9996 )'CGBSVX',
     $                                 FACT, TRANS, N, KL, KU, IMAT, 7,
     $                                 RESULT( 7 )
                                 END IF
                                 NFAIL = NFAIL + 1
                                 NRUN = NRUN + 1
                              END IF
                           END IF

*                    --- Test CGBSVXX ---

*                    Restore the matrices A and B.

c                     write(*,*) 'begin cgbsvxx testing'

                     CALL CLACPY( 'Full', KL+KU+1, N, ASAV, LDA, A,
     $                          LDA )
                     CALL CLACPY( 'Full', N, NRHS, BSAV, LDB, B, LDB )

                     IF( .NOT.PREFAC )
     $                  CALL CLASET( 'Full', 2*KL+KU+1, N,
     $                               CMPLX( ZERO ), CMPLX( ZERO ),
     $                               AFB, LDAFB )
                     CALL CLASET( 'Full', N, NRHS,
     $                            CMPLX( ZERO ), CMPLX( ZERO ),
     $                               X, LDB )
                     IF( IEQUED.GT.1 .AND. N.GT.0 ) THEN
*
*                       Equilibrate the matrix if FACT = 'F' and
*                       EQUED = 'R', 'C', or 'B'.
*
                        CALL CLAQGB( N, N, KL, KU, A, LDA, S,
     $                       S( N+1 ), ROWCND, COLCND, AMAX, EQUED )
                     END IF
*
*                    Solve the system and compute the condition number
*                    and error bounds using CGBSVXX.
*
                     SRNAMT = 'CGBSVXX'
                     N_ERR_BNDS = 3
                     CALL CGBSVXX( FACT, TRANS, N, KL, KU, NRHS, A, LDA,
     $                    AFB, LDAFB, IWORK, EQUED, S, S( N+1 ), B, LDB,
     $                    X, LDB, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS,
     $                    ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK,
     $                    RWORK, INFO )
*
*                    Check the error code from CGBSVXX.
*
                     IF( INFO.EQ.N+1 ) GOTO 90
                     IF( INFO.NE.IZERO ) THEN
                        CALL ALAERH( PATH, 'CGBSVXX', INFO, IZERO,
     $                               FACT // TRANS, N, N, -1, -1, NRHS,
     $                               IMAT, NFAIL, NERRS, NOUT )
                        GOTO 90
                     END IF
*
*                    Compare rpvgrw_svxx from CGESVXX with the computed
*                    reciprocal pivot growth factor RPVGRW
*

                     IF ( INFO .GT. 0 .AND. INFO .LT. N+1 ) THEN
                        RPVGRW = CLA_GBRPVGRW(N, KL, KU, INFO, A, LDA,
     $                       AFB, LDAFB)
                     ELSE
                        RPVGRW = CLA_GBRPVGRW(N, KL, KU, N, A, LDA,
     $                       AFB, LDAFB)
                     ENDIF

                     RESULT( 7 ) = ABS( RPVGRW-rpvgrw_svxx ) /
     $                             MAX( rpvgrw_svxx, RPVGRW ) /
     $                             SLAMCH( 'E' )
*
                     IF( .NOT.PREFAC ) THEN
*
*                       Reconstruct matrix from factors and compute
*                       residual.
*
                        CALL CGBT01( N, N, KL, KU, A, LDA, AFB, LDAFB,
     $                       IWORK, WORK( 2*NRHS+1 ), RESULT( 1 ) )
                        K1 = 1
                     ELSE
                        K1 = 2
                     END IF
*
                     IF( INFO.EQ.0 ) THEN
                        TRFCON = .FALSE.
*
*                       Compute residual of the computed solution.
*
                        CALL CLACPY( 'Full', N, NRHS, BSAV, LDB, WORK,
     $                               LDB )
                        CALL CGBT02( TRANS, N, N, KL, KU, NRHS, ASAV,
     $                               LDA, X, LDB, WORK, LDB, RWORK,
     $                               RESULT( 2 ) )
*
*                       Check solution from generated exact solution.
*
                        IF( NOFACT .OR. ( PREFAC .AND. LSAME( EQUED,
     $                      'N' ) ) ) THEN
                           CALL CGET04( N, NRHS, X, LDB, XACT, LDB,
     $                                  RCONDC, RESULT( 3 ) )
                        ELSE
                           IF( ITRAN.EQ.1 ) THEN
                              ROLDC = ROLDO
                           ELSE
                              ROLDC = ROLDI
                           END IF
                           CALL CGET04( N, NRHS, X, LDB, XACT, LDB,
     $                                  ROLDC, RESULT( 3 ) )
                        END IF
                     ELSE
                        TRFCON = .TRUE.
                     END IF
*
*                    Compare RCOND from CGBSVXX with the computed value
*                    in RCONDC.
*
                     RESULT( 6 ) = SGET06( RCOND, RCONDC )
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     IF( .NOT.TRFCON ) THEN
                        DO 45 K = K1, NTESTS
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                           CALL ALADHD( NOUT, PATH )
                              IF( PREFAC ) THEN
                                 WRITE( NOUT, FMT = 9995 )'CGBSVXX',
     $                                FACT, TRANS, N, KL, KU, EQUED,
     $                                IMAT, K, RESULT( K )
                              ELSE
                                 WRITE( NOUT, FMT = 9996 )'CGBSVXX',
     $                                FACT, TRANS, N, KL, KU, IMAT, K,
     $                                RESULT( K )
                              END IF
                              NFAIL = NFAIL + 1
                           END IF
 45                     CONTINUE
                        NRUN = NRUN + 7 - K1
                     ELSE
                        IF( RESULT( 1 ).GE.THRESH .AND. .NOT.PREFAC )
     $                       THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT,
     $                             TRANS, N, KL, KU, EQUED, IMAT, 1,
     $                             RESULT( 1 )
                           ELSE
                              WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT,
     $                             TRANS, N, KL, KU, IMAT, 1,
     $                             RESULT( 1 )
                           END IF
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        END IF
                        IF( RESULT( 6 ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT,
     $                             TRANS, N, KL, KU, EQUED, IMAT, 6,
     $                             RESULT( 6 )
                           ELSE
                              WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT,
     $                             TRANS, N, KL, KU, IMAT, 6,
     $                             RESULT( 6 )
                           END IF
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        END IF
                        IF( RESULT( 7 ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALADHD( NOUT, PATH )
                           IF( PREFAC ) THEN
                              WRITE( NOUT, FMT = 9995 )'CGBSVXX', FACT,
     $                             TRANS, N, KL, KU, EQUED, IMAT, 7,
     $                             RESULT( 7 )
                           ELSE
                              WRITE( NOUT, FMT = 9996 )'CGBSVXX', FACT,
     $                             TRANS, N, KL, KU, IMAT, 7,
     $                             RESULT( 7 )
                           END IF
                           NFAIL = NFAIL + 1
                           NRUN = NRUN + 1
                        END IF
*
                     END IF
*
   90                   CONTINUE
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )
*

*     Test Error Bounds from CGBSVXX

      CALL CEBCHVXX(THRESH, PATH)

 9999 FORMAT( ' *** In CDRVGB, LA=', I5, ' is too small for N=', I5,
     $      ', KU=', I5, ', KL=', I5, / ' ==> Increase LA to at least ',
     $      I5 )
 9998 FORMAT( ' *** In CDRVGB, LAFB=', I5, ' is too small for N=', I5,
     $      ', KU=', I5, ', KL=', I5, /
     $      ' ==> Increase LAFB to at least ', I5 )
 9997 FORMAT( 1X, A, ', N=', I5, ', KL=', I5, ', KU=', I5, ', type ',
     $      I1, ', test(', I1, ')=', G12.5 )
 9996 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',',
     $      I5, ',...), type ', I1, ', test(', I1, ')=', G12.5 )
 9995 FORMAT( 1X, A, '( ''', A1, ''',''', A1, ''',', I5, ',', I5, ',',
     $      I5, ',...), EQUED=''', A1, ''', type ', I1, ', test(', I1,
     $      ')=', G12.5 )
*
      RETURN
*
*     End of CDRVGBX
*
      END
