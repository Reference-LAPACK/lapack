      SUBROUTINE CCHKPS( DOTYPE, NN, NVAL, NNB, NBVAL, NRANK, RANKVAL,
     $                   THRESH, TSTERR, NMAX, A, AFAC, PERM, PIV, WORK,
     $                   RWORK, NOUT )
*
*  -- LAPACK test routine (version 3.2.1) --
*     Craig Lucas, University of Manchester / NAG Ltd.
*     -- April 2009
*
*     .. Scalar Arguments ..
      REAL               THRESH
      INTEGER            NMAX, NN, NNB, NOUT, NRANK
      LOGICAL            TSTERR
*     ..
*     .. Array Arguments ..
      COMPLEX            A( * ), AFAC( * ), PERM( * ), WORK( * )
      REAL               RWORK( * )
      INTEGER            NBVAL( * ), NVAL( * ), PIV( * ), RANKVAL( * )
      LOGICAL            DOTYPE( * )
*     ..
*
*  Purpose
*  =======
*
*  CCHKPS tests CPSTRF.
*
*  Arguments
*  =========
*
*  DOTYPE  (input) LOGICAL array, dimension (NTYPES)
*          The matrix types to be used for testing.  Matrices of type j
*          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
*          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix dimension N.
*
*  NNB     (input) INTEGER
*          The number of values of NB contained in the vector NBVAL.
*
*  NBVAL   (input) INTEGER array, dimension (NBVAL)
*          The values of the block size NB.
*
*  NRANK   (input) INTEGER
*          The number of values of RANK contained in the vector RANKVAL.
*
*  RANKVAL (input) INTEGER array, dimension (NBVAL)
*          The values of the block size NB.
*
*  THRESH  (input) REAL
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULT >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  TSTERR  (input) LOGICAL
*          Flag that indicates whether error exits are to be tested.
*
*  NMAX    (input) INTEGER
*          The maximum value permitted for N, used in dimensioning the
*          work arrays.
*
*  A       (workspace) COMPLEX array, dimension (NMAX*NMAX)
*
*  AFAC    (workspace) COMPLEX array, dimension (NMAX*NMAX)
*
*  PERM    (workspace) COMPLEX array, dimension (NMAX*NMAX)
*
*  PIV     (workspace) INTEGER array, dimension (NMAX)
*
*  WORK    (workspace) COMPLEX array, dimension (NMAX*3)
*
*  RWORK   (workspace) REAL array, dimension (NMAX)
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 9 )
*     ..
*     .. Local Scalars ..
      REAL               ANORM, CNDNUM, RESULT, TOL
      INTEGER            COMPRANK, I, IMAT, IN, INB, INFO, IRANK, IUPLO,
     $                   IZERO, KL, KU, LDA, MODE, N, NB, NERRS, NFAIL,
     $                   NIMAT, NRUN, RANK, RANKDIFF
      CHARACTER          DIST, TYPE, UPLO
      CHARACTER*3        PATH
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      CHARACTER          UPLOS( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, CERRPS, CLACPY,
     $                   CLATB5, CLATMT, CPST01, CPSTRF, XLAENV
*     ..
*     .. Scalars in Common ..
      INTEGER            INFOT, NUNIT
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, REAL, CEILING
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Complex Precision'
      PATH( 2: 3 ) = 'PS'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 100 I = 1, 4
         ISEED( I ) = ISEEDY( I )
  100 CONTINUE
*
*     Test the error exits
*
      IF( TSTERR )
     $   CALL CERRPS( PATH, NOUT )
      INFOT = 0
*
*     Do for each value of N in NVAL
*
      DO 150 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         IZERO = 0
         DO 140 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 140
*
*              Do for each value of RANK in RANKVAL
*
            DO 130 IRANK = 1, NRANK
*
*              Only repeat test 3 to 5 for different ranks
*              Other tests use full rank
*
               IF( ( IMAT.LT.3 .OR. IMAT.GT.5 ) .AND. IRANK.GT.1 )
     $            GO TO 130
*
               RANK = CEILING( ( N * REAL( RANKVAL( IRANK ) ) )
     $              / 100.E+0 )
*
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
               DO 120 IUPLO = 1, 2
                  UPLO = UPLOS( IUPLO )
*
*              Set up parameters with CLATB5 and generate a test matrix
*              with CLATMT.
*
                  CALL CLATB5( PATH, IMAT, N, TYPE, KL, KU, ANORM,
     $                         MODE, CNDNUM, DIST )
*
                  SRNAMT = 'CLATMT'
                  CALL CLATMT( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                         CNDNUM, ANORM, RANK, KL, KU, UPLO, A,
     $                         LDA, WORK, INFO )
*
*              Check error code from CLATMT.
*
                  IF( INFO.NE.0 ) THEN
                    CALL ALAERH( PATH, 'CLATMT', INFO, 0, UPLO, N,
     $                           N, -1, -1, -1, IMAT, NFAIL, NERRS,
     $                           NOUT )
                     GO TO 120
                  END IF
*
*              Do for each value of NB in NBVAL
*
                  DO 110 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
*
*                 Compute the pivoted L*L' or U'*U factorization
*                 of the matrix.
*
                     CALL CLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     SRNAMT = 'CPSTRF'
*
*                 Use default tolerance
*
                     TOL = -ONE
                     CALL CPSTRF( UPLO, N, AFAC, LDA, PIV, COMPRANK,
     $                            TOL, RWORK, INFO )
*
*                 Check error code from CPSTRF.
*
                     IF( (INFO.LT.IZERO)
     $                    .OR.(INFO.NE.IZERO.AND.RANK.EQ.N)
     $                    .OR.(INFO.LE.IZERO.AND.RANK.LT.N) ) THEN
                        CALL ALAERH( PATH, 'CPSTRF', INFO, IZERO,
     $                               UPLO, N, N, -1, -1, NB, IMAT,
     $                               NFAIL, NERRS, NOUT )
                        GO TO 110
                     END IF
*
*                 Skip the test if INFO is not 0.
*
                     IF( INFO.NE.0 )
     $                  GO TO 110
*
*                 Reconstruct matrix from factors and compute residual.
*
*                 PERM holds permuted L*L^T or U^T*U
*
                     CALL CPST01( UPLO, N, A, LDA, AFAC, LDA, PERM, LDA,
     $                            PIV, RWORK, RESULT, COMPRANK )
*
*                 Print information about the tests that did not pass
*                 the threshold or where computed rank was not RANK.
*
                     IF( N.EQ.0 )
     $                  COMPRANK = 0
                     RANKDIFF = RANK - COMPRANK
                     IF( RESULT.GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 )UPLO, N, RANK,
     $                     RANKDIFF, NB, IMAT, RESULT
                        NFAIL = NFAIL + 1
                     END IF
                     NRUN = NRUN + 1
  110             CONTINUE
*
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', RANK =', I3,
     $      ', Diff =', I5, ', NB =', I4, ', type ', I2, ', Ratio =',
     $      G12.5 )
      RETURN
*
*     End of CCHKPS
*
      END
