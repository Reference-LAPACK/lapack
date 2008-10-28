      SUBROUTINE DCKLSE( NN, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
     $                   NMAX, A, AF, B, BF, X, WORK, RWORK, NIN, NOUT,
     $                   INFO )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, NIN, NMATS, NMAX, NN, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * )
      DOUBLE PRECISION   A( * ), AF( * ), B( * ), BF( * ), RWORK( * ),
     $                   WORK( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DCKLSE tests DGGLSE - a subroutine for solving linear equality
*  constrained least square problem (LSE).
*
*  Arguments
*  =========
*
*  NN      (input) INTEGER
*          The number of values of (M,P,N) contained in the vectors
*          (MVAL, PVAL, NVAL).
*
*  MVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix row(column) dimension M.
*
*  PVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix row(column) dimension P.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column(row) dimension N.
*
*  NMATS   (input) INTEGER
*          The number of matrix types to be tested for each combination
*          of matrix dimensions.  If NMATS >= NTYPES (the maximum
*          number of matrix types), then all the different types are
*          generated for testing.  If NMATS < NTYPES, another input line
*          is read to get the numbers of the matrix types to be used.
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator.  The array
*          elements should be between 0 and 4095, otherwise they will be
*          reduced mod 4096, and ISEED(4) must be odd.
*          On exit, the next seed in the random number sequence after
*          all the test matrices have been generated.
*
*  THRESH  (input) DOUBLE PRECISION
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULT >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  NMAX    (input) INTEGER
*          The maximum value permitted for M or N, used in dimensioning
*          the work arrays.
*
*  A       (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  AF      (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  B       (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  BF      (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  X       (workspace) DOUBLE PRECISION array, dimension (5*NMAX)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (NMAX)
*
*  NIN     (input) INTEGER
*          The unit number for input.
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  INFO    (output) INTEGER
*          = 0 :  successful exit
*          > 0 :  If DLATMS returns an error code, the absolute value
*                 of it is returned.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 7 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRSTT
      CHARACTER          DISTA, DISTB, TYPE
      CHARACTER*3        PATH
      INTEGER            I, IINFO, IK, IMAT, KLA, KLB, KUA, KUB, LDA,
     $                   LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN,
     $                   NT, P
      DOUBLE PRECISION   ANORM, BNORM, CNDNMA, CNDNMB
*     ..
*     .. Local Arrays ..
      LOGICAL            DOTYPE( NTYPES )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHDG, ALAREQ, ALASUM, DLARHS, DLATB9, DLATMS,
     $                   DLSETS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 3 ) = 'LSE'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
      LDA = NMAX
      LDB = NMAX
      LWORK = NMAX*NMAX
*
*     Check for valid input values.
*
      DO 10 IK = 1, NN
         M = MVAL( IK )
         P = PVAL( IK )
         N = NVAL( IK )
         IF( P.GT.N .OR. N.GT.M+P ) THEN
            IF( FIRSTT ) THEN
               WRITE( NOUT, FMT = * )
               FIRSTT = .FALSE.
            END IF
            WRITE( NOUT, FMT = 9997 )M, P, N
         END IF
   10 CONTINUE
      FIRSTT = .TRUE.
*
*     Do for each value of M in MVAL.
*
      DO 40 IK = 1, NN
         M = MVAL( IK )
         P = PVAL( IK )
         N = NVAL( IK )
         IF( P.GT.N .OR. N.GT.M+P )
     $      GO TO 40
*
         DO 30 IMAT = 1, NTYPES
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 30
*
*           Set up parameters with DLATB9 and generate test
*           matrices A and B with DLATMS.
*
            CALL DLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB,
     $                   ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB,
     $                   DISTA, DISTB )
*
            CALL DLATMS( M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA,
     $                   ANORM, KLA, KUA, 'No packing', A, LDA, WORK,
     $                   IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 30
            END IF
*
            CALL DLATMS( P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB,
     $                   BNORM, KLB, KUB, 'No packing', B, LDB, WORK,
     $                   IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 30
            END IF
*
*           Generate the right-hand sides C and D for the LSE.
*
            CALL DLARHS( 'DGE', 'New solution', 'Upper', 'N', M, N,
     $                   MAX( M-1, 0 ), MAX( N-1, 0 ), 1, A, LDA,
     $                   X( 4*NMAX+1 ), MAX( N, 1 ), X, MAX( M, 1 ),
     $                   ISEED, IINFO )
*
            CALL DLARHS( 'DGE', 'Computed', 'Upper', 'N', P, N,
     $                   MAX( P-1, 0 ), MAX( N-1, 0 ), 1, B, LDB,
     $                   X( 4*NMAX+1 ), MAX( N, 1 ), X( 2*NMAX+1 ),
     $                   MAX( P, 1 ), ISEED, IINFO )
*
            NT = 2
*
            CALL DLSETS( M, P, N, A, AF, LDA, B, BF, LDB, X,
     $                   X( NMAX+1 ), X( 2*NMAX+1 ), X( 3*NMAX+1 ),
     $                   X( 4*NMAX+1 ), WORK, LWORK, RWORK,
     $                   RESULT( 1 ) )
*
*           Print information about the tests that did not
*           pass the threshold.
*
            DO 20 I = 1, NT
               IF( RESULT( I ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                     FIRSTT = .FALSE.
                     CALL ALAHDG( NOUT, PATH )
                  END IF
                  WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I,
     $               RESULT( I )
                  NFAIL = NFAIL + 1
               END IF
   20       CONTINUE
            NRUN = NRUN + NT
*
   30    CONTINUE
   40 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
*
 9999 FORMAT( ' DLATMS in DCKLSE   INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' *** Invalid input  for LSE:  M = ', I6, ', P = ', I6,
     $      ', N = ', I6, ';', / '     must satisfy P <= N <= P+M  ',
     $      '(this set of values will be skipped)' )
      RETURN
*
*     End of DCKLSE
*
      END
