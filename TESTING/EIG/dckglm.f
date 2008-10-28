      SUBROUTINE DCKGLM( NN, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
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
*  DCKGLM tests DGGGLM - subroutine for solving generalized linear
*                        model problem.
*
*  Arguments
*  =========
*
*  NN      (input) INTEGER
*          The number of values of N, M and P contained in the vectors
*          NVAL, MVAL and PVAL.
*
*  MVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension M.
*
*  PVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension P.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix row dimension N.
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
*          included in the output file if RESID >= THRESH.  To have
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
*  X       (workspace) DOUBLE PRECISION array, dimension (4*NMAX)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (NMAX)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
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
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRSTT
      CHARACTER          DISTA, DISTB, TYPE
      CHARACTER*3        PATH
      INTEGER            I, IINFO, IK, IMAT, KLA, KLB, KUA, KUB, LDA,
     $                   LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, P
      DOUBLE PRECISION   ANORM, BNORM, CNDNMA, CNDNMB, RESID
*     ..
*     .. Local Arrays ..
      LOGICAL            DOTYPE( NTYPES )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLARND
      EXTERNAL           DLARND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHDG, ALAREQ, ALASUM, DGLMTS, DLATB9, DLATMS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Initialize constants.
*
      PATH( 1: 3 ) = 'GLM'
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
         IF( M.GT.N .OR. N.GT.M+P ) THEN
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
         IF( M.GT.N .OR. N.GT.M+P )
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
            CALL DLATMS( N, M, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA,
     $                   ANORM, KLA, KUA, 'No packing', A, LDA, WORK,
     $                   IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 30
            END IF
*
            CALL DLATMS( N, P, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB,
     $                   BNORM, KLB, KUB, 'No packing', B, LDB, WORK,
     $                   IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 30
            END IF
*
*           Generate random left hand side vector of GLM
*
            DO 20 I = 1, N
               X( I ) = DLARND( 2, ISEED )
   20       CONTINUE
*
            CALL DGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, X,
     $                   X( NMAX+1 ), X( 2*NMAX+1 ), X( 3*NMAX+1 ),
     $                   WORK, LWORK, RWORK, RESID )
*
*           Print information about the tests that did not
*           pass the threshold.
*
            IF( RESID.GE.THRESH ) THEN
               IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                  FIRSTT = .FALSE.
                  CALL ALAHDG( NOUT, PATH )
               END IF
               WRITE( NOUT, FMT = 9998 )N, M, P, IMAT, 1, RESID
               NFAIL = NFAIL + 1
            END IF
            NRUN = NRUN + 1
*
   30    CONTINUE
   40 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
*
 9999 FORMAT( ' DLATMS in DCKGLM INFO = ', I5 )
 9998 FORMAT( ' N=', I4, ' M=', I4, ', P=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' *** Invalid input  for GLM:  M = ', I6, ', P = ', I6,
     $      ', N = ', I6, ';', / '     must satisfy M <= N <= M+P  ',
     $      '(this set of values will be skipped)' )
      RETURN
*
*     End of DCKGLM
*
      END
