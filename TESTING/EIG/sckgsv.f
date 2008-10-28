      SUBROUTINE SCKGSV( NM, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
     $                   NMAX, A, AF, B, BF, U, V, Q, ALPHA, BETA, R,
     $                   IWORK, WORK, RWORK, NIN, NOUT, INFO )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, NIN, NM, NMATS, NMAX, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * ),
     $                   PVAL( * )
      REAL               A( * ), AF( * ), ALPHA( * ), B( * ), BETA( * ),
     $                   BF( * ), Q( * ), R( * ), RWORK( * ), U( * ),
     $                   V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SCKGSV tests SGGSVD:
*         the GSVD for M-by-N matrix A and P-by-N matrix B.
*
*  Arguments
*  =========
*
*  NM      (input) INTEGER
*          The number of values of M contained in the vector MVAL.
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix row dimension M.
*
*  PVAL    (input) INTEGER array, dimension (NP)
*          The values of the matrix row dimension P.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension N.
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
*  THRESH  (input) REAL
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULT >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  NMAX    (input) INTEGER
*          The maximum value permitted for M or N, used in dimensioning
*          the work arrays.
*
*  A       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  AF      (workspace) REAL array, dimension (NMAX*NMAX)
*
*  B       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  BF      (workspace) REAL array, dimension (NMAX*NMAX)
*
*  U       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  V       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  Q       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  ALPHA   (workspace) REAL array, dimension (NMAX)
*
*  BETA    (workspace) REAL array, dimension (NMAX)
*
*  R       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  WORK    (workspace) REAL array, dimension (NMAX*NMAX)
*
*  RWORK   (workspace) REAL array, dimension (NMAX)
*
*  NIN     (input) INTEGER
*          The unit number for input.
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  INFO    (output) INTEGER
*          = 0 :  successful exit
*          > 0 :  If SLATMS returns an error code, the absolute value
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
      INTEGER            I, IINFO, IM, IMAT, KLA, KLB, KUA, KUB, LDA,
     $                   LDB, LDQ, LDR, LDU, LDV, LWORK, M, MODEA,
     $                   MODEB, N, NFAIL, NRUN, NT, P
      REAL               ANORM, BNORM, CNDNMA, CNDNMB
*     ..
*     .. Local Arrays ..
      LOGICAL            DOTYPE( NTYPES )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHDG, ALAREQ, ALASUM, SGSVTS, SLATB9, SLATMS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 3 ) = 'GSV'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
      LDA = NMAX
      LDB = NMAX
      LDU = NMAX
      LDV = NMAX
      LDQ = NMAX
      LDR = NMAX
      LWORK = NMAX*NMAX
*
*     Do for each value of M in MVAL.
*
      DO 30 IM = 1, NM
         M = MVAL( IM )
         P = PVAL( IM )
         N = NVAL( IM )
*
         DO 20 IMAT = 1, NTYPES
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 20
*
*           Set up parameters with SLATB9 and generate test
*           matrices A and B with SLATMS.
*
            CALL SLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB,
     $                   ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB,
     $                   DISTA, DISTB )
*
*           Generate M by N matrix A
*
            CALL SLATMS( M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA,
     $                   ANORM, KLA, KUA, 'No packing', A, LDA, WORK,
     $                   IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 20
            END IF
*
            CALL SLATMS( P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB,
     $                   BNORM, KLB, KUB, 'No packing', B, LDB, WORK,
     $                   IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 20
            END IF
*
            NT = 6
*
            CALL SGSVTS( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V,
     $                   LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK,
     $                   LWORK, RWORK, RESULT )
*
*           Print information about the tests that did not
*           pass the threshold.
*
            DO 10 I = 1, NT
               IF( RESULT( I ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                     FIRSTT = .FALSE.
                     CALL ALAHDG( NOUT, PATH )
                  END IF
                  WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I,
     $               RESULT( I )
                  NFAIL = NFAIL + 1
               END IF
   10       CONTINUE
            NRUN = NRUN + NT
   20    CONTINUE
   30 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
*
 9999 FORMAT( ' SLATMS in SCKGSV   INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
      RETURN
*
*     End of SCKGSV
*
      END
