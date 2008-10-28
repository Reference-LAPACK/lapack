      SUBROUTINE STIMGE( LINE, NM, MVAL, NNS, NSVAL, NNB, NBVAL, NLDA,
     $                   LDAVAL, TIMMIN, A, B, WORK, IWORK, RESLTS,
     $                   LDR1, LDR2, LDR3, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, LDR3, NLDA, NM, NNB, NNS, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), LDAVAL( * ), MVAL( * ), NBVAL( * ),
     $                   NSVAL( * )
      REAL               A( * ), B( * ), RESLTS( LDR1, LDR2, LDR3, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  STIMGE times SGETRF, -TRS, and -TRI.
*
*  Arguments
*  =========
*
*  LINE    (input) CHARACTER*80
*          The input line that requested this routine.  The first six
*          characters contain either the name of a subroutine or a
*          generic path name.  The remaining characters may be used to
*          specify the individual routines to be timed.  See ATIMIN for
*          a full description of the format of the input line.
*
*  NM      (input) INTEGER
*          The number of values of M contained in the vector MVAL.
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix size M.
*
*  NNS     (input) INTEGER
*          The number of values of NRHS contained in the vector NSVAL.
*
*  NSVAL   (input) INTEGER array, dimension (NNS)
*          The values of the number of right hand sides NRHS.
*
*  NNB     (input) INTEGER
*          The number of values of NB contained in the vector NBVAL.
*
*  NBVAL   (input) INTEGER array, dimension (NNB)
*          The values of the blocksize NB.
*
*  NLDA    (input) INTEGER
*          The number of values of LDA contained in the vector LDAVAL.
*
*  LDAVAL  (input) INTEGER array, dimension (NLDA)
*          The values of the leading dimension of the array A.
*
*  TIMMIN  (input) REAL
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) REAL array, dimension (LDAMAX*NMAX)
*          where LDAMAX and NMAX are the maximum values permitted
*          for LDA and N.
*
*  B       (workspace) REAL array, dimension (LDAMAX*NMAX)
*
*  WORK    (workspace) REAL array, dimension (LDAMAX*NBMAX)
*          where NBMAX is the maximum value of the block size NB.
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  RESLTS  (output) REAL array, dimension
*                   (LDR1,LDR2,LDR3,NSUBS)
*          The timing results for each subroutine over the relevant
*          values of N and NB.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(4,NNB).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NM).
*
*  LDR3    (input) INTEGER
*          The third dimension of RESLTS.  LDR3 >= max(1,NLDA).
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 3 )
*     ..
*     .. Local Scalars ..
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, ILDA, IM, INB, INFO, ISUB, LDA,
     $                   LDB, M, N, NB, NRHS
      REAL               OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER(32)      SUBNAM( NSUBS )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      REAL               SECOND, SMFLOP, SOPLA
      EXTERNAL           SECOND, SMFLOP, SOPLA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, SGETRF, SGETRI, SGETRS, SLACPY,
     $                   SPRTBL, STIMMG, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'SGETRF', 'SGETRS', 'SGETRI' /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'GE'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 130
*
*     Check that N <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 2, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 130
      END IF
*
*     Do for each value of M:
*
      DO 100 IM = 1, NM
*
         M = MVAL( IM )
         N = M
*
*        Do for each value of LDA:
*
         DO 90 ILDA = 1, NLDA
            LDA = LDAVAL( ILDA )
*
*           Do for each value of NB in NBVAL.  Only the blocked
*           routines are timed in this loop since the other routines
*           are independent of NB.
*
            DO 50 INB = 1, NNB
               NB = NBVAL( INB )
               CALL XLAENV( 1, NB )
*
*              Time SGETRF
*
               IF( TIMSUB( 1 ) ) THEN
                  CALL STIMMG( 1, M, N, A, LDA, 0, 0 )
                  IC = 0
                  S1 = SECOND( )
   10             CONTINUE
                  CALL SGETRF( M, N, A, LDA, IWORK, INFO )
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL STIMMG( 1, M, N, A, LDA, 0, 0 )
                     GO TO 10
                  END IF
*
*                 Subtract the time used in STIMMG.
*
                  ICL = 1
                  S1 = SECOND( )
   20             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL STIMMG( 1, M, N, A, LDA, 0, 0 )
                     GO TO 20
                  END IF
*
                  TIME = ( TIME-UNTIME ) / REAL( IC )
                  OPS = SOPLA( 'SGETRF', M, N, 0, 0, NB )
                  RESLTS( INB, IM, ILDA, 1 ) = SMFLOP( OPS, TIME, INFO )
*
               ELSE
                  IC = 0
                  CALL STIMMG( 1, M, N, A, LDA, 0, 0 )
               END IF
*
*              Generate another matrix and factor it using SGETRF so
*              that the factored form can be used in timing the other
*              routines.
*
               IF( IC.NE.1 )
     $            CALL SGETRF( M, N, A, LDA, IWORK, INFO )
*
*              Time SGETRI
*
               IF( TIMSUB( 3 ) ) THEN
                  CALL SLACPY( 'Full', M, M, A, LDA, B, LDA )
                  IC = 0
                  S1 = SECOND( )
   30             CONTINUE
                  CALL SGETRI( M, B, LDA, IWORK, WORK, LDA*NB, INFO )
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL SLACPY( 'Full', M, M, A, LDA, B, LDA )
                     GO TO 30
                  END IF
*
*                 Subtract the time used in SLACPY.
*
                  ICL = 1
                  S1 = SECOND( )
   40             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL SLACPY( 'Full', M, M, A, LDA, B, LDA )
                     GO TO 40
                  END IF
*
                  TIME = ( TIME-UNTIME ) / REAL( IC )
                  OPS = SOPLA( 'SGETRI', M, M, 0, 0, NB )
                  RESLTS( INB, IM, ILDA, 3 ) = SMFLOP( OPS, TIME, INFO )
               END IF
   50       CONTINUE
*
*           Time SGETRS
*
            IF( TIMSUB( 2 ) ) THEN
               DO 80 I = 1, NNS
                  NRHS = NSVAL( I )
                  LDB = LDA
                  CALL STIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                  IC = 0
                  S1 = SECOND( )
   60             CONTINUE
                  CALL SGETRS( 'No transpose', M, NRHS, A, LDA, IWORK,
     $                         B, LDB, INFO )
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL STIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 60
                  END IF
*
*                 Subtract the time used in STIMMG.
*
                  ICL = 1
                  S1 = SECOND( )
   70             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL STIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 70
                  END IF
*
                  TIME = ( TIME-UNTIME ) / REAL( IC )
                  OPS = SOPLA( 'SGETRS', M, NRHS, 0, 0, 0 )
                  RESLTS( I, IM, ILDA, 2 ) = SMFLOP( OPS, TIME, INFO )
   80          CONTINUE
            END IF
   90    CONTINUE
  100 CONTINUE
*
*     Print a table of results for each timed routine.
*
      DO 120 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 120
         WRITE( NOUT, FMT = 9998 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
         IF( NLDA.GT.1 ) THEN
            DO 110 I = 1, NLDA
               WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
  110       CONTINUE
         END IF
         WRITE( NOUT, FMT = * )
         IF( ISUB.EQ.1 ) THEN
            CALL SPRTBL( 'NB', 'N', NNB, NBVAL, NM, MVAL, NLDA, RESLTS,
     $                   LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.2 ) THEN
            CALL SPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 2 ), LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.3 ) THEN
            CALL SPRTBL( 'NB', 'N', NNB, NBVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 3 ), LDR1, LDR2, NOUT )
         END IF
  120 CONTINUE
*
  130 CONTINUE
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
      RETURN
*
*     End of STIMGE
*
      END
