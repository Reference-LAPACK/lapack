      SUBROUTINE ZTIMGE( LINE, NM, MVAL, NNS, NSVAL, NNB, NBVAL, NLDA,
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
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), LDAVAL( * ), MVAL( * ), NBVAL( * ),
     $                   NSVAL( * )
      DOUBLE PRECISION   RESLTS( LDR1, LDR2, LDR3, * )
      COMPLEX*16         A( * ), B( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZTIMGE times ZGETRF, -TRS, and -TRI.
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
*  TIMMIN  (input) DOUBLE PRECISION
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) COMPLEX*16 array, dimension (LDAMAX*NMAX)
*          where LDAMAX and NMAX are the maximum values permitted
*          for LDA and N.
*
*  B       (workspace) COMPLEX*16 array, dimension (LDAMAX*NMAX)
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LDAMAX*NBMAX)
*          where NBMAX is the maximum value of the block size NB.
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  RESLTS  (output) DOUBLE PRECISION array, dimension
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
      DOUBLE PRECISION   OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER(32)      SUBNAM( NSUBS )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      DOUBLE PRECISION   DMFLOP, DOPLA, DSECND
      EXTERNAL           DMFLOP, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, DPRTBL, XLAENV, ZGETRF, ZGETRI,
     $                   ZGETRS, ZLACPY, ZTIMMG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'ZGETRF', 'ZGETRS', 'ZGETRI' /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Zomplex precision'
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
*              Time ZGETRF
*
               IF( TIMSUB( 1 ) ) THEN
                  CALL ZTIMMG( 1, M, N, A, LDA, 0, 0 )
                  IC = 0
                  S1 = DSECND( )
   10             CONTINUE
                  CALL ZGETRF( M, N, A, LDA, IWORK, INFO )
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL ZTIMMG( 1, M, N, A, LDA, 0, 0 )
                     GO TO 10
                  END IF
*
*                 Subtract the time used in ZTIMMG.
*
                  ICL = 1
                  S1 = DSECND( )
   20             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL ZTIMMG( 1, M, N, A, LDA, 0, 0 )
                     GO TO 20
                  END IF
*
                  TIME = ( TIME-UNTIME ) / DBLE( IC )
                  OPS = DOPLA( 'ZGETRF', M, N, 0, 0, NB )
                  RESLTS( INB, IM, ILDA, 1 ) = DMFLOP( OPS, TIME, INFO )
*
               ELSE
                  IC = 0
                  CALL ZTIMMG( 1, M, N, A, LDA, 0, 0 )
               END IF
*
*              Generate another matrix and factor it using ZGETRF so
*              that the factored form can be used in timing the other
*              routines.
*
               IF( IC.NE.1 )
     $            CALL ZGETRF( M, N, A, LDA, IWORK, INFO )
*
*              Time ZGETRI
*
               IF( TIMSUB( 3 ) ) THEN
                  CALL ZLACPY( 'Full', M, M, A, LDA, B, LDA )
                  IC = 0
                  S1 = DSECND( )
   30             CONTINUE
                  CALL ZGETRI( M, B, LDA, IWORK, WORK, LDA*NB, INFO )
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL ZLACPY( 'Full', M, M, A, LDA, B, LDA )
                     GO TO 30
                  END IF
*
*                 Subtract the time used in ZLACPY.
*
                  ICL = 1
                  S1 = DSECND( )
   40             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL ZLACPY( 'Full', M, M, A, LDA, B, LDA )
                     GO TO 40
                  END IF
*
                  TIME = ( TIME-UNTIME ) / DBLE( IC )
                  OPS = DOPLA( 'ZGETRI', M, M, 0, 0, NB )
                  RESLTS( INB, IM, ILDA, 3 ) = DMFLOP( OPS, TIME, INFO )
               END IF
   50       CONTINUE
*
*           Time ZGETRS
*
            IF( TIMSUB( 2 ) ) THEN
               DO 80 I = 1, NNS
                  NRHS = NSVAL( I )
                  LDB = LDA
                  CALL ZTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                  IC = 0
                  S1 = DSECND( )
   60             CONTINUE
                  CALL ZGETRS( 'No transpose', M, NRHS, A, LDA, IWORK,
     $                         B, LDB, INFO )
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL ZTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 60
                  END IF
*
*                 Subtract the time used in ZTIMMG.
*
                  ICL = 1
                  S1 = DSECND( )
   70             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL ZTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 70
                  END IF
*
                  TIME = ( TIME-UNTIME ) / DBLE( IC )
                  OPS = DOPLA( 'ZGETRS', M, NRHS, 0, 0, 0 )
                  RESLTS( I, IM, ILDA, 2 ) = DMFLOP( OPS, TIME, INFO )
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
            CALL DPRTBL( 'NB', 'N', NNB, NBVAL, NM, MVAL, NLDA, RESLTS,
     $                   LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.2 ) THEN
            CALL DPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 2 ), LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.3 ) THEN
            CALL DPRTBL( 'NB', 'N', NNB, NBVAL, NM, MVAL, NLDA,
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
*     End of ZTIMGE
*
      END
