      SUBROUTINE ZTIMSY( LINE, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NLDA,
     $                   LDAVAL, TIMMIN, A, B, WORK, IWORK, RESLTS,
     $                   LDR1, LDR2, LDR3, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, LDR3, NLDA, NN, NNB, NNS, NOUT
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), LDAVAL( * ), NBVAL( * ),
     $                   NSVAL( * ), NVAL( * )
      DOUBLE PRECISION   RESLTS( LDR1, LDR2, LDR3, * )
      COMPLEX*16         A( * ), B( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZTIMSY times ZSYTRF, -TRS, and -TRI.
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
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix size N.
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
*  A       (workspace) DOUBLE PRECISION array, dimension (LDAMAX*NMAX)
*          where LDAMAX and NMAX are the maximum values permitted
*          for LDA and N.
*
*  B       (workspace) DOUBLE PRECISION array, dimension (LDAMAX*NMAX)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (NMAX)
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  RESLTS  (output) DOUBLE PRECISION array, dimension
*                   (LDR1,LDR2,LDR3,NSUBS)
*          The timing results for each subroutine over the relevant
*          values of N, NB, and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(4,NNB).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NN).
*
*  LDR3    (input) INTEGER
*          The third dimension of RESLTS.  LDR3 >= max(1,2*NLDA).
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
      CHARACTER          UPLO
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, I3, IC, ICL, ILDA, IN, INB, INFO, ISUB,
     $                   IUPLO, LDA, LDB, LWORK, MAT, N, NB, NRHS
      DOUBLE PRECISION   OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER          UPLOS( 2 )
      CHARACTER(32)      SUBNAM( NSUBS )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      LOGICAL            LSAME
      DOUBLE PRECISION   DMFLOP, DOPLA, DSECND
      EXTERNAL           LSAME, DMFLOP, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, DPRTBL, XLAENV, ZLACPY, ZSYTRF,
     $                   ZSYTRI, ZSYTRS, ZTIMMG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Data statements ..
      DATA               UPLOS / 'U', 'L' /
      DATA               SUBNAM / 'ZSYTRF', 'ZSYTRS', 'ZSYTRI' /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'SY'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 150
*
*     Check that N <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 2, CNAME, NN, NVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 150
      END IF
*
*     Do first for UPLO = 'U', then for UPLO = 'L'
*
      DO 110 IUPLO = 1, 2
         UPLO = UPLOS( IUPLO )
         IF( LSAME( UPLO, 'U' ) ) THEN
            MAT = 8
         ELSE
            MAT = -8
         END IF
*
*        Do for each value of N in NVAL.
*
         DO 100 IN = 1, NN
            N = NVAL( IN )
*
*           Do for each value of LDA:
*
            DO 90 ILDA = 1, NLDA
               LDA = LDAVAL( ILDA )
               I3 = ( IUPLO-1 )*NLDA + ILDA
*
*              Do for each value of NB in NBVAL.  Only the blocked
*              routines are timed in this loop since the other routines
*              are independent of NB.
*
               IF( TIMSUB( 1 ) ) THEN
                  DO 30 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
                     LWORK = MAX( 2*N, NB*N )
*
*                    Time ZSYTRF
*
                     CALL ZTIMMG( MAT, N, N, A, LDA, 0, 0 )
                     IC = 0
                     S1 = DSECND( )
   10                CONTINUE
                     CALL ZSYTRF( UPLO, N, A, LDA, IWORK, B, LWORK,
     $                            INFO )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL ZTIMMG( MAT, N, N, A, LDA, 0, 0 )
                        GO TO 10
                     END IF
*
*                    Subtract the time used in ZTIMMG.
*
                     ICL = 1
                     S1 = DSECND( )
   20                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL ZTIMMG( MAT, N, N, B, LDA, 0, 0 )
                        GO TO 20
                     END IF
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
                     OPS = DOPLA( 'ZSYTRF', N, N, 0, 0, NB )
                     RESLTS( INB, IN, I3, 1 ) = DMFLOP( OPS, TIME,
     $                  INFO )
*
   30             CONTINUE
               ELSE
*
*                 If ZSYTRF was not timed, generate a matrix and
*                 factor it using ZSYTRF anyway so that the factored
*                 form of the matrix can be used in timing the other
*                 routines.
*
                  CALL ZTIMMG( MAT, N, N, A, LDA, 0, 0 )
                  NB = 1
                  CALL XLAENV( 1, NB )
                  CALL ZSYTRF( UPLO, N, A, LDA, IWORK, B, LWORK, INFO )
               END IF
*
*              Time ZSYTRI
*
               IF( TIMSUB( 3 ) ) THEN
                  CALL ZLACPY( UPLO, N, N, A, LDA, B, LDA )
                  IC = 0
                  S1 = DSECND( )
   40             CONTINUE
                  CALL ZSYTRI( UPLO, N, B, LDA, IWORK, WORK, INFO )
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL ZLACPY( UPLO, N, N, A, LDA, B, LDA )
                     GO TO 40
                  END IF
*
*                 Subtract the time used in ZLACPY.
*
                  ICL = 1
                  S1 = DSECND( )
   50             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL ZLACPY( UPLO, N, N, A, LDA, B, LDA )
                     GO TO 50
                  END IF
*
                  TIME = ( TIME-UNTIME ) / DBLE( IC )
                  OPS = DOPLA( 'ZSYTRI', N, N, 0, 0, 0 )
                  RESLTS( 1, IN, I3, 3 ) = DMFLOP( OPS, TIME, INFO )
               END IF
*
*              Time ZSYTRS
*
               IF( TIMSUB( 2 ) ) THEN
                  DO 80 I = 1, NNS
                     NRHS = NSVAL( I )
                     LDB = LDA
                     CALL ZTIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                     IC = 0
                     S1 = DSECND( )
   60                CONTINUE
                     CALL ZSYTRS( UPLO, N, NRHS, A, LDA, IWORK, B, LDB,
     $                            INFO )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL ZTIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                        GO TO 60
                     END IF
*
*                    Subtract the time used in ZTIMMG.
*
                     ICL = 1
                     S1 = DSECND( )
   70                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL ZTIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                        GO TO 70
                     END IF
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
                     OPS = DOPLA( 'ZSYTRS', N, NRHS, 0, 0, 0 )
                     RESLTS( I, IN, I3, 2 ) = DMFLOP( OPS, TIME, INFO )
   80             CONTINUE
               END IF
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
*
*     Print tables of results for each timed routine.
*
      DO 140 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 140
         WRITE( NOUT, FMT = 9998 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
         IF( NLDA.GT.1 ) THEN
            DO 120 I = 1, NLDA
               WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
  120       CONTINUE
         END IF
         WRITE( NOUT, FMT = * )
         DO 130 IUPLO = 1, 2
            WRITE( NOUT, FMT = 9996 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) )),
     $     UPLOS( IUPLO )
            I3 = ( IUPLO-1 )*NLDA + 1
            IF( ISUB.EQ.1 ) THEN
               CALL DPRTBL( 'NB', 'N', NNB, NBVAL, NN, NVAL, NLDA,
     $                      RESLTS( 1, 1, I3, 1 ), LDR1, LDR2, NOUT )
            ELSE IF( ISUB.EQ.2 ) THEN
               CALL DPRTBL( 'NRHS', 'N', NNS, NSVAL, NN, NVAL, NLDA,
     $                      RESLTS( 1, 1, I3, 2 ), LDR1, LDR2, NOUT )
            ELSE IF( ISUB.EQ.3 ) THEN
               CALL DPRTBL( ' ', 'N', 1, NBVAL, NN, NVAL, NLDA,
     $                      RESLTS( 1, 1, I3, 3 ), LDR1, LDR2, NOUT )
            END IF
  130    CONTINUE
  140 CONTINUE
*
  150 CONTINUE
 9999 FORMAT( 1X, A, ' timing run not attempted' )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
 9996 FORMAT( 5X, A, ' with UPLO = ''', A1, '''', / )
      RETURN
*
*     End of ZTIMSY
*
      END
