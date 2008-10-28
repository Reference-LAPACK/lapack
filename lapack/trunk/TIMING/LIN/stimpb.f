      SUBROUTINE STIMPB( LINE, NN, NVAL, NK, KVAL, NNS, NSVAL, NNB,
     $                   NBVAL, NLDA, LDAVAL, TIMMIN, A, B, IWORK,
     $                   RESLTS, LDR1, LDR2, LDR3, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, LDR3, NK, NLDA, NN, NNB, NNS, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), KVAL( * ), LDAVAL( * ), NBVAL( * ),
     $                   NSVAL( * ), NVAL( * )
      REAL               A( * ), B( * ), RESLTS( LDR1, LDR2, LDR3, * )
*     ..
*
*  Purpose
*  =======
*
*  STIMPB times SPBTRF and -TRS.
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
*  NK      (input) INTEGER
*          The number of values of K contained in the vector KVAL.
*
*  KVAL    (input) INTEGER array, dimension (NK)
*          The values of the band width K.
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
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  RESLTS  (output) REAL array, dimension
*                   (LDR1,LDR2,LDR3,NSUBS)
*          The timing results for each subroutine over the relevant
*          values of N, K, NB, and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(4,NNB).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NK).
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
      PARAMETER          ( NSUBS = 2 )
*     ..
*     .. Local Scalars ..
      CHARACTER          UPLO
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, I3, IC, ICL, IK, ILDA, IN, INB, INFO, ISUB,
     $                   IUPLO, K, LDA, LDB, MAT, N, NB, NRHS
      REAL               OPS, S1, S2, TIME, UNTIME
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
      REAL               SECOND, SMFLOP, SOPLA
      EXTERNAL           LSAME, SECOND, SMFLOP, SOPLA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, SPBTRF, SPBTRS, SPRTBL, STIMMG,
     $                   XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Data statements ..
      DATA               UPLOS / 'U', 'L' /
      DATA               SUBNAM / 'SPBTRF', 'SPBTRS' /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'PB'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 140
*
*     Check that K+1 <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 0, CNAME, NK, KVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 140
      END IF
*
*     Do for each value of the matrix size N:
*
      DO 130 IN = 1, NN
         N = NVAL( IN )
*
*        Do first for UPLO = 'U', then for UPLO = 'L'
*
         DO 90 IUPLO = 1, 2
            UPLO = UPLOS( IUPLO )
            IF( LSAME( UPLO, 'U' ) ) THEN
               MAT = 5
            ELSE
               MAT = -5
            END IF
*
*           Do for each value of LDA:
*
            DO 80 ILDA = 1, NLDA
               LDA = LDAVAL( ILDA )
               I3 = ( IUPLO-1 )*NLDA + ILDA
*
*              Do for each value of the band width K:
*
               DO 70 IK = 1, NK
                  K = KVAL( IK )
                  K = MAX( 0, MIN( K, N-1 ) )
*
*                 Time SPBTRF
*
                  IF( TIMSUB( 1 ) ) THEN
*
*                    Do for each value of NB in NBVAL.  Only SPBTRF is
*                    timed in this loop since the other routines are
*                    independent of NB.
*
                     DO 30 INB = 1, NNB
                        NB = NBVAL( INB )
                        CALL XLAENV( 1, NB )
                        CALL STIMMG( MAT, N, N, A, LDA, K, K )
                        IC = 0
                        S1 = SECOND( )
   10                   CONTINUE
                        CALL SPBTRF( UPLO, N, K, A, LDA, INFO )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL STIMMG( MAT, N, N, A, LDA, K, K )
                           GO TO 10
                        END IF
*
*                       Subtract the time used in STIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
   20                   CONTINUE
                        CALL STIMMG( MAT, N, N, A, LDA, K, K )
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC )
     $                     GO TO 20
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPLA( 'SPBTRF', N, N, K, K, NB )
                        RESLTS( INB, IK, I3, 1 ) = SMFLOP( OPS, TIME,
     $                     INFO )
   30                CONTINUE
                  ELSE
                     IC = 0
                     CALL STIMMG( MAT, N, N, A, LDA, K, K )
                  END IF
*
*                 Generate another matrix and factor it using SPBTRF so
*                 that the factored form can be used in timing the other
*                 routines.
*
                  NB = 1
                  CALL XLAENV( 1, NB )
                  IF( IC.NE.1 )
     $               CALL SPBTRF( UPLO, N, K, A, LDA, INFO )
*
*                 Time SPBTRS
*
                  IF( TIMSUB( 2 ) ) THEN
                     DO 60 I = 1, NNS
                        NRHS = NSVAL( I )
                        LDB = N
                        CALL STIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
   40                   CONTINUE
                        CALL SPBTRS( UPLO, N, K, NRHS, A, LDA, B, LDB,
     $                               INFO )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL STIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                           GO TO 40
                        END IF
*
*                       Subtract the time used in STIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
   50                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL STIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                           GO TO 50
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPLA( 'SPBTRS', N, NRHS, K, K, 0 )
                        RESLTS( I, IK, I3, 2 ) = SMFLOP( OPS, TIME,
     $                     INFO )
   60                CONTINUE
                  END IF
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
*
*        Print tables of results for each timed routine.
*
         DO 120 ISUB = 1, NSUBS
            IF( .NOT.TIMSUB( ISUB ) )
     $         GO TO 120
*
*           Print header for routine names.
*
            IF( IN.EQ.1 .OR. CNAME.EQ.'SPB   ' ) THEN
               WRITE( NOUT, FMT = 9998 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
               IF( NLDA.GT.1 ) THEN
                  DO 100 I = 1, NLDA
                     WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
  100             CONTINUE
               END IF
            END IF
            WRITE( NOUT, FMT = * )
            DO 110 IUPLO = 1, 2
               WRITE( NOUT, FMT = 9996 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) )), N,
     $            UPLOS( IUPLO )
               I3 = ( IUPLO-1 )*NLDA + 1
               IF( ISUB.EQ.1 ) THEN
                  CALL SPRTBL( 'NB', 'K', NNB, NBVAL, NK, KVAL, NLDA,
     $                         RESLTS( 1, 1, I3, 1 ), LDR1, LDR2, NOUT )
               ELSE IF( ISUB.EQ.2 ) THEN
                  CALL SPRTBL( 'NRHS', 'K', NNS, NSVAL, NK, KVAL, NLDA,
     $                         RESLTS( 1, 1, I3, 2 ), LDR1, LDR2, NOUT )
               END IF
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
*
  140 CONTINUE
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
 9996 FORMAT( 5X, A, ' with M =', I6, ', UPLO = ''', A1, '''', / )
      RETURN
*
*     End of STIMPB
*
      END
