      SUBROUTINE DTIMB2( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NINC,
     $                   INCVAL, NLDA, LDAVAL, LA, TIMMIN, A, X, Y,
     $                   RESLTS, LDR1, LDR2, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    LINE
      INTEGER            LA, LDR1, LDR2, NINC, NK, NLDA, NM, NN, NOUT
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            INCVAL( * ), KVAL( * ), LDAVAL( * ), MVAL( * ),
     $                   NVAL( * )
      DOUBLE PRECISION   A( * ), RESLTS( LDR1, LDR2, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DTIMB2 times the BLAS 2 routines.
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
*          The values of the matrix row dimension M.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension N.
*
*  NK      (input) INTEGER
*          The number of values of K contained in the vector KVAL.
*
*  KVAL    (input) INTEGER array, dimension (NK)
*          The values of the band width K.
*
*  NINC    (input) INTEGER
*          The number of values of INCX contained in the vector INCVAL.
*
*  INCVAL  (input) INTEGER array, dimension (NINC)
*          The values of INCX, the increment between successive values
*          of the vector X.
*
*  NLDA    (input) INTEGER
*          The number of values of LDA contained in the vector LDAVAL.
*
*  LDAVAL  (input) INTEGER array, dimension (NLDA)
*          The values of the leading dimension of the array A.
*
*  LA      (input) INTEGER
*          The size of the array A.
*
*  TIMMIN  (input) DOUBLE PRECISION
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) DOUBLE PRECISION array, dimension (LA)
*
*  X       (workspace) DOUBLE PRECISION array, dimension (NMAX*INCMAX)
*             where NMAX and INCMAX are the maximum values permitted
*             for N and INCX.
*
*  Y       (workspace) DOUBLE PRECISION array, dimension (NMAX*INCMAX)
*             where NMAX and INCMAX are the maximum values permitted
*             for N and INCX.
*
*  RESLTS  (output) DOUBLE PRECISION array, dimension (LDR1,LDR2,p),
*             where p = NLDA*NINC.
*          The timing results for each subroutine over the relevant
*          values of M, N, K, INCX, and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(1,NM,NK).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NN).
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 16 )
      INTEGER            NTRANS, NUPLOS
      PARAMETER          ( NTRANS = 2, NUPLOS = 2 )
      DOUBLE PRECISION   ALPHA, BETA
      PARAMETER          ( ALPHA = 1.0D0, BETA = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            IXANDY
      CHARACTER          TRANSA, UPLO
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, I3, IC, ICL, IINC, IK, ILDA, IM, IMAT, IN,
     $                   INCX, INFO, ISUB, ITA, IUPLO, J, K, LDA, M, N,
     $                   NX, NY
      DOUBLE PRECISION   OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER          TRANS( NTRANS ), UPLOS( NUPLOS )
      CHARACTER(32)      NAMES( NSUBS )
      INTEGER            LAVAL( 1 )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      DOUBLE PRECISION   DMFLOP, DOPBL2, DSECND
      EXTERNAL           DMFLOP, DOPBL2, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, DGBMV, DGEMV, DGER, DPRTBL,
     $                   DSBMV, DSPMV, DSPR, DSPR2, DSYMV, DSYR, DSYR2,
     $                   DTBMV, DTBSV, DTIMMG, DTPMV, DTPSV, DTRMV,
     $                   DTRSV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Data statements ..
      DATA               TRANS / 'N', 'T' /
      DATA               UPLOS / 'U', 'L' /
      DATA               NAMES / 'DGEMV ', 'DGBMV ', 'DSYMV ', 'DSBMV ',
     $                   'DSPMV ', 'DTRMV ', 'DTBMV ', 'DTPMV ',
     $                   'DTRSV ', 'DTBSV ', 'DTPSV ', 'DGER  ',
     $                   'DSYR  ', 'DSPR  ', 'DSYR2 ', 'DSPR2 ' /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'B2'
      CALL ATIMIN( PATH, LINE, NSUBS, NAMES, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 1070
*
*     Time each routine
*
      DO 1060 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 1060
*
*        Check the input values.  The conditions are
*           M <= LDA for general storage
*           K <= LDA for banded storage
*           N*(N+1)/2 <= LA  for packed storage
*
         CNAME = NAMES( ISUB )
         IF( CNAME( 2: 3 ).EQ.'GE' ) THEN
            CALL ATIMCK( 1, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
         ELSE IF( CNAME( 3: 3 ).EQ.'B' ) THEN
            CALL ATIMCK( 0, CNAME, NK, KVAL, NLDA, LDAVAL, NOUT, INFO )
         ELSE IF( CNAME( 3: 3 ).EQ.'P' ) THEN
            LAVAL( 1 ) = LA
            CALL ATIMCK( 4, CNAME, NN, NVAL, 1, LAVAL, NOUT, INFO )
         ELSE
            CALL ATIMCK( 2, CNAME, NN, NVAL, NLDA, LDAVAL, NOUT, INFO )
         END IF
         IF( INFO.GT.0 ) THEN
            WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
            GO TO 1060
         END IF
*
*        Print header.
*
         WRITE( NOUT, FMT = 9998 )CNAME(1:ILA_LEN_TRIM(CNAME))
         IXANDY = ISUB.LE.5 .OR. ISUB.EQ.12 .OR. ISUB.EQ.15 .OR.
     $            ISUB.EQ.16
         IF( CNAME( 3: 3 ).NE.'P' ) THEN
            IF( NLDA*NINC.EQ.1 ) THEN
               IF( IXANDY ) THEN
                  WRITE( NOUT, FMT = 9997 )LDAVAL( 1 ), INCVAL( 1 )
               ELSE
                  WRITE( NOUT, FMT = 9996 )LDAVAL( 1 ), INCVAL( 1 )
               END IF
            ELSE
               DO 20 I = 1, NLDA
                  DO 10 J = 1, NINC
                     IF( IXANDY ) THEN
                        WRITE( NOUT, FMT = 9993 )( I-1 )*NINC + J,
     $                     LDAVAL( I ), INCVAL( J )
                     ELSE
                        WRITE( NOUT, FMT = 9992 )( I-1 )*NINC + J,
     $                     LDAVAL( I ), INCVAL( J )
                     END IF
   10             CONTINUE
   20          CONTINUE
            END IF
         ELSE
            IF( NINC.EQ.1 ) THEN
               IF( IXANDY ) THEN
                  WRITE( NOUT, FMT = 9995 )INCVAL( 1 )
               ELSE
                  WRITE( NOUT, FMT = 9994 )INCVAL( 1 )
               END IF
            ELSE
               DO 30 J = 1, NINC
                  IF( IXANDY ) THEN
                     WRITE( NOUT, FMT = 9991 )J, INCVAL( J )
                  ELSE
                     WRITE( NOUT, FMT = 9990 )J, INCVAL( J )
                  END IF
   30          CONTINUE
            END IF
         END IF
*
*        Time DGEMV
*
         IF( CNAME.EQ.'DGEMV ' ) THEN
            DO 100 ITA = 1, NTRANS
               TRANSA = TRANS( ITA )
               I3 = 0
               DO 90 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 80 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 70 IM = 1, NM
                        M = MVAL( IM )
                        DO 60 IN = 1, NN
                           N = NVAL( IN )
                           IF( TRANSA.EQ.'N' ) THEN
                              NX = N
                              NY = M
                           ELSE
                              NX = M
                              NY = N
                           END IF
                           CALL DTIMMG( 1, M, N, A, LDA, 0, 0 )
                           CALL DTIMMG( 0, 1, NX, X, INCX, 0, 0 )
                           CALL DTIMMG( 0, 1, NY, Y, INCX, 0, 0 )
                           IC = 0
                           S1 = DSECND( )
   40                      CONTINUE
                           CALL DGEMV( TRANSA, M, N, ALPHA, A, LDA, X,
     $                                 INCX, BETA, Y, INCX )
                           S2 = DSECND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL DTIMMG( 0, 1, NY, Y, INCX, 0, 0 )
                              GO TO 40
                           END IF
*
*                          Subtract the time used in DTIMMG.
*
                           ICL = 1
                           S1 = DSECND( )
   50                      CONTINUE
                           S2 = DSECND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL DTIMMG( 0, 1, NY, Y, INCX, 0, 0 )
                              GO TO 50
                           END IF
*
                           TIME = ( TIME-UNTIME ) / DBLE( IC )
                           OPS = DOPBL2( CNAME, M, N, 0, 0 )
                           RESLTS( IM, IN, I3 ) = DMFLOP( OPS, TIME, 0 )
   60                   CONTINUE
   70                CONTINUE
   80             CONTINUE
   90          CONTINUE
               WRITE( NOUT, FMT = 9989 )TRANSA
               CALL DPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  100       CONTINUE
*
         ELSE IF( CNAME.EQ.'DGBMV ' ) THEN
*
*           Time DGBMV
*
            DO 170 ITA = 1, NTRANS
               TRANSA = TRANS( ITA )
               I3 = 0
               DO 160 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 150 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 140 IK = 1, NK
                        K = KVAL( IK )
                        DO 130 IN = 1, NN
                           N = NVAL( IN )
                           M = N
                           CALL DTIMMG( -2, M, N, A, LDA, K, K )
                           CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           CALL DTIMMG( 0, 1, M, Y, INCX, 0, 0 )
                           IC = 0
                           S1 = DSECND( )
  110                      CONTINUE
                           CALL DGBMV( TRANSA, M, N, K, K, ALPHA, A,
     $                                 LDA, X, INCX, BETA, Y, INCX )
                           S2 = DSECND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL DTIMMG( 0, 1, M, Y, INCX, 0, 0 )
                              GO TO 110
                           END IF
*
*                          Subtract the time used in DTIMMG.
*
                           ICL = 1
                           S1 = DSECND( )
  120                      CONTINUE
                           S2 = DSECND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL DTIMMG( 0, 1, M, Y, INCX, 0, 0 )
                              GO TO 120
                           END IF
*
                           TIME = ( TIME-UNTIME ) / DBLE( IC )
                           OPS = DOPBL2( CNAME, M, N, K, K )
                           RESLTS( IK, IN, I3 ) = DMFLOP( OPS, TIME, 0 )
  130                   CONTINUE
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
               WRITE( NOUT, FMT = 9988 )TRANSA
               CALL DPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  170       CONTINUE
*
         ELSE IF( CNAME.EQ.'DSYMV ' ) THEN
*
*           Time DSYMV
*
            DO 230 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 6
               IF( UPLO.EQ.'L' )
     $            IMAT = -6
               I3 = 0
               DO 220 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 210 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 200 IN = 1, NN
                        N = NVAL( IN )
                        CALL DTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                        CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        IC = 0
                        S1 = DSECND( )
  180                   CONTINUE
                        CALL DSYMV( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                              BETA, Y, INCX )
                        S2 = DSECND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                           GO TO 180
                        END IF
*
*                       Subtract the time used in DTIMMG.
*
                        ICL = 1
                        S1 = DSECND( )
  190                   CONTINUE
                        S2 = DSECND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                           GO TO 190
                        END IF
*
                        TIME = ( TIME-UNTIME ) / DBLE( IC )
                        OPS = DOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, I3 ) = DMFLOP( OPS, TIME, 0 )
  200                CONTINUE
  210             CONTINUE
  220          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  230       CONTINUE
*
         ELSE IF( CNAME.EQ.'DSBMV ' ) THEN
*
*           Time DSBMV
*
            DO 300 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 8
               IF( UPLO.EQ.'L' )
     $            IMAT = -8
               I3 = 0
               DO 290 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 280 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 270 IK = 1, NK
                        K = KVAL( IK )
                        DO 260 IN = 1, NN
                           N = NVAL( IN )
                           CALL DTIMMG( IMAT, N, N, A, LDA, K, K )
                           CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                           IC = 0
                           S1 = DSECND( )
  240                      CONTINUE
                           CALL DSBMV( UPLO, N, K, ALPHA, A, LDA, X,
     $                                 INCX, BETA, Y, INCX )
                           S2 = DSECND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                              GO TO 240
                           END IF
*
*                          Subtract the time used in DTIMMG.
*
                           ICL = 1
                           S1 = DSECND( )
  250                      CONTINUE
                           S2 = DSECND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                              GO TO 250
                           END IF
*
                           TIME = ( TIME-UNTIME ) / DBLE( IC )
                           OPS = DOPBL2( CNAME, N, N, K, K )
                           RESLTS( IK, IN, I3 ) = DMFLOP( OPS, TIME, 0 )
  260                   CONTINUE
  270                CONTINUE
  280             CONTINUE
  290          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL DPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  300       CONTINUE
*
         ELSE IF( CNAME.EQ.'DSPMV ' ) THEN
*
*           Time DSPMV
*
            DO 350 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 7
               IF( UPLO.EQ.'L' )
     $            IMAT = -7
               ILDA = 1
               LDA = LDAVAL( ILDA )
               DO 340 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  DO 330 IN = 1, NN
                     N = NVAL( IN )
                     CALL DTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0, 0 )
                     CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                     CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                     IC = 0
                     S1 = DSECND( )
  310                CONTINUE
                     CALL DSPMV( UPLO, N, ALPHA, A, X, INCX, BETA, Y,
     $                           INCX )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        GO TO 310
                     END IF
*
*                    Subtract the time used in DTIMMG.
*
                     ICL = 1
                     S1 = DSECND( )
  320                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        GO TO 320
                     END IF
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
                     OPS = DOPBL2( CNAME, N, N, 0, 0 )
                     RESLTS( 1, IN, IINC ) = DMFLOP( OPS, TIME, 0 )
  330             CONTINUE
  340          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC, RESLTS,
     $                      LDR1, LDR2, NOUT )
  350       CONTINUE
*
         ELSE IF( CNAME.EQ.'DTRMV ' ) THEN
*
*           Time DTRMV
*
            DO 420 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 9
               IF( UPLO.EQ.'L' )
     $            IMAT = -9
               DO 410 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  I3 = 0
                  DO 400 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 390 IINC = 1, NINC
                        INCX = INCVAL( IINC )
                        I3 = I3 + 1
                        DO 380 IN = 1, NN
                           N = NVAL( IN )
                           CALL DTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           IC = 0
                           S1 = DSECND( )
  360                      CONTINUE
                           CALL DTRMV( UPLO, TRANSA, 'Non-unit', N, A,
     $                                 LDA, X, INCX )
                           S2 = DSECND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              GO TO 360
                           END IF
*
*                          Subtract the time used in DTIMMG.
*
                           ICL = 1
                           S1 = DSECND( )
  370                      CONTINUE
                           S2 = DSECND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              GO TO 370
                           END IF
*
                           TIME = ( TIME-UNTIME ) / DBLE( IC )
                           OPS = DOPBL2( CNAME, N, N, 0, 0 )
                           RESLTS( 1, IN, I3 ) = DMFLOP( OPS, TIME, 0 )
  380                   CONTINUE
  390                CONTINUE
  400             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  410          CONTINUE
  420       CONTINUE
*
         ELSE IF( CNAME.EQ.'DTRSV ' ) THEN
*
*           Time DTRSV
*
            DO 490 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 9
               IF( UPLO.EQ.'L' )
     $            IMAT = -9
               DO 480 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  I3 = 0
                  DO 470 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 460 IINC = 1, NINC
                        INCX = INCVAL( IINC )
                        I3 = I3 + 1
                        DO 450 IN = 1, NN
                           N = NVAL( IN )
                           CALL DTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           IC = 0
                           S1 = DSECND( )
  430                      CONTINUE
                           CALL DTRSV( UPLO, TRANSA, 'Non-unit', N, A,
     $                                 LDA, X, INCX )
                           S2 = DSECND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              GO TO 430
                           END IF
*
*                          Subtract the time used in DTIMMG.
*
                           ICL = 1
                           S1 = DSECND( )
  440                      CONTINUE
                           S2 = DSECND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              GO TO 440
                           END IF
*
                           TIME = ( TIME-UNTIME ) / DBLE( IC )
                           OPS = DOPBL2( CNAME, N, N, 0, 0 )
                           RESLTS( 1, IN, I3 ) = DMFLOP( OPS, TIME, 0 )
  450                   CONTINUE
  460                CONTINUE
  470             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  480          CONTINUE
  490       CONTINUE
*
         ELSE IF( CNAME.EQ.'DTBMV ' ) THEN
*
*           Time DTBMV
*
            DO 570 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 11
               IF( UPLO.EQ.'L' )
     $            IMAT = -11
               DO 560 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  I3 = 0
                  DO 550 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 540 IINC = 1, NINC
                        INCX = INCVAL( IINC )
                        I3 = I3 + 1
                        DO 530 IK = 1, NK
                           K = KVAL( IK )
                           DO 520 IN = 1, NN
                              N = NVAL( IN )
                              CALL DTIMMG( IMAT, N, N, A, LDA, K, K )
                              CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
  500                         CONTINUE
                              CALL DTBMV( UPLO, TRANSA, 'Non-unit', N,
     $                                    K, A, LDA, X, INCX )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                                 GO TO 500
                              END IF
*
*                             Subtract the time used in DTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
  510                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                                 GO TO 510
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
                              OPS = DOPBL2( CNAME, N, N, K, K )
                              RESLTS( IK, IN, I3 ) = DMFLOP( OPS, TIME,
     $                           0 )
  520                      CONTINUE
  530                   CONTINUE
  540                CONTINUE
  550             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL DPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NINC*NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  560          CONTINUE
  570       CONTINUE
*
         ELSE IF( CNAME.EQ.'DTBSV ' ) THEN
*
*           Time DTBSV
*
            DO 650 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 11
               IF( UPLO.EQ.'L' )
     $            IMAT = -11
               DO 640 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  I3 = 0
                  DO 630 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 620 IINC = 1, NINC
                        INCX = INCVAL( IINC )
                        I3 = I3 + 1
                        DO 610 IK = 1, NK
                           K = KVAL( IK )
                           DO 600 IN = 1, NN
                              N = NVAL( IN )
                              CALL DTIMMG( IMAT, N, N, A, LDA, K, K )
                              CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
  580                         CONTINUE
                              CALL DTBSV( UPLO, TRANSA, 'Non-unit', N,
     $                                    K, A, LDA, X, INCX )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                                 GO TO 580
                              END IF
*
*                             Subtract the time used in DTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
  590                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                                 GO TO 590
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
                              OPS = DOPBL2( CNAME, N, N, K, K )
                              RESLTS( IK, IN, I3 ) = DMFLOP( OPS, TIME,
     $                           0 )
  600                      CONTINUE
  610                   CONTINUE
  620                CONTINUE
  630             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL DPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NINC*NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  640          CONTINUE
  650       CONTINUE
*
         ELSE IF( CNAME.EQ.'DTPMV ' ) THEN
*
*           Time DTPMV
*
            DO 710 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 10
               IF( UPLO.EQ.'L' )
     $            IMAT = -10
               DO 700 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  ILDA = 1
                  LDA = LDAVAL( ILDA )
                  DO 690 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     DO 680 IN = 1, NN
                        N = NVAL( IN )
                        CALL DTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        IC = 0
                        S1 = DSECND( )
  660                   CONTINUE
                        CALL DTPMV( UPLO, TRANSA, 'Non-unit', N, A, X,
     $                              INCX )
                        S2 = DSECND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           GO TO 660
                        END IF
*
*                       Subtract the time used in DTIMMG.
*
                        ICL = 1
                        S1 = DSECND( )
  670                   CONTINUE
                        S2 = DSECND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           GO TO 670
                        END IF
*
                        TIME = ( TIME-UNTIME ) / DBLE( IC )
                        OPS = DOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, IINC ) = DMFLOP( OPS, TIME, 0 )
  680                CONTINUE
  690             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC,
     $                         RESLTS, LDR1, LDR2, NOUT )
  700          CONTINUE
  710       CONTINUE
*
         ELSE IF( CNAME.EQ.'DTPSV ' ) THEN
*
*           Time DTPSV
*
            DO 770 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 10
               IF( UPLO.EQ.'L' )
     $            IMAT = -10
               DO 760 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  ILDA = 1
                  LDA = LDAVAL( ILDA )
                  DO 750 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     DO 740 IN = 1, NN
                        N = NVAL( IN )
                        CALL DTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        IC = 0
                        S1 = DSECND( )
  720                   CONTINUE
                        CALL DTPSV( UPLO, TRANSA, 'Non-unit', N, A, X,
     $                              INCX )
                        S2 = DSECND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           GO TO 720
                        END IF
*
*                       Subtract the time used in DTIMMG.
*
                        ICL = 1
                        S1 = DSECND( )
  730                   CONTINUE
                        S2 = DSECND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           GO TO 730
                        END IF
*
                        TIME = ( TIME-UNTIME ) / DBLE( IC )
                        OPS = DOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, IINC ) = DMFLOP( OPS, TIME, 0 )
  740                CONTINUE
  750             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC,
     $                         RESLTS, LDR1, LDR2, NOUT )
  760          CONTINUE
  770       CONTINUE
*
         ELSE IF( CNAME.EQ.'DGER  ' ) THEN
*
*           Time DGER
*
            I3 = 0
            DO 830 ILDA = 1, NLDA
               LDA = LDAVAL( ILDA )
               DO 820 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  I3 = I3 + 1
                  DO 810 IM = 1, NM
                     M = MVAL( IM )
                     DO 800 IN = 1, NN
                        N = NVAL( IN )
                        CALL DTIMMG( 0, 1, M, X, INCX, 0, 0 )
                        CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        CALL DTIMMG( 1, M, N, A, LDA, 0, 0 )
                        IC = 0
                        S1 = DSECND( )
  780                   CONTINUE
                        CALL DGER( M, N, ALPHA, X, INCX, Y, INCX, A,
     $                             LDA )
                        S2 = DSECND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL DTIMMG( 1, M, N, A, LDA, 0, 0 )
                           GO TO 780
                        END IF
*
*                       Subtract the time used in DTIMMG.
*
                        ICL = 1
                        S1 = DSECND( )
  790                   CONTINUE
                        S2 = DSECND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL DTIMMG( 1, M, N, A, LDA, 0, 0 )
                           GO TO 790
                        END IF
*
                        TIME = ( TIME-UNTIME ) / DBLE( IC )
                        OPS = DOPBL2( CNAME, M, N, 0, 0 )
                        RESLTS( IM, IN, I3 ) = DMFLOP( OPS, TIME, 0 )
  800                CONTINUE
  810             CONTINUE
  820          CONTINUE
  830       CONTINUE
            WRITE( NOUT, FMT = 9985 )
            CALL DPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NINC*NLDA,
     $                   RESLTS, LDR1, LDR2, NOUT )
*
         ELSE IF( CNAME.EQ.'DSYR  ' ) THEN
*
*           Time DSYR
*
            DO 890 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 6
               IF( UPLO.EQ.'L' )
     $            IMAT = -6
               I3 = 0
               DO 880 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 870 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 860 IN = 1, NN
                        N = NVAL( IN )
                        CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        CALL DTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                        IC = 0
                        S1 = DSECND( )
  840                   CONTINUE
                        CALL DSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
                        S2 = DSECND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL DTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 840
                        END IF
*
*                       Subtract the time used in DTIMMG.
*
                        ICL = 1
                        S1 = DSECND( )
  850                   CONTINUE
                        S2 = DSECND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL DTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 850
                        END IF
*
                        TIME = ( TIME-UNTIME ) / DBLE( IC )
                        OPS = DOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, I3 ) = DMFLOP( OPS, TIME, 0 )
  860                CONTINUE
  870             CONTINUE
  880          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  890       CONTINUE
*
         ELSE IF( CNAME.EQ.'DSYR2 ' ) THEN
*
*           Time DSYR2
*
            DO 950 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 6
               IF( UPLO.EQ.'L' )
     $            IMAT = -6
               I3 = 0
               DO 940 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 930 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 920 IN = 1, NN
                        N = NVAL( IN )
                        CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        CALL DTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                        IC = 0
                        S1 = DSECND( )
  900                   CONTINUE
                        CALL DSYR2( UPLO, N, ALPHA, X, INCX, Y, INCX, A,
     $                              LDA )
                        S2 = DSECND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL DTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 900
                        END IF
*
*                       Subtract the time used in DTIMMG.
*
                        ICL = 1
                        S1 = DSECND( )
  910                   CONTINUE
                        S2 = DSECND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL DTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 910
                        END IF
*
                        TIME = ( TIME-UNTIME ) / DBLE( IC )
                        OPS = DOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, I3 ) = DMFLOP( OPS, TIME, 0 )
  920                CONTINUE
  930             CONTINUE
  940          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  950       CONTINUE
*
         ELSE IF( CNAME.EQ.'DSPR  ' ) THEN
*
*           Time DSPR
*
            DO 1000 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 7
               IF( UPLO.EQ.'L' )
     $            IMAT = -7
               ILDA = 1
               LDA = LDAVAL( ILDA )
               DO 990 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  DO 980 IN = 1, NN
                     N = NVAL( IN )
                     CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                     CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                     CALL DTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0, 0 )
                     IC = 0
                     S1 = DSECND( )
  960                CONTINUE
                     CALL DSPR( UPLO, N, ALPHA, X, INCX, A )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL DTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 960
                     END IF
*
*                    Subtract the time used in DTIMMG.
*
                     ICL = 1
                     S1 = DSECND( )
  970                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL DTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 970
                     END IF
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
                     OPS = DOPBL2( CNAME, N, N, 0, 0 )
                     RESLTS( 1, IN, IINC ) = DMFLOP( OPS, TIME, 0 )
  980             CONTINUE
  990          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC, RESLTS,
     $                      LDR1, LDR2, NOUT )
 1000       CONTINUE
*
         ELSE IF( CNAME.EQ.'DSPR2 ' ) THEN
*
*           Time DSPR2
*
            DO 1050 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 7
               IF( UPLO.EQ.'L' )
     $            IMAT = -7
               ILDA = 1
               LDA = LDAVAL( ILDA )
               DO 1040 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  DO 1030 IN = 1, NN
                     N = NVAL( IN )
                     CALL DTIMMG( 0, 1, N, X, INCX, 0, 0 )
                     CALL DTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                     CALL DTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0, 0 )
                     IC = 0
                     S1 = DSECND( )
 1010                CONTINUE
                     CALL DSPR2( UPLO, N, ALPHA, X, INCX, Y, INCX, A )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL DTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 1010
                     END IF
*
*                    Subtract the time used in DTIMMG.
*
                     ICL = 1
                     S1 = DSECND( )
 1020                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL DTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 1020
                     END IF
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
                     OPS = DOPBL2( CNAME, N, N, 0, 0 )
                     RESLTS( 1, IN, IINC ) = DMFLOP( OPS, TIME, 0 )
 1030             CONTINUE
 1040          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL DPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC, RESLTS,
     $                      LDR1, LDR2, NOUT )
 1050       CONTINUE
         END IF
         WRITE( NOUT, FMT = 9984 )
 1060 CONTINUE
 1070 CONTINUE
*
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9997 FORMAT( 5X, 'with LDA = ', I5, ' and INCX = INCY = ', I5 )
 9996 FORMAT( 5X, 'with LDA = ', I5, ' and INCX = ', I5 )
 9995 FORMAT( 5X, 'with INCX = INCY = ', I5 )
 9994 FORMAT( 5X, 'with INCX = ', I5 )
 9993 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5,
     $      ' and INCX = INCY = ', I5 )
 9992 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5, ' and INCX = ', I5 )
 9991 FORMAT( 5X, 'line ', I2, ' with INCX = INCY = ', I5 )
 9990 FORMAT( 5X, 'line ', I2, ' with INCX = ', I5 )
 9989 FORMAT( / 1X, 'DGEMV  with TRANS = ''', A1, '''', / )
 9988 FORMAT( / 1X, 'DGBMV  with TRANS = ''', A1,
     $      ''', M = N and KL = K', 'U ', '= K', / )
 9987 FORMAT( / 1X, A, ' with UPLO = ''', A1, ''', TRANS = ''', A1,
     $      '''', / )
 9986 FORMAT( / 1X, A, ' with UPLO = ''', A1, '''', / )
 9985 FORMAT( / 1X, 'DGER', / )
 9984 FORMAT( / / / / / )
      RETURN
*
*     End of DTIMB2
*
      END
