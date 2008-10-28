      SUBROUTINE CTIMB2( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NINC,
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
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            INCVAL( * ), KVAL( * ), LDAVAL( * ), MVAL( * ),
     $                   NVAL( * )
      REAL               RESLTS( LDR1, LDR2, * )
      COMPLEX            A( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  CTIMB2 times the BLAS 2 routines.
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
*  TIMMIN  (input) REAL
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) COMPLEX array, dimension (LA)
*
*  X       (workspace) COMPLEX array, dimension (NMAX*INCMAX)
*             where NMAX and INCMAX are the maximum values permitted
*             for N and INCX.
*
*  Y       (workspace) COMPLEX array, dimension (NMAX*INCMAX)
*             where NMAX and INCMAX are the maximum values permitted
*             for N and INCX.
*
*  RESLTS  (output) REAL array, dimension (LDR1,LDR2,p),
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
      INTEGER            NSUBS, NTRANS, NUPLOS
      PARAMETER          ( NSUBS = 21, NTRANS = 3, NUPLOS = 2 )
      REAL               RALPHA
      PARAMETER          ( RALPHA = 1.0E0 )
      COMPLEX            ALPHA, BETA
      PARAMETER          ( ALPHA = ( 1.0E0, 0.0E0 ),
     $                   BETA = ( 1.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            IXANDY
      CHARACTER          TRANSA, UPLO
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, I3, IC, ICL, IINC, IK, ILDA, IM, IMAT, IN,
     $                   INCX, INFO, ISUB, ITA, IUPLO, J, K, LDA, M, N,
     $                   NX, NY
      REAL               OPS, S1, S2, TIME, UNTIME
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
      REAL               SECOND, SMFLOP, SOPBL2
      EXTERNAL           SECOND, SMFLOP, SOPBL2
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, CGBMV, CGEMV, CGERC, CGERU,
     $                   CHBMV, CHEMV, CHER, CHER2, CHPMV, CHPR, CHPR2,
     $                   CSPMV, CSPR, CSYMV, CSYR, CTBMV, CTBSV, CTIMMG,
     $                   CTPMV, CTPSV, CTRMV, CTRSV, SPRTBL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Data statements ..
      DATA               TRANS / 'N', 'T', 'C' /
      DATA               UPLOS / 'U', 'L' /
      DATA               NAMES / 'CGEMV ', 'CGBMV ', 'CHEMV ', 'CHBMV ',
     $                   'CHPMV ', 'CTRMV ', 'CTBMV ', 'CTPMV ',
     $                   'CTRSV ', 'CTBSV ', 'CTPSV ', 'CGERU ',
     $                   'CGERC ', 'CHER  ', 'CHPR  ', 'CHER2 ',
     $                   'CHPR2 ', 'CSYMV ', 'CSYR  ', 'CSPMV ',
     $                   'CSPR  ' /
*     ..
*     .. Executable Statements ..
*
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'B2'
      CALL ATIMIN( PATH, LINE, NSUBS, NAMES, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 1350
*
*     Time each routine
*
      DO 1340 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 1340
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
            GO TO 1340
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
*        Time CGEMV
*
         IF( CNAME.EQ.'CGEMV ' ) THEN
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
                           CALL CTIMMG( 1, M, N, A, LDA, 0, 0 )
                           CALL CTIMMG( 0, 1, NX, X, INCX, 0, 0 )
                           CALL CTIMMG( 0, 1, NY, Y, INCX, 0, 0 )
                           IC = 0
                           S1 = SECOND( )
   40                      CONTINUE
                           CALL CGEMV( TRANSA, M, N, ALPHA, A, LDA, X,
     $                                 INCX, BETA, Y, INCX )
                           S2 = SECOND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL CTIMMG( 0, 1, NY, Y, INCX, 0, 0 )
                              GO TO 40
                           END IF
*
*                          Subtract the time used in CTIMMG.
*
                           ICL = 1
                           S1 = SECOND( )
   50                      CONTINUE
                           S2 = SECOND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL CTIMMG( 0, 1, NY, Y, INCX, 0, 0 )
                              GO TO 50
                           END IF
*
                           TIME = ( TIME-UNTIME ) / REAL( IC )
                           OPS = SOPBL2( CNAME, M, N, 0, 0 )
                           RESLTS( IM, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
   60                   CONTINUE
   70                CONTINUE
   80             CONTINUE
   90          CONTINUE
               WRITE( NOUT, FMT = 9989 )TRANSA
               CALL SPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  100       CONTINUE
*
*        Time CGBMV
*
         ELSE IF( CNAME.EQ.'CGBMV ' ) THEN
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
                           CALL CTIMMG( -2, M, N, A, LDA, K, K )
                           CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           CALL CTIMMG( 0, 1, M, Y, INCX, 0, 0 )
                           IC = 0
                           S1 = SECOND( )
  110                      CONTINUE
                           CALL CGBMV( TRANSA, M, N, K, K, ALPHA, A,
     $                                 LDA, X, INCX, BETA, Y, INCX )
                           S2 = SECOND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL CTIMMG( 0, 1, M, Y, INCX, 0, 0 )
                              GO TO 110
                           END IF
*
*                          Subtract the time used in CTIMMG.
*
                           ICL = 1
                           S1 = SECOND( )
  120                      CONTINUE
                           S2 = SECOND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL CTIMMG( 0, 1, M, Y, INCX, 0, 0 )
                              GO TO 120
                           END IF
*
                           TIME = ( TIME-UNTIME ) / REAL( IC )
                           OPS = SOPBL2( CNAME, M, N, K, K )
                           RESLTS( IK, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
  130                   CONTINUE
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
               WRITE( NOUT, FMT = 9988 )TRANSA
               CALL SPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  170       CONTINUE
*
*        Time CHEMV
*
         ELSE IF( CNAME.EQ.'CHEMV ' ) THEN
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
                        CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                        CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
  180                   CONTINUE
                        CALL CHEMV( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                              BETA, Y, INCX )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                           GO TO 180
                        END IF
*
*                       Subtract the time used in CTIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
  190                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                           GO TO 190
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
  200                CONTINUE
  210             CONTINUE
  220          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  230       CONTINUE
*
*        Time CSYMV
*
         ELSE IF( CNAME.EQ.'CSYMV ' ) THEN
            DO 290 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 8
               IF( UPLO.EQ.'L' )
     $            IMAT = -8
               I3 = 0
               DO 280 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 270 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 260 IN = 1, NN
                        N = NVAL( IN )
                        CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                        CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
  240                   CONTINUE
                        CALL CSYMV( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                              BETA, Y, INCX )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                           GO TO 240
                        END IF
*
*                       Subtract the time used in CTIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
  250                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                           GO TO 250
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
  260                CONTINUE
  270             CONTINUE
  280          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  290       CONTINUE
*
*        Time CHBMV
*
         ELSE IF( CNAME.EQ.'CHBMV ' ) THEN
            DO 360 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 5
               IF( UPLO.EQ.'L' )
     $            IMAT = -5
               I3 = 0
               DO 350 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 340 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 330 IK = 1, NK
                        K = KVAL( IK )
                        DO 320 IN = 1, NN
                           N = NVAL( IN )
                           CALL CTIMMG( IMAT, N, N, A, LDA, K, K )
                           CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                           IC = 0
                           S1 = SECOND( )
  300                      CONTINUE
                           CALL CHBMV( UPLO, N, K, ALPHA, A, LDA, X,
     $                                 INCX, BETA, Y, INCX )
                           S2 = SECOND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                              GO TO 300
                           END IF
*
*                          Subtract the time used in CTIMMG.
*
                           ICL = 1
                           S1 = SECOND( )
  310                      CONTINUE
                           S2 = SECOND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                              GO TO 310
                           END IF
*
                           TIME = ( TIME-UNTIME ) / REAL( IC )
                           OPS = SOPBL2( CNAME, N, N, K, K )
                           RESLTS( IK, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
  320                   CONTINUE
  330                CONTINUE
  340             CONTINUE
  350          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
  360       CONTINUE
*
*        Time CHPMV
*
         ELSE IF( CNAME.EQ.'CHPMV ' ) THEN
            DO 410 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 7
               IF( UPLO.EQ.'L' )
     $            IMAT = -7
               ILDA = 1
               LDA = LDAVAL( ILDA )
               DO 400 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  DO 390 IN = 1, NN
                     N = NVAL( IN )
                     CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0, 0 )
                     CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                     CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                     IC = 0
                     S1 = SECOND( )
  370                CONTINUE
                     CALL CHPMV( UPLO, N, ALPHA, A, X, INCX, BETA, Y,
     $                           INCX )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        GO TO 370
                     END IF
*
*                    Subtract the time used in CTIMMG.
*
                     ICL = 1
                     S1 = SECOND( )
  380                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        GO TO 380
                     END IF
*
                     TIME = ( TIME-UNTIME ) / REAL( IC )
                     OPS = SOPBL2( CNAME, N, N, 0, 0 )
                     RESLTS( 1, IN, IINC ) = SMFLOP( OPS, TIME, 0 )
  390             CONTINUE
  400          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC, RESLTS,
     $                      LDR1, LDR2, NOUT )
  410       CONTINUE
*
*        Time CSPMV
*
         ELSE IF( CNAME.EQ.'CSPMV ' ) THEN
            DO 460 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 9
               IF( UPLO.EQ.'L' )
     $            IMAT = -9
               ILDA = 1
               LDA = LDAVAL( ILDA )
               DO 450 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  DO 440 IN = 1, NN
                     N = NVAL( IN )
                     CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0, 0 )
                     CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                     CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                     IC = 0
                     S1 = SECOND( )
  420                CONTINUE
                     CALL CSPMV( UPLO, N, ALPHA, A, X, INCX, BETA, Y,
     $                           INCX )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        GO TO 420
                     END IF
*
*                    Subtract the time used in CTIMMG.
*
                     ICL = 1
                     S1 = SECOND( )
  430                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        GO TO 430
                     END IF
*
                     TIME = ( TIME-UNTIME ) / REAL( IC )
                     OPS = SOPBL2( CNAME, N, N, 0, 0 )
                     RESLTS( 1, IN, IINC ) = SMFLOP( OPS, TIME, 0 )
  440             CONTINUE
  450          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC, RESLTS,
     $                      LDR1, LDR2, NOUT )
  460       CONTINUE
*
*        Time CTRMV
*
         ELSE IF( CNAME.EQ.'CTRMV ' ) THEN
            DO 530 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 11
               IF( UPLO.EQ.'L' )
     $            IMAT = -11
               DO 520 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  I3 = 0
                  DO 510 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 500 IINC = 1, NINC
                        INCX = INCVAL( IINC )
                        I3 = I3 + 1
                        DO 490 IN = 1, NN
                           N = NVAL( IN )
                           CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           IC = 0
                           S1 = SECOND( )
  470                      CONTINUE
                           CALL CTRMV( UPLO, TRANSA, 'Non-unit', N, A,
     $                                 LDA, X, INCX )
                           S2 = SECOND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              GO TO 470
                           END IF
*
*                          Subtract the time used in CTIMMG.
*
                           ICL = 1
                           S1 = SECOND( )
  480                      CONTINUE
                           S2 = SECOND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              GO TO 480
                           END IF
*
                           TIME = ( TIME-UNTIME ) / REAL( IC )
                           OPS = SOPBL2( CNAME, N, N, 0, 0 )
                           RESLTS( 1, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
  490                   CONTINUE
  500                CONTINUE
  510             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  520          CONTINUE
  530       CONTINUE
*
*        Time CTRSV
*
         ELSE IF( CNAME.EQ.'CTRSV ' ) THEN
            DO 600 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 11
               IF( UPLO.EQ.'L' )
     $            IMAT = -11
               DO 590 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  I3 = 0
                  DO 580 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 570 IINC = 1, NINC
                        INCX = INCVAL( IINC )
                        I3 = I3 + 1
                        DO 560 IN = 1, NN
                           N = NVAL( IN )
                           CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           IC = 0
                           S1 = SECOND( )
  540                      CONTINUE
                           CALL CTRSV( UPLO, TRANSA, 'Non-unit', N, A,
     $                                 LDA, X, INCX )
                           S2 = SECOND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              GO TO 540
                           END IF
*
*                          Subtract the time used in CTIMMG.
*
                           ICL = 1
                           S1 = SECOND( )
  550                      CONTINUE
                           S2 = SECOND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              GO TO 550
                           END IF
*
                           TIME = ( TIME-UNTIME ) / REAL( IC )
                           OPS = SOPBL2( CNAME, N, N, 0, 0 )
                           RESLTS( 1, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
  560                   CONTINUE
  570                CONTINUE
  580             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  590          CONTINUE
  600       CONTINUE
*
*        Time CTBMV
*
         ELSE IF( CNAME.EQ.'CTBMV ' ) THEN
            DO 680 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 13
               IF( UPLO.EQ.'L' )
     $            IMAT = -13
               DO 670 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  I3 = 0
                  DO 660 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 650 IINC = 1, NINC
                        INCX = INCVAL( IINC )
                        I3 = I3 + 1
                        DO 640 IK = 1, NK
                           K = KVAL( IK )
                           DO 630 IN = 1, NN
                              N = NVAL( IN )
                              CALL CTIMMG( IMAT, N, N, A, LDA, K, K )
                              CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              IC = 0
                              S1 = SECOND( )
  610                         CONTINUE
                              CALL CTBMV( UPLO, TRANSA, 'Non-unit', N,
     $                                    K, A, LDA, X, INCX )
                              S2 = SECOND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                                 GO TO 610
                              END IF
*
*                             Subtract the time used in CTIMMG.
*
                              ICL = 1
                              S1 = SECOND( )
  620                         CONTINUE
                              S2 = SECOND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                                 GO TO 620
                              END IF
*
                              TIME = ( TIME-UNTIME ) / REAL( IC )
                              OPS = SOPBL2( CNAME, N, N, K, K )
                              RESLTS( IK, IN, I3 ) = SMFLOP( OPS, TIME,
     $                           0 )
  630                      CONTINUE
  640                   CONTINUE
  650                CONTINUE
  660             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL SPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NINC*NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  670          CONTINUE
  680       CONTINUE
*
*        Time CTBSV
*
         ELSE IF( CNAME.EQ.'CTBSV ' ) THEN
            DO 760 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 13
               IF( UPLO.EQ.'L' )
     $            IMAT = -13
               DO 750 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  I3 = 0
                  DO 740 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 730 IINC = 1, NINC
                        INCX = INCVAL( IINC )
                        I3 = I3 + 1
                        DO 720 IK = 1, NK
                           K = KVAL( IK )
                           DO 710 IN = 1, NN
                              N = NVAL( IN )
                              CALL CTIMMG( IMAT, N, N, A, LDA, K, K )
                              CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                              IC = 0
                              S1 = SECOND( )
  690                         CONTINUE
                              CALL CTBSV( UPLO, TRANSA, 'Non-unit', N,
     $                                    K, A, LDA, X, INCX )
                              S2 = SECOND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                                 GO TO 690
                              END IF
*
*                             Subtract the time used in CTIMMG.
*
                              ICL = 1
                              S1 = SECOND( )
  700                         CONTINUE
                              S2 = SECOND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                                 GO TO 700
                              END IF
*
                              TIME = ( TIME-UNTIME ) / REAL( IC )
                              OPS = SOPBL2( CNAME, N, N, K, K )
                              RESLTS( IK, IN, I3 ) = SMFLOP( OPS, TIME,
     $                           0 )
  710                      CONTINUE
  720                   CONTINUE
  730                CONTINUE
  740             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL SPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NINC*NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  750          CONTINUE
  760       CONTINUE
*
*        Time CTPMV
*
         ELSE IF( CNAME.EQ.'CTPMV ' ) THEN
            DO 820 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 12
               IF( UPLO.EQ.'L' )
     $            IMAT = -12
               DO 810 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  ILDA = 1
                  LDA = LDAVAL( ILDA )
                  DO 800 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     DO 790 IN = 1, NN
                        N = NVAL( IN )
                        CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
  770                   CONTINUE
                        CALL CTPMV( UPLO, TRANSA, 'Non-unit', N, A, X,
     $                              INCX )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           GO TO 770
                        END IF
*
*                       Subtract the time used in CTIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
  780                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           GO TO 780
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, IINC ) = SMFLOP( OPS, TIME, 0 )
  790                CONTINUE
  800             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC,
     $                         RESLTS, LDR1, LDR2, NOUT )
  810          CONTINUE
  820       CONTINUE
*
*        Time CTPSV
*
         ELSE IF( CNAME.EQ.'CTPSV ' ) THEN
            DO 880 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 12
               IF( UPLO.EQ.'L' )
     $            IMAT = -12
               DO 870 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  ILDA = 1
                  LDA = LDAVAL( ILDA )
                  DO 860 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     DO 850 IN = 1, NN
                        N = NVAL( IN )
                        CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
  830                   CONTINUE
                        CALL CTPSV( UPLO, TRANSA, 'Non-unit', N, A, X,
     $                              INCX )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           GO TO 830
                        END IF
*
*                       Subtract the time used in CTIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
  840                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                           GO TO 840
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, IINC ) = SMFLOP( OPS, TIME, 0 )
  850                CONTINUE
  860             CONTINUE
                  WRITE( NOUT, FMT = 9987 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC,
     $                         RESLTS, LDR1, LDR2, NOUT )
  870          CONTINUE
  880       CONTINUE
*
*        Time CGERU
*
         ELSE IF( CNAME.EQ.'CGERU ' ) THEN
            I3 = 0
            DO 940 ILDA = 1, NLDA
               LDA = LDAVAL( ILDA )
               DO 930 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  I3 = I3 + 1
                  DO 920 IM = 1, NM
                     M = MVAL( IM )
                     DO 910 IN = 1, NN
                        N = NVAL( IN )
                        CALL CTIMMG( 0, 1, M, X, INCX, 0, 0 )
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        CALL CTIMMG( 1, M, N, A, LDA, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
  890                   CONTINUE
                        CALL CGERU( M, N, ALPHA, X, INCX, Y, INCX, A,
     $                              LDA )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL CTIMMG( 1, M, N, A, LDA, 0, 0 )
                           GO TO 890
                        END IF
*
*                       Subtract the time used in CTIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
  900                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL CTIMMG( 1, M, N, A, LDA, 0, 0 )
                           GO TO 900
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPBL2( CNAME, M, N, 0, 0 )
                        RESLTS( IM, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
  910                CONTINUE
  920             CONTINUE
  930          CONTINUE
  940       CONTINUE
            WRITE( NOUT, FMT = 9985 )CNAME(1:ILA_LEN_TRIM(CNAME))
            CALL SPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NINC*NLDA,
     $                   RESLTS, LDR1, LDR2, NOUT )
*
*        Time CGERC
*
         ELSE IF( CNAME.EQ.'CGERC ' ) THEN
            I3 = 0
            DO 1000 ILDA = 1, NLDA
               LDA = LDAVAL( ILDA )
               DO 990 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  I3 = I3 + 1
                  DO 980 IM = 1, NM
                     M = MVAL( IM )
                     DO 970 IN = 1, NN
                        N = NVAL( IN )
                        CALL CTIMMG( 0, 1, M, X, INCX, 0, 0 )
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        CALL CTIMMG( 1, M, N, A, LDA, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
  950                   CONTINUE
                        CALL CGERC( M, N, ALPHA, X, INCX, Y, INCX, A,
     $                              LDA )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL CTIMMG( 1, M, N, A, LDA, 0, 0 )
                           GO TO 950
                        END IF
*
*                       Subtract the time used in CTIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
  960                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL CTIMMG( 1, M, N, A, LDA, 0, 0 )
                           GO TO 960
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPBL2( CNAME, M, N, 0, 0 )
                        RESLTS( IM, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
  970                CONTINUE
  980             CONTINUE
  990          CONTINUE
 1000       CONTINUE
            WRITE( NOUT, FMT = 9985 )CNAME(1:ILA_LEN_TRIM(CNAME))
            CALL SPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NINC*NLDA,
     $                   RESLTS, LDR1, LDR2, NOUT )
*
*        Time CHER
*
         ELSE IF( CNAME.EQ.'CHER  ' ) THEN
            DO 1060 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 6
               IF( UPLO.EQ.'L' )
     $            IMAT = -6
               I3 = 0
               DO 1050 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 1040 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 1030 IN = 1, NN
                        N = NVAL( IN )
                        CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
 1010                   CONTINUE
                        CALL CHER( UPLO, N, RALPHA, X, INCX, A, LDA )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 1010
                        END IF
*
*                       Subtract the time used in CTIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
 1020                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 1020
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
 1030                CONTINUE
 1040             CONTINUE
 1050          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
 1060       CONTINUE
*
*        Time CSYR
*
         ELSE IF( CNAME.EQ.'CSYR  ' ) THEN
            DO 1120 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 8
               IF( UPLO.EQ.'L' )
     $            IMAT = -8
               I3 = 0
               DO 1110 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 1100 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 1090 IN = 1, NN
                        N = NVAL( IN )
                        CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
 1070                   CONTINUE
                        CALL CSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 1070
                        END IF
*
*                       Subtract the time used in CTIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
 1080                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 1080
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
 1090                CONTINUE
 1100             CONTINUE
 1110          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
 1120       CONTINUE
*
*        Time CHER2
*
         ELSE IF( CNAME.EQ.'CHER2 ' ) THEN
            DO 1180 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 6
               IF( UPLO.EQ.'L' )
     $            IMAT = -6
               I3 = 0
               DO 1170 ILDA = 1, NLDA
                  LDA = LDAVAL( ILDA )
                  DO 1160 IINC = 1, NINC
                     INCX = INCVAL( IINC )
                     I3 = I3 + 1
                     DO 1150 IN = 1, NN
                        N = NVAL( IN )
                        CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                        CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                        CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                        IC = 0
                        S1 = SECOND( )
 1130                   CONTINUE
                        CALL CHER2( UPLO, N, ALPHA, X, INCX, Y, INCX, A,
     $                              LDA )
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN ) THEN
                           CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 1130
                        END IF
*
*                       Subtract the time used in CTIMMG.
*
                        ICL = 1
                        S1 = SECOND( )
 1140                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
                        ICL = ICL + 1
                        IF( ICL.LE.IC ) THEN
                           CALL CTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           GO TO 1140
                        END IF
*
                        TIME = ( TIME-UNTIME ) / REAL( IC )
                        OPS = SOPBL2( CNAME, N, N, 0, 0 )
                        RESLTS( 1, IN, I3 ) = SMFLOP( OPS, TIME, 0 )
 1150                CONTINUE
 1160             CONTINUE
 1170          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC*NLDA,
     $                      RESLTS, LDR1, LDR2, NOUT )
 1180       CONTINUE
*
*        Time CHPR
*
         ELSE IF( CNAME.EQ.'CHPR  ' ) THEN
            DO 1230 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 7
               IF( UPLO.EQ.'L' )
     $            IMAT = -7
               ILDA = 1
               LDA = LDAVAL( ILDA )
               DO 1220 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  DO 1210 IN = 1, NN
                     N = NVAL( IN )
                     CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                     CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                     CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0, 0 )
                     IC = 0
                     S1 = SECOND( )
 1190                CONTINUE
                     CALL CHPR( UPLO, N, RALPHA, X, INCX, A )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 1190
                     END IF
*
*                    Subtract the time used in CTIMMG.
*
                     ICL = 1
                     S1 = SECOND( )
 1200                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 1200
                     END IF
*
                     TIME = ( TIME-UNTIME ) / REAL( IC )
                     OPS = SOPBL2( CNAME, N, N, 0, 0 )
                     RESLTS( 1, IN, IINC ) = SMFLOP( OPS, TIME, 0 )
 1210             CONTINUE
 1220          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC, RESLTS,
     $                      LDR1, LDR2, NOUT )
 1230       CONTINUE
*
*        Time CSPR
*
         ELSE IF( CNAME.EQ.'CSPR  ' ) THEN
            DO 1280 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 9
               IF( UPLO.EQ.'L' )
     $            IMAT = -9
               ILDA = 1
               LDA = LDAVAL( ILDA )
               DO 1270 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  DO 1260 IN = 1, NN
                     N = NVAL( IN )
                     CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                     CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                     CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0, 0 )
                     IC = 0
                     S1 = SECOND( )
 1240                CONTINUE
                     CALL CSPR( UPLO, N, ALPHA, X, INCX, A )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 1240
                     END IF
*
*                    Subtract the time used in CTIMMG.
*
                     ICL = 1
                     S1 = SECOND( )
 1250                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 1250
                     END IF
*
                     TIME = ( TIME-UNTIME ) / REAL( IC )
                     OPS = SOPBL2( CNAME, N, N, 0, 0 )
                     RESLTS( 1, IN, IINC ) = SMFLOP( OPS, TIME, 0 )
 1260             CONTINUE
 1270          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC, RESLTS,
     $                      LDR1, LDR2, NOUT )
 1280       CONTINUE
*
*        Time CHPR2
*
         ELSE IF( CNAME.EQ.'CHPR2 ' ) THEN
            DO 1330 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IMAT = 7
               IF( UPLO.EQ.'L' )
     $            IMAT = -7
               ILDA = 1
               LDA = LDAVAL( ILDA )
               DO 1320 IINC = 1, NINC
                  INCX = INCVAL( IINC )
                  DO 1310 IN = 1, NN
                     N = NVAL( IN )
                     CALL CTIMMG( 0, 1, N, X, INCX, 0, 0 )
                     CALL CTIMMG( 0, 1, N, Y, INCX, 0, 0 )
                     CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0, 0 )
                     IC = 0
                     S1 = SECOND( )
 1290                CONTINUE
                     CALL CHPR2( UPLO, N, ALPHA, X, INCX, Y, INCX, A )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 1290
                     END IF
*
*                    Subtract the time used in CTIMMG.
*
                     ICL = 1
                     S1 = SECOND( )
 1300                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL CTIMMG( IMAT, N, N, A, N*( N+1 ) / 2, 0,
     $                               0 )
                        GO TO 1300
                     END IF
*
                     TIME = ( TIME-UNTIME ) / REAL( IC )
                     OPS = SOPBL2( CNAME, N, N, 0, 0 )
                     RESLTS( 1, IN, IINC ) = SMFLOP( OPS, TIME, 0 )
 1310             CONTINUE
 1320          CONTINUE
               WRITE( NOUT, FMT = 9986 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO
               CALL SPRTBL( ' ', 'N', 1, NVAL, NN, NVAL, NINC, RESLTS,
     $                      LDR1, LDR2, NOUT )
 1330       CONTINUE
         END IF
         WRITE( NOUT, FMT = 9984 )
 1340 CONTINUE
 1350 CONTINUE
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
 9989 FORMAT( / 1X, 'CGEMV  with TRANS = ''', A1, '''', / )
 9988 FORMAT( / 1X, 'CGBMV  with TRANS = ''', A1,
     $      ''', M = N and KL = K', 'U ', '= K', / )
 9987 FORMAT( / 1X, A, ' with UPLO = ''', A1, ''', TRANS = ''', A1,
     $      '''', / )
 9986 FORMAT( / 1X, A, ' with UPLO = ''', A1, '''', / )
 9985 FORMAT( / 1X, A, / )
 9984 FORMAT( / / / / / )
      RETURN
*
*     End of CTIMB2
*
      END
