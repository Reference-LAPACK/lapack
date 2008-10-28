      SUBROUTINE STIMB3( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NLDA,
     $                   LDAVAL, TIMMIN, A, B, C, RESLTS, LDR1, LDR2,
     $                   NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    LINE
      INTEGER            LDR1, LDR2, NK, NLDA, NM, NN, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            KVAL( * ), LDAVAL( * ), MVAL( * ), NVAL( * )
      REAL               A( * ), B( * ), C( * ), RESLTS( LDR1, LDR2, * )
*     ..
*
*  Purpose
*  =======
*
*  STIMB3 times the Level 3 BLAS routines.
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
*          The values of K.  K is used as the intermediate matrix
*          dimension for SGEMM (the product of an M x K matrix and a
*          K x N matrix) and as the dimension of the rank-K update in
*          SSYRK and SSYR2K.
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
*             where LDAMAX and NMAX are the maximum values permitted
*             for LDA and N.
*
*  B       (workspace) REAL array, dimension (LDAMAX*NMAX)
*
*  C       (workspace) REAL array, dimension (LDAMAX*NMAX)
*
*  RESLTS  (output) REAL array, dimension (LDR1,LDR2,NLDA)
*          The timing results for each subroutine over the relevant
*          values of M, N, K, and LDA.
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
      PARAMETER          ( NSUBS = 6 )
      INTEGER            NTRANS, NSIDES, NUPLOS
      PARAMETER          ( NTRANS = 2, NSIDES = 2, NUPLOS = 2 )
      REAL               ALPHA, BETA
      PARAMETER          ( ALPHA = 1.0E0, BETA = 1.0E0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          SIDE, TRANSA, TRANSB, UPLO
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, IK, ILDA, IM, IMAT, IN, INFO,
     $                   ISIDE, ISUB, ITA, ITB, IUPLO, K, LDA, M, N
      REAL               OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER          SIDES( NSIDES ), TRANS( NTRANS ),
     $                   UPLOS( NUPLOS )
      CHARACTER(32)      NAMES( NSUBS )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      LOGICAL            LSAME
      REAL               SECOND, SMFLOP, SOPBL3
      EXTERNAL           LSAME, SECOND, SMFLOP, SOPBL3
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, SGEMM, SPRTBL, SSYMM, SSYR2K,
     $                   SSYRK, STIMMG, STRMM, STRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Data statements ..
      DATA               NAMES / 'SGEMM ', 'SSYMM ', 'SSYRK ', 'SSYR2K',
     $                   'STRMM ', 'STRSM ' /
      DATA               TRANS / 'N', 'T' /
      DATA               SIDES / 'L', 'R' /
      DATA               UPLOS / 'U', 'L' /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'B3'
      CALL ATIMIN( PATH, LINE, NSUBS, NAMES, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 480
*
*     Check that M <= LDA.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 1, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 480
      END IF
*
*     Time each routine.
*
      DO 470 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 470
*
*        Print header.
*
         CNAME = NAMES( ISUB )
         WRITE( NOUT, FMT = 9998 )CNAME(1:ILA_LEN_TRIM(CNAME))
         IF( NLDA.EQ.1 ) THEN
            WRITE( NOUT, FMT = 9997 )LDAVAL( 1 )
         ELSE
            DO 10 I = 1, NLDA
               WRITE( NOUT, FMT = 9996 )I, LDAVAL( I )
   10       CONTINUE
         END IF
*
*        Time SGEMM
*
         IF( CNAME.EQ.'SGEMM ' ) THEN
            DO 90 ITA = 1, NTRANS
               TRANSA = TRANS( ITA )
               DO 80 ITB = 1, NTRANS
                  TRANSB = TRANS( ITB )
                  DO 70 IK = 1, NK
                     K = KVAL( IK )
                     DO 60 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 50 IM = 1, NM
                           M = MVAL( IM )
                           DO 40 IN = 1, NN
                              N = NVAL( IN )
                              IF( TRANSA.EQ.'N' ) THEN
                                 CALL STIMMG( 1, M, K, A, LDA, 0, 0 )
                              ELSE
                                 CALL STIMMG( 1, K, M, A, LDA, 0, 0 )
                              END IF
                              IF( TRANSB.EQ.'N' ) THEN
                                 CALL STIMMG( 0, K, N, B, LDA, 0, 0 )
                              ELSE
                                 CALL STIMMG( 0, N, K, B, LDA, 0, 0 )
                              END IF
                              CALL STIMMG( 1, M, N, C, LDA, 0, 0 )
                              IC = 0
                              S1 = SECOND( )
   20                         CONTINUE
                              CALL SGEMM( TRANSA, TRANSB, M, N, K,
     $                                    ALPHA, A, LDA, B, LDA, BETA,
     $                                    C, LDA )
                              S2 = SECOND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL STIMMG( 1, M, N, C, LDA, 0, 0 )
                                 GO TO 20
                              END IF
*
*                             Subtract the time used in STIMMG.
*
                              ICL = 1
                              S1 = SECOND( )
   30                         CONTINUE
                              S2 = SECOND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL STIMMG( 1, M, N, C, LDA, 0, 0 )
                                 GO TO 30
                              END IF
*
                              TIME = ( TIME-UNTIME ) / REAL( IC )
                              OPS = SOPBL3( CNAME, M, N, K )
                              RESLTS( IM, IN, ILDA ) = SMFLOP( OPS,
     $                           TIME, 0 )
   40                      CONTINUE
   50                   CONTINUE
   60                CONTINUE
                     IF( IK.EQ.1 )
     $                  WRITE( NOUT, FMT = 9995 )TRANSA, TRANSB
                     WRITE( NOUT, FMT = 9994 )KVAL( IK )
                     CALL SPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
*
         ELSE IF( CNAME.EQ.'SSYMM ' ) THEN
*
*           Time SSYMM
*
            DO 160 ISIDE = 1, NSIDES
               SIDE = SIDES( ISIDE )
               DO 150 IUPLO = 1, NUPLOS
                  UPLO = UPLOS( IUPLO )
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     IMAT = 6
                  ELSE
                     IMAT = -6
                  END IF
                  DO 140 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 130 IM = 1, NM
                        M = MVAL( IM )
                        DO 120 IN = 1, NN
                           N = NVAL( IN )
                           IF( ISIDE.EQ.1 ) THEN
                              CALL STIMMG( IMAT, M, M, A, LDA, 0, 0 )
                              CALL STIMMG( 0, M, N, B, LDA, 0, 0 )
                           ELSE
                              CALL STIMMG( 0, M, N, B, LDA, 0, 0 )
                              CALL STIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           END IF
                           CALL STIMMG( 1, M, N, C, LDA, 0, 0 )
                           IC = 0
                           S1 = SECOND( )
  100                      CONTINUE
                           CALL SSYMM( SIDE, UPLO, M, N, ALPHA, A, LDA,
     $                                 B, LDA, BETA, C, LDA )
                           S2 = SECOND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL STIMMG( 1, M, N, C, LDA, 0, 0 )
                              GO TO 100
                           END IF
*
*                          Subtract the time used in STIMMG.
*
                           ICL = 1
                           S1 = SECOND( )
  110                      CONTINUE
                           S2 = SECOND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL STIMMG( 1, M, N, C, LDA, 0, 0 )
                              GO TO 110
                           END IF
*
                           TIME = ( TIME-UNTIME ) / REAL( IC )
                           OPS = SOPBL3( CNAME, M, N, ISIDE-1 )
                           RESLTS( IM, IN, ILDA ) = SMFLOP( OPS, TIME,
     $                        0 )
  120                   CONTINUE
  130                CONTINUE
  140             CONTINUE
                  WRITE( NOUT, FMT = 9993 )SIDE, UPLO
                  CALL SPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  150          CONTINUE
  160       CONTINUE
*
         ELSE IF( CNAME.EQ.'SSYRK ' ) THEN
*
*           Time SSYRK
*
            DO 230 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IMAT = 6
               ELSE
                  IMAT = -6
               END IF
               DO 220 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  DO 210 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 200 IK = 1, NK
                        K = KVAL( IK )
                        IF( TRANSA.EQ.'N' ) THEN
                           CALL STIMMG( 1, N, K, A, LDA, 0, 0 )
                        ELSE
                           CALL STIMMG( 1, K, N, A, LDA, 0, 0 )
                        END IF
                        DO 190 IN = 1, NN
                           N = NVAL( IN )
                           CALL STIMMG( IMAT, N, N, C, LDA, 0, 0 )
                           IC = 0
                           S1 = SECOND( )
  170                      CONTINUE
                           CALL SSYRK( UPLO, TRANSA, N, K, ALPHA, A,
     $                                 LDA, BETA, C, LDA )
                           S2 = SECOND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL STIMMG( IMAT, N, N, C, LDA, 0, 0 )
                              GO TO 170
                           END IF
*
*                          Subtract the time used in STIMMG.
*
                           ICL = 1
                           S1 = SECOND( )
  180                      CONTINUE
                           S2 = SECOND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL STIMMG( IMAT, N, N, C, LDA, 0, 0 )
                              GO TO 180
                           END IF
*
                           TIME = ( TIME-UNTIME ) / REAL( IC )
                           OPS = SOPBL3( CNAME, N, N, K )
                           RESLTS( IK, IN, ILDA ) = SMFLOP( OPS, TIME,
     $                        0 )
  190                   CONTINUE
  200                CONTINUE
  210             CONTINUE
                  WRITE( NOUT, FMT = 9992 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                  CALL SPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  220          CONTINUE
  230       CONTINUE
*
         ELSE IF( CNAME.EQ.'SSYR2K' ) THEN
*
*           Time SSYR2K
*
            DO 300 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IMAT = 6
               ELSE
                  IMAT = -6
               END IF
               DO 290 ITB = 1, NTRANS
                  TRANSB = TRANS( ITB )
                  DO 280 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 270 IK = 1, NK
                        K = KVAL( IK )
                        IF( TRANSB.EQ.'N' ) THEN
                           CALL STIMMG( 1, N, K, A, LDA, 0, 0 )
                           CALL STIMMG( 0, N, K, B, LDA, 0, 0 )
                        ELSE
                           CALL STIMMG( 1, K, N, A, LDA, 0, 0 )
                           CALL STIMMG( 0, K, N, B, LDA, 0, 0 )
                        END IF
                        DO 260 IN = 1, NN
                           N = NVAL( IN )
                           CALL STIMMG( IMAT, N, N, C, LDA, 0, 0 )
                           IC = 0
                           S1 = SECOND( )
  240                      CONTINUE
                           CALL SSYR2K( UPLO, TRANSB, N, K, ALPHA, A,
     $                                  LDA, B, LDA, BETA, C, LDA )
                           S2 = SECOND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL STIMMG( IMAT, N, N, C, LDA, 0, 0 )
                              GO TO 240
                           END IF
*
*                          Subtract the time used in STIMMG.
*
                           ICL = 1
                           S1 = SECOND( )
  250                      CONTINUE
                           S2 = SECOND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL STIMMG( IMAT, N, N, C, LDA, 0, 0 )
                              GO TO 250
                           END IF
*
                           TIME = ( TIME-UNTIME ) / REAL( IC )
                           OPS = SOPBL3( CNAME, N, N, K )
                           RESLTS( IK, IN, ILDA ) = SMFLOP( OPS, TIME,
     $                        0 )
  260                   CONTINUE
  270                CONTINUE
  280             CONTINUE
                  WRITE( NOUT, FMT = 9992 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSB
                  CALL SPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  290          CONTINUE
  300       CONTINUE
*
         ELSE IF( CNAME.EQ.'STRMM ' ) THEN
*
*           Time STRMM
*
            DO 380 ISIDE = 1, NSIDES
               SIDE = SIDES( ISIDE )
               DO 370 IUPLO = 1, NUPLOS
                  UPLO = UPLOS( IUPLO )
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     IMAT = 9
                  ELSE
                     IMAT = -9
                  END IF
                  DO 360 ITA = 1, NTRANS
                     TRANSA = TRANS( ITA )
                     DO 350 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 340 IM = 1, NM
                           M = MVAL( IM )
                           DO 330 IN = 1, NN
                              N = NVAL( IN )
                              IF( ISIDE.EQ.1 ) THEN
                                 CALL STIMMG( IMAT, M, M, A, LDA, 0, 0 )
                              ELSE
                                 CALL STIMMG( IMAT, N, N, A, LDA, 0, 0 )
                              END IF
                              CALL STIMMG( 0, M, N, B, LDA, 0, 0 )
                              IC = 0
                              S1 = SECOND( )
  310                         CONTINUE
                              CALL STRMM( SIDE, UPLO, TRANSA,
     $                                    'Non-unit', M, N, ALPHA, A,
     $                                    LDA, B, LDA )
                              S2 = SECOND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL STIMMG( 0, M, N, B, LDA, 0, 0 )
                                 GO TO 310
                              END IF
*
*                             Subtract the time used in STIMMG.
*
                              ICL = 1
                              S1 = SECOND( )
  320                         CONTINUE
                              S2 = SECOND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL STIMMG( 0, M, N, B, LDA, 0, 0 )
                                 GO TO 320
                              END IF
*
                              TIME = ( TIME-UNTIME ) / REAL( IC )
                              OPS = SOPBL3( CNAME, M, N, ISIDE-1 )
                              RESLTS( IM, IN, ILDA ) = SMFLOP( OPS,
     $                           TIME, 0 )
  330                      CONTINUE
  340                   CONTINUE
  350                CONTINUE
                     WRITE( NOUT, FMT = 9991 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), SIDE, UPLO, TRANSA
                     CALL SPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
  360             CONTINUE
  370          CONTINUE
  380       CONTINUE
*
         ELSE IF( CNAME.EQ.'STRSM ' ) THEN
*
*           Time STRSM
*
            DO 460 ISIDE = 1, NSIDES
               SIDE = SIDES( ISIDE )
               DO 450 IUPLO = 1, NUPLOS
                  UPLO = UPLOS( IUPLO )
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     IMAT = 9
                  ELSE
                     IMAT = -9
                  END IF
                  DO 440 ITA = 1, NTRANS
                     TRANSA = TRANS( ITA )
                     DO 430 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 420 IM = 1, NM
                           M = MVAL( IM )
                           DO 410 IN = 1, NN
                              N = NVAL( IN )
                              IF( ISIDE.EQ.1 ) THEN
                                 CALL STIMMG( IMAT, M, M, A, LDA, 0, 0 )
                              ELSE
                                 CALL STIMMG( IMAT, N, N, A, LDA, 0, 0 )
                              END IF
                              CALL STIMMG( 0, M, N, B, LDA, 0, 0 )
                              IC = 0
                              S1 = SECOND( )
  390                         CONTINUE
                              CALL STRSM( SIDE, UPLO, TRANSA,
     $                                    'Non-unit', M, N, ALPHA, A,
     $                                    LDA, B, LDA )
                              S2 = SECOND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL STIMMG( 0, M, N, B, LDA, 0, 0 )
                                 GO TO 390
                              END IF
*
*                             Subtract the time used in STIMMG.
*
                              ICL = 1
                              S1 = SECOND( )
  400                         CONTINUE
                              S2 = SECOND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL STIMMG( 0, M, N, B, LDA, 0, 0 )
                                 GO TO 400
                              END IF
*
                              TIME = ( TIME-UNTIME ) / REAL( IC )
                              OPS = SOPBL3( CNAME, M, N, ISIDE-1 )
                              RESLTS( IM, IN, ILDA ) = SMFLOP( OPS,
     $                           TIME, 0 )
  410                      CONTINUE
  420                   CONTINUE
  430                CONTINUE
                     WRITE( NOUT, FMT = 9991 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), SIDE, UPLO, TRANSA
                     CALL SPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
  440             CONTINUE
  450          CONTINUE
  460       CONTINUE
         END IF
         WRITE( NOUT, FMT = 9990 )
  470 CONTINUE
  480 CONTINUE
*
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9997 FORMAT( 5X, 'with LDA = ', I5 )
 9996 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
 9995 FORMAT( / 1X, 'SGEMM  with TRANSA = ''', A1, ''', TRANSB = ''',
     $      A1, '''' )
 9994 FORMAT( / 1X, 'K = ', I4, / )
 9993 FORMAT( / 1X, 'SSYMM  with SIDE = ''', A1, ''', UPLO = ''', A1,
     $      '''', / )
 9992 FORMAT( / 1X, A, ' with UPLO = ''', A1, ''', TRANS = ''', A1,
     $      '''', / )
 9991 FORMAT( / 1X, A, ' with SIDE = ''', A1, ''', UPLO = ''', A1,
     $      ''',', ' TRANS = ''', A1, '''', / )
 9990 FORMAT( / / / / / )
      RETURN
*
*     End of STIMB3
*
      END
