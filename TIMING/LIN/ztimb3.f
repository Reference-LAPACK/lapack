      SUBROUTINE ZTIMB3( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NLDA,
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
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            KVAL( * ), LDAVAL( * ), MVAL( * ), NVAL( * )
      DOUBLE PRECISION   RESLTS( LDR1, LDR2, * )
      COMPLEX*16         A( * ), B( * ), C( * )
*     ..
*
*  Purpose
*  =======
*
*  ZTIMB3 times the Level 3 BLAS routines.
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
*          dimension for ZGEMM (the product of an M x K matrix and a
*          K x N matrix) and as the dimension of the rank-K update in
*          ZHERK and ZSYRK.
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
*             where LDAMAX and NMAX are the maximum values permitted
*             for LDA and N.
*
*  B       (workspace) COMPLEX*16 array, dimension (LDAMAX*NMAX)
*
*  C       (workspace) COMPLEX*16 array, dimension (LDAMAX*NMAX)
*
*  RESLTS  (output) DOUBLE PRECISION array, dimension (LDR1,LDR2,NLDA)
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
      INTEGER            NSUBS, NTRANS, NSIDES, NUPLOS
      PARAMETER          ( NSUBS = 9, NTRANS = 3, NSIDES = 2,
     $                   NUPLOS = 2 )
      DOUBLE PRECISION   RALPHA, RBETA
      PARAMETER          ( RALPHA = 1.0D0, RBETA = 1.0D0 )
      COMPLEX*16         ALPHA, BETA
      PARAMETER          ( ALPHA = ( 1.0D0, 0.0D0 ),
     $                   BETA = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          SIDE, TRANSA, TRANSB, UPLO
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, IK, ILDA, IM, IMAT, IN, INFO,
     $                   ISIDE, ISUB, ITA, ITB, IUPLO, K, LDA, M, N
      DOUBLE PRECISION   OPS, S1, S2, TIME, UNTIME
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
      DOUBLE PRECISION   DMFLOP, DOPBL3, DSECND
      EXTERNAL           LSAME, DMFLOP, DOPBL3, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, DPRTBL, ZGEMM, ZHEMM, ZHER2K,
     $                   ZHERK, ZSYMM, ZSYR2K, ZSYRK, ZTIMMG, ZTRMM,
     $                   ZTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Data statements ..
      DATA               NAMES / 'ZGEMM ', 'ZHEMM ', 'ZSYMM ', 'ZHERK ',
     $                   'ZHER2K', 'ZSYRK ', 'ZSYR2K', 'ZTRMM ',
     $                   'ZTRSM ' /
      DATA               TRANS / 'N', 'T', 'C' /
      DATA               SIDES / 'L', 'R' /
      DATA               UPLOS / 'U', 'L' /
*     ..
*     .. Executable Statements ..
*
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'B3'
      CALL ATIMIN( PATH, LINE, NSUBS, NAMES, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 690
*
*     Check that M <= LDA.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 1, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 690
      END IF
*
*     Time each routine.
*
      DO 680 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 680
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
*        Time ZGEMM
*
         IF( CNAME.EQ.'ZGEMM ' ) THEN
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
                                 CALL ZTIMMG( 1, M, K, A, LDA, 0, 0 )
                              ELSE
                                 CALL ZTIMMG( 1, K, M, A, LDA, 0, 0 )
                              END IF
                              IF( TRANSB.EQ.'N' ) THEN
                                 CALL ZTIMMG( 0, K, N, B, LDA, 0, 0 )
                              ELSE
                                 CALL ZTIMMG( 0, N, K, B, LDA, 0, 0 )
                              END IF
                              CALL ZTIMMG( 1, M, N, C, LDA, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
   20                         CONTINUE
                              CALL ZGEMM( TRANSA, TRANSB, M, N, K,
     $                                    ALPHA, A, LDA, B, LDA, BETA,
     $                                    C, LDA )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL ZTIMMG( 1, M, N, C, LDA, 0, 0 )
                                 GO TO 20
                              END IF
*
*                             Subtract the time used in ZTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
   30                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL ZTIMMG( 1, M, N, C, LDA, 0, 0 )
                                 GO TO 30
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
                              OPS = DOPBL3( CNAME, M, N, K )
                              RESLTS( IM, IN, ILDA ) = DMFLOP( OPS,
     $                           TIME, 0 )
   40                      CONTINUE
   50                   CONTINUE
   60                CONTINUE
                     IF( IK.EQ.1 )
     $                  WRITE( NOUT, FMT = 9995 )TRANSA, TRANSB
                     WRITE( NOUT, FMT = 9994 )KVAL( IK )
                     CALL DPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
*
*        Time ZHEMM
*
         ELSE IF( CNAME.EQ.'ZHEMM ' ) THEN
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
                              CALL ZTIMMG( IMAT, M, M, A, LDA, 0, 0 )
                              CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                           ELSE
                              CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                              CALL ZTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           END IF
                           CALL ZTIMMG( 1, M, N, C, LDA, 0, 0 )
                           IC = 0
                           S1 = DSECND( )
  100                      CONTINUE
                           CALL ZHEMM( SIDE, UPLO, M, N, ALPHA, A, LDA,
     $                                 B, LDA, BETA, C, LDA )
                           S2 = DSECND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL ZTIMMG( 1, M, N, C, LDA, 0, 0 )
                              GO TO 100
                           END IF
*
*                          Subtract the time used in ZTIMMG.
*
                           ICL = 1
                           S1 = DSECND( )
  110                      CONTINUE
                           S2 = DSECND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL ZTIMMG( 1, M, N, C, LDA, 0, 0 )
                              GO TO 110
                           END IF
*
                           TIME = ( TIME-UNTIME ) / DBLE( IC )
                           OPS = DOPBL3( CNAME, M, N, ISIDE-1 )
                           RESLTS( IM, IN, ILDA ) = DMFLOP( OPS, TIME,
     $                        0 )
  120                   CONTINUE
  130                CONTINUE
  140             CONTINUE
                  WRITE( NOUT, FMT = 9993 )'ZHEMM ', SIDE, UPLO
                  CALL DPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  150          CONTINUE
  160       CONTINUE
*
*        Time ZSYMM
*
         ELSE IF( CNAME.EQ.'ZSYMM ' ) THEN
            DO 230 ISIDE = 1, NSIDES
               SIDE = SIDES( ISIDE )
               DO 220 IUPLO = 1, NUPLOS
                  UPLO = UPLOS( IUPLO )
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     IMAT = 8
                  ELSE
                     IMAT = -8
                  END IF
                  DO 210 ILDA = 1, NLDA
                     LDA = LDAVAL( ILDA )
                     DO 200 IM = 1, NM
                        M = MVAL( IM )
                        DO 190 IN = 1, NN
                           N = NVAL( IN )
                           IF( ISIDE.EQ.1 ) THEN
                              CALL ZTIMMG( IMAT, M, M, A, LDA, 0, 0 )
                              CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                           ELSE
                              CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                              CALL ZTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                           END IF
                           CALL ZTIMMG( 1, M, N, C, LDA, 0, 0 )
                           IC = 0
                           S1 = DSECND( )
  170                      CONTINUE
                           CALL ZSYMM( SIDE, UPLO, M, N, ALPHA, A, LDA,
     $                                 B, LDA, BETA, C, LDA )
                           S2 = DSECND( )
                           TIME = S2 - S1
                           IC = IC + 1
                           IF( TIME.LT.TIMMIN ) THEN
                              CALL ZTIMMG( 1, M, N, C, LDA, 0, 0 )
                              GO TO 170
                           END IF
*
*                          Subtract the time used in ZTIMMG.
*
                           ICL = 1
                           S1 = DSECND( )
  180                      CONTINUE
                           S2 = DSECND( )
                           UNTIME = S2 - S1
                           ICL = ICL + 1
                           IF( ICL.LE.IC ) THEN
                              CALL ZTIMMG( 1, M, N, C, LDA, 0, 0 )
                              GO TO 180
                           END IF
*
                           TIME = ( TIME-UNTIME ) / DBLE( IC )
                           OPS = DOPBL3( CNAME, M, N, ISIDE-1 )
                           RESLTS( IM, IN, ILDA ) = DMFLOP( OPS, TIME,
     $                        0 )
  190                   CONTINUE
  200                CONTINUE
  210             CONTINUE
                  WRITE( NOUT, FMT = 9993 )'ZSYMM ', SIDE, UPLO
                  CALL DPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NLDA,
     $                         RESLTS, LDR1, LDR2, NOUT )
  220          CONTINUE
  230       CONTINUE
*
*        Time ZHERK
*
         ELSE IF( CNAME.EQ.'ZHERK ' ) THEN
            DO 300 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IMAT = 6
               ELSE
                  IMAT = -6
               END IF
               DO 290 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  IF( TRANSA.NE.'T' ) THEN
                     DO 280 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 270 IK = 1, NK
                           K = KVAL( IK )
                           IF( TRANSA.EQ.'N' ) THEN
                              CALL ZTIMMG( 1, N, K, A, LDA, 0, 0 )
                           ELSE
                              CALL ZTIMMG( 1, K, N, A, LDA, 0, 0 )
                           END IF
                           DO 260 IN = 1, NN
                              N = NVAL( IN )
                              CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
  240                         CONTINUE
                              CALL ZHERK( UPLO, TRANSA, N, K, RALPHA, A,
     $                                    LDA, RBETA, C, LDA )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                                 GO TO 240
                              END IF
*
*                             Subtract the time used in ZTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
  250                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                                 GO TO 250
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
                              OPS = DOPBL3( CNAME, N, N, K )
                              RESLTS( IK, IN, ILDA ) = DMFLOP( OPS,
     $                           TIME, 0 )
  260                      CONTINUE
  270                   CONTINUE
  280                CONTINUE
                     WRITE( NOUT, FMT = 9992 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                     CALL DPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
                  END IF
  290          CONTINUE
  300       CONTINUE
*
*        Time ZHER2K
*
         ELSE IF( CNAME.EQ.'ZHER2K' ) THEN
            DO 370 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IMAT = 6
               ELSE
                  IMAT = -6
               END IF
               DO 360 ITB = 1, NTRANS
                  TRANSB = TRANS( ITB )
                  IF( TRANSB.NE.'T' ) THEN
                     DO 350 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 340 IK = 1, NK
                           K = KVAL( IK )
                           IF( TRANSB.EQ.'N' ) THEN
                              CALL ZTIMMG( 1, N, K, A, LDA, 0, 0 )
                              CALL ZTIMMG( 0, N, K, B, LDA, 0, 0 )
                           ELSE
                              CALL ZTIMMG( 1, K, N, A, LDA, 0, 0 )
                              CALL ZTIMMG( 0, K, N, B, LDA, 0, 0 )
                           END IF
                           DO 330 IN = 1, NN
                              N = NVAL( IN )
                              CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
  310                         CONTINUE
                              CALL ZHER2K( UPLO, TRANSB, N, K, ALPHA, A,
     $                                     LDA, B, LDA, RBETA, C, LDA )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                                 GO TO 310
                              END IF
*
*                             Subtract the time used in ZTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
  320                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                                 GO TO 320
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
                              OPS = DOPBL3( CNAME, N, N, K )
                              RESLTS( IK, IN, ILDA ) = DMFLOP( OPS,
     $                           TIME, 0 )
  330                      CONTINUE
  340                   CONTINUE
  350                CONTINUE
                     WRITE( NOUT, FMT = 9992 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSB
                     CALL DPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
                  END IF
  360          CONTINUE
  370       CONTINUE
*
*        Time ZSYRK
*
         ELSE IF( CNAME.EQ.'ZSYRK ' ) THEN
            DO 440 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IMAT = 8
               ELSE
                  IMAT = -8
               END IF
               DO 430 ITA = 1, NTRANS
                  TRANSA = TRANS( ITA )
                  IF( TRANSA.NE.'C' ) THEN
                     DO 420 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 410 IK = 1, NK
                           K = KVAL( IK )
                           IF( TRANSA.EQ.'N' ) THEN
                              CALL ZTIMMG( 1, N, K, A, LDA, 0, 0 )
                           ELSE
                              CALL ZTIMMG( 1, K, N, A, LDA, 0, 0 )
                           END IF
                           DO 400 IN = 1, NN
                              N = NVAL( IN )
                              CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
  380                         CONTINUE
                              CALL ZSYRK( UPLO, TRANSA, N, K, ALPHA, A,
     $                                    LDA, BETA, C, LDA )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                                 GO TO 380
                              END IF
*
*                             Subtract the time used in ZTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
  390                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                                 GO TO 390
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
                              OPS = DOPBL3( CNAME, N, N, K )
                              RESLTS( IK, IN, ILDA ) = DMFLOP( OPS,
     $                           TIME, 0 )
  400                      CONTINUE
  410                   CONTINUE
  420                CONTINUE
                     WRITE( NOUT, FMT = 9992 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSA
                     CALL DPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
                  END IF
  430          CONTINUE
  440       CONTINUE
*
*        Time ZSYR2K
*
         ELSE IF( CNAME.EQ.'ZSYR2K' ) THEN
            DO 510 IUPLO = 1, NUPLOS
               UPLO = UPLOS( IUPLO )
               IF( LSAME( UPLO, 'U' ) ) THEN
                  IMAT = 8
               ELSE
                  IMAT = -8
               END IF
               DO 500 ITB = 1, NTRANS
                  TRANSB = TRANS( ITB )
                  IF( TRANSB.NE.'C' ) THEN
                     DO 490 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 480 IK = 1, NK
                           K = KVAL( IK )
                           IF( TRANSB.EQ.'N' ) THEN
                              CALL ZTIMMG( 1, N, K, A, LDA, 0, 0 )
                              CALL ZTIMMG( 0, N, K, B, LDA, 0, 0 )
                           ELSE
                              CALL ZTIMMG( 1, K, N, A, LDA, 0, 0 )
                              CALL ZTIMMG( 0, K, N, B, LDA, 0, 0 )
                           END IF
                           DO 470 IN = 1, NN
                              N = NVAL( IN )
                              CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
  450                         CONTINUE
                              CALL ZSYR2K( UPLO, TRANSB, N, K, ALPHA, A,
     $                                     LDA, B, LDA, BETA, C, LDA )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                                 GO TO 450
                              END IF
*
*                             Subtract the time used in ZTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
  460                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL ZTIMMG( IMAT, N, N, C, LDA, 0, 0 )
                                 GO TO 460
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
                              OPS = DOPBL3( CNAME, N, N, K )
                              RESLTS( IK, IN, ILDA ) = DMFLOP( OPS,
     $                           TIME, 0 )
  470                      CONTINUE
  480                   CONTINUE
  490                CONTINUE
                     WRITE( NOUT, FMT = 9992 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), UPLO, TRANSB
                     CALL DPRTBL( 'K', 'N', NK, KVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
                  END IF
  500          CONTINUE
  510       CONTINUE
*
*        Time ZTRMM
*
         ELSE IF( CNAME.EQ.'ZTRMM ' ) THEN
            DO 590 ISIDE = 1, NSIDES
               SIDE = SIDES( ISIDE )
               DO 580 IUPLO = 1, NUPLOS
                  UPLO = UPLOS( IUPLO )
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     IMAT = 11
                  ELSE
                     IMAT = -11
                  END IF
                  DO 570 ITA = 1, NTRANS
                     TRANSA = TRANS( ITA )
                     DO 560 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 550 IM = 1, NM
                           M = MVAL( IM )
                           DO 540 IN = 1, NN
                              N = NVAL( IN )
                              IF( ISIDE.EQ.1 ) THEN
                                 CALL ZTIMMG( IMAT, M, M, A, LDA, 0, 0 )
                              ELSE
                                 CALL ZTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                              END IF
                              CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
  520                         CONTINUE
                              CALL ZTRMM( SIDE, UPLO, TRANSA,
     $                                    'Non-unit', M, N, ALPHA, A,
     $                                    LDA, B, LDA )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                                 GO TO 520
                              END IF
*
*                             Subtract the time used in ZTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
  530                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                                 GO TO 530
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
                              OPS = DOPBL3( CNAME, M, N, ISIDE-1 )
                              RESLTS( IM, IN, ILDA ) = DMFLOP( OPS,
     $                           TIME, 0 )
  540                      CONTINUE
  550                   CONTINUE
  560                CONTINUE
                     WRITE( NOUT, FMT = 9991 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), SIDE, UPLO, TRANSA
                     CALL DPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
  570             CONTINUE
  580          CONTINUE
  590       CONTINUE
*
*        Time ZTRSM
*
         ELSE IF( CNAME.EQ.'ZTRSM ' ) THEN
            DO 670 ISIDE = 1, NSIDES
               SIDE = SIDES( ISIDE )
               DO 660 IUPLO = 1, NUPLOS
                  UPLO = UPLOS( IUPLO )
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     IMAT = 11
                  ELSE
                     IMAT = -11
                  END IF
                  DO 650 ITA = 1, NTRANS
                     TRANSA = TRANS( ITA )
                     DO 640 ILDA = 1, NLDA
                        LDA = LDAVAL( ILDA )
                        DO 630 IM = 1, NM
                           M = MVAL( IM )
                           DO 620 IN = 1, NN
                              N = NVAL( IN )
                              IF( ISIDE.EQ.1 ) THEN
                                 CALL ZTIMMG( IMAT, M, M, A, LDA, 0, 0 )
                              ELSE
                                 CALL ZTIMMG( IMAT, N, N, A, LDA, 0, 0 )
                              END IF
                              CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
  600                         CONTINUE
                              CALL ZTRSM( SIDE, UPLO, TRANSA,
     $                                    'Non-unit', M, N, ALPHA, A,
     $                                    LDA, B, LDA )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                                 GO TO 600
                              END IF
*
*                             Subtract the time used in ZTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
  610                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL ZTIMMG( 0, M, N, B, LDA, 0, 0 )
                                 GO TO 610
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
                              OPS = DOPBL3( CNAME, M, N, ISIDE-1 )
                              RESLTS( IM, IN, ILDA ) = DMFLOP( OPS,
     $                           TIME, 0 )
  620                      CONTINUE
  630                   CONTINUE
  640                CONTINUE
                     WRITE( NOUT, FMT = 9991 )
     $     CNAME(1:ILA_LEN_TRIM(CNAME)), SIDE, UPLO, TRANSA
                     CALL DPRTBL( 'M', 'N', NM, MVAL, NN, NVAL, NLDA,
     $                            RESLTS, LDR1, LDR2, NOUT )
  650             CONTINUE
  660          CONTINUE
  670       CONTINUE
         END IF
         WRITE( NOUT, FMT = 9990 )
  680 CONTINUE
  690 CONTINUE
*
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9997 FORMAT( 5X, 'with LDA = ', I5 )
 9996 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
 9995 FORMAT( / 1X, 'ZGEMM  with TRANSA = ''', A1, ''', TRANSB = ''',
     $      A1, '''' )
 9994 FORMAT( / 1X, 'K = ', I4, / )
 9993 FORMAT( / 1X, A, ' with SIDE = ''', A1, ''', UPLO = ''', A1,
     $      '''', / )
 9992 FORMAT( / 1X, A, ' with UPLO = ''', A1, ''', TRANS = ''', A1,
     $      '''', / )
 9991 FORMAT( / 1X, A, ' with SIDE = ''', A1, ''', UPLO = ''', A1,
     $      ''',', ' TRANS = ''', A1, '''', / )
 9990 FORMAT( / / / / / )
      RETURN
*
*     End of ZTIMB3
*
      END
