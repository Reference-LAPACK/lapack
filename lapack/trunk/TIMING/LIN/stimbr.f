      SUBROUTINE STIMBR( LINE, NM, MVAL, NVAL, NK, KVAL, NNB, NBVAL,
     $                   NXVAL, NLDA, LDAVAL, TIMMIN, A, B, D, TAU,
     $                   WORK, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, LDR3, NK, NLDA, NM, NNB, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            KVAL( * ), LDAVAL( * ), MVAL( * ), NBVAL( * ),
     $                   NVAL( * ), NXVAL( * )
      REAL               A( * ), B( * ), D( * ),
     $                   RESLTS( LDR1, LDR2, LDR3, * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  STIMBR times SGEBRD, SORGBR, and SORMBR.
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
*          The number of values of M and N contained in the vectors
*          MVAL and NVAL.  The matrix sizes are used in pairs (M,N).
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix row dimension M.
*
*  NVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix column dimension N.
*
*  NK      (input) INTEGER
*          The number of values of K contained in the vector KVAL.
*
*  KVAL    (input) INTEGER array, dimension (NK)
*          The values of the matrix dimension K.
*
*  NNB     (input) INTEGER
*          The number of values of NB and NX contained in the
*          vectors NBVAL and NXVAL.  The blocking parameters are used
*          in pairs (NB,NX).
*
*  NBVAL   (input) INTEGER array, dimension (NNB)
*          The values of the blocksize NB.
*
*  NXVAL   (input) INTEGER array, dimension (NNB)
*          The values of the crossover point NX.
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
*          where LDAMAX and NMAX are the maximum values of LDA and N.
*
*  B       (workspace) REAL array, dimension (LDAMAX*NMAX)
*
*  D       (workspace) REAL array, dimension
*                      (2*max(min(M,N))-1)
*
*  TAU     (workspace) REAL array, dimension
*                      (2*max(min(M,N)))
*
*  WORK    (workspace) REAL array, dimension (LDAMAX*NBMAX)
*          where NBMAX is the maximum value of NB.
*
*  RESLTS  (output) REAL array, dimension (LDR1,LDR2,LDR3,6)
*          The timing results for each subroutine over the relevant
*          values of (M,N), (NB,NX), and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(1,NNB).
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
*  Internal Parameters
*  ===================
*
*  MODE    INTEGER
*          The matrix type.  MODE = 3 is a geometric distribution of
*          eigenvalues.  See CLATMS for further details.
*
*  COND    REAL
*          The condition number of the matrix.  The singular values are
*          set to values from DMAX to DMAX/COND.
*
*  DMAX    REAL
*          The magnitude of the largest singular value.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 3 )
      INTEGER            MODE
      REAL               COND, DMAX
      PARAMETER          ( MODE = 3, COND = 100.0E0, DMAX = 1.0E0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          LABK, LABM, LABN, SIDE, TRANS, VECT
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, I3, I4, IC, ICL, IK, ILDA, IM, INB, INFO,
     $                   INFO2, ISIDE, ISUB, ITOFF, ITRAN, IVECT, K, K1,
     $                   LDA, LW, M, M1, MINMN, N, N1, NB, NQ, NX
      REAL               OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER          SIDES( 2 ), TRANSS( 2 ), VECTS( 2 )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            ISEED( 4 ), RESEED( 4 )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      REAL               SECOND, SMFLOP, SOPLA
      EXTERNAL           SECOND, SMFLOP, SOPLA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, ICOPY, SGEBRD, SLACPY, SLATMS,
     $                   SORGBR, SORMBR, SPRTB4, SPRTB5, STIMMG, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'SGEBRD', 'SORGBR', 'SORMBR' / ,
     $                   SIDES / 'L', 'R' / , VECTS / 'Q', 'P' / ,
     $                   TRANSS / 'N', 'T' /
      DATA               ISEED / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'BR'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 220
*
*     Check that M <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 1, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 220
      END IF
*
*     Check that N <= LDA and K <= LDA for SORMBR
*
      IF( TIMSUB( 3 ) ) THEN
         CALL ATIMCK( 2, CNAME, NM, NVAL, NLDA, LDAVAL, NOUT, INFO )
         CALL ATIMCK( 3, CNAME, NK, KVAL, NLDA, LDAVAL, NOUT, INFO2 )
         IF( INFO.GT.0 .OR. INFO2.GT.0 ) THEN
            WRITE( NOUT, FMT = 9999 )
     $     SUBNAM( 3 )(1:ILA_LEN_TRIM( SUBNAM( 3 ) ))
            TIMSUB( 3 ) = .FALSE.
         END IF
      END IF
*
*     Do for each pair of values (M,N):
*
      DO 140 IM = 1, NM
         M = MVAL( IM )
         N = NVAL( IM )
         MINMN = MIN( M, N )
         CALL ICOPY( 4, ISEED, 1, RESEED, 1 )
*
*        Do for each value of LDA:
*
         DO 130 ILDA = 1, NLDA
            LDA = LDAVAL( ILDA )
*
*           Do for each pair of values (NB, NX) in NBVAL and NXVAL.
*
            DO 120 INB = 1, NNB
               NB = NBVAL( INB )
               CALL XLAENV( 1, NB )
               NX = NXVAL( INB )
               CALL XLAENV( 3, NX )
               LW = MAX( M+N, MAX( 1, NB )*( M+N ) )
*
*              Generate a test matrix of size M by N.
*
               CALL ICOPY( 4, RESEED, 1, ISEED, 1 )
               CALL SLATMS( M, N, 'Uniform', ISEED, 'Nonsym', TAU, MODE,
     $                      COND, DMAX, M, N, 'No packing', B, LDA,
     $                      WORK, INFO )
*
               IF( TIMSUB( 1 ) ) THEN
*
*                 SGEBRD:  Block reduction to bidiagonal form
*
                  CALL SLACPY( 'Full', M, N, B, LDA, A, LDA )
                  IC = 0
                  S1 = SECOND( )
   10             CONTINUE
                  CALL SGEBRD( M, N, A, LDA, D, D( MINMN ), TAU,
     $                         TAU( MINMN+1 ), WORK, LW, INFO )
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL SLACPY( 'Full', M, N, B, LDA, A, LDA )
                     GO TO 10
                  END IF
*
*                 Subtract the time used in SLACPY.
*
                  ICL = 1
                  S1 = SECOND( )
   20             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL SLACPY( 'Full', M, N, A, LDA, B, LDA )
                     GO TO 20
                  END IF
*
                  TIME = ( TIME-UNTIME ) / REAL( IC )
                  OPS = SOPLA( 'SGEBRD', M, N, 0, 0, NB )
                  RESLTS( INB, IM, ILDA, 1 ) = SMFLOP( OPS, TIME, INFO )
               ELSE
*
*                 If SGEBRD was not timed, generate a matrix and reduce
*                 it using SGEBRD anyway so that the orthogonal
*                 transformations may be used in timing the other
*                 routines.
*
                  CALL SLACPY( 'Full', M, N, B, LDA, A, LDA )
                  CALL SGEBRD( M, N, A, LDA, D, D( MINMN ), TAU,
     $                         TAU( MINMN+1 ), WORK, LW, INFO )
*
               END IF
*
               IF( TIMSUB( 2 ) ) THEN
*
*                 SORGBR:  Generate one of the orthogonal matrices Q or
*                 P' from the reduction to bidiagonal form
*                 A = Q * B * P'.
*
                  DO 50 IVECT = 1, 2
                     IF( IVECT.EQ.1 ) THEN
                        VECT = 'Q'
                        M1 = M
                        N1 = MIN( M, N )
                        K1 = N
                     ELSE
                        VECT = 'P'
                        M1 = MIN( M, N )
                        N1 = N
                        K1 = M
                     END IF
                     I3 = ( IVECT-1 )*NLDA
                     LW = MAX( 1, MAX( 1, NB )*MIN( M, N ) )
                     CALL SLACPY( 'Full', M, N, A, LDA, B, LDA )
                     IC = 0
                     S1 = SECOND( )
   30                CONTINUE
                     CALL SORGBR( VECT, M1, N1, K1, B, LDA, TAU, WORK,
     $                            LW, INFO )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL SLACPY( 'Full', M, N, A, LDA, B, LDA )
                        GO TO 30
                     END IF
*
*                    Subtract the time used in SLACPY.
*
                     ICL = 1
                     S1 = SECOND( )
   40                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL SLACPY( 'Full', M, N, A, LDA, B, LDA )
                        GO TO 40
                     END IF
*
                     TIME = ( TIME-UNTIME ) / REAL( IC )
*
*                    Op count for SORGBR:
*
                     IF( IVECT.EQ.1 ) THEN
                        IF( M1.GE.K1 ) THEN
                           OPS = SOPLA( 'SORGQR', M1, N1, K1, -1, NB )
                        ELSE
                           OPS = SOPLA( 'SORGQR', M1-1, M1-1, M1-1, -1,
     $                           NB )
                        END IF
                     ELSE
                        IF( K1.LT.N1 ) THEN
                           OPS = SOPLA( 'SORGLQ', M1, N1, K1, -1, NB )
                        ELSE
                           OPS = SOPLA( 'SORGLQ', N1-1, N1-1, N1-1, -1,
     $                           NB )
                        END IF
                     END IF
*
                     RESLTS( INB, IM, I3+ILDA, 2 ) = SMFLOP( OPS, TIME,
     $                  INFO )
   50             CONTINUE
               END IF
*
               IF( TIMSUB( 3 ) ) THEN
*
*                 SORMBR:  Multiply an m by n matrix B by one of the
*                 orthogonal matrices Q or P' from the reduction to
*                 bidiagonal form A = Q * B * P'.
*
                  DO 110 IVECT = 1, 2
                     IF( IVECT.EQ.1 ) THEN
                        VECT = 'Q'
                        K1 = N
                        NQ = M
                     ELSE
                        VECT = 'P'
                        K1 = M
                        NQ = N
                     END IF
                     I3 = ( IVECT-1 )*NLDA
                     I4 = 2
                     DO 100 ISIDE = 1, 2
                        SIDE = SIDES( ISIDE )
                        DO 90 IK = 1, NK
                           K = KVAL( IK )
                           IF( ISIDE.EQ.1 ) THEN
                              M1 = NQ
                              N1 = K
                              LW = MAX( 1, MAX( 1, NB )*N1 )
                           ELSE
                              M1 = K
                              N1 = NQ
                              LW = MAX( 1, MAX( 1, NB )*M1 )
                           END IF
                           ITOFF = 0
                           DO 80 ITRAN = 1, 2
                              TRANS = TRANSS( ITRAN )
                              CALL STIMMG( 0, M1, N1, B, LDA, 0, 0 )
                              IC = 0
                              S1 = SECOND( )
   60                         CONTINUE
                              CALL SORMBR( VECT, SIDE, TRANS, M1, N1,
     $                                     K1, A, LDA, TAU, B, LDA,
     $                                     WORK, LW, INFO )
                              S2 = SECOND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL STIMMG( 0, M1, N1, B, LDA, 0, 0 )
                                 GO TO 60
                              END IF
*
*                             Subtract the time used in STIMMG.
*
                              ICL = 1
                              S1 = SECOND( )
   70                         CONTINUE
                              S2 = SECOND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL STIMMG( 0, M1, N1, B, LDA, 0, 0 )
                                 GO TO 70
                              END IF
*
                              TIME = ( TIME-UNTIME ) / REAL( IC )
                              IF( IVECT.EQ.1 ) THEN
*
*                                Op count for SORMBR, VECT = 'Q':
*
                                 IF( NQ.GE.K1 ) THEN
                                    OPS = SOPLA( 'SORMQR', M1, N1, K1,
     $                                    ISIDE-1, NB )
                                 ELSE IF( ISIDE.EQ.1 ) THEN
                                    OPS = SOPLA( 'SORMQR', M1-1, N1,
     $                                    NQ-1, ISIDE-1, NB )
                                 ELSE
                                    OPS = SOPLA( 'SORMQR', M1, N1-1,
     $                                    NQ-1, ISIDE-1, NB )
                                 END IF
                              ELSE
*
*                                Op count for SORMBR, VECT = 'P':
*
                                 IF( NQ.GT.K1 ) THEN
                                    OPS = SOPLA( 'SORMLQ', M1, N1, K1,
     $                                    ISIDE-1, NB )
                                 ELSE IF( ISIDE.EQ.1 ) THEN
                                    OPS = SOPLA( 'SORMLQ', M1-1, N1,
     $                                    NQ-1, ISIDE-1, NB )
                                 ELSE
                                    OPS = SOPLA( 'SORMLQ', M1, N1-1,
     $                                    NQ-1, ISIDE-1, NB )
                                 END IF
                              END IF
*
                              RESLTS( INB, IM, I3+ILDA,
     $                           I4+ITOFF+IK ) = SMFLOP( OPS, TIME,
     $                           INFO )
                              ITOFF = NK
   80                      CONTINUE
   90                   CONTINUE
                        I4 = 2*NK + 2
  100                CONTINUE
  110             CONTINUE
               END IF
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
*
*     Print a table of results for each timed routine.
*
      DO 210 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 210
         WRITE( NOUT, FMT = 9998 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
         IF( NLDA.GT.1 ) THEN
            DO 150 I = 1, NLDA
               WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
  150       CONTINUE
         END IF
         IF( ISUB.EQ.1 ) THEN
            WRITE( NOUT, FMT = * )
            CALL SPRTB4( '(  NB,  NX)', 'M', 'N', NNB, NBVAL, NXVAL, NM,
     $                   MVAL, NVAL, NLDA, RESLTS( 1, 1, 1, ISUB ),
     $                   LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.2 ) THEN
            DO 160 IVECT = 1, 2
               I3 = ( IVECT-1 )*NLDA + 1
               IF( IVECT.EQ.1 ) THEN
                  LABK = 'N'
                  LABM = 'M'
                  LABN = 'K'
               ELSE
                  LABK = 'M'
                  LABM = 'K'
                  LABN = 'N'
               END IF
               WRITE( NOUT, FMT = 9996 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) )),
     $     VECTS( IVECT ),
     $            LABK, LABM, LABN
               CALL SPRTB4( '(  NB,  NX)', LABM, LABN, NNB, NBVAL,
     $                      NXVAL, NM, MVAL, NVAL, NLDA,
     $                      RESLTS( 1, 1, I3, ISUB ), LDR1, LDR2, NOUT )
  160       CONTINUE
         ELSE IF( ISUB.EQ.3 ) THEN
            DO 200 IVECT = 1, 2
               I3 = ( IVECT-1 )*NLDA + 1
               I4 = 3
               DO 190 ISIDE = 1, 2
                  IF( ISIDE.EQ.1 ) THEN
                     IF( IVECT.EQ.1 ) THEN
                        LABM = 'M'
                        LABN = 'K'
                     ELSE
                        LABM = 'K'
                        LABN = 'M'
                     END IF
                     LABK = 'N'
                  ELSE
                     IF( IVECT.EQ.1 ) THEN
                        LABM = 'N'
                        LABN = 'K'
                     ELSE
                        LABM = 'K'
                        LABN = 'N'
                     END IF
                     LABK = 'M'
                  END IF
                  DO 180 ITRAN = 1, 2
                     DO 170 IK = 1, NK
                        WRITE( NOUT, FMT = 9995 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) )),
     $                     VECTS( IVECT ), SIDES( ISIDE ),
     $                     TRANSS( ITRAN ), LABK, KVAL( IK )
                        CALL SPRTB5( 'NB', LABM, LABN, NNB, NBVAL, NM,
     $                               MVAL, NVAL, NLDA,
     $                               RESLTS( 1, 1, I3, I4 ), LDR1, LDR2,
     $                               NOUT )
                        I4 = I4 + 1
  170                CONTINUE
  180             CONTINUE
  190          CONTINUE
  200       CONTINUE
         END IF
  210 CONTINUE
  220 CONTINUE
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
 9996 FORMAT( / 5X, A, ' with VECT = ''', A1, ''', ', A1, ' = MIN(',
     $      A1, ',', A1, ')', / )
 9995 FORMAT( / 5X, A, ' with VECT = ''', A1, ''', SIDE = ''', A1,
     $      ''', TRANS = ''', A1, ''', ', A1, ' =', I6, / )
      RETURN
*
*     End of STIMBR
*
      END
