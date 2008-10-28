      SUBROUTINE DTIM21( LINE, NSIZES, NN, NTYPES, DOTYPE, NPARMS, NNB,
     $                   NSHFTS, MAXBS, LDAS, TIMMIN, NOUT, ISEED, A, H,
     $                   Z, W, WORK, LWORK, LLWORK, IWORK, TIMES, LDT1,
     $                   LDT2, LDT3, OPCNTS, LDO1, LDO2, LDO3, INFO )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            INFO, LDO1, LDO2, LDO3, LDT1, LDT2, LDT3,
     $                   LWORK, NOUT, NPARMS, NSIZES, NTYPES
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * ), LLWORK( * )
      INTEGER            ISEED( * ), IWORK( * ), LDAS( * ), MAXBS( * ),
     $                   NN( * ), NNB( * ), NSHFTS( * )
      DOUBLE PRECISION   A( * ), H( * ), OPCNTS( LDO1, LDO2, LDO3, * ),
     $                   TIMES( LDT1, LDT2, LDT3, * ), W( * ),
     $                   WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*     DTIM21 times the LAPACK routines for the DOUBLE PRECISION
*     non-symmetric eigenvalue problem.
*
*     For each N value in NN(1:NSIZES) and .TRUE. value in
*     DOTYPE(1:NTYPES), a matrix will be generated and used to test the
*     selected routines.  Thus, NSIZES*(number of .TRUE. values in
*     DOTYPE) matrices will be generated.
*
*  Arguments
*  =========
*
*  LINE    (input) CHARACTER*80
*          On entry, LINE contains the input line which requested
*          this routine.  This line may contain a subroutine name,
*          such as DGEHRD, indicating that only routine SGEHRD will
*          be timed, or it may contain a generic name, such as DHS.
*          In this case, the rest of the line is scanned for the
*          first 12 non-blank characters, corresponding to the twelve
*          combinations of subroutine and options:
*          LAPACK:
*          1: DGEHRD
*          2: DHSEQR(JOB='E')
*          3: DHSEQR(JOB='S')
*          4: DHSEQR(JOB='I')
*          5: DTREVC(JOB='L')
*          6: DTREVC(JOB='R')
*          7: DHSEIN(JOB='L')
*          8: DHSEIN(JOB='R')
*          EISPACK:
*           9: ORTHES (compare with DGEHRD)
*          10: HQR    (compare w/ DHSEQR -- JOB='E')
*          11: HQR2   (compare w/ DHSEQR(JOB='I') plus DTREVC(JOB='R'))
*          12: INVIT  (compare with DHSEIN)
*          If a character is 'T' or 't', the corresponding routine in
*          this path is timed.  If the entire line is blank, all the
*          routines in the path are timed.
*
*  NSIZES  (input) INTEGER
*          The number of values of N contained in the vector NN.
*
*  NN      (input) INTEGER array, dimension( NSIZES )
*          The values of the matrix size N to be tested.  For each
*          N value in the array NN, and each .TRUE. value in DOTYPE,
*          a matrix A will be generated and used to test the routines.
*
*  NTYPES  (input) INTEGER
*          The number of types in DOTYPE.  Only the first MAXTYP
*          elements will be examined.  Exception: if NSIZES=1 and
*          NTYPES=MAXTYP+1, and DOTYPE=MAXTYP*f,t, then the input
*          value of A will be used.
*
*  DOTYPE  (input) LOGICAL
*          If DOTYPE(j) is .TRUE., then a matrix of type j will be
*          generated.  The matrix A has the form X**(-1) T X, where
*          X is orthogonal (for j=1--4) or has condition sqrt(ULP)
*          (for j=5--8), and T has random O(1) entries in the upper
*          triangle and:
*          (j=1,5) evenly spaced entries 1, ..., ULP with random signs
*          (j=2,6) geometrically spaced entries 1, ..., ULP with random
*                  signs
*          (j=3,7) "clustered" entries 1, ULP,..., ULP with random
*                  signs
*          (j=4,8) real or complex conjugate paired eigenvalues
*                  randomly chosen from ( ULP, 1 )
*          on the diagonal.
*
*  NPARMS  (input) INTEGER
*          The number of values in each of the arrays NNB, NSHFTS,
*          MAXBS, and LDAS.  For each matrix A generated according to
*          NN and DOTYPE, tests will be run with (NB,NSHIFT,MAXB,LDA)=
*          (NNB(1), NSHFTS(1), MAXBS(1), LDAS(1)),...,
*          (NNB(NPARMS), NSHFTS(NPARMS), MAXBS(NPARMS), LDAS(NPARMS))
*
*  NNB     (input) INTEGER array, dimension( NPARMS )
*          The values of the blocksize ("NB") to be tested.
*
*  NSHFTS  (input) INTEGER array, dimension( NPARMS )
*          The values of the number of shifts ("NSHIFT") to be tested.
*
*  MAXBS   (input) INTEGER array, dimension( NPARMS )
*          The values of "MAXB", the size of largest submatrix to be
*          processed by DLAHQR (EISPACK method), to be tested.
*
*  LDAS    (input) INTEGER array, dimension( NPARMS )
*          The values of LDA, the leading dimension of all matrices,
*          to be tested.
*
*  TIMMIN  (input) DOUBLE PRECISION
*          The minimum time a subroutine will be timed.
*
*  NOUT    (input) INTEGER
*          If NOUT > 0 then NOUT specifies the unit number
*          on which the output will be printed.  If NOUT <= 0, no
*          output is printed.
*
*  ISEED   (input/output) INTEGER array, dimension( 4 )
*          The random seed used by the random number generator, used
*          by the test matrix generator.  It is used and updated on
*          each call to DTIM21
*
*  A       (workspace) DOUBLE PRECISION array,
*                      dimension( max(NN)*max(LDAS) )
*          (a) During the testing of DGEHRD, the original matrix to
*              be tested.
*          (b) Later, the Schur form of the original matrix.
*
*  H       (workspace) DOUBLE PRECISION array,
*                      dimension( max(NN)*max(LDAS) )
*          The Hessenberg form of the original matrix.
*
*  Z       (workspace) DOUBLE PRECISION array,
*                      dimension( max(NN)*max(LDAS) )
*          Various output arrays: from DGEHRD and DHSEQR, the
*          orthogonal reduction matrices; from DTREVC and DHSEIN,
*          the eigenvector matrices.
*
*  W       (workspace) DOUBLE PRECISION array,
*                      dimension( 2*max(LDAS) )
*          Treated as an LDA x 2 matrix whose 1st column holds WR, the
*          real parts of the eigenvalues, and whose 2nd column holds
*          WI, the imaginary parts of the eigenvalues of A.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension( LWORK )
*
*  LWORK   (input) INTEGER
*          Number of elements in WORK.  It must be at least
*          (a)  max(NN)*( 3*max(NNB) + 2 )
*          (b)  max(NN)*( max(NNB+NSHFTS) + 1 )
*          (c)  max(NSHFTS)*( max(NSHFTS) + max(NN) )
*          (d)  max(MAXBS)*( max(MAXBS) + max(NN) )
*          (e)  ( max(NN) + 2 )**2  +  max(NN)
*          (f)  NSIZES*NTYPES*NPARMS
*
*  LLWORK  (workspace) LOGICAL array, dimension( max( max(NN), NPARMS ))
*
*  IWORK   (workspace) INTEGER array, dimension( 2*max(NN) )
*          Workspace needed for parameters IFAILL and IFAILR in call
*          to DHSEIN.
*
*  TIMES   (output) DOUBLE PRECISION array,
*                   dimension (LDT1,LDT2,LDT3,NSUBS)
*          TIMES(i,j,k,l) will be set to the run time (in seconds) for
*          subroutine l, with N=NN(k), matrix type j, and LDA=LDAS(i),
*          MAXB=MAXBS(i), NBLOCK=NNB(i), and NSHIFT=NSHFTS(i).
*
*  LDT1    (input) INTEGER
*          The first dimension of TIMES.  LDT1 >= min( 1, NPARMS ).
*
*  LDT2    (input) INTEGER
*          The second dimension of TIMES.  LDT2 >= min( 1, NTYPES ).
*
*  LDT3    (input) INTEGER
*          The third dimension of TIMES.  LDT3 >= min( 1, NSIZES ).
*
*  OPCNTS  (output) DOUBLE PRECISION array,
*                   dimension (LDO1,LDO2,LDO3,NSUBS)
*          OPCNTS(i,j,k,l) will be set to the number of floating-point
*          operations executed by subroutine l, with N=NN(k), matrix
*          type j, and LDA=LDAS(i), MAXB=MAXBS(i), NBLOCK=NNB(i), and
*          NSHIFT=NSHFTS(i).
*
*  LDO1    (input) INTEGER
*          The first dimension of OPCNTS.  LDO1 >= min( 1, NPARMS ).
*
*  LDO2    (input) INTEGER
*          The second dimension of OPCNTS.  LDO2 >= min( 1, NTYPES ).
*
*  LDO3    (input) INTEGER
*          The third dimension of OPCNTS.  LDO3 >= min( 1, NSIZES ).
*
*  INFO    (output) INTEGER
*          Error flag.  It will be set to zero if no error occurred.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXTYP, NSUBS
      PARAMETER          ( MAXTYP = 8, NSUBS = 12 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            RUNHQR, RUNHRD, RUNORT, RUNQRE, RUNQRS
      INTEGER            IC, ICONDS, IINFO, IMODE, IN, IPAR, ISUB,
     $                   ITEMP, ITYPE, J, J1, J2, J3, J4, JC, JR, LASTL,
     $                   LASTNL, LDA, LDAMIN, LDH, LDT, LDW, MAXB,
     $                   MBMAX, MTYPES, N, NB, NBMAX, NMAX, NSBMAX,
     $                   NSHIFT, NSMAX
      DOUBLE PRECISION   CONDS, RTULP, RTULPI, S1, S2, TIME, ULP,
     $                   ULPINV, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER          ADUMMA( 1 )
      CHARACTER*4        PNAMES( 4 )
      CHARACTER*9        SUBNAM( NSUBS )
      INTEGER            INPARM( NSUBS ), IOLDSD( 4 ), KCONDS( MAXTYP ),
     $                   KMODE( MAXTYP )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DOPLA, DSECND
      EXTERNAL           DLAMCH, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMIN, DGEHRD, DHSEIN, DHSEQR, DLACPY, DLASET,
     $                   DLATME, DPRTBE, DTREVC, HQR, HQR2, INVIT,
     $                   ORTHES, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Scalars in Common ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'DGEHRD', 'DHSEQR(E)', 'DHSEQR(S)',
     $                   'DHSEQR(V)', 'DTREVC(L)', 'DTREVC(R)',
     $                   'DHSEIN(L)', 'DHSEIN(R)', 'ORTHES', 'HQR',
     $                   'HQR2', 'INVIT' /
      DATA               INPARM / 2, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1 /
      DATA               PNAMES / 'LDA', 'NB', 'NS', 'MAXB' /
      DATA               KMODE / 4, 3, 1, 5, 4, 3, 1, 5 /
      DATA               KCONDS / 4*1, 4*2 /
*     ..
*     .. Executable Statements ..
*
*     Quick Return
*
      INFO = 0
      IF( NSIZES.LE.0 .OR. NTYPES.LE.0 .OR. NPARMS.LE.0 )
     $   RETURN
*
*     Extract the timing request from the input line.
*
      CALL ATIMIN( 'DHS', LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   RETURN
*
*     Compute Maximum Values
*
      NMAX = 0
      DO 10 J1 = 1, NSIZES
         NMAX = MAX( NMAX, NN( J1 ) )
   10 CONTINUE
*
      LDAMIN = 2*MAX( 1, NMAX )
      NBMAX = 0
      NSMAX = 0
      MBMAX = 0
      NSBMAX = 0
      DO 20 J1 = 1, NPARMS
         LDAMIN = MIN( LDAMIN, LDAS( J1 ) )
         NBMAX = MAX( NBMAX, NNB( J1 ) )
         NSMAX = MAX( NSMAX, NSHFTS( J1 ) )
         MBMAX = MAX( MBMAX, MAXBS( J1 ) )
         NSBMAX = MAX( NSBMAX, NNB( J1 )+NSHFTS( J1 ) )
   20 CONTINUE
*
*     Check that N <= LDA for the input values.
*
      IF( NMAX.GT.LDAMIN ) THEN
         INFO = -10
         WRITE( NOUT, FMT = 9999 )LINE( 1: 6 )
 9999    FORMAT( 1X, A, ' timing run not attempted -- N > LDA', / )
         RETURN
      END IF
*
*     Check LWORK
*
      IF( LWORK.LT.MAX( NMAX*MAX( 3*NBMAX+2, NSBMAX+1 ),
     $    NSMAX*( NSMAX+NMAX ), MBMAX*( MBMAX+NMAX ),
     $    ( NMAX+1 )*( NMAX+4 ), NSIZES*NTYPES*NPARMS ) ) THEN
         INFO = -19
         WRITE( NOUT, FMT = 9998 )LINE( 1: 6 )
 9998    FORMAT( 1X, A, ' timing run not attempted -- LWORK too small.',
     $         / )
         RETURN
      END IF
*
*     Check to see whether DGEHRD or DHSEQR must be run.
*
*     RUNQRE -- if DHSEQR must be run to get eigenvalues.
*     RUNQRS -- if DHSEQR must be run to get Schur form.
*     RUNHRD -- if DGEHRD must be run.
*
      RUNQRS = .FALSE.
      RUNQRE = .FALSE.
      RUNHRD = .FALSE.
      IF( TIMSUB( 5 ) .OR. TIMSUB( 6 ) )
     $   RUNQRS = .TRUE.
      IF( ( TIMSUB( 7 ) .OR. TIMSUB( 8 ) ) )
     $   RUNQRE = .TRUE.
      IF( TIMSUB( 2 ) .OR. TIMSUB( 3 ) .OR. TIMSUB( 4 ) .OR. RUNQRS .OR.
     $    RUNQRE )RUNHRD = .TRUE.
      IF( TIMSUB( 3 ) .OR. TIMSUB( 4 ) .OR. RUNQRS )
     $   RUNQRE = .FALSE.
      IF( TIMSUB( 4 ) )
     $   RUNQRS = .FALSE.
*
*     Check to see whether ORTHES or HQR must be run.
*
*     RUNHQR -- if HQR must be run to get eigenvalues.
*     RUNORT -- if ORTHES must be run.
*
      RUNHQR = .FALSE.
      RUNORT = .FALSE.
      IF( TIMSUB( 12 ) )
     $   RUNHQR = .TRUE.
      IF( TIMSUB( 10 ) .OR. TIMSUB( 11 ) .OR. RUNHQR )
     $   RUNORT = .TRUE.
      IF( TIMSUB( 10 ) .OR. TIMSUB( 11 ) )
     $   RUNHQR = .FALSE.
      IF( TIMSUB( 9 ) )
     $   RUNORT = .FALSE.
*
*     Various Constants
*
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTULP = SQRT( ULP )
      RTULPI = ONE / RTULP
*
*     Zero out OPCNTS, TIMES
*
      DO 60 J4 = 1, NSUBS
         DO 50 J3 = 1, NSIZES
            DO 40 J2 = 1, NTYPES
               DO 30 J1 = 1, NPARMS
                  OPCNTS( J1, J2, J3, J4 ) = ZERO
                  TIMES( J1, J2, J3, J4 ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
*
*     Do for each value of N:
*
      DO 620 IN = 1, NSIZES
*
         N = NN( IN )
*
*        Do for each .TRUE. value in DOTYPE:
*
         MTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.EQ.MAXTYP+1 .AND. NSIZES.EQ.1 )
     $      MTYPES = NTYPES
         DO 610 ITYPE = 1, MTYPES
            IF( .NOT.DOTYPE( ITYPE ) )
     $         GO TO 610
*
*           Save random number seed for error messages
*
            DO 70 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   70       CONTINUE
*
*-----------------------------------------------------------------------
*
*           Time the LAPACK Routines
*
*           Generate A
*
            IF( ITYPE.LE.MAXTYP ) THEN
               IMODE = KMODE( ITYPE )
               ICONDS = KCONDS( ITYPE )
               IF( ICONDS.EQ.1 ) THEN
                  CONDS = ONE
               ELSE
                  CONDS = RTULPI
               END IF
               ADUMMA( 1 ) = ' '
               CALL DLATME( N, 'S', ISEED, WORK, IMODE, ULPINV, ONE,
     $                      ADUMMA, 'T', 'T', 'T', WORK( N+1 ), 4,
     $                      CONDS, N, N, ONE, A, N, WORK( 2*N+1 ),
     $                      IINFO )
            END IF
*
*           Time DGEHRD for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 1 ) ) THEN
               DO 110 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
*
*                 If this combination of (NB,LDA) has occurred before,
*                 just use that value.
*
                  LASTNL = 0
                  DO 80 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) .AND. NB.EQ.
     $                   MIN( N, NNB( J ) ) )LASTNL = J
   80             CONTINUE
*
                  IF( LASTNL.EQ.0 ) THEN
                     CALL XLAENV( 1, NB )
                     CALL XLAENV( 2, 2 )
                     CALL XLAENV( 3, NB )
*
*                    Time DGEHRD
*
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
   90                CONTINUE
                     CALL DLACPY( 'Full', N, N, A, N, H, LDA )
*
                     CALL DGEHRD( N, 1, N, H, LDA, WORK, WORK( N+1 ),
     $                            LWORK-N, IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 1 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 610
                     END IF
*
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 90
*
*                    Subtract the time used in DLACPY.
*
                     S1 = DSECND( )
                     DO 100 J = 1, IC
                        CALL DLACPY( 'Full', N, N, A, N, Z, LDA )
  100                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 1 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 1 ) = DOPLA( 'DGEHRD', N,
     $                  1, N, 0, NB )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 1 ) = OPCNTS( LASTNL,
     $                  ITYPE, IN, 1 )
                     TIMES( IPAR, ITYPE, IN, 1 ) = TIMES( LASTNL, ITYPE,
     $                  IN, 1 )
                  END IF
  110          CONTINUE
               LDH = LDA
            ELSE
               IF( RUNHRD ) THEN
                  CALL DLACPY( 'Full', N, N, A, N, H, N )
*
                  CALL DGEHRD( N, 1, N, H, N, WORK, WORK( N+1 ),
     $                         LWORK-N, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 1 ), IINFO, N,
     $                  ITYPE, 0, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 610
                  END IF
                  LDH = N
               END IF
            END IF
*
*           Time DHSEQR with JOB='E' for each 4-tuple
*           NNB(j), NSHFTS(j), MAXBS(j), LDAS(j)
*
            IF( TIMSUB( 2 ) ) THEN
               DO 140 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = 1
                  NSHIFT = NSHFTS( IPAR )
                  MAXB = MAXBS( IPAR )
                  CALL XLAENV( 4, NSHIFT )
                  CALL XLAENV( 8, MAXB )
*
*                 Time DHSEQR with JOB='E'
*
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
  120             CONTINUE
                  CALL DLACPY( 'Full', N, N, H, LDH, A, LDA )
*
                  CALL DHSEQR( 'E', 'N', N, 1, N, A, LDA, W, W( LDA+1 ),
     $                         Z, LDA, WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 2 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 610
                  END IF
*
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 120
*
*                 Subtract the time used in DLACPY.
*
                  S1 = DSECND( )
                  DO 130 J = 1, IC
                     CALL DLACPY( 'Full', N, N, H, LDH, Z, LDA )
  130             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 2 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 2 ) = OPS / DBLE( IC )
  140          CONTINUE
               LDT = 0
               LDW = LDA
            ELSE
               IF( RUNQRE ) THEN
                  CALL DLACPY( 'Full', N, N, H, LDH, A, N )
*
                  CALL DHSEQR( 'E', 'N', N, 1, N, A, N, W, W( N+1 ), Z,
     $                         N, WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 2 ), IINFO, N,
     $                  ITYPE, 0, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 610
                  END IF
                  LDT = 0
                  LDW = N
               END IF
            END IF
*
*           Time DHSEQR with JOB='S' for each 4-tuple
*           NNB(j), NSHFTS(j), MAXBS(j), LDAS(j)
*
            IF( TIMSUB( 3 ) ) THEN
               DO 170 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NSHIFT = NSHFTS( IPAR )
                  MAXB = MAXBS( IPAR )
                  NB = 1
                  CALL XLAENV( 4, NSHIFT )
                  CALL XLAENV( 8, MAXB )
*
*                 Time DHSEQR with JOB='S'
*
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
  150             CONTINUE
                  CALL DLACPY( 'Full', N, N, H, LDH, A, LDA )
*
                  CALL DHSEQR( 'S', 'N', N, 1, N, A, LDA, W, W( LDA+1 ),
     $                         Z, LDA, WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 3 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 610
                  END IF
*
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 150
*
*                 Subtract the time used in DLACPY.
*
                  S1 = DSECND( )
                  DO 160 J = 1, IC
                     CALL DLACPY( 'Full', N, N, H, LDH, Z, LDA )
  160             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 3 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 3 ) = OPS / DBLE( IC )
  170          CONTINUE
               LDT = LDA
               LDW = LDA
            ELSE
               IF( RUNQRS ) THEN
                  CALL DLACPY( 'Full', N, N, H, LDH, A, N )
*
                  CALL DHSEQR( 'S', 'N', N, 1, N, A, N, W, W( N+1 ), Z,
     $                         N, WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 3 ), IINFO, N,
     $                  ITYPE, 0, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 610
                  END IF
                  LDT = N
                  LDW = N
               END IF
            END IF
*
*           Time DHSEQR with JOB='I' for each 4-tuple
*           NNB(j), NSHFTS(j), MAXBS(j), LDAS(j)
*
            IF( TIMSUB( 4 ) ) THEN
               DO 200 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NSHIFT = NSHFTS( IPAR )
                  MAXB = MAXBS( IPAR )
                  NB = 1
                  CALL XLAENV( 4, NSHIFT )
                  CALL XLAENV( 8, MAXB )
*
*                 Time DHSEQR with JOB='I'
*
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
  180             CONTINUE
                  CALL DLACPY( 'Full', N, N, H, LDH, A, LDA )
*
                  CALL DHSEQR( 'S', 'I', N, 1, N, A, LDA, W, W( LDA+1 ),
     $                         Z, LDA, WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 4 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 610
                  END IF
*
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 180
*
*                 Subtract the time used in DLACPY.
*
                  S1 = DSECND( )
                  DO 190 J = 1, IC
                     CALL DLACPY( 'Full', N, N, H, LDH, Z, LDA )
  190             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 4 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 4 ) = OPS / DBLE( IC )
  200          CONTINUE
               LDT = LDA
               LDW = LDA
            END IF
*
*           Time DTREVC and DHSEIN with various values of LDA
*
*           Select All Eigenvectors
*
            DO 210 J = 1, N
               LLWORK( J ) = .TRUE.
  210       CONTINUE
*
            DO 370 IPAR = 1, NPARMS
               LDA = LDAS( IPAR )
*
*              If this value of LDA has come up before, just use
*              the value previously computed.
*
               LASTL = 0
               DO 220 J = 1, IPAR - 1
                  IF( LDA.EQ.LDAS( J ) )
     $               LASTL = J
  220          CONTINUE
*
*              Time DTREVC
*
               IF( ( TIMSUB( 5 ) .OR. TIMSUB( 6 ) ) .AND. LASTL.EQ.0 )
     $              THEN
*
*                 Copy T (which is in A) if necessary to get right LDA.
*
                  IF( LDA.GT.LDT ) THEN
                     DO 240 JC = N, 1, -1
                        DO 230 JR = N, 1, -1
                           A( JR+( JC-1 )*LDA ) = A( JR+( JC-1 )*LDT )
  230                   CONTINUE
  240                CONTINUE
                  ELSE IF( LDA.LT.LDT ) THEN
                     DO 260 JC = 1, N
                        DO 250 JR = 1, N
                           A( JR+( JC-1 )*LDA ) = A( JR+( JC-1 )*LDT )
  250                   CONTINUE
  260                CONTINUE
                  END IF
                  LDT = LDA
*
*                 Time DTREVC for Left Eigenvectors
*
                  IF( TIMSUB( 5 ) ) THEN
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  270                CONTINUE
*
                     CALL DTREVC( 'L', 'A', LLWORK, N, A, LDA, Z, LDA,
     $                            Z, LDA, N, ITEMP, WORK, IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 5 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 610
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 270
*
                     TIMES( IPAR, ITYPE, IN, 5 ) = TIME / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 5 ) = OPS / DBLE( IC )
                  END IF
*
*                 Time DTREVC for Right Eigenvectors
*
                  IF( TIMSUB( 6 ) ) THEN
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  280                CONTINUE
                     CALL DTREVC( 'R', 'A', LLWORK, N, A, LDA, Z, LDA,
     $                            Z, LDA, N, ITEMP, WORK, IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 6 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 610
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 280
*
                     TIMES( IPAR, ITYPE, IN, 6 ) = TIME / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 6 ) = OPS / DBLE( IC )
                  END IF
               ELSE
                  IF( TIMSUB( 5 ) ) THEN
                     OPCNTS( IPAR, ITYPE, IN, 5 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 5 )
                     TIMES( IPAR, ITYPE, IN, 5 ) = TIMES( LASTL, ITYPE,
     $                  IN, 5 )
                  END IF
                  IF( TIMSUB( 6 ) ) THEN
                     OPCNTS( IPAR, ITYPE, IN, 6 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 6 )
                     TIMES( IPAR, ITYPE, IN, 6 ) = TIMES( LASTL, ITYPE,
     $                  IN, 6 )
                  END IF
               END IF
*
*              Time DHSEIN
*
               IF( ( TIMSUB( 7 ) .OR. TIMSUB( 8 ) ) .AND. LASTL.EQ.0 )
     $              THEN
*
*                 Copy H if necessary to get right LDA.
*
                  IF( LDA.GT.LDH ) THEN
                     DO 300 JC = N, 1, -1
                        DO 290 JR = N, 1, -1
                           H( JR+( JC-1 )*LDA ) = H( JR+( JC-1 )*LDH )
  290                   CONTINUE
                        W( JC+LDA ) = W( JC+LDH )
  300                CONTINUE
                  ELSE IF( LDA.LT.LDH ) THEN
                     DO 320 JC = 1, N
                        DO 310 JR = 1, N
                           H( JR+( JC-1 )*LDA ) = H( JR+( JC-1 )*LDH )
  310                   CONTINUE
                        W( JC+LDA ) = W( JC+LDH )
  320                CONTINUE
                  END IF
                  LDH = LDA
*
*                 Copy W if necessary to get right LDA.
*
                  IF( LDA.GT.LDW ) THEN
                     DO 330 J = N, 1, -1
                        W( J+LDA ) = W( J+LDW )
  330                CONTINUE
                  ELSE IF( LDA.LT.LDW ) THEN
                     DO 340 J = 1, N
                        W( J+LDA ) = W( J+LDW )
  340                CONTINUE
                  END IF
                  LDW = LDA
*
*                 Time DHSEIN for Left Eigenvectors
*
                  IF( TIMSUB( 7 ) ) THEN
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  350                CONTINUE
*
                     CALL DHSEIN( 'L', 'Q', 'N', LLWORK, N, H, LDA, W,
     $                            W( LDA+1 ), Z, LDA, Z, LDA, N, ITEMP,
     $                            WORK, IWORK, IWORK( N+1 ), IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 7 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 610
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 350
*
                     TIMES( IPAR, ITYPE, IN, 7 ) = TIME / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 7 ) = OPS / DBLE( IC )
                  END IF
*
*                 Time DHSEIN for Right Eigenvectors
*
                  IF( TIMSUB( 8 ) ) THEN
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  360                CONTINUE
*
                     CALL DHSEIN( 'R', 'Q', 'N', LLWORK, N, H, LDA, W,
     $                            W( LDA+1 ), Z, LDA, Z, LDA, N, ITEMP,
     $                            WORK, IWORK, IWORK( N+1 ), IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 8 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 610
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 360
*
                     TIMES( IPAR, ITYPE, IN, 8 ) = TIME / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 8 ) = OPS / DBLE( IC )
                  END IF
               ELSE
                  IF( TIMSUB( 7 ) ) THEN
                     OPCNTS( IPAR, ITYPE, IN, 7 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 7 )
                     TIMES( IPAR, ITYPE, IN, 7 ) = TIMES( LASTL, ITYPE,
     $                  IN, 7 )
                  END IF
                  IF( TIMSUB( 8 ) ) THEN
                     OPCNTS( IPAR, ITYPE, IN, 8 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 8 )
                     TIMES( IPAR, ITYPE, IN, 8 ) = TIMES( LASTL, ITYPE,
     $                  IN, 8 )
                  END IF
               END IF
  370       CONTINUE
*
*-----------------------------------------------------------------------
*
*           Time the EISPACK Routines
*
*           Restore random number seed
*
            DO 380 J = 1, 4
               ISEED( J ) = IOLDSD( J )
  380       CONTINUE
*
*           Re-generate A
*
            IF( ITYPE.LE.MAXTYP ) THEN
               IMODE = KMODE( ITYPE )
               IF( ICONDS.EQ.1 ) THEN
                  CONDS = ONE
               ELSE
                  CONDS = RTULPI
               END IF
               CALL DLATME( N, 'S', ISEED, WORK, IMODE, ULPINV, ONE,
     $                      ADUMMA, 'T', 'T', 'T', WORK( N+1 ), 4,
     $                      CONDS, N, N, ONE, A, N, WORK( 2*N+1 ),
     $                      IINFO )
            END IF
*
*           Time ORTHES for each LDAS(j)
*
            IF( TIMSUB( 9 ) ) THEN
               DO 420 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 390 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  390             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time ORTHES
*
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
*
  400                CONTINUE
                     CALL DLACPY( 'Full', N, N, A, N, H, LDA )
*
                     CALL ORTHES( LDA, N, 1, N, H, WORK )
*
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 400
*
*                    Subtract the time used in DLACPY.
*
                     S1 = DSECND( )
                     DO 410 J = 1, IC
                        CALL DLACPY( 'Full', N, N, A, N, Z, LDA )
  410                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
*
*                     OPS1 = ( 20*N**3 - 3*N**2 - 23*N ) / 6 - 17
*
                     TIMES( IPAR, ITYPE, IN, 9 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 9 ) = OPS / DBLE( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 9 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 9 )
                     TIMES( IPAR, ITYPE, IN, 9 ) = TIMES( LASTL, ITYPE,
     $                  IN, 9 )
                  END IF
                  LDH = LDA
  420          CONTINUE
            ELSE
               IF( RUNORT ) THEN
                  CALL DLACPY( 'Full', N, N, A, N, H, N )
*
                  CALL ORTHES( N, N, 1, N, H, WORK )
*
                  LDH = N
               END IF
            END IF
*
*           Time HQR for each LDAS(j)
*
            IF( TIMSUB( 10 ) ) THEN
               DO 460 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 430 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  430             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time HQR
*
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  440                CONTINUE
                     CALL DLACPY( 'Full', N, N, H, LDH, A, LDA )
*
                     CALL HQR( LDA, N, 1, N, A, W, W( LDA+1 ), IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 10 ), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 610
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 440
*
*                    Subtract the time used in DLACPY.
*
                     S1 = DSECND( )
                     DO 450 J = 1, IC
                        CALL DLACPY( 'Full', N, N, H, LDH, Z, LDA )
  450                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 10 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 10 ) = OPS / DBLE( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 10 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 10 )
                     TIMES( IPAR, ITYPE, IN, 10 ) = TIMES( LASTL, ITYPE,
     $                  IN, 10 )
                  END IF
                  LDW = LDA
  460          CONTINUE
            ELSE
               IF( RUNHQR ) THEN
                  CALL DLACPY( 'Full', N, N, A, N, H, N )
*
                  CALL HQR( N, N, 1, N, A, W, W( N+1 ), IINFO )
*
                  LDW = N
               END IF
            END IF
*
*           Time HQR2 for each LDAS(j)
*
            IF( TIMSUB( 11 ) ) THEN
               DO 500 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 470 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  470             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time HQR2
*
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  480                CONTINUE
                     CALL DLACPY( 'Full', N, N, H, LDH, A, LDA )
                     CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDA )
*
                     CALL HQR2( LDA, N, 1, N, A, W, W( LDA+1 ), Z,
     $                          IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 11 ), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 610
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 480
*
*                    Subtract the time used in DLACPY.
*
                     S1 = DSECND( )
                     DO 490 J = 1, IC
                        CALL DLACPY( 'Full', N, N, H, LDH, Z, LDA )
  490                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 11 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 11 ) = OPS / DBLE( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 11 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 11 )
                     TIMES( IPAR, ITYPE, IN, 11 ) = TIMES( LASTL, ITYPE,
     $                  IN, 11 )
                  END IF
                  LDW = LDA
  500          CONTINUE
            END IF
*
*           Time INVIT for each LDAS(j)
*
*           Select All Eigenvectors
*
            DO 510 J = 1, N
               LLWORK( J ) = .TRUE.
  510       CONTINUE
*
            IF( TIMSUB( 12 ) ) THEN
               DO 600 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 520 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  520             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Copy H if necessary to get right LDA.
*
                     IF( LDA.GT.LDH ) THEN
                        DO 540 JC = N, 1, -1
                           DO 530 JR = N, 1, -1
                              H( JR+( JC-1 )*LDA ) = H( JR+( JC-1 )*
     $                           LDH )
  530                      CONTINUE
  540                   CONTINUE
                     ELSE IF( LDA.LT.LDH ) THEN
                        DO 560 JC = 1, N
                           DO 550 JR = 1, N
                              H( JR+( JC-1 )*LDA ) = H( JR+( JC-1 )*
     $                           LDH )
  550                      CONTINUE
  560                   CONTINUE
                     END IF
                     LDH = LDA
*
*                    Copy W if necessary to get right LDA.
*
                     IF( LDA.GT.LDW ) THEN
                        DO 570 J = N, 1, -1
                           W( J+LDA ) = W( J+LDW )
  570                   CONTINUE
                     ELSE IF( LDA.LT.LDW ) THEN
                        DO 580 J = 1, N
                           W( J+LDA ) = W( J+LDW )
  580                   CONTINUE
                     END IF
                     LDW = LDA
*
*                    Time INVIT for right eigenvectors.
*
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  590                CONTINUE
*
                     CALL INVIT( LDA, N, H, W, W( LDA+1 ), LLWORK, N,
     $                           ITEMP, Z, IINFO, WORK( 2*N+1 ), WORK,
     $                           WORK( N+1 ) )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 12 ), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 610
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 590
*
*                    TIME = TIME / DOUBLE PRECISION( IC )
*                    OPS1 = OPS / DOUBLE PRECISION( IC )
*                    OPCNTS( IPAR, ITYPE, IN, 12 ) = OPS1
*                    TIMES( IPAR, ITYPE, IN, 12 ) = DMFLOP( OPS1, TIME,
*     $                  IINFO )
*
                     TIMES( IPAR, ITYPE, IN, 12 ) = TIME / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 12 ) = OPS / DBLE( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 12 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 12 )
                     TIMES( IPAR, ITYPE, IN, 12 ) = TIMES( LASTL, ITYPE,
     $                  IN, 12 )
                  END IF
  600          CONTINUE
            END IF
*
  610    CONTINUE
  620 CONTINUE
*
*-----------------------------------------------------------------------
*
*     Print a table of results for each timed routine.
*
      ISUB = 1
      IF( TIMSUB( ISUB ) ) THEN
         CALL DPRTBE( SUBNAM( ISUB ), MTYPES, DOTYPE, NSIZES, NN,
     $                INPARM( ISUB ), PNAMES, NPARMS, LDAS, NNB, NSHFTS,
     $                MAXBS, OPCNTS( 1, 1, 1, ISUB ), LDO1, LDO2,
     $                TIMES( 1, 1, 1, ISUB ), LDT1, LDT2, WORK, LLWORK,
     $                NOUT )
      END IF
*
      DO 630 IN = 1, NPARMS
         NNB( IN ) = 1
  630 CONTINUE
*
      DO 640 ISUB = 2, NSUBS
         IF( TIMSUB( ISUB ) ) THEN
            CALL DPRTBE( SUBNAM( ISUB ), MTYPES, DOTYPE, NSIZES, NN,
     $                   INPARM( ISUB ), PNAMES, NPARMS, LDAS, NNB,
     $                   NSHFTS, MAXBS, OPCNTS( 1, 1, 1, ISUB ), LDO1,
     $                   LDO2, TIMES( 1, 1, 1, ISUB ), LDT1, LDT2, WORK,
     $                   LLWORK, NOUT )
         END IF
  640 CONTINUE
*
      RETURN
*
*     End of DTIM21
*
 9997 FORMAT( ' DTIM21: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', ITYPE=', I6, ', IPAR=', I6, ', ISEED=(',
     $      3( I5, ',' ), I5, ')' )
*
      END
