      SUBROUTINE CTIM21( LINE, NSIZES, NN, NTYPES, DOTYPE, NPARMS, NNB,
     $                   NSHFTS, MAXBS, LDAS, TIMMIN, NOUT, ISEED, A,
     $                   ARE, AIM, H, HRE, HIM, Z, ZRE, ZIM, W, WRE,
     $                   WIM, WORK, WORKRE, WORKIM, LWORK, RWORK,
     $                   LLWORK, IWORK, TIMES, LDT1, LDT2, LDT3, OPCNTS,
     $                   LDO1, LDO2, LDO3, INFO )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            INFO, LDO1, LDO2, LDO3, LDT1, LDT2, LDT3,
     $                   LWORK, NOUT, NPARMS, NSIZES, NTYPES
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * ), LLWORK( * )
      INTEGER            ISEED( * ), IWORK( * ), LDAS( * ), MAXBS( * ),
     $                   NN( * ), NNB( * ), NSHFTS( * )
      REAL               AIM( * ), ARE( * ), HIM( * ), HRE( * ),
     $                   OPCNTS( LDO1, LDO2, LDO3, * ), RWORK( * ),
     $                   TIMES( LDT1, LDT2, LDT3, * ), WIM( * ),
     $                   WORKIM( * ), WORKRE( * ), WRE( * ), ZIM( * ),
     $                   ZRE( * )
      COMPLEX            A( * ), H( * ), W( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*     CTIM21 times the LAPACK routines for the COMPLEX non-symmetric
*     eigenvalue problem.
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
*          such as CGEHRD, indicating that only routine CGEHRD will
*          be timed, or it may contain a generic name, such as CHS.
*          In this case, the rest of the line is scanned for the
*          first 12 non-blank characters, corresponding to the twelve
*          combinations of subroutine and options:
*          LAPACK:
*          1: CGEHRD
*          2: CHSEQR(JOB='E')
*          3: CHSEQR(JOB='S')
*          4: CHSEQR(JOB='I')
*          5: CTREVC(JOB='L')
*          6: CTREVC(JOB='R')
*          7: CHSEIN(JOB='L')
*          8: CHSEIN(JOB='R')
*          EISPACK:
*           9: CORTH  (compare with CGEHRD)
*          10: COMQR  (compare w/ CHSEQR -- JOB='E')
*          11: COMQR2 (compare w/ CHSEQR(JOB='I') plus CTREVC(JOB='R'))
*          12: CINVIT (compare with CHSEIN)
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
*          X is unitary (for j=1--4) or has condition sqrt(ULP)
*          (for j=5--8), and T has random O(1) entries in the upper
*          triangle and:
*          (j=1,5) evenly spaced entries 1, ..., ULP with random
*                  arguments
*          (j=2,6) geometrically spaced entries 1, ..., ULP with random
*                  arguments
*          (j=3,7) "clustered" entries 1, ULP,..., ULP with random
*                  arguments
*          (j=4,8) eigenvalues randomly chosen from ( ULP, 1 ) with
*                  random arguments
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
*          processed by CLAHQR (EISPACK method), to be tested.
*
*  LDAS    (input) INTEGER array, dimension( NPARMS )
*          The values of LDA, the leading dimension of all matrices,
*          to be tested.
*
*  TIMMIN  (input) REAL
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
*          each call to CTIM21
*
*  A       (workspace) COMPLEX array,
*                      dimension( max(NN)*max(LDAS) )
*          (a) During the testing of CGEHRD, the original matrix to
*              be tested.
*          (b) Later, the Schur form of the original matrix.
*
*  ARE     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          (a) During the testing of CORTH, the real part of the
*          (b) Later, the Schur form of the original matrix.
*          May be equivalenced with first half of A in calling routine.
*
*  AIM     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          (a) During the testing of CORTH, the imaginary part of the
*              original matrix to be tested.
*          (b) Later, the Schur form of the original matrix.
*          May be equivalenced with second half of A in calling
*          routine.
*
*  H       (workspace) COMPLEX array,
*                      dimension( max(NN)*max(LDAS) )
*          The Hessenberg form of the original matrix.
*
*  HRE     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          The real part of the Hessenberg form of the original matrix.
*          May be equivalenced with first half of H in calling routine.
*          Used for testing EISPACK routines.
*
*  HIM     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          The imaginary part of the Hessenberg form of the original
*          matrix. May be equivalenced with second half of H in calling
*          routine. Used for testing EISPACK routines.
*
*  Z       (workspace) COMPLEX array,
*                      dimension( max(NN)*max(LDAS) )
*          Various output arrays: from CGEHRD and CHSEQR, the
*          unitary reduction matrices; from CTREVC and CHSEIN,
*          the eigenvector matrices.
*
*  ZRE     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          Various output arrays in testing EISPACK routines.
*          May be equivalenced with first half of Z in calling routine.
*
*  ZIM     (workspace) REAL array,
*                      dimension( max(NN)*max(LDAS) )
*          Various output arrays in testing EISPACK routines.
*          May be equivalenced with second half of Z in calling
*          routine.
*
*  W       (workspace) COMPLEX array, dimension( 2*max(LDAS) )
*          Holds computed eigenvalues.
*
*  WRE     (workspace) REAL array,
*                      dimension( 2*max(LDAS) )
*          Holds real parts of computed eigenvalues. Used for testing
*          EISPACK routines. May be equivalenced with first half of W
*          in calling routine.
*
*  WIM     (workspace) REAL array,
*                      dimension( 2*max(LDAS) )
*          Holds imaginary parts of computed eigenvalues. Used for
*          testing EISPACK routines. May be equivalenced with second
*          half of W in calling routine.
*
*  WORK    (workspace) COMPLEX array, dimension( LWORK )
*
*  WORKRE  (workspace) REAL array, dimension( LWORK )
*          May be equivalenced with first half of WORK in calling
*          routine.
*
*  WORKIM  (workspace) REAL array, dimension( LWORK )
*          May be equivalenced with second half of WORK in calling
*          routine.
*
*  LWORK   (input) INTEGER
*          Number of elements in WORK.  It must be at least:
*          (a)  max(NN)*( 3*max(NNB) + 2 )
*          (b)  max(NN)*( max(NNB+NSHFTS) + 1 )
*          (c)  max(NSHFTS)*( max(NSHFTS) + max(NN) )
*          (d)  max(MAXBS)*( max(MAXBS) + max(NN) )
*          (e)  max(NN)**2  +  max(NN)
*          (f)  4*max(NN)
*
*  RWORK   (workspace) REAL array, dimension
*                   ( max(max(NN),NSIZES*NTYPES*NPARMS) )
*          This should *not* be EQUIVALENCEd with any part of WORK.
*
*  LLWORK  (workspace) LOGICAL array, dimension( max( max(NN), NPARMS ))
*
*  IWORK   (workspace) INTEGER array, dimension( 2*max(NN) )
*          Workspace needed for parameters IFAILL and IFAILR in call
*          to CHSEIN.
*
*  TIMES   (output) REAL array,
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
*  OPCNTS  (output) REAL array,
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
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            RUNHQR, RUNHRD, RUNORT, RUNQRE, RUNQRS
      INTEGER            IC, ICONDS, IINFO, IMODE, IN, IPAR, ISUB,
     $                   ITEMP, ITYPE, J, J1, J2, J3, J4, JC, JR, LASTL,
     $                   LASTNL, LDA, LDAMIN, LDH, LDT, MAXB, MBMAX,
     $                   MTYPES, N, NB, NBMAX, NMAX, NSBMAX,
     $                   NSHIFT, NSMAX
      REAL               CONDS, RTULP, RTULPI, S1, S2, TIME, ULP,
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
      REAL               SLAMCH, SECOND, SOPLA
      EXTERNAL           SLAMCH, SECOND, SOPLA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMIN, CGEHRD, CHSEIN, CHSEQR, CINVIT, CLACPY,
     $                   CLATME, COMQR, COMQR2, CORTH, CTREVC, SLACPY,
     $                   SLASET, SPRTBE, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, MAX, MIN, REAL, SQRT
*     ..
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      REAL               ITCNT, OPS
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'CGEHRD', 'CHSEQR(E)', 'CHSEQR(S)',
     $                   'CHSEQR(V)', 'CTREVC(L)', 'CTREVC(R)',
     $                   'CHSEIN(L)', 'CHSEIN(R)', 'CORTH', 'COMQR',
     $                   'COMQR2', 'CINVIT' /
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
*
*     Extract the timing request from the input line.
*
      CALL ATIMIN( 'CHS', LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
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
 9999    FORMAT( 1X, A, ' timing run not attempted -- N < LDA', / )
         RETURN
      END IF
*
*     Check LWORK
*
      IF( LWORK.LT.MAX( NMAX*MAX( 4, 3*NBMAX+2, NSBMAX+1 ),
     $    NSMAX*( NSMAX+NMAX ), MBMAX*( MBMAX+NMAX ),
     $    ( NMAX+1 )*NMAX ) ) THEN
         INFO = -29
         WRITE( NOUT, FMT = 9998 )LINE( 1: 6 )
 9998    FORMAT( 1X, A, ' timing run not attempted -- LWORK too small.',
     $         / )
         RETURN
      END IF
*
*     Check to see whether CGEHRD or CHSEQR must be run.
*
*     RUNQRE -- if CHSEQR must be run to get eigenvalues.
*     RUNQRS -- if CHSEQR must be run to get Schur form.
*     RUNHRD -- if CGEHRD must be run.
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
*     Check to see whether CORTH or COMQR must be run.
*
*     RUNHQR -- if COMQR must be run to get eigenvalues.
*     RUNORT -- if CORTH must be run.
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
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
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
      DO 550 IN = 1, NSIZES
*
         N = NN( IN )
*
*        Do for each .TRUE. value in DOTYPE:
*
         MTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.EQ.MAXTYP+1 .AND. NSIZES.EQ.1 )
     $      MTYPES = NTYPES
         DO 540 ITYPE = 1, MTYPES
            IF( .NOT.DOTYPE( ITYPE ) )
     $         GO TO 540
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
               CALL CLATME( N, 'S', ISEED, WORK, IMODE, ULPINV, CONE,
     $                      ADUMMA, 'T', 'T', 'T', RWORK, 4, CONDS, N,
     $                      N, ONE, A, N, WORK( N+1 ), IINFO )
            END IF
*
*           Time CGEHRD for each pair NNB(j), LDAS(j)
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
*                    Time CGEHRD
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
   90                CONTINUE
                     CALL CLACPY( 'Full', N, N, A, N, H, LDA )
*
                     CALL CGEHRD( N, 1, N, H, LDA, WORK,
     $                            WORK( N+1 ), LWORK-N, IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 1 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 540
                     END IF
*
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 90
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 100 J = 1, IC
                        CALL CLACPY( 'Full', N, N, A, N, Z, LDA )
  100                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 1 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 1 ) = SOPLA( 'CGEHRD', N,
     $                       1, N, 0, NB )
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
                  CALL CLACPY( 'Full', N, N, A, N, H, N )
*
                  CALL CGEHRD( N, 1, N, H, N, WORK, WORK( N+1 ),
     $                         LWORK-N, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 1 ), IINFO, N,
     $                  ITYPE, 0, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 540
                  END IF
                  LDH = N
               END IF
            END IF
*
*           Time CHSEQR with JOB='E' for each 4-tuple
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
*                 Time CHSEQR with JOB='E'
*
                  IC = 0
                  OPS = ZERO
                  S1 = SECOND( )
  120             CONTINUE
                  CALL CLACPY( 'Full', N, N, H, LDH, A, LDA )
*
                  CALL CHSEQR( 'E', 'N', N, 1, N, A, LDA, W, Z, LDA,
     $                         WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 2 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 540
                  END IF
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 120
*
*                 Subtract the time used in CLACPY.
*
                  S1 = SECOND( )
                  DO 130 J = 1, IC
                     CALL CLACPY( 'Full', N, N, H, LDH, Z, LDA )
  130             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 2 ) = MAX( TIME-UNTIME,
     $               ZERO ) / REAL( IC )
                  OPCNTS( IPAR, ITYPE, IN, 2 ) = OPS / REAL( IC )
  140          CONTINUE
               LDT = 0
            ELSE
               IF( RUNQRE ) THEN
                  CALL CLACPY( 'Full', N, N, H, LDH, A, N )
*
                  CALL CHSEQR( 'E', 'N', N, 1, N, A, N, W, Z, N,
     $                         WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 2 ), IINFO, N,
     $                  ITYPE, 0, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 540
                  END IF
                  LDT = 0
               END IF
            END IF
*
*           Time CHSEQR with JOB='S' for each 4-tuple
*           NNB(j), NSHFTS(j), MAXBS(j), LDAS(j)
*
            IF( TIMSUB( 3 ) ) THEN
               DO 170 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = 1
                  NSHIFT = NSHFTS( IPAR )
                  MAXB = MAXBS( IPAR )
                  CALL XLAENV( 4, NSHIFT )
                  CALL XLAENV( 8, MAXB )
*
*                 Time CHSEQR with JOB='S'
*
                  IC = 0
                  OPS = ZERO
                  S1 = SECOND( )
  150             CONTINUE
                  CALL CLACPY( 'Full', N, N, H, LDH, A, LDA )
*
                  CALL CHSEQR( 'S', 'N', N, 1, N, A, LDA, W, Z, LDA,
     $                         WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 3 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 540
                  END IF
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 150
*
*                 Subtract the time used in CLACPY.
*
                  S1 = SECOND( )
                  DO 160 J = 1, IC
                     CALL CLACPY( 'Full', N, N, H, LDH, Z, LDA )
  160             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 3 ) = MAX( TIME-UNTIME,
     $               ZERO ) / REAL( IC )
                  OPCNTS( IPAR, ITYPE, IN, 3 ) = OPS / REAL( IC )
  170          CONTINUE
               LDT = LDA
            ELSE
               IF( RUNQRS ) THEN
                  CALL CLACPY( 'Full', N, N, H, LDH, A, N )
*
                  CALL CHSEQR( 'S', 'N', N, 1, N, A, N, W, Z, N,
     $                         WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 3 ), IINFO, N,
     $                  ITYPE, 0, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 540
                  END IF
                  LDT = N
               END IF
            END IF
*
*           Time CHSEQR with JOB='I' for each 4-tuple
*           NNB(j), NSHFTS(j), MAXBS(j), LDAS(j)
*
            IF( TIMSUB( 4 ) ) THEN
               DO 200 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = 1
                  NSHIFT = NSHFTS( IPAR )
                  MAXB = MAXBS( IPAR )
                  CALL XLAENV( 4, NSHIFT )
                  CALL XLAENV( 8, MAXB )
*
*                 Time CHSEQR with JOB='I'
*
                  IC = 0
                  OPS = ZERO
                  S1 = SECOND( )
  180             CONTINUE
                  CALL CLACPY( 'Full', N, N, H, LDH, A, LDA )
*
                  CALL CHSEQR( 'S', 'I', N, 1, N, A, LDA, W, Z, LDA,
     $                         WORK, LWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 4 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 540
                  END IF
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 180
*
*                 Subtract the time used in CLACPY.
*
                  S1 = SECOND( )
                  DO 190 J = 1, IC
                     CALL CLACPY( 'Full', N, N, H, LDH, Z, LDA )
  190             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 4 ) = MAX( TIME-UNTIME,
     $               ZERO ) / REAL( IC )
                  OPCNTS( IPAR, ITYPE, IN, 4 ) = OPS / REAL( IC )
  200          CONTINUE
               LDT = LDA
            END IF
*
*           Time CTREVC and CHSEIN with various values of LDA
*
*           Select All Eigenvectors
*
            DO 210 J = 1, N
               LLWORK( J ) = .TRUE.
  210       CONTINUE
*
            DO 350 IPAR = 1, NPARMS
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
*              Time CTREVC
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
*                 Time CTREVC for Left Eigenvectors
*
                  IF( TIMSUB( 5 ) ) THEN
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  270                CONTINUE
*
                     CALL CTREVC( 'L', 'A', LLWORK, N, A, LDA, Z, LDA,
     $                            Z, LDA, N, ITEMP, WORK, RWORK, IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 5 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 540
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 270
*
                     TIMES( IPAR, ITYPE, IN, 5 ) = TIME / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 5 ) = OPS / REAL( IC )
                  END IF
*
*                 Time CTREVC for Right Eigenvectors
*
                  IF( TIMSUB( 6 ) ) THEN
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  280                CONTINUE
                     CALL CTREVC( 'R', 'A', LLWORK, N, A, LDA, Z, LDA,
     $                            Z, LDA, N, ITEMP, WORK, RWORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 6 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 540
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 280
*
                     TIMES( IPAR, ITYPE, IN, 6 ) = TIME / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 6 ) = OPS / REAL( IC )
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
*              Time CHSEIN
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
  300                CONTINUE
                  ELSE IF( LDA.LT.LDH ) THEN
                     DO 320 JC = 1, N
                        DO 310 JR = 1, N
                           H( JR+( JC-1 )*LDA ) = H( JR+( JC-1 )*LDH )
  310                   CONTINUE
  320                CONTINUE
                  END IF
                  LDH = LDA
*
*                 Time CHSEIN for Left Eigenvectors
*
                  IF( TIMSUB( 7 ) ) THEN
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  330                CONTINUE
*
                     CALL CHSEIN( 'L', 'Q', 'N', LLWORK, N, H, LDA, W,
     $                            Z, LDA, Z, LDA, N, ITEMP, WORK, RWORK,
     $                            IWORK, IWORK( N+1 ), IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 7 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 540
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 330
*
                     TIMES( IPAR, ITYPE, IN, 7 ) = TIME / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 7 ) = OPS / REAL( IC )
                  END IF
*
*                 Time CHSEIN for Right Eigenvectors
*
                  IF( TIMSUB( 8 ) ) THEN
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  340                CONTINUE
*
                     CALL CHSEIN( 'R', 'Q', 'N', LLWORK, N, H, LDA, W,
     $                            Z, LDA, Z, LDA, N, ITEMP, WORK, RWORK,
     $                            IWORK, IWORK( N+1 ), IINFO )
*
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 8 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 540
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 340
*
                     TIMES( IPAR, ITYPE, IN, 8 ) = TIME / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 8 ) = OPS / REAL( IC )
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
  350       CONTINUE
*
*-----------------------------------------------------------------------
*
*           Time the EISPACK Routines
*
*           Restore random number seed
*
            DO 360 J = 1, 4
               ISEED( J ) = IOLDSD( J )
  360       CONTINUE
*
*           Re-generate A, copy to ARE and AIM
*
            IF( ITYPE.LE.MAXTYP ) THEN
               IMODE = KMODE( ITYPE )
               IF( ICONDS.EQ.1 ) THEN
                  CONDS = ONE
               ELSE
                  CONDS = RTULPI
               END IF
               CALL CLATME( N, 'S', ISEED, WORK, IMODE, ULPINV, CONE,
     $                      ADUMMA, 'T', 'T', 'T', RWORK, 4, CONDS, N,
     $                      N, ONE, H, N, WORK( N+1 ), IINFO )
               DO 370 J = 1, N*N
                  ARE( J ) = REAL( H( J ) )
                  AIM( J ) = AIMAG( H( J ) )
  370          CONTINUE
            END IF
*
*           Time CORTH for each LDAS(j)
*
            IF( TIMSUB( 9 ) ) THEN
               DO 410 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 380 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  380             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CORTH
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  390                CONTINUE
                     CALL SLACPY( 'Full', N, N, ARE, N, HRE, LDA )
                     CALL SLACPY( 'Full', N, N, AIM, N, HIM, LDA )
                     CALL CORTH( LDA, N, 1, N, HRE, HIM, WORKRE,
     $                           WORKIM )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 390
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 400 J = 1, IC
                        CALL SLACPY( 'Full', N, N, ARE, N, ZRE, LDA )
                        CALL SLACPY( 'Full', N, N, AIM, N, ZIM, LDA )
  400                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 9 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 9 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 9 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 9 )
                     TIMES( IPAR, ITYPE, IN, 9 ) = TIMES( LASTL, ITYPE,
     $                  IN, 9 )
                  END IF
                  LDH = LDA
  410          CONTINUE
            ELSE
               IF( RUNORT ) THEN
                  CALL SLACPY( 'Full', N, N, ARE, N, HRE, N )
                  CALL SLACPY( 'Full', N, N, AIM, N, HIM, N )
                  CALL CORTH( N, N, 1, N, HRE, HIM, WORKRE, WORKIM )
                  LDH = N
               END IF
            END IF
*
*           Time COMQR for each LDAS(j)
*
            IF( TIMSUB( 10 ) ) THEN
               DO 450 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 420 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  420             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time COMQR
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  430                CONTINUE
                     CALL SLACPY( 'Full', N, N, HRE, LDH, ARE, LDA )
                     CALL SLACPY( 'Full', N, N, HIM, LDH, AIM, LDA )
                     CALL COMQR( LDA, N, 1, N, ARE, AIM, WRE, WIM,
     $                           IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 10 ), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 540
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 430
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 440 J = 1, IC
                        CALL SLACPY( 'Full', N, N, HRE, LDH, ZRE, LDA )
                        CALL SLACPY( 'Full', N, N, HIM, LDH, ZIM, LDA )
  440                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 10 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 10 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 10 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 10 )
                     TIMES( IPAR, ITYPE, IN, 10 ) = TIMES( LASTL, ITYPE,
     $                  IN, 10 )
                  END IF
  450          CONTINUE
            ELSE
               IF( RUNHQR ) THEN
                  CALL SLACPY( 'Full', N, N, HRE, LDH, ARE, N )
                  CALL SLACPY( 'Full', N, N, HIM, LDH, AIM, N )
                  CALL COMQR( N, N, 1, N, ARE, AIM, WRE, WIM, IINFO )
               END IF
            END IF
*
*           Time COMQR2 for each LDAS(j)
*
            IF( TIMSUB( 11 ) ) THEN
               DO 490 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 460 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  460             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time COMQR2
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  470                CONTINUE
                     CALL SLACPY( 'Full', N, N, HRE, LDH, ARE, LDA )
                     CALL SLACPY( 'Full', N, N, HIM, LDH, AIM, LDA )
                     CALL SLASET( 'Full', N, N, ZERO, ONE, ZRE, LDA )
                     CALL SLASET( 'Full', N, N, ZERO, ZERO, ZIM, LDA )
                     CALL SLASET( 'Full', 1, N, ZERO, ZERO, WORKRE, 1 )
                     CALL SLASET( 'Full', 1, N, ZERO, ZERO, WORKIM, 1 )
                     CALL COMQR2( LDA, N, 1, N, WORKRE, WORKIM, ARE,
     $                            AIM, WRE, WIM, ZRE, ZIM, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 11 ), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 540
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 470
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 480 J = 1, IC
                        CALL SLACPY( 'Full', N, N, HRE, LDH, ZRE, LDA )
                        CALL SLACPY( 'Full', N, N, HIM, LDH, ZIM, LDA )
  480                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 11 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 11 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 11 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 11 )
                     TIMES( IPAR, ITYPE, IN, 11 ) = TIMES( LASTL, ITYPE,
     $                  IN, 11 )
                  END IF
  490          CONTINUE
            END IF
*
*           Time CINVIT for each LDAS(j)
*
*           Select All Eigenvectors
*
            DO 500 J = 1, N
               LLWORK( J ) = .TRUE.
  500       CONTINUE
*
            IF( TIMSUB( 12 ) ) THEN
               DO 530 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 510 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  510             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Copy H if necessary to get right LDA.
*
                     IF( LDA.NE.LDH ) THEN
                        CALL SLACPY( 'Full', N, N, HRE, LDH, ZRE, LDA )
                        CALL SLACPY( 'Full', N, N, HIM, LDH, ZIM, LDA )
                        CALL SLACPY( 'Full', N, N, ZRE, LDA, HRE, LDA )
                        CALL SLACPY( 'Full', N, N, ZIM, LDA, HIM, LDA )
                     END IF
                     LDH = LDA
*
*                    Time CINVIT for right eigenvectors.
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  520                CONTINUE
                     CALL CINVIT( LDA, N, HRE, HIM, WRE, WIM, LLWORK, N,
     $                            ITEMP, ZRE, ZIM, IINFO, WORKRE( N+1 ),
     $                            WORKIM( N+1 ), WORKRE, WORKIM )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 12 ), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 540
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 520
*
*                    TIME = TIME / REAL( IC )
*                    OPS1 = OPS / REAL( IC )
*                    OPCNTS( IPAR, ITYPE, IN, 12 ) = OPS1
*                    TIMES( IPAR, ITYPE, IN, 12 ) = SMFLOP( OPS1, TIME,
*     $                  IINFO )
                     TIMES( IPAR, ITYPE, IN, 12 ) = TIME / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 12 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 12 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 12 )
                     TIMES( IPAR, ITYPE, IN, 12 ) = TIMES( LASTL, ITYPE,
     $                  IN, 12 )
                  END IF
  530          CONTINUE
            END IF
*
  540    CONTINUE
  550 CONTINUE
*
*-----------------------------------------------------------------------
*
*     Print a table of results for each timed routine.
*
      ISUB = 1
      IF( TIMSUB( ISUB ) ) THEN
         CALL SPRTBE( SUBNAM( ISUB ), MTYPES, DOTYPE, NSIZES, NN,
     $                 INPARM( ISUB ), PNAMES, NPARMS, LDAS, NNB,
     $                 NSHFTS, MAXBS, OPCNTS( 1, 1, 1, ISUB ), LDO1,
     $                 LDO2, TIMES( 1, 1, 1, ISUB ), LDT1, LDT2,
     $                 RWORK, LLWORK, NOUT )
      END IF
*
      DO 555 IN = 1, NPARMS
         NNB( IN ) = 1
  555 CONTINUE
*
      DO 560 ISUB = 2, NSUBS
         IF( TIMSUB( ISUB ) ) THEN
            CALL SPRTBE( SUBNAM( ISUB ), MTYPES, DOTYPE, NSIZES, NN,
     $                   INPARM( ISUB ), PNAMES, NPARMS, LDAS, NNB,
     $                   NSHFTS, MAXBS, OPCNTS( 1, 1, 1, ISUB ), LDO1,
     $                   LDO2, TIMES( 1, 1, 1, ISUB ), LDT1, LDT2,
     $                   RWORK, LLWORK, NOUT )
         END IF
  560 CONTINUE
*
 9997 FORMAT( ' CTIM21: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', ITYPE=', I6, ', IPAR=', I6, ', ISEED=(',
     $      3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of CTIM21
*
      END
