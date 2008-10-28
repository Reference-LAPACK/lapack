      SUBROUTINE CTIM51( LINE, NSIZES, NN, NTYPES, DOTYPE, NPARMS, NNB,
     $                   NSHFTS, NEISPS, MINNBS, MINBKS, LDAS, TIMMIN,
     $                   NOUT, ISEED, A, AR, AI, B, BR, BI, H, HR, HI,
     $                   T, TR, TI, Q, QR, QI, Z, ZR, ZI, W, WR, WORK,
     $                   LWORK, RWORK, LLWORK, TIMES, LDT1, LDT2, LDT3,
     $                   OPCNTS, LDO1, LDO2, LDO3, INFO )
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
      INTEGER            ISEED( * ), LDAS( * ), MINBKS( * ),
     $                   MINNBS( * ), NEISPS( * ), NN( * ), NNB( * ),
     $                   NSHFTS( * )
      REAL               AI( * ), AR( * ), BI( * ), BR( * ), HI( * ),
     $                   HR( * ), OPCNTS( LDO1, LDO2, LDO3, * ),
     $                   QI( * ), QR( * ), RWORK( * ), TI( * ),
     $                   TIMES( LDT1, LDT2, LDT3, * ), TR( * ), WR( * ),
     $                   ZI( * ), ZR( * )
      COMPLEX            A( * ), B( * ), H( * ), Q( * ), T( * ), W( * ),
     $                   WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  CTIM51 times the LAPACK routines for the COMPLEX non-symmetric
*  generalized eigenvalue problem   A x = w B x.
*
*  For each N value in NN(1:NSIZES) and .TRUE. value in
*  DOTYPE(1:NTYPES), a pair of matrices will be generated and used to
*  test the selected routines.  Thus, NSIZES*(number of .TRUE. values
*  in DOTYPE) matrices will be generated.
*
*  Arguments
*  =========
*
*  LINE    (input) CHARACTER*80
*          The input line which requested this routine.  This line may
*          contain a subroutine name, such as CGGHRD, indicating that
*          only routine CGGHRD will be timed, or it may contain a
*          generic name, such as CHG.  In this case, the rest of the
*          line is scanned for the first 18 non-blank characters,
*          corresponding to the eighteen combinations of subroutine and
*          options:
*          LAPACK:                                     Table Heading:
*           1: CGGHRD(no Q, no Z) (+CGEQRF, etc.)      'CGGHRD(N)'
*           2: CGGHRD(Q only)     (+CGEQRF, etc.)      'CGGHRD(Q)'
*           3: CGGHRD(Z only)     (+CGEQRF, etc.)      'CGGHRD(Z)'
*           4: CGGHRD(Q and Z)    (+CGEQRF, etc.)      'CGGHRD(Q,Z)'
*           5: CHGEQZ(Eigenvalues only)                'CHGEQZ(E)'
*           6: CHGEQZ(Schur form only)                 'CHGEQZ(S)'
*           7: CHGEQZ(Schur form and Q)                'CHGEQZ(Q)'
*           8: CHGEQZ(Schur form and Z)                'CHGEQZ(Z)'
*           9: CHGEQZ(Schur form, Q and Z)             'CHGEQZ(Q,Z)'
*          10: CTGEVC(SIDE='L', HOWMNY='A')            'CTGEVC(L,A)'
*          11: CTGEVC(SIDE='L', HOWMNY='B')            'CTGEVC(L,B)'
*          12: CTGEVC(SIDE='R', HOWMNY='A')            'CTGEVC(R,A)'
*          13: CTGEVC(SIDE='R', HOWMNY='B')            'CTGEVC(R,B)'
*          EISPACK:                       Compare w/:  Table Heading:
*          14: CQZHES w/ matz=.false.            1      'CQZHES(F)'
*          15: CQZHES w/ matz=.true.             3      'CQZHES(T)'
*          16: CQZVAL w/ matz=.false.            5      'CQZVAL(F)'
*          17: CQZVAL w/ matz=.true.             8      'CQZVAL(T)'
*          18: CQZVEC                           13      'CQZVEC'
*          If a character is 'T' or 't', the corresponding routine in
*          this path is timed.  If the entire line is blank, all the
*          routines in the path are timed.
*
*          Note that since QZHES does more than SGGHRD, the
*          "SGGHRD" timing also includes the time for the calls
*          to SGEQRF, SORMQR, and (if Q is computed) SORGQR
*          which are necessary to get the same functionality
*          as QZHES.
*
*  NSIZES  (input) INTEGER
*          The number of values of N contained in the vector NN.
*
*  NN      (input) INTEGER array, dimension (NSIZES)
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
*          If DOTYPE(j) is .TRUE., then a pair of matrices (A,B) of
*          type j will be generated.  A and B have the form  U T1 V
*          and  U T2 V , resp., where U and V are orthogonal, T1 and
*          T2 are upper triangular.  T2 has random O(1) entries in the
*          strict upper triangle and ( 0, 1, 0, 1, 1, ..., 1, 0 ) on
*          the diagonal, while T1 has random O(1) entries in the strict
*          upper triangle, its diagonal will have the values:
*          (j=1)   0, 0, 1, 1, ULP,..., ULP, 0.
*          (j=2)   0, 0, 1, 1, 1-d, 1-2*d, ..., 1-(N-5)*d=ULP, 0.
*
*                                  2        N-5
*          (j=3)   0, 0, 1, 1, a, a , ..., a   =ULP, 0.
*          (j=4)   0, 0, 1, r1, r2, ..., r(N-4), 0, where r1, etc.
*                  are random numbers in (ULP,1).
*
*  NPARMS  (input) INTEGER
*          The number of values in each of the arrays NNB, NSHFTS,
*          NEISPS, and LDAS.  For each pair of matrices A,B generated
*          according to NN and DOTYPE, tests will be run with
*          (NB,NSHIFT,NEISP,LDA)= (NNB(1), NSHFTS(1), NEISPS(1),
*          LDAS(1)),..., (NNB(NPARMS), NSHFTS(NPARMS), NEISPS(NPARMS),
*          LDAS(NPARMS))
*
*  NNB     (input) INTEGER array, dimension (NPARMS)
*          The values of the blocksize ("NB") to be tested.  They must
*          be at least 1.  Currently, this is only used by CGEQRF,
*          etc., in the timing of CGGHRD.
*
*  NSHFTS  (input) INTEGER array, dimension (NPARMS)
*          The values of the number of shifts ("NSHIFT") to be tested.
*          (Currently not used.)
*
*  NEISPS  (input) INTEGER array, dimension (NPARMS)
*          The values of "NEISP", the size of largest submatrix to be
*          processed by CLAEQZ (EISPACK method), to be tested.
*          (Currently not used.)
*
*  MINNBS  (input) INTEGER array, dimension (NPARMS)
*          The values of "MINNB", the minimum size of a product of
*          transformations which may be applied as a blocked
*          transformation, to be tested.  (Currently not used.)
*
*  MINBKS  (input) INTEGER array, dimension (NPARMS)
*          The values of "MINBK", the minimum number of rows/columns
*          to be updated with a blocked transformation, to be tested.
*          (Currently not used.)
*
*  LDAS    (input) INTEGER array, dimension (NPARMS)
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
*  ISEED   (input/output) INTEGER array, dimension (4)
*          The random seed used by the random number generator, used
*          by the test matrix generator.  It is used and updated on
*          each call to CTIM51.
*
*  A       (workspace) COMPLEX array, dimension (max(NN)*max(LDAS))
*          (a) During the testing of CGGHRD, "A", the original
*              left-hand-side matrix to be tested.
*          (b) Later, "S", the Schur form of the original "A" matrix.
*
*  AR, AI  (workspace) REAL arrays, dimension
*                      (max(NN)*max(LDAS))
*          The real and imaginary parts of A, stored separately for
*          the benefit of CQZHES, CQZVAL, and CQZVEC.  These may be
*          equivalenced with A by the calling routine.
*
*  B       (workspace) COMPLEX array, dimension (max(NN)*max(LDAS))
*          (a) During the testing of CGGHRD, "B", the original
*              right-hand-side matrix to be tested.
*          (b) Later, "P", the Schur form of the original "B" matrix.
*
*  BR, BI  (workspace) REAL arrays, dimension
*                      (max(NN)*max(LDAS))
*          The real and imaginary parts of B, stored separately for
*          the benefit of CQZHES, CQZVAL, and CQZVEC.  These may be
*          equivalenced with B by the calling routine.
*
*  H       (workspace) COMPLEX array, dimension (max(NN)*max(LDAS))
*          (a) During the testing of CGGHRD and CHGEQZ, "H", the
*              Hessenberg form of the original "A" matrix.
*          (b) During the testing of CTGEVC, "L", the matrix of left
*              eigenvectors.
*
*  HR, HI  (workspace) REAL arrays, dimension
*                      (max(NN)*max(LDAS))
*          The real and imaginary parts of H, stored separately for
*          the benefit of CQZHES, CQZVAL, and CQZVEC.  These may be
*          equivalenced with H by the calling routine.
*
*  T       (workspace) COMPLEX array, dimension (max(NN)*max(LDAS))
*          (a) During the testing of CGGHRD and CHGEQZ, "T", the
*              triangular form of the original "B" matrix.
*          (b) During the testing of CTGEVC, "R", the matrix of right
*              eigenvectors.
*
*  TR, TI  (workspace) REAL arrays, dimension
*                      (max(NN)*max(LDAS))
*          The real and imaginary parts of T, stored separately for
*          the benefit of CQZHES, CQZVAL, and CQZVEC.  These may be
*          equivalenced with T by the calling routine.
*
*  Q       (workspace) COMPLEX array, dimension (max(NN)*max(LDAS))
*          The orthogonal matrix on the left generated by CGGHRD.  If
*          CHGEQZ computes only Q or Z, then that matrix is stored here.
*          If both Q and Z are computed, the Q matrix goes here.
*
*  QR, QI  (workspace) REAL arrays, dimension
*                      (max(NN)*max(LDAS))
*          The real and imaginary parts of Q, stored separately for
*          the benefit of CQZVAL.  These may be equivalenced with Q by
*          the calling routine.
*
*  Z       (workspace) COMPLEX array, dimension (max(NN)*max(LDAS))
*          The orthogonal matrix on the right generated by CGGHRD.
*          If CHGEQZ computes both Q and Z, the Z matrix is stored here.
*          Also used as scratch space for timing the CLACPY calls.
*
*  ZR, ZI  (workspace) REAL arrays, dimension
*                      (max(NN)*max(LDAS))
*          The real and imaginary parts of Z, stored separately for
*          the benefit of CQZHES, CQZVAL, and CQZVEC.  These may be
*          equivalenced with Z by the calling routine.
*
*  W       (workspace) COMPLEX array, dimension (2*max(LDAS))
*          Treated as an LDA x 2 matrix whose 1st column holds
*          ALPHA, the diagonal entries of "S", and whose 2nd column
*          holds BETA, the diagonal entries of "P".
*
*  WR      (workspace) REAL array, dimension (3*max(LDAS))
*          Treated as an LDA x 3 matrix whose 1st and 2nd columns hold
*          the real and imaginary parts of ALPHA (see the description
*          of W), and whose 3rd column holds BETA (real part only.)
*          This may be equivalenced to W by the calling routine.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER [the following formulae are certainly wrong]
*          Number of elements in WORK.  LWORK >= 4*max(NN).
*
*  RWORK   (workspace) REAL array, dimension
*                      (max( 2*max(NN), NSIZES*NTYPES*NPARMS ))
*
*  LLWORK  (workspace) LOGICAL array, dimension (max( max(NN), NPARMS ))
*
*  TIMES   (output) REAL array, dimension
*                   (LDT1,LDT2,LDT3,NSUBS)
*          TIMES(i,j,k,l) will be set to the run time (in seconds) for
*          subroutine l, with N=NN(k), matrix type j, and LDA=LDAS(i),
*          NEISP=NEISPS(i), NBLOCK=NNB(i), NSHIFT=NSHFTS(i),
*          MINNB=MINNBS(i), and MINBLK=MINBKS(i).
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
*  OPCNTS  (output) REAL array, dimension
*                   (LDO1,LDO2,LDO3,NSUBS)
*          OPCNTS(i,j,k,l) will be set to the number of floating-point
*          operations executed by subroutine l, with N=NN(k), matrix
*          type j, and LDA=LDAS(i), NEISP=NEISPS(i), NBLOCK=NNB(i),
*          NSHIFT=NSHFTS(i), MINNB=MINNBS(i), and MINBLK=MINBKS(i).
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
      PARAMETER          ( MAXTYP = 4, NSUBS = 18 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            RUNEQ, RUNES, RUNHES, RUNHRD, RUNQZ
      INTEGER            IC, IINFO, IN, IPAR, ISUB, ITEMP, ITYPE, J, J1,
     $                   J2, J3, J4, JC, JR, LASTL, LDA, LDAMIN, LDH,
     $                   LDQ, LDS, LDW, MINBLK, MINNB, MTYPES, N, N1,
     $                   NB, NBSMAX, NEISP, NMAX, NSHIFT
      REAL               S1, S2, TIME, ULP, UNTIME
      COMPLEX            CTEMP
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER(32)      PNAMES( 6 )
      CHARACTER*11       SUBNAM( NSUBS )
      INTEGER            INPARM( NSUBS ), IOLDSD( 4 ), KATYPE( MAXTYP )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      REAL               SECOND, SLAMCH, SOPLA
      COMPLEX            CLARND
      EXTERNAL           SECOND, SLAMCH, SOPLA, CLARND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMIN, CHGEQZ, CLACPY, CLAQZH, CLARFG, CLATM4,
     $                   CQZHES, CQZVAL, CQZVEC, CTGEVC, CUNM2R, SLACPY,
     $                   SLASET, SPRTBG, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CONJG, MAX, MIN, REAL, SIGN
*     ..
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      REAL               ITCNT, OPS
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'CGGHRD(N)', 'CGGHRD(Q)', 'CGGHRD(Z)',
     $                   'CGGHRD(Q,Z)', 'CHGEQZ(E)', 'CHGEQZ(S)',
     $                   'CHGEQZ(Q)', 'CHGEQZ(Z)', 'CHGEQZ(Q,Z)',
     $                   'CTGEVC(L,A)', 'CTGEVC(L,B)', 'CTGEVC(R,A)',
     $                   'CTGEVC(R,B)', 'CQZHES(F)', 'CQZHES(T)',
     $                   'CQZVAL(F)', 'CQZVAL(T)', 'CQZVEC' /
      DATA               INPARM / 4*2, 5*1, 4*1, 5*1 /
      DATA               PNAMES / '   LDA', '    NB', '    NS',
     $                   ' NEISP', ' MINNB', 'MINBLK' /
      DATA               KATYPE / 5, 8, 7, 9 /
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
      CALL ATIMIN( 'CHG', LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
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
      NBSMAX = 0
      DO 20 J1 = 1, NPARMS
         LDAMIN = MIN( LDAMIN, LDAS( J1 ) )
         NBSMAX = MAX( NBSMAX, NNB( J1 )+NSHFTS( J1 ) )
   20 CONTINUE
*
*     Check that N <= LDA for the input values.
*
      IF( NMAX.GT.LDAMIN ) THEN
         INFO = -12
         WRITE( NOUT, FMT = 9999 )LINE( 1: 6 )
 9999    FORMAT( 1X, A, ' timing run not attempted -- N > LDA', / )
         RETURN
      END IF
*
*     Check LWORK
*
      IF( LWORK.LT.4*NMAX )
     $     THEN
         INFO = -24
         WRITE( NOUT, FMT = 9998 )LINE( 1: 6 )
 9998    FORMAT( 1X, A, ' timing run not attempted -- LWORK too small.',
     $         / )
         RETURN
      END IF
*
*     Check to see whether CGGHRD or CHGEQZ must be run.
*        RUNHRD -- if CGGHRD must be run.
*        RUNES  -- if CHGEQZ must be run to get Schur form.
*        RUNEQ  -- if CHGEQZ must be run to get Schur form and Q.
*
      RUNHRD = .FALSE.
      RUNES = .FALSE.
      RUNEQ = .FALSE.
*
      IF( TIMSUB( 10 ) .OR. TIMSUB( 12 ) )
     $   RUNES = .TRUE.
      IF( TIMSUB( 11 ) .OR. TIMSUB( 13 ) )
     $   RUNEQ = .TRUE.
      IF( TIMSUB( 5 ) .OR. TIMSUB( 6 ) .OR. TIMSUB( 7 ) .OR.
     $    TIMSUB( 8 ) .OR. TIMSUB( 9 ) .OR. RUNES .OR. RUNEQ )
     $    RUNHRD = .TRUE.
*
      IF( TIMSUB( 6 ) .OR. TIMSUB( 7 ) .OR. TIMSUB( 8 ) .OR.
     $    TIMSUB( 9 ) .OR. RUNEQ )RUNES = .FALSE.
      IF( TIMSUB( 7 ) .OR. TIMSUB( 8 ) .OR. TIMSUB( 9 ) )
     $   RUNEQ = .FALSE.
      IF( TIMSUB( 1 ) .OR. TIMSUB( 2 ) .OR. TIMSUB( 3 ) .OR.
     $    TIMSUB( 4 ) )RUNHRD = .FALSE.
*
*     Check to see whether CQZHES or CQZVAL must be run.
*
*     RUNHES -- if CQZHES must be run.
*     RUNQZ  -- if CQZVAL must be run (w/ MATZ=.TRUE.).
*
      RUNHES = .FALSE.
      RUNQZ = .FALSE.
*
      IF( TIMSUB( 18 ) )
     $   RUNQZ = .TRUE.
      IF( TIMSUB( 16 ) .OR. TIMSUB( 17 ) .OR. RUNQZ )
     $   RUNHES = .TRUE.
      IF( TIMSUB( 17 ) )
     $   RUNQZ = .FALSE.
      IF( TIMSUB( 14 ) .OR. TIMSUB( 15 ) )
     $   RUNHES = .FALSE.
*
*     Various Constants
*
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
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
      DO 940 IN = 1, NSIZES
*
         N = NN( IN )
         N1 = MAX( 1, N )
*
*        Do for each .TRUE. value in DOTYPE:
*
         MTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.EQ.MAXTYP+1 .AND. NSIZES.EQ.1 )
     $      MTYPES = NTYPES
         DO 930 ITYPE = 1, MTYPES
            IF( .NOT.DOTYPE( ITYPE ) )
     $         GO TO 930
*
*           Save random number seed for error messages
*
            DO 70 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   70       CONTINUE
*
*           Time the LAPACK Routines
*
*           Generate A and B
*
            IF( ITYPE.LE.MAXTYP ) THEN
*
*              Generate A (w/o rotation)
*
               CALL CLATM4( KATYPE( ITYPE ), N, 3, 1, .TRUE., ONE, ULP,
     $                      ONE, 2, ISEED, A, N1 )
               IF( 3.LE.N )
     $            A( 3+2*N1 ) = CONE
*
*              Generate B (w/o rotation)
*
               CALL CLATM4( 8, N, 3, 1, .FALSE., ONE, ONE, ONE, 2,
     $                      ISEED, B, N1 )
               IF( 2.LE.N )
     $            B( 2+N1 ) = CONE
*
               IF( N.GT.0 ) THEN
*
*                 Include rotations
*
*                 Generate U, V as Householder transformations times a
*                 diagonal matrix.  (Note that CLARFG makes Q(jc+ic)
*                 and Z(jc+ic) real.)
*
                  DO 90 JC = 1, N - 1
                     IC = ( JC-1 )*N1
                     DO 80 JR = JC, N
                        Q( JR+IC ) = CLARND( 3, ISEED )
                        Z( JR+IC ) = CLARND( 3, ISEED )
   80                CONTINUE
                     CALL CLARFG( N+1-JC, Q( JC+IC ), Q( JC+1+IC ), 1,
     $                            WORK( JC ) )
                     WORK( 2*N+JC ) = SIGN( ONE, REAL( Q( JC+IC ) ) )
                     Q( JC+IC ) = CONE
                     CALL CLARFG( N+1-JC, Z( JC+IC ), Z( JC+1+IC ), 1,
     $                            WORK( N+JC ) )
                     WORK( 3*N+JC ) = SIGN( ONE, REAL( Z( JC+IC ) ) )
                     Z( JC+IC ) = CONE
   90             CONTINUE
                  IC = ( N-1 )*N1
                  CTEMP = CLARND( 3, ISEED )
                  Q( N+IC ) = CONE
                  WORK( N ) = CZERO
                  WORK( 3*N ) = CTEMP / ABS( CTEMP )
                  CTEMP = CLARND( 3, ISEED )
                  Z( N+IC ) = CONE
                  WORK( 2*N ) = CZERO
                  WORK( 4*N ) = CTEMP / ABS( CTEMP )
*
*                 Apply the diagonal matrices
*
                  DO 110 JC = 1, N
                     DO 100 JR = 1, N
                        A( JR+IC ) = WORK( 2*N+JR )*
     $                               CONJG( WORK( 3*N+JC ) )*A( JR+IC )
                        B( JR+IC ) = WORK( 2*N+JR )*
     $                               CONJG( WORK( 3*N+JC ) )*B( JR+IC )
  100                CONTINUE
  110             CONTINUE
                  CALL CUNM2R( 'L', 'N', N, N, N-1, Q, N1, WORK, A, N1,
     $                         WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 120
                  CALL CUNM2R( 'R', 'C', N, N, N-1, Z, N1, WORK( N+1 ),
     $                         A, N1, WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 120
                  CALL CUNM2R( 'L', 'N', N, N, N-1, Q, N1, WORK, B, N1,
     $                         WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 120
                  CALL CUNM2R( 'R', 'C', N, N, N-1, Z, N1, WORK( N+1 ),
     $                         B, N1, WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 120
               END IF
  120          CONTINUE
            END IF
*
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*           Time CGGHRD
*
*           Time CGEQRF+CGGHRD('N','N',...) for each pair
*           (LDAS(j),NNB(j))
*
            IF( TIMSUB( 1 ) ) THEN
               DO 160 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = NNB( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 1 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 1 ) = ZERO
                     GO TO 160
                  END IF
*
*                 If this value of (NB,LDA) has occurred before,
*                 just use that value.
*
                  LASTL = 0
                  DO 130 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) .AND. NB.EQ.NNB( J ) )
     $                  LASTL = J
  130             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CGGHRD, computing neither Q nor Z
*                    (Actually, time CGEQRF + CUNMQR + CGGHRD.)
*
                     CALL XLAENV( 1, NB )
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  140                CONTINUE
                     CALL CLACPY( 'Full', N, N, A, N1, H, LDA )
                     CALL CLACPY( 'Full', N, N, B, N1, T, LDA )
                     CALL CLAQZH( .FALSE., .FALSE., N, 1, N, H, LDA, T,
     $                            LDA, Q, LDA, Z, LDA, WORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 1 )(1:ILA_LEN_TRIM( SUBNAM( 1 ) )), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
*
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 140
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 150 J = 1, IC
                        CALL CLACPY( 'Full', N, N, A, N1, Z, LDA )
                        CALL CLACPY( 'Full', N, N, B, N1, Z, LDA )
  150                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 1 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 1 ) = OPS / REAL( IC ) +
     $                  SOPLA( 'CGEQRF', N, N, 0, 0, NB ) +
     $                  SOPLA( 'CUNMQR', N, N, 0, 0, NB )
                     LDH = LDA
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 1 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 1 )
                     TIMES( IPAR, ITYPE, IN, 1 ) = TIMES( LASTL, ITYPE,
     $                  IN, 1 )
                  END IF
  160          CONTINUE
            ELSE IF( RUNHRD ) THEN
               CALL CLACPY( 'Full', N, N, A, N1, H, N1 )
               CALL CLACPY( 'Full', N, N, B, N1, T, N1 )
               CALL CLAQZH( .FALSE., .FALSE., N, 1, N, H, N1, T, N1, Q,
     $                      N1, Z, N1, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 1 )(1:ILA_LEN_TRIM( SUBNAM( 1 ) )), IINFO, N,
     $               ITYPE, 0, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 930
               END IF
               LDH = N
            END IF
*
*           Time CGGHRD('I','N',...) for each pair (LDAS(j),NNB(j))
*
            IF( TIMSUB( 2 ) ) THEN
               DO 200 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = NNB( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 2 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 2 ) = ZERO
                     GO TO 200
                  END IF
*
*                 If this value of (NB,LDA) has occurred before,
*                 just use that value.
*
                  LASTL = 0
                  DO 170 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) .AND. NB.EQ.NNB( J ) )
     $                  LASTL = J
  170             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CGGHRD, computing Q but not Z
*                    (Actually, CGEQRF + CUNMQR + CUNGQR + CGGHRD.)
*
                     CALL XLAENV( 1, NB )
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  180                CONTINUE
                     CALL CLACPY( 'Full', N, N, A, N1, H, LDA )
                     CALL CLACPY( 'Full', N, N, B, N1, T, LDA )
                     CALL CLAQZH( .TRUE., .FALSE., N, 1, N, H, LDA, T,
     $                            LDA, Q, LDA, Z, LDA, WORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 2 )(1:ILA_LEN_TRIM( SUBNAM( 2 ) )), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
*
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 180
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 190 J = 1, IC
                        CALL CLACPY( 'Full', N, N, A, N1, Z, LDA )
                        CALL CLACPY( 'Full', N, N, B, N1, Z, LDA )
  190                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 2 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 2 ) = OPS / REAL( IC ) +
     $                  SOPLA( 'CGEQRF', N, N, 0, 0, NB ) +
     $                  SOPLA( 'CUNMQR', N, N, 0, 0, NB ) +
     $                  SOPLA( 'CUNGQR', N, N, 0, 0, NB )
                     LDH = LDA
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 2 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 2 )
                     TIMES( IPAR, ITYPE, IN, 2 ) = TIMES( LASTL, ITYPE,
     $                  IN, 2 )
                  END IF
  200          CONTINUE
            END IF
*
*           Time CGGHRD('N','I',...) for each pair (LDAS(j),NNB(j))
*
            IF( TIMSUB( 3 ) ) THEN
               DO 240 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = NNB( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 3 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 3 ) = ZERO
                     GO TO 240
                  END IF
*
*                 If this value of (NB,LDA) has occurred before,
*                 just use that value.
*
                  LASTL = 0
                  DO 210 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) .AND. NB.EQ.NNB( J ) )
     $                  LASTL = J
  210             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CGGHRD, computing Z but not Q
*                    (Actually, CGEQRF + CUNMQR + CGGHRD.)
*
                     CALL XLAENV( 1, NB )
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  220                CONTINUE
                     CALL CLACPY( 'Full', N, N, A, N1, H, LDA )
                     CALL CLACPY( 'Full', N, N, B, N1, T, LDA )
                     CALL CLAQZH( .FALSE., .TRUE., N, 1, N, H, LDA, T,
     $                            LDA, Q, LDA, Z, LDA, WORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 3 )(1:ILA_LEN_TRIM( SUBNAM( 3 ) )), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
*
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 220
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 230 J = 1, IC
                        CALL CLACPY( 'Full', N, N, A, N1, Z, LDA )
                        CALL CLACPY( 'Full', N, N, B, N1, Z, LDA )
  230                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 3 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 3 ) = OPS / REAL( IC ) +
     $                  SOPLA( 'CGEQRF', N, N, 0, 0, NB ) +
     $                  SOPLA( 'CUNMQR', N, N, 0, 0, NB )
                     LDH = LDA
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 3 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 3 )
                     TIMES( IPAR, ITYPE, IN, 3 ) = TIMES( LASTL, ITYPE,
     $                  IN, 3 )
                  END IF
  240          CONTINUE
            END IF
*
*           Time CGGHRD('I','I',...) for each pair (LDAS(j),NNB(j))
*
            IF( TIMSUB( 4 ) ) THEN
               DO 280 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = NNB( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 4 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 4 ) = ZERO
                     GO TO 280
                  END IF
*
*                 If this value of (NB,LDA) has occurred before,
*                 just use that value.
*
                  LASTL = 0
                  DO 250 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) .AND. NB.EQ.NNB( J ) )
     $                  LASTL = J
  250             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CGGHRD, computing Q and Z
*                    (Actually, CGEQRF + CUNMQR + CUNGQR + CGGHRD.)
*
                     CALL XLAENV( 1, NB )
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  260                CONTINUE
                     CALL CLACPY( 'Full', N, N, A, N1, H, LDA )
                     CALL CLACPY( 'Full', N, N, B, N1, T, LDA )
                     CALL CLAQZH( .TRUE., .TRUE., N, 1, N, H, LDA, T,
     $                            LDA, Q, LDA, Z, LDA, WORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 4 )(1:ILA_LEN_TRIM( SUBNAM( 4 ) )), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
*
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 260
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 270 J = 1, IC
                        CALL CLACPY( 'Full', N, N, A, N1, Z, LDA )
                        CALL CLACPY( 'Full', N, N, B, N1, Z, LDA )
  270                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 4 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 4 ) = OPS / REAL( IC ) +
     $                  SOPLA( 'CGEQRF', N, N, 0, 0, NB ) +
     $                  SOPLA( 'CUNMQR', N, N, 0, 0, NB ) +
     $                  SOPLA( 'CUNGQR', N, N, 0, 0, NB )
                     LDH = LDA
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 4 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 4 )
                     TIMES( IPAR, ITYPE, IN, 4 ) = TIMES( LASTL, ITYPE,
     $                  IN, 4 )
                  END IF
  280          CONTINUE
            END IF
*
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*           Time CHGEQZ
*
*           Time CHGEQZ with JOB='E' for each value of LDAS(j)
*
            IF( TIMSUB( 5 ) ) THEN
               DO 320 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 5 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 5 ) = ZERO
                     GO TO 320
                  END IF
*
*                 If this value of LDA has occurred before,
*                 just use that value.
*
                  LASTL = 0
                  DO 290 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  290             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CHGEQZ with JOB='E'
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  300                CONTINUE
                     CALL CLACPY( 'Full', N, N, H, LDH, A, LDA )
                     CALL CLACPY( 'Full', N, N, T, LDH, B, LDA )
                     CALL CHGEQZ( 'E', 'N', 'N', N, 1, N, A, LDA, B,
     $                            LDA, W, W( LDA+1 ), Q, LDA, Z, LDA,
     $                            WORK, LWORK, RWORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 5 )(1:ILA_LEN_TRIM( SUBNAM( 5 ) )), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 300
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 310 J = 1, IC
                        CALL CLACPY( 'Full', N, N, H, LDH, Z, LDA )
                        CALL CLACPY( 'Full', N, N, T, LDH, Z, LDA )
  310                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 5 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 5 ) = OPS / REAL( IC )
                     LDS = 0
                     LDQ = 0
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 5 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 5 )
                     TIMES( IPAR, ITYPE, IN, 5 ) = TIMES( LASTL, ITYPE,
     $                  IN, 5 )
                  END IF
  320          CONTINUE
            END IF
*
*           Time CHGEQZ with JOB='S', COMPQ=COMPZ='N' for each value
*           of LDAS(j)
*
            IF( TIMSUB( 6 ) ) THEN
               DO 360 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 6 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 6 ) = ZERO
                     GO TO 360
                  END IF
*
*                 If this value of LDA has occurred before,
*                 just use that value.
*
                  LASTL = 0
                  DO 330 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  330             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                 Time CHGEQZ with JOB='S', COMPQ=COMPZ='N'
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  340                CONTINUE
                     CALL CLACPY( 'Full', N, N, H, LDH, A, LDA )
                     CALL CLACPY( 'Full', N, N, T, LDH, B, LDA )
                     CALL CHGEQZ( 'S', 'N', 'N', N, 1, N, A, LDA, B,
     $                            LDA, W, W( LDA+1 ), Q, LDA, Z, LDA,
     $                            WORK, LWORK, RWORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 6 )(1:ILA_LEN_TRIM( SUBNAM( 6 ) )), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 340
*
*                 Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 350 J = 1, IC
                        CALL CLACPY( 'Full', N, N, H, LDH, Z, LDA )
                        CALL CLACPY( 'Full', N, N, T, LDH, Z, LDA )
  350                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 6 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 6 ) = OPS / REAL( IC )
                     LDS = LDA
                     LDQ = 0
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 6 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 6 )
                     TIMES( IPAR, ITYPE, IN, 6 ) = TIMES( LASTL, ITYPE,
     $                  IN, 6 )
                  END IF
  360          CONTINUE
            ELSE IF( RUNES ) THEN
               CALL CLACPY( 'Full', N, N, H, LDH, A, N1 )
               CALL CLACPY( 'Full', N, N, T, LDH, B, N1 )
               CALL CHGEQZ( 'S', 'N', 'N', N, 1, N, A, N1, B, N1, W,
     $                      W( N1+1 ), Q, N1, Z, N1, WORK,
     $                      LWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 6 )(1:ILA_LEN_TRIM( SUBNAM( 6 ) )), IINFO, N,
     $               ITYPE, 0, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 930
               END IF
               LDS = N1
               LDQ = 0
            END IF
*
*           Time CHGEQZ with JOB='S', COMPQ='I', COMPZ='N' for each
*           value of LDAS(j)
*
            IF( TIMSUB( 7 ) ) THEN
               DO 400 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 7 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 7 ) = ZERO
                     GO TO 400
                  END IF
*
*                 If this value of LDA has occurred before,
*                 just use that value.
*
                  LASTL = 0
                  DO 370 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  370             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                 Time CHGEQZ with JOB='S', COMPQ='I', COMPZ='N'
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  380                CONTINUE
                     CALL CLACPY( 'Full', N, N, H, LDH, A, LDA )
                     CALL CLACPY( 'Full', N, N, T, LDH, B, LDA )
                     CALL CHGEQZ( 'S', 'I', 'N', N, 1, N, A, LDA, B,
     $                            LDA, W, W( LDA+1 ), Q, LDA, Z, LDA,
     $                            WORK, LWORK, RWORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 7 )(1:ILA_LEN_TRIM( SUBNAM( 7 ) )), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 380
*
*                 Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 390 J = 1, IC
                        CALL CLACPY( 'Full', N, N, H, LDH, Z, LDA )
                        CALL CLACPY( 'Full', N, N, T, LDH, Z, LDA )
  390                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 7 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 7 ) = OPS / REAL( IC )
                     LDS = LDA
                     LDQ = LDA
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 7 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 7 )
                     TIMES( IPAR, ITYPE, IN, 7 ) = TIMES( LASTL, ITYPE,
     $                  IN, 7 )
                  END IF
  400          CONTINUE
            ELSE IF( RUNEQ ) THEN
               CALL CLACPY( 'Full', N, N, H, LDH, A, N1 )
               CALL CLACPY( 'Full', N, N, T, LDH, B, N1 )
               CALL CHGEQZ( 'S', 'I', 'N', N, 1, N, A, N1, B, N1, W,
     $                      W( N1+1 ), Q, N1, Z, N1, WORK,
     $                      LWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 7 )(1:ILA_LEN_TRIM( SUBNAM( 7 ) )), IINFO, N,
     $               ITYPE, 0, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 930
               END IF
               LDS = N1
               LDQ = N1
            END IF
*
*           Time CHGEQZ with JOB='S', COMPQ='N', COMPZ='I' for each
*           value of LDAS(j)
*
            IF( TIMSUB( 8 ) ) THEN
               DO 440 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 8 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 8 ) = ZERO
                     GO TO 440
                  END IF
*
*                 If this value of LDA has occurred before,
*                 just use that value.
*
                  LASTL = 0
                  DO 410 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  410             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
                     NB = MIN( N, NNB( IPAR ) )
                     NSHIFT = NSHFTS( IPAR )
                     NEISP = NEISPS( IPAR )
                     MINNB = MINNBS( IPAR )
                     MINBLK = MINBKS( IPAR )
                     CALL XLAENV( 1, NB )
                     CALL XLAENV( 2, MINNB )
                     CALL XLAENV( 8, NEISP )
                     CALL XLAENV( 4, NSHIFT )
                     CALL XLAENV( 5, MINBLK )
*
*                 Time CHGEQZ with JOB='S', COMPQ='N', COMPZ='I'
*                 (Note that the "Z" matrix is stored in the array Q)
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  420                CONTINUE
                     CALL CLACPY( 'Full', N, N, H, LDH, A, LDA )
                     CALL CLACPY( 'Full', N, N, T, LDH, B, LDA )
                     CALL CHGEQZ( 'S', 'N', 'I', N, 1, N, A, LDA, B,
     $                            LDA, W, W( LDA+1 ), Z, LDA, Q, LDA,
     $                            WORK, LWORK, RWORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 8 )(1:ILA_LEN_TRIM( SUBNAM( 8 ) )), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 420
*
*                 Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 430 J = 1, IC
                        CALL CLACPY( 'Full', N, N, H, LDH, Z, LDA )
                        CALL CLACPY( 'Full', N, N, T, LDH, Z, LDA )
  430                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 8 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 8 ) = OPS / REAL( IC )
                     LDS = LDA
                     LDQ = LDA
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 8 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 8 )
                     TIMES( IPAR, ITYPE, IN, 8 ) = TIMES( LASTL, ITYPE,
     $                  IN, 8 )
                  END IF
  440          CONTINUE
            END IF
*
*           Time CHGEQZ with JOB='S', COMPQ='I', COMPZ='I' for each
*           value of LDAS(j)
*
            IF( TIMSUB( 9 ) ) THEN
               DO 480 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 9 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 9 ) = ZERO
                     GO TO 480
                  END IF
*
*                 If this value of LDA has occurred before,
*                 just use that value.
*
                  LASTL = 0
                  DO 450 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  450             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                 Time CHGEQZ with JOB='S', COMPQ='I', COMPZ='I'
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  460                CONTINUE
                     CALL CLACPY( 'Full', N, N, H, LDH, A, LDA )
                     CALL CLACPY( 'Full', N, N, T, LDH, B, LDA )
                     CALL CHGEQZ( 'S', 'I', 'I', N, 1, N, A, LDA, B,
     $                            LDA, W, W( LDA+1 ), Q, LDA, Z, LDA,
     $                            WORK, LWORK, RWORK, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 9 )(1:ILA_LEN_TRIM( SUBNAM( 9 ) )), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 460
*
*                 Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 470 J = 1, IC
                        CALL CLACPY( 'Full', N, N, H, LDH, Z, LDA )
                        CALL CLACPY( 'Full', N, N, T, LDH, Z, LDA )
  470                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 9 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 9 ) = OPS / REAL( IC )
                     LDS = LDA
                     LDQ = LDA
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 9 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 9 )
                     TIMES( IPAR, ITYPE, IN, 9 ) = TIMES( LASTL, ITYPE,
     $                  IN, 9 )
                  END IF
  480          CONTINUE
            END IF
*
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*           Time CTGEVC
*
            IF( TIMSUB( 10 ) .OR. TIMSUB( 11 ) .OR. TIMSUB( 12 ) .OR.
     $          TIMSUB( 13 ) ) THEN
               DO 610 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     DO 490 J = 10, 13
                        IF( TIMSUB( J ) ) THEN
                           TIMES( IPAR, ITYPE, IN, J ) = ZERO
                           OPCNTS( IPAR, ITYPE, IN, J ) = ZERO
                        END IF
  490                CONTINUE
                     GO TO 610
                  END IF
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 500 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  500             CONTINUE
*
*                 Time CTGEVC if this is a new value of LDA
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Copy S (which is in A) and P (which is in B)
*                    if necessary to get right LDA.
*
                     IF( LDA.GT.LDS ) THEN
                        DO 520 JC = N, 1, -1
                           DO 510 JR = N, 1, -1
                              A( JR+( JC-1 )*LDA ) = A( JR+( JC-1 )*
     $                           LDS )
                              B( JR+( JC-1 )*LDA ) = B( JR+( JC-1 )*
     $                           LDS )
  510                      CONTINUE
  520                   CONTINUE
                     ELSE IF( LDA.LT.LDS ) THEN
                        DO 540 JC = 1, N
                           DO 530 JR = 1, N
                              A( JR+( JC-1 )*LDA ) = A( JR+( JC-1 )*
     $                           LDS )
                              B( JR+( JC-1 )*LDA ) = B( JR+( JC-1 )*
     $                           LDS )
  530                      CONTINUE
  540                   CONTINUE
                     END IF
                     LDS = LDA
*
*                    Time CTGEVC for Left Eigenvectors only,
*                    without back transforming
*
                     IF( TIMSUB( 10 ) ) THEN
                        IC = 0
                        OPS = ZERO
                        S1 = SECOND( )
  550                   CONTINUE
                        CALL CTGEVC( 'L', 'A', LLWORK, N, A, LDA, B,
     $                               LDA, H, LDA, T, LDA, N, ITEMP,
     $                               WORK, RWORK, IINFO )
                        IF( IINFO.NE.0 ) THEN
                           WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 10 )(1:ILA_LEN_TRIM( SUBNAM( 10 ) )),
     $                        IINFO, N, ITYPE, IPAR, IOLDSD
                           INFO = ABS( IINFO )
                           GO TO 930
                        END IF
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN )
     $                     GO TO 550
*
                        TIMES( IPAR, ITYPE, IN, 10 ) = TIME / REAL( IC )
                        OPCNTS( IPAR, ITYPE, IN, 10 ) = OPS / REAL( IC )
                     END IF
*
*                    Time CTGEVC for Left Eigenvectors only,
*                    with back transforming
*
                     IF( TIMSUB( 11 ) ) THEN
                        IC = 0
                        OPS = ZERO
                        S1 = SECOND( )
  560                   CONTINUE
                        CALL CLACPY( 'Full', N, N, Q, LDQ, H, LDA )
                        CALL CTGEVC( 'L', 'B', LLWORK, N, A, LDA, B,
     $                               LDA, H, LDA, T, LDA, N, ITEMP,
     $                               WORK, RWORK, IINFO )
                        IF( IINFO.NE.0 ) THEN
                           WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 11 )(1:ILA_LEN_TRIM( SUBNAM( 11 ) )),
     $                        IINFO, N, ITYPE, IPAR, IOLDSD
                           INFO = ABS( IINFO )
                           GO TO 930
                        END IF
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN )
     $                     GO TO 560
*
*                       Subtract the time used in CLACPY.
*
                        S1 = SECOND( )
                        DO 570 J = 1, IC
                           CALL CLACPY( 'Full', N, N, Q, LDQ, H, LDA )
  570                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
*
                        TIMES( IPAR, ITYPE, IN, 11 ) = MAX( TIME-UNTIME,
     $                     ZERO ) / REAL( IC )
                        OPCNTS( IPAR, ITYPE, IN, 11 ) = OPS / REAL( IC )
                     END IF
*
*                    Time CTGEVC for Right Eigenvectors only,
*                    without back transforming
*
                     IF( TIMSUB( 12 ) ) THEN
                        IC = 0
                        OPS = ZERO
                        S1 = SECOND( )
  580                   CONTINUE
                        CALL CTGEVC( 'R', 'A', LLWORK, N, A, LDA, B,
     $                               LDA, H, LDA, T, LDA, N, ITEMP,
     $                               WORK, RWORK, IINFO )
                        IF( IINFO.NE.0 ) THEN
                           WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 12 )(1:ILA_LEN_TRIM( SUBNAM( 12 ) )),
     $                        IINFO, N, ITYPE, IPAR, IOLDSD
                           INFO = ABS( IINFO )
                           GO TO 930
                        END IF
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN )
     $                     GO TO 580
*
                        TIMES( IPAR, ITYPE, IN, 12 ) = TIME / REAL( IC )
                        OPCNTS( IPAR, ITYPE, IN, 12 ) = OPS / REAL( IC )
                     END IF
*
*                    Time CTGEVC for Right Eigenvectors only,
*                    with back transforming
*
                     IF( TIMSUB( 13 ) ) THEN
                        IC = 0
                        OPS = ZERO
                        S1 = SECOND( )
  590                   CONTINUE
                        CALL CLACPY( 'Full', N, N, Q, LDQ, T, LDA )
                        CALL CTGEVC( 'R', 'B', LLWORK, N, A, LDA, B,
     $                               LDA, H, LDA, T, LDA, N, ITEMP,
     $                               WORK, RWORK, IINFO )
                        IF( IINFO.NE.0 ) THEN
                           WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 13 )(1:ILA_LEN_TRIM( SUBNAM( 13 ) )),
     $                        IINFO, N, ITYPE, IPAR, IOLDSD
                           INFO = ABS( IINFO )
                           GO TO 930
                        END IF
                        S2 = SECOND( )
                        TIME = S2 - S1
                        IC = IC + 1
                        IF( TIME.LT.TIMMIN )
     $                     GO TO 590
*
*                       Subtract the time used in CLACPY.
*
                        S1 = SECOND( )
                        DO 600 J = 1, IC
                           CALL CLACPY( 'Full', N, N, Q, LDQ, T, LDA )
  600                   CONTINUE
                        S2 = SECOND( )
                        UNTIME = S2 - S1
*
                        TIMES( IPAR, ITYPE, IN, 13 ) = MAX( TIME-UNTIME,
     $                     ZERO ) / REAL( IC )
                        OPCNTS( IPAR, ITYPE, IN, 13 ) = OPS / REAL( IC )
                     END IF
*
                  ELSE
*
*                    If this LDA has previously appeared, use the
*                    previously computed value(s).
*
                     IF( TIMSUB( 10 ) ) THEN
                        OPCNTS( IPAR, ITYPE, IN, 10 ) = OPCNTS( LASTL,
     $                     ITYPE, IN, 10 )
                        TIMES( IPAR, ITYPE, IN, 10 ) = TIMES( LASTL,
     $                     ITYPE, IN, 10 )
                     END IF
                     IF( TIMSUB( 11 ) ) THEN
                        OPCNTS( IPAR, ITYPE, IN, 11 ) = OPCNTS( LASTL,
     $                     ITYPE, IN, 11 )
                        TIMES( IPAR, ITYPE, IN, 11 ) = TIMES( LASTL,
     $                     ITYPE, IN, 11 )
                     END IF
                     IF( TIMSUB( 12 ) ) THEN
                        OPCNTS( IPAR, ITYPE, IN, 12 ) = OPCNTS( LASTL,
     $                     ITYPE, IN, 12 )
                        TIMES( IPAR, ITYPE, IN, 12 ) = TIMES( LASTL,
     $                     ITYPE, IN, 12 )
                     END IF
                     IF( TIMSUB( 13 ) ) THEN
                        OPCNTS( IPAR, ITYPE, IN, 13 ) = OPCNTS( LASTL,
     $                     ITYPE, IN, 13 )
                        TIMES( IPAR, ITYPE, IN, 13 ) = TIMES( LASTL,
     $                     ITYPE, IN, 13 )
                     END IF
                  END IF
  610          CONTINUE
            END IF
*
*           Time the EISPACK Routines
*
*           Restore random number seed
*
            DO 620 J = 1, 4
               ISEED( J ) = IOLDSD( J )
  620       CONTINUE
*
*           Re-generate A
*
            IF( ITYPE.LE.MAXTYP ) THEN
*
*              Generate A (w/o rotation)
*
               CALL CLATM4( KATYPE( ITYPE ), N, 3, 1, .TRUE., ONE, ULP,
     $                      ONE, 2, ISEED, H, N1 )
               IF( N.GE.3 )
     $            H( 3+2*N1 ) = ONE
*
*              Generate B (w/o rotation)
*
               CALL CLATM4( 8, N, 3, 1, .FALSE., ONE, ONE, ONE, 2,
     $                      ISEED, T, N1 )
               IF( N.GE.2 )
     $            T( 2+N1 ) = ONE
*
               IF( N.GT.0 ) THEN
*
*                 Include rotations
*
*                 Generate U, V as Householder transformations times a
*                 diagonal matrix.  (Note that CLARFG makes Q(jc+ic)
*                 and Z(jc+ic) real.)
*
                  DO 640 JC = 1, N - 1
                     IC = ( JC-1 )*N1
                     DO 630 JR = JC, N
                        Q( JR+IC ) = CLARND( 3, ISEED )
                        Z( JR+IC ) = CLARND( 3, ISEED )
  630                CONTINUE
                     CALL CLARFG( N+1-JC, Q( JC+IC ), Q( JC+1+IC ), 1,
     $                            WORK( JC ) )
                     WORK( 2*N+JC ) = SIGN( ONE, REAL( Q( JC+IC ) ) )
                     Q( JC+IC ) = ONE
                     CALL CLARFG( N+1-JC, Z( JC+IC ), Z( JC+1+IC ), 1,
     $                            WORK( N+JC ) )
                     WORK( 3*N+JC ) = SIGN( ONE, REAL( Z( JC+IC ) ) )
                     Z( JC+IC ) = ONE
  640             CONTINUE
                  IC = ( N-1 )*N1
                  CTEMP = CLARND( 3, ISEED )
                  Q( N+IC ) = CONE
                  WORK( N ) = CZERO
                  WORK( 3*N ) = CTEMP / ABS( CTEMP )
                  CTEMP = CLARND( 3, ISEED )
                  Z( N+IC ) = CONE
                  WORK( 2*N ) = CZERO
                  WORK( 4*N ) = CTEMP / ABS( CTEMP )
*
*                 Apply the diagonal matrices
*
                  DO 660 JC = 1, N
                     DO 650 JR = 1, N
                        H( JR+IC ) = WORK( 2*N+JR )*
     $                               CONJG( WORK( 3*N+JC ) )*H( JR+IC )
                        T( JR+IC ) = WORK( 2*N+JR )*WORK( 3*N+JC )*
     $                               CONJG( WORK( 3*N+JC ) )*T( JR+IC )
  650                CONTINUE
  660             CONTINUE
                  CALL CUNM2R( 'L', 'N', N, N, N-1, Q, N1, WORK, H, N1,
     $                         WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 670
                  CALL CUNM2R( 'R', 'C', N, N, N-1, Z, N1, WORK( N+1 ),
     $                         H, N1, WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 670
                  CALL CUNM2R( 'L', 'N', N, N, N-1, Q, N1, WORK, T, N1,
     $                         WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 670
                  CALL CUNM2R( 'R', 'C', N, N, N-1, Z, N1, WORK( N+1 ),
     $                         T, N1, WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 670
               END IF
  670          CONTINUE
*
*              Copy real and imaginary parts into separate arrays
*
               DO 680 J = 1, N*N
                  AR( J ) = REAL( H( J ) )
                  AI( J ) = AIMAG( H( J ) )
                  BR( J ) = REAL( T( J ) )
                  BI( J ) = AIMAG( T( J ) )
  680          CONTINUE
            END IF
*
*           Time CQZHES w/ MATZ=.FALSE. for each LDAS(j)
*
            IF( TIMSUB( 14 ) ) THEN
               DO 720 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 14 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 14 ) = ZERO
                     GO TO 720
                  END IF
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 690 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  690             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CQZHES( ...,.FALSE.,..)
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  700                CONTINUE
                     CALL SLACPY( 'Full', N, N, AR, N1, HR, LDA )
                     CALL SLACPY( 'Full', N, N, AI, N1, HI, LDA )
                     CALL SLACPY( 'Full', N, N, BR, N1, TR, LDA )
                     CALL SLACPY( 'Full', N, N, BI, N1, TI, LDA )
                     CALL CQZHES( LDA, N, HR, HI, TR, TI, .FALSE., QR,
     $                            QI )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 700
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 710 J = 1, IC
                        CALL SLACPY( 'Full', N, N, AR, N1, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, AI, N1, ZI, LDA )
                        CALL SLACPY( 'Full', N, N, BR, N1, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, BI, N1, ZI, LDA )
  710                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 14 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 14 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 14 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 14 )
                     TIMES( IPAR, ITYPE, IN, 14 ) = TIMES( LASTL, ITYPE,
     $                  IN, 14 )
                  END IF
                  LDH = LDA
  720          CONTINUE
            ELSE IF( RUNHES ) THEN
               CALL SLACPY( 'Full', N, N, AR, N1, HR, N1 )
               CALL SLACPY( 'Full', N, N, AI, N1, HI, N1 )
               CALL SLACPY( 'Full', N, N, BR, N1, TR, N1 )
               CALL SLACPY( 'Full', N, N, BI, N1, TI, N1 )
               CALL CQZHES( N1, N, HR, HI, TR, TI, .FALSE., QR, QI )
               LDH = N1
            END IF
*
*           Time CQZHES w/ MATZ=.TRUE. for each LDAS(j)
*
            IF( TIMSUB( 15 ) ) THEN
               DO 760 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 15 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 15 ) = ZERO
                     GO TO 760
                  END IF
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 730 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  730             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CQZHES( ...,.TRUE.,..)
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  740                CONTINUE
                     CALL SLACPY( 'Full', N, N, AR, N1, HR, LDA )
                     CALL SLACPY( 'Full', N, N, AI, N1, HI, LDA )
                     CALL SLACPY( 'Full', N, N, BR, N1, TR, LDA )
                     CALL SLACPY( 'Full', N, N, BI, N1, TI, LDA )
                     CALL CQZHES( LDA, N, HR, HI, TR, TI, .TRUE., QR,
     $                            QI )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 740
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 750 J = 1, IC
                        CALL SLACPY( 'Full', N, N, AR, N1, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, AI, N1, ZI, LDA )
                        CALL SLACPY( 'Full', N, N, BR, N1, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, BI, N1, ZI, LDA )
  750                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 15 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 15 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 15 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 15 )
                     TIMES( IPAR, ITYPE, IN, 15 ) = TIMES( LASTL, ITYPE,
     $                  IN, 15 )
                  END IF
                  LDH = LDA
  760          CONTINUE
            END IF
*
*           Time CQZVAL w/ MATZ=.FALSE. for each LDAS(j)
*
            IF( TIMSUB( 16 ) ) THEN
               DO 800 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 16 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 16 ) = ZERO
                     GO TO 800
                  END IF
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 770 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  770             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CQZVAL with MATZ=.FALSE.
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  780                CONTINUE
                     CALL SLACPY( 'Full', N, N, HR, LDH, AR, LDA )
                     CALL SLACPY( 'Full', N, N, HI, LDH, AI, LDA )
                     CALL SLACPY( 'Full', N, N, TR, LDH, BR, LDA )
                     CALL SLACPY( 'Full', N, N, TI, LDH, BI, LDA )
                     CALL CQZVAL( LDA, N, AR, AI, BR, BI, ZERO, WR,
     $                            WR( LDA+1 ), WR( 2*LDA+1 ), .FALSE.,
     $                            QR, QI, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 16 )(1:ILA_LEN_TRIM( SUBNAM( 16 ) )), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
*
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 780
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 790 J = 1, IC
                        CALL SLACPY( 'Full', N, N, HR, LDH, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, HI, LDH, ZI, LDA )
                        CALL SLACPY( 'Full', N, N, TR, LDH, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, TI, LDH, ZI, LDA )
  790                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 16 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 16 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 16 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 16 )
                     TIMES( IPAR, ITYPE, IN, 16 ) = TIMES( LASTL, ITYPE,
     $                  IN, 16 )
                  END IF
                  LDS = 0
                  LDW = LDA
  800          CONTINUE
            END IF
*
*           Time CQZVAL w/ MATZ=.TRUE. for each LDAS(j)
*
            IF( TIMSUB( 17 ) ) THEN
               DO 840 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 17 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 17 ) = ZERO
                     GO TO 840
                  END IF
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 810 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  810             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time CQZVAL with MATZ=.TRUE.
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  820                CONTINUE
                     CALL SLACPY( 'Full', N, N, HR, LDH, AR, LDA )
                     CALL SLACPY( 'Full', N, N, HI, LDH, AI, LDA )
                     CALL SLACPY( 'Full', N, N, TR, LDH, BR, LDA )
                     CALL SLACPY( 'Full', N, N, TI, LDH, BI, LDA )
                     CALL SLASET( 'Full', N, N, ZERO, ONE, QR, LDA )
                     CALL SLASET( 'Full', N, N, ZERO, ONE, QI, LDA )
                     CALL CQZVAL( LDA, N, AR, AI, BR, BI, ZERO, WR,
     $                            WR( LDA+1 ), WR( 2*LDA+1 ), .TRUE.,
     $                            QR, QI, IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 17 )(1:ILA_LEN_TRIM( SUBNAM( 17 ) )), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 930
                     END IF
*
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 820
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 830 J = 1, IC
                        CALL SLACPY( 'Full', N, N, HR, LDH, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, HI, LDH, ZI, LDA )
                        CALL SLACPY( 'Full', N, N, TR, LDH, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, TI, LDH, ZI, LDA )
                        CALL SLASET( 'Full', N, N, ZERO, ONE, ZR, LDA )
                        CALL SLASET( 'Full', N, N, ZERO, ONE, ZI, LDA )
  830                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 17 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 17 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 17 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 17 )
                     TIMES( IPAR, ITYPE, IN, 17 ) = TIMES( LASTL, ITYPE,
     $                  IN, 17 )
                  END IF
                  LDS = LDA
                  LDW = LDA
  840          CONTINUE
            ELSE IF( RUNQZ ) THEN
               CALL SLACPY( 'Full', N, N, HR, LDH, AR, N1 )
               CALL SLACPY( 'Full', N, N, HI, LDH, AI, N1 )
               CALL SLACPY( 'Full', N, N, TR, LDH, BR, N1 )
               CALL SLACPY( 'Full', N, N, TI, LDH, BI, N1 )
               CALL SLASET( 'Full', N, N, ZERO, ONE, QR, N1 )
               CALL SLASET( 'Full', N, N, ZERO, ONE, QI, N1 )
               CALL CQZVAL( N1, N, AR, AI, BR, BI, ZERO, WR, WR( N1+1 ),
     $                      WR( 2*N1+1 ), .TRUE., QR, QI, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( 17 )(1:ILA_LEN_TRIM( SUBNAM( 17 ) )), IINFO, N,
     $               ITYPE, IPAR, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 930
               END IF
*
               LDS = N1
               LDW = N1
            END IF
*
*           Time CQZVEC for each LDAS(j)
*
            IF( TIMSUB( 18 ) ) THEN
               DO 920 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  IF( LDA.LT.N1 ) THEN
                     TIMES( IPAR, ITYPE, IN, 18 ) = ZERO
                     OPCNTS( IPAR, ITYPE, IN, 18 ) = ZERO
                     GO TO 920
                  END IF
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 850 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  850             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Copy W if necessary to get right LDA.
*
                     IF( LDA.GT.LDW ) THEN
                        DO 870 JC = 3, 1, -1
                           DO 860 JR = N, 1, -1
                              WR( JR+( JC-1 )*LDA ) = WR( JR+( JC-1 )*
     $                           LDW )
  860                      CONTINUE
  870                   CONTINUE
                     ELSE IF( LDA.LT.LDW ) THEN
                        DO 890 JC = 1, 3
                           DO 880 JR = 1, N
                              WR( JR+( JC-1 )*LDA ) = WR( JR+( JC-1 )*
     $                           LDW )
  880                      CONTINUE
  890                   CONTINUE
                     END IF
                     LDW = LDA
*
*                    Time CQZVEC
*
                     IC = 0
                     OPS = ZERO
                     S1 = SECOND( )
  900                CONTINUE
                     CALL SLACPY( 'Full', N, N, AR, LDS, HR, LDA )
                     CALL SLACPY( 'Full', N, N, AI, LDS, HI, LDA )
                     CALL SLACPY( 'Full', N, N, BR, LDS, TR, LDA )
                     CALL SLACPY( 'Full', N, N, BI, LDS, TI, LDA )
                     CALL SLACPY( 'Full', N, N, QR, LDS, ZR, LDA )
                     CALL SLACPY( 'Full', N, N, QI, LDS, ZI, LDA )
                     CALL CQZVEC( LDA, N, HR, HI, TR, TI, WR,
     $                            WR( LDA+1 ), WR( 2*LDA+1 ), ZR, ZI )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 900
*
*                    Subtract the time used in CLACPY.
*
                     S1 = SECOND( )
                     DO 910 J = 1, IC
                        CALL SLACPY( 'Full', N, N, AR, LDS, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, AI, LDS, ZI, LDA )
                        CALL SLACPY( 'Full', N, N, BR, LDS, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, BI, LDS, ZI, LDA )
                        CALL SLACPY( 'Full', N, N, QR, LDS, ZR, LDA )
                        CALL SLACPY( 'Full', N, N, QI, LDS, ZI, LDA )
  910                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 18 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / REAL( IC )
                     OPCNTS( IPAR, ITYPE, IN, 18 ) = OPS / REAL( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 18 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 18 )
                     TIMES( IPAR, ITYPE, IN, 18 ) = TIMES( LASTL, ITYPE,
     $                  IN, 18 )
                  END IF
  920          CONTINUE
            END IF
*
  930    CONTINUE
  940 CONTINUE
*
*     Print a table of results for each timed routine.
*
      DO 950 ISUB = 1, NSUBS
         IF( TIMSUB( ISUB ) ) THEN
            CALL SPRTBG( SUBNAM( ISUB ), MTYPES, DOTYPE, NSIZES, NN,
     $                   INPARM( ISUB ), PNAMES, NPARMS, LDAS, NNB,
     $                   NSHFTS, NEISPS, MINNBS, MINBKS,
     $                   OPCNTS( 1, 1, 1, ISUB ), LDO1, LDO2,
     $                   TIMES( 1, 1, 1, ISUB ), LDT1, LDT2, RWORK,
     $                   LLWORK, NOUT )
         END IF
  950 CONTINUE
*
      RETURN
*
*     End of CTIM51
*
 9997 FORMAT( ' CTIM51: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', ITYPE=', I6, ', IPAR=', I6, ', ISEED=(',
     $      3( I5, ',' ), I5, ')' )
*
      END
