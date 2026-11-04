*> \brief \b SSTEDS computes all eigenvalues and selected eigenvectors of
*> a symmetric tridiagonal matrix by boundary-row divide and conquer
*> and top-down selected-vector reconstruction.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SSTEDS( N, D, E, M, ISEL, Z, LDZ, WORK, LWORK,
*                          IWORK, LIWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDZ, LIWORK, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            ISEL( * ), IWORK( * )
*       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SSTEDS computes all eigenvalues of a real symmetric tridiagonal
*> matrix T = X * diag(LAMBDA) * X**T and the M eigenvectors
*> X(:,ISEL(1:M)).  Full eigenvector matrices are formed only at the
*> bounded-size leaves; above them the second pass propagates blocks
*> with M columns.  If M=N, the requested output is itself N-by-N.
*>
*> The first pass is the values-only boundary-row divide and conquer
*> method of SSTEBR: every subproblem passes only its eigenvalues and
*> the first and last rows of its eigenvector matrix to its parent,
*> and every merge builds its secular eigenvectors from Lowner-
*> corrected weights (SLAED8_BR, SLAED9_BR).  In addition each merge is
*> recorded in O(n_v) numbers (SLAES1), O(N log N) in total, and the
*> leaf eigenvector matrices are kept.
*>
*> The second pass (SLAES2) pushes the N-by-M selection matrix from the
*> root of the tree to the leaves; every intermediate is at most N-by-M.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is REAL array, dimension (N)
*>          On entry, the diagonal elements of the tridiagonal matrix.
*>          On exit, if INFO = 0, all eigenvalues in ascending order.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is REAL array, dimension (N-1)
*>          On entry, the subdiagonal elements of the tridiagonal
*>          matrix.  On exit, E has been destroyed.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of wanted eigenvectors.  0 <= M <= N.
*> \endverbatim
*>
*> \param[in] ISEL
*> \verbatim
*>          ISEL is INTEGER array, dimension (M)
*>          The indices (1 = smallest eigenvalue) of the wanted
*>          eigenvectors, strictly increasing, 1 <= ISEL(i) <= N.
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is REAL array, dimension (LDZ, M)
*>          On exit, if INFO = 0, Z(:,i) is the eigenvector belonging to
*>          the eigenvalue D(ISEL(i)).
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  With SMLSIZ = ILAENV(9,
*>          'SSTEDS',...), L the number of merge levels for N and
*>          NB = 64:
*>          LWORK >= N*SMLSIZ + 7*N*L + 19*N + 2*N*M + N*NB.
*>          For N <= 1, LWORK >= 1.
*>          If LWORK = -1, then a workspace query is assumed.
*>          A query does not reference D, E, ISEL, or Z.  If the
*>          required size cannot be represented by INTEGER, INFO=-9
*>          is returned, also for a query.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*> \endverbatim
*>
*> \param[in] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          The dimension of the array IWORK.
*>          LIWORK >= 3*N*L + 32*N + 3 + 2*M.
*>          For N <= 1, LIWORK >= 1.
*>          If LIWORK = -1, then a workspace query is assumed.
*>          Workspace sizes must fit INTEGER; see LWORK above.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          = 1:  An eigenvalue computation failed, but its submatrix
*>                location cannot be encoded in INTEGER.
*>          > 1:  The algorithm failed to compute an eigenvalue while
*>                working on a submatrix lying in rows and columns
*>                INFO/(N+1) through MOD(INFO,N+1).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup stebr
*
*> \par Contributors:
*  ==================
*>
*> Ruixin Jiang and Shaoshuai Zhang, University of Electronic Science
*> and Technology of China
*
*  =====================================================================
      SUBROUTINE SSTEDS( N, D, E, M, ISEL, Z, LDZ, WORK, LWORK,
     $                   IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine --
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, LIWORK, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEL( * ), IWORK( * )
      REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NB
      PARAMETER          ( NB = 64 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            B, BASE, DPTR, FINISH, G, I, IBCNT, IBL,
     $                   IBLKOF, IBNT, IBOFF, IBR, IBSIZ, IBST, IBTLV,
     $                   IC, IDREC, IDX, IINDXQ, IIREC, IIW0, ILC, ILT,
     $                   INTB, IP, IPAN, IPT, IPTR, IRUN, ITMP, IV, IVT,
     $                   IWW, J, KB, LDREC, LIREC, LIWMIN, LNT, LWMIN,
     $                   N1, N2, NBLOCK, NNEW, NRUN,
     $                   NTPTR, R, SMLSIZ, START, SZ, TLV, TLVN,
     $                   LIMIT, COEF, IMAX
      REAL               EPS, ORGNRM, TINY
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      REAL               SLAMCH, SLANST, SROUNDUP_LWORK
      EXTERNAL           ILAENV, SLAMCH, SLANST, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAES0, SLAES2, SLAMRG, SLASCL, SLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, HUGE, MAX, MOD, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.EQ.0 .AND. .NOT.LQUERY ) THEN
         DO 10 I = 1, M
            IF( ISEL( I ).LT.1 .OR. ISEL( I ).GT.N )
     $         INFO = -5
            IF( I.GT.1 ) THEN
               IF( ISEL( I ).LE.ISEL( I-1 ) )
     $            INFO = -5
            END IF
   10    CONTINUE
      END IF
*
*     Workspace.
*
      IF( INFO.EQ.0 ) THEN
         LWMIN = 1
         LIWMIN = 1
         LDREC = 0
         LIREC = 0
         LNT = 0
         SMLSIZ = 0
         IF( N.GT.1 ) THEN
            SMLSIZ = ILAENV( 9, 'SSTEDS', ' ', 0, 0, 0, 0 )
            TLVN = 0
            SZ = N
   20       CONTINUE
            IF( SZ.GT.SMLSIZ ) THEN
               SZ = SZ / 2 + MOD( SZ, 2 )
               TLVN = TLVN + 1
               GO TO 20
            END IF
*
*           Reject unrepresentable workspace sizes, including queries,
*           before any integer product can overflow.  After this test
*           the smaller record and offset expressions also fit.
*
            IMAX = HUGE( 0 )
            LIMIT = IMAX / N
            COEF = 7*TLVN + 19 + NB
            IF( SMLSIZ.GT.LIMIT ) THEN
               INFO = -9
               GO TO 25
            END IF
            IF( COEF.GT.LIMIT-SMLSIZ ) THEN
               INFO = -9
               GO TO 25
            END IF
            COEF = COEF + SMLSIZ
            IF( M.GT.( LIMIT-COEF ) / 2 ) THEN
               INFO = -9
               GO TO 25
            END IF
            LWMIN = N*( COEF+2*M )
            LDREC = N*SMLSIZ + 7*N*TLVN + 2*N
            LIREC = 3*N*TLVN + 2*N
            LNT = 14*N
            LIWMIN = N*( 3*TLVN+32 ) + 3 + 2*M
         END IF
         WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
         IWORK( 1 ) = LIWMIN
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -9
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         END IF
      END IF
   25 CONTINUE
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSTEDS', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
      IF( M.GT.0 )
     $   CALL SLASET( 'A', N, M, ZERO, ZERO, Z, LDZ )
      IF( N.EQ.1 ) THEN
         IF( M.EQ.1 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Real workspace layout.
*
      IDREC = 1
      IBL = IDREC + LDREC
      IBR = IBL + N
      IWW = IBR + N
      IV = IWW + 13*N
      IVT = IV + N
      IC = IVT + N
      ITMP = IC + N*M
      IPAN = ITMP + N*M
*
*     Integer workspace layout.
*
      IIREC = 1
      INTB = IIREC + LIREC
      IINDXQ = INTB + LNT
      IBLKOF = IINDXQ + N
      IBST = IBLKOF + N
      IBSIZ = IBST + N
      IBTLV = IBSIZ + N
      IBNT = IBTLV + N
      IBCNT = IBNT + N
      IBOFF = IBCNT + N
      IIW0 = IBOFF + N
      IP = IIW0 + 4*N + 2
      IPT = IP + N
      IDX = IPT + N
      IRUN = IDX + N
      ILT = IRUN + N + 1
      ILC = ILT + M
*
*     Zero matrix: identity eigenvectors.
*
      ORGNRM = SLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO ) THEN
         DO 30 I = 1, M
            Z( ISEL( I ), I ) = ONE
   30    CONTINUE
         RETURN
      END IF
*
*     First pass on each unreduced block (splitting test of SSTEDC).
*
      EPS = SLAMCH( 'Epsilon' )
      DPTR = 1
      IPTR = 1
      NTPTR = INTB
      NBLOCK = 0
      START = 1
   40 CONTINUE
      IF( START.LE.N ) THEN
         FINISH = START
   50    CONTINUE
         IF( FINISH.LT.N ) THEN
            TINY = EPS*SQRT( ABS( D( FINISH ) ) )*
     $             SQRT( ABS( D( FINISH+1 ) ) )
            IF( ABS( E( FINISH ) ).GT.TINY ) THEN
               FINISH = FINISH + 1
               GO TO 50
            END IF
         END IF
         SZ = FINISH - START + 1
         NBLOCK = NBLOCK + 1
         IWORK( IBST+NBLOCK-1 ) = START
         IWORK( IBSIZ+NBLOCK-1 ) = SZ
         IWORK( IBNT+NBLOCK-1 ) = NTPTR
         DO 60 I = START, FINISH
            IWORK( IBLKOF+I-1 ) = NBLOCK
   60    CONTINUE
         IF( SZ.EQ.1 ) THEN
            IWORK( IBTLV+NBLOCK-1 ) = 0
            IWORK( IINDXQ+START-1 ) = 1
         ELSE
            ORGNRM = SLANST( 'M', SZ, D( START ), E( START ) )
            CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, SZ, 1, D( START ), SZ,
     $                   INFO )
            CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, SZ-1, 1, E( START ),
     $                   SZ-1, INFO )
            CALL SLAES0( SZ, D( START ), E( START ), SMLSIZ, TLV,
     $                   IWORK( NTPTR ), 7, IWORK( IINDXQ+START-1 ),
     $                   WORK( IBL ), WORK( IBR ), WORK( IDREC ), LDREC,
     $                   DPTR, IWORK( IIREC ), LIREC, IPTR, WORK( IWW ),
     $                   IWORK( IIW0 ), INFO )
            IF( INFO.NE.0 ) THEN
               IF( INFO.GT.0 ) THEN
                  IF( INFO.EQ.1 ) RETURN
                  I = INFO / ( SZ+1 ) + START - 1
                  J = MOD( INFO, ( SZ+1 ) ) + START - 1
                  IF( I.GT.( HUGE( 0 )-J ) / ( N+1 ) ) THEN
                     INFO = 1
                  ELSE
                     INFO = I*( N+1 ) + J
                  END IF
               ELSE
                  INFO = -9
               END IF
               RETURN
            END IF
            CALL SLASCL( 'G', 0, 0, ONE, ORGNRM, SZ, 1, D( START ), SZ,
     $                   INFO )
            IWORK( IBTLV+NBLOCK-1 ) = TLV
            NTPTR = NTPTR + 7*( 2**( TLV+1 )-1 )
         END IF
         START = FINISH + 1
         GO TO 40
      END IF
*
*     Global ascending order: merge the sorted runs of the blocks.
*
      DO 80 B = 1, NBLOCK
         START = IWORK( IBST+B-1 )
         SZ = IWORK( IBSIZ+B-1 )
         IWORK( IRUN+B-1 ) = START
         DO 70 I = 1, SZ
            G = START - 1 + IWORK( IINDXQ+START+I-2 )
            WORK( IV+START+I-2 ) = D( G )
            IWORK( IP+START+I-2 ) = G
   70    CONTINUE
   80 CONTINUE
      NRUN = NBLOCK
      IWORK( IRUN+NRUN ) = N + 1
   90 CONTINUE
      IF( NRUN.GT.1 ) THEN
         NNEW = 0
         DO 120 R = 1, NRUN, 2
            NNEW = NNEW + 1
            IWORK( IRUN+NNEW-1 ) = IWORK( IRUN+R-1 )
            IF( R.LT.NRUN ) THEN
               BASE = IWORK( IRUN+R-1 )
               N1 = IWORK( IRUN+R ) - BASE
               N2 = IWORK( IRUN+R+1 ) - IWORK( IRUN+R )
               CALL SLAMRG( N1, N2, WORK( IV+BASE-1 ), 1, 1,
     $                      IWORK( IDX ) )
               DO 100 I = 1, N1 + N2
                  J = BASE - 1 + IWORK( IDX+I-1 )
                  WORK( IVT+I-1 ) = WORK( IV+J-1 )
                  IWORK( IPT+I-1 ) = IWORK( IP+J-1 )
  100          CONTINUE
               DO 110 I = 1, N1 + N2
                  WORK( IV+BASE+I-2 ) = WORK( IVT+I-1 )
                  IWORK( IP+BASE+I-2 ) = IWORK( IPT+I-1 )
  110          CONTINUE
            END IF
  120    CONTINUE
         NRUN = NNEW
         IWORK( IRUN+NRUN ) = N + 1
         GO TO 90
      END IF
      DO 130 I = 1, N
         D( I ) = WORK( IV+I-1 )
  130 CONTINUE
      IF( M.EQ.0 )
     $   GO TO 200
*
*     Distribute the wanted columns to their blocks.
*
      DO 140 B = 1, NBLOCK
         IWORK( IBCNT+B-1 ) = 0
  140 CONTINUE
      DO 150 I = 1, M
         B = IWORK( IBLKOF+IWORK( IP+ISEL( I )-1 )-1 )
         IWORK( IBCNT+B-1 ) = IWORK( IBCNT+B-1 ) + 1
  150 CONTINUE
      J = 0
      DO 160 B = 1, NBLOCK
         IWORK( IBOFF+B-1 ) = J
         J = J + IWORK( IBCNT+B-1 )
         IWORK( IBCNT+B-1 ) = 0
  160 CONTINUE
      DO 170 I = 1, M
         G = IWORK( IP+ISEL( I )-1 )
         B = IWORK( IBLKOF+G-1 )
         J = IWORK( IBOFF+B-1 ) + IWORK( IBCNT+B-1 )
         IWORK( ILT+J ) = I
         IWORK( ILC+J ) = G - IWORK( IBST+B-1 ) + 1
         IWORK( IBCNT+B-1 ) = IWORK( IBCNT+B-1 ) + 1
  170 CONTINUE
*
*     Second pass on each block: reconstruct every selected column.
*
      DO 180 B = 1, NBLOCK
         KB = IWORK( IBCNT+B-1 )
         IF( KB.EQ.0 )
     $      GO TO 180
         START = IWORK( IBST+B-1 )
         SZ = IWORK( IBSIZ+B-1 )
         J = IWORK( IBOFF+B-1 )
         IF( SZ.EQ.1 ) THEN
            Z( START, IWORK( ILT+J ) ) = ONE
         ELSE
            TLV = IWORK( IBTLV+B-1 )
            CALL SLAES2( SZ, TLV, IWORK( IWORK( IBNT+B-1 ) ), 7,
     $                   WORK( IDREC ), IWORK( IIREC ), KB,
     $                   IWORK( ILC+J ), IWORK( ILT+J ), Z( START, 1 ),
     $                   LDZ, WORK( IC ), WORK( ITMP ),
     $                   WORK( IPAN ), NB, INFO )
            IF( INFO.NE.0 ) THEN
               J = START + SZ - 1
               IF( START.GT.( HUGE( 0 )-J ) / ( N+1 ) ) THEN
                  INFO = 1
               ELSE
                  INFO = START*( N+1 ) + J
               END IF
               RETURN
            END IF
         END IF
  180 CONTINUE
*
  200 CONTINUE
      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of SSTEDS
*
      END
