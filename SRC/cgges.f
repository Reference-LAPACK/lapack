*> \brief <b> CGGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CGGES + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgges.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgges.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgges.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB,
*                         SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK,
*                         LWORK, RWORK, BWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVSL, JOBVSR, SORT
*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
*       ..
*       .. Array Arguments ..
*       LOGICAL            BWORK( * )
*       REAL               RWORK( * )
*       COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ),
*      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
*      $                   WORK( * )
*       ..
*       .. Function Arguments ..
*       LOGICAL            SELCTG
*       EXTERNAL           SELCTG
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CGGES computes for a pair of N-by-N complex nonsymmetric matrices
*> (A,B), the generalized eigenvalues, the generalized complex Schur
*> form (S, T), and optionally left and/or right Schur vectors (VSL
*> and VSR). This gives the generalized Schur factorization
*>
*>         (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
*>
*> where (VSR)**H is the conjugate-transpose of VSR.
*>
*> Optionally, it also orders the eigenvalues so that a selected cluster
*> of eigenvalues appears in the leading diagonal blocks of the upper
*> triangular matrix S and the upper triangular matrix T. The leading
*> columns of VSL and VSR then form an unitary basis for the
*> corresponding left and right eigenspaces (deflating subspaces).
*>
*> (If only the generalized eigenvalues are needed, use the driver
*> CGGEV instead, which is faster.)
*>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
*> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
*> usually represented as the pair (alpha,beta), as there is a
*> reasonable interpretation for beta=0, and even for both being zero.
*>
*> A pair of matrices (S,T) is in generalized complex Schur form if S
*> and T are upper triangular and, in addition, the diagonal elements
*> of T are non-negative real numbers.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBVSL
*> \verbatim
*>          JOBVSL is CHARACTER*1
*>          = 'N':  do not compute the left Schur vectors;
*>          = 'V':  compute the left Schur vectors.
*> \endverbatim
*>
*> \param[in] JOBVSR
*> \verbatim
*>          JOBVSR is CHARACTER*1
*>          = 'N':  do not compute the right Schur vectors;
*>          = 'V':  compute the right Schur vectors.
*> \endverbatim
*>
*> \param[in] SORT
*> \verbatim
*>          SORT is CHARACTER*1
*>          Specifies whether or not to order the eigenvalues on the
*>          diagonal of the generalized Schur form.
*>          = 'N':  Eigenvalues are not ordered;
*>          = 'S':  Eigenvalues are ordered (see SELCTG).
*> \endverbatim
*>
*> \param[in] SELCTG
*> \verbatim
*>          SELCTG is a LOGICAL FUNCTION of two COMPLEX arguments
*>          SELCTG must be declared EXTERNAL in the calling subroutine.
*>          If SORT = 'N', SELCTG is not referenced.
*>          If SORT = 'S', SELCTG is used to select eigenvalues to sort
*>          to the top left of the Schur form.
*>          An eigenvalue ALPHA(j)/BETA(j) is selected if
*>          SELCTG(ALPHA(j),BETA(j)) is true.
*>
*>          Note that a selected complex eigenvalue may no longer satisfy
*>          SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since
*>          ordering may change the value of complex eigenvalues
*>          (especially if the eigenvalue is ill-conditioned), in this
*>          case INFO is set to N+2 (See INFO below).
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VSL, and VSR.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA, N)
*>          On entry, the first of the pair of matrices.
*>          On exit, A has been overwritten by its generalized Schur
*>          form S.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB, N)
*>          On entry, the second of the pair of matrices.
*>          On exit, B has been overwritten by its generalized Schur
*>          form T.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] SDIM
*> \verbatim
*>          SDIM is INTEGER
*>          If SORT = 'N', SDIM = 0.
*>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
*>          for which SELCTG is true.
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is COMPLEX array, dimension (N)
*>          On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the
*>          generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j),
*>          j=1,...,N  are the diagonals of the complex Schur form (A,B)
*>          output by CGGES. The  BETA(j) will be non-negative real.
*>
*>          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
*>          underflow, and BETA(j) may even be zero.  Thus, the user
*>          should avoid naively computing the ratio alpha/beta.
*>          However, ALPHA will be always less than and usually
*>          comparable with norm(A) in magnitude, and BETA always less
*>          than and usually comparable with norm(B).
*> \endverbatim
*>
*> \param[out] VSL
*> \verbatim
*>          VSL is COMPLEX array, dimension (LDVSL,N)
*>          If JOBVSL = 'V', VSL will contain the left Schur vectors.
*>          Not referenced if JOBVSL = 'N'.
*> \endverbatim
*>
*> \param[in] LDVSL
*> \verbatim
*>          LDVSL is INTEGER
*>          The leading dimension of the matrix VSL. LDVSL >= 1, and
*>          if JOBVSL = 'V', LDVSL >= N.
*> \endverbatim
*>
*> \param[out] VSR
*> \verbatim
*>          VSR is COMPLEX array, dimension (LDVSR,N)
*>          If JOBVSR = 'V', VSR will contain the right Schur vectors.
*>          Not referenced if JOBVSR = 'N'.
*> \endverbatim
*>
*> \param[in] LDVSR
*> \verbatim
*>          LDVSR is INTEGER
*>          The leading dimension of the matrix VSR. LDVSR >= 1, and
*>          if JOBVSR = 'V', LDVSR >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,2*N).
*>          For good performance, LWORK must generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (8*N)
*> \endverbatim
*>
*> \param[out] BWORK
*> \verbatim
*>          BWORK is LOGICAL array, dimension (N)
*>          Not referenced if SORT = 'N'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          =1,...,N:
*>                The QZ iteration failed.  (A,B) are not in Schur
*>                form, but ALPHA(j) and BETA(j) should be correct for
*>                j=INFO+1,...,N.
*>          > N:  =N+1: other than QZ iteration failed in CHGEQZ
*>                =N+2: after reordering, roundoff changed values of
*>                      some complex eigenvalues so that leading
*>                      eigenvalues in the Generalized Schur form no
*>                      longer satisfy SELCTG=.TRUE.  This could also
*>                      be caused due to scaling.
*>                =N+3: reordering failed in CTGSEN.
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
*> \ingroup gges
*
*  =====================================================================
      SUBROUTINE CGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B,
     $                  LDB,
     $                  SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK,
     $                  LWORK, RWORK, BWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVSL, JOBVSR, SORT
      INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
     $                   WORK( * )
*     ..
*     .. Function Arguments ..
      LOGICAL            SELCTG
      EXTERNAL           SELCTG
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ),
     $                   CONE = ( 1.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL,
     $                   LQUERY, WANTST
      INTEGER            I, ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT,
     $                   ILO, IRIGHT, IROWS, IRWRK, ITAU, IWRK, LWKMIN,
     $                   LWKOPT
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PVSL,
     $                   PVSR, SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      REAL               DIF( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEQRF, CGGBAK, CGGBAL, CGGHRD, CHGEQZ,
     $                   CLACPY,
     $                   CLASCL, CLASET, CTGSEN, CUNGQR, CUNMQR, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               CLANGE, SLAMCH, SROUNDUP_LWORK
      EXTERNAL           LSAME, ILAENV, CLANGE, SLAMCH,
     $                   SROUNDUP_LWORK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVSL, 'N' ) ) THEN
         IJOBVL = 1
         ILVSL = .FALSE.
      ELSE IF( LSAME( JOBVSL, 'V' ) ) THEN
         IJOBVL = 2
         ILVSL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVSL = .FALSE.
      END IF
*
      IF( LSAME( JOBVSR, 'N' ) ) THEN
         IJOBVR = 1
         ILVSR = .FALSE.
      ELSE IF( LSAME( JOBVSR, 'V' ) ) THEN
         IJOBVR = 2
         ILVSR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVSR = .FALSE.
      END IF
*
      WANTST = LSAME( SORT, 'S' )
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTST ) .AND.
     $         ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDVSL.LT.1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) THEN
         INFO = -16
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.)
*
      IF( INFO.EQ.0 ) THEN
         LWKMIN = MAX( 1, 2*N )
         LWKOPT = MAX( 1, N + N*ILAENV( 1, 'CGEQRF', ' ', N, 1, N,
     $                 0 ) )
         LWKOPT = MAX( LWKOPT, N +
     $                 N*ILAENV( 1, 'CUNMQR', ' ', N, 1, N, -1 ) )
         IF( ILVSL ) THEN
            LWKOPT = MAX( LWKOPT, N +
     $                    N*ILAENV( 1, 'CUNGQR', ' ', N, 1, N, -1 ) )
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY )
     $      INFO = -18
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGGES ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
*
      IF( ILASCL )
     $   CALL CLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      BNRM = CLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
*
      IF( ILBSCL )
     $   CALL CLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
*
*     Permute the matrix to make it more nearly triangular
*     (Real Workspace: need 6*N)
*
      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      CALL CGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ),
     $             RWORK( IRIGHT ), RWORK( IRWRK ), IERR )
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Complex Workspace: need N, prefer N*NB)
*
      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL CGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWRK ), LWORK+1-IWRK, IERR )
*
*     Apply the orthogonal transformation to matrix A
*     (Complex Workspace: need N, prefer N*NB)
*
      CALL CUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ),
     $             LWORK+1-IWRK, IERR )
*
*     Initialize VSL
*     (Complex Workspace: need N, prefer N*NB)
*
      IF( ILVSL ) THEN
         CALL CLASET( 'Full', N, N, CZERO, CONE, VSL, LDVSL )
         IF( IROWS.GT.1 ) THEN
            CALL CLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                   VSL( ILO+1, ILO ), LDVSL )
         END IF
         CALL CUNGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL,
     $                WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
*     Initialize VSR
*
      IF( ILVSR )
     $   CALL CLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR )
*
*     Reduce to generalized Hessenberg form
*     (Workspace: none needed)
*
      CALL CGGHRD( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL,
     $             LDVSL, VSR, LDVSR, IERR )
*
      SDIM = 0
*
*     Perform QZ algorithm, computing Schur vectors if desired
*     (Complex Workspace: need N)
*     (Real Workspace: need N)
*
      IWRK = ITAU
      CALL CHGEQZ( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ),
     $             LWORK+1-IWRK, RWORK( IRWRK ), IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 30
      END IF
*
*     Sort eigenvalues ALPHA/BETA if desired
*     (Workspace: none needed)
*
      IF( WANTST ) THEN
*
*        Undo scaling on eigenvalues before selecting
*
         IF( ILASCL )
     $      CALL CLASCL( 'G', 0, 0, ANRM, ANRMTO, N, 1, ALPHA, N,
     $                   IERR )
         IF( ILBSCL )
     $      CALL CLASCL( 'G', 0, 0, BNRM, BNRMTO, N, 1, BETA, N,
     $                   IERR )
*
*        Select eigenvalues
*
         DO 10 I = 1, N
            BWORK( I ) = SELCTG( ALPHA( I ), BETA( I ) )
   10    CONTINUE
*
         CALL CTGSEN( 0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB,
     $                ALPHA,
     $                BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR,
     $                DIF, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, IERR )
         IF( IERR.EQ.1 )
     $      INFO = N + 3
*
      END IF
*
*     Apply back-permutation to VSL and VSR
*     (Workspace: none needed)
*
      IF( ILVSL )
     $   CALL CGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ),
     $                RWORK( IRIGHT ), N, VSL, LDVSL, IERR )
      IF( ILVSR )
     $   CALL CGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ),
     $                RWORK( IRIGHT ), N, VSR, LDVSR, IERR )
*
*     Undo scaling
*
      IF( ILASCL ) THEN
         CALL CLASCL( 'U', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR )
         CALL CLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
      END IF
*
      IF( ILBSCL ) THEN
         CALL CLASCL( 'U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR )
         CALL CLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      END IF
*
      IF( WANTST ) THEN
*
*        Check if reordering is correct
*
         LASTSL = .TRUE.
         SDIM = 0
         DO 20 I = 1, N
            CURSL = SELCTG( ALPHA( I ), BETA( I ) )
            IF( CURSL )
     $         SDIM = SDIM + 1
            IF( CURSL .AND. .NOT.LASTSL )
     $         INFO = N + 2
            LASTSL = CURSL
   20    CONTINUE
*
      END IF
*
   30 CONTINUE
*
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
      RETURN
*
*     End of CGGES
*
      END
