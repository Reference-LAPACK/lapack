*> \brief \b DLAES1 merges two subproblems of the boundary-row DC tree
*> and records the merge transform (used by DSTEDS).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAES1( N, CUTPNT, D, INDXQ, RHO, BLO, BHI, K, NGIV,
*                          IREC, LIREC, LIU, DREC, LDREC, LDU, WORK,
*                          IWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            CUTPNT, INFO, K, LDREC, LDU, LIREC, LIU, N,
*      $                   NGIV
*       DOUBLE PRECISION   RHO
*       ..
*       .. Array Arguments ..
*       INTEGER            INDXQ( * ), IREC( * ), IWORK( * )
*       DOUBLE PRECISION   BHI( * ), BLO( * ), D( * ), DREC( * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAES1 performs one merge of the values-only boundary-row divide
*> and conquer method of DSTEBR (it calls DLAED8_BR and DLAED9_BR in
*> the same way as DLAED7_BR does with two boundary rows) and, in
*> addition, records the merge transform W_v so that DSTEDS can later
*> apply it to selected eigenvector coefficients (see DLAES3):
*>
*>    X_v = diag( X_L, X_R ) * W_v,
*>    W_v = G_1 * ... * G_NGIV * PI * diag( Y, I_{N-K} ),
*>
*> where G_i are the deflation rotations of DLAED8_BR, PI(PERM(j),j)=1
*> and Y is the K-by-K secular eigenvector matrix built by DLAED9_BR
*> from the Lowner-corrected weights.
*>
*> Unlike DLAED7_BR, the root merge is also processed with two
*> boundary rows, so that the corrected weights and the root
*> representation (TAU, NEAR) are available for the root as well.
*>
*> The record is
*>    IREC:  PERM(N), GIVCOL(2,NGIV)
*>    DREC:  GIVNUM(2,NGIV), followed by
*>           K >= 3: DLAMBDA(K), W(K) (corrected), TAU(K), NEAR(2,K)
*>           K <= 2: DLAMBDA(K), W(K) (input to DLAED4), RHO
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*>          The order of the merged problem.
*> \param[in] CUTPNT
*>          The order of the left subproblem.
*> \param[in,out] D
*>          On entry the eigenvalues of the two subproblems, each in
*>          its internal order; on exit those of the merged problem in
*>          internal order (active roots, then deflated values).
*> \param[in,out] INDXQ
*>          On entry the sorting permutations of the two subproblems
*>          (local numbering); on exit the one of the merged problem.
*> \param[in] RHO
*>          The off-diagonal element at the cut.
*> \param[in,out] BLO
*>          First boundary rows of the children / of the merged node.
*> \param[in,out] BHI
*>          Last boundary rows of the children / of the merged node.
*> \param[out] K
*>          The number of non-deflated eigenvalues.
*> \param[out] NGIV
*>          The number of deflation rotations.
*> \param[out] IREC
*>          Integer record, dimension (LIREC), LIREC >= 3*N.
*> \param[in] LIREC
*> \param[out] LIU
*>          The number of entries of IREC used.
*> \param[out] DREC
*>          Real record, dimension (LDREC), LDREC >= 7*N.
*> \param[in] LDREC
*> \param[out] LDU
*>          The number of entries of DREC used.
*> \param[out] WORK
*>          Dimension (13*N).
*> \param[out] IWORK
*>          Dimension (2*N).
*> \param[out] INFO
*>          = 0: success; > 0: DLAED4 failed; -1: record too small.
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
      SUBROUTINE DLAES1( N, CUTPNT, D, INDXQ, RHO, BLO, BHI, K, NGIV,
     $                   IREC, LIREC, LIU, DREC, LDREC, LDU, WORK,
     $                   IWORK, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, INFO, K, LDREC, LDU, LIREC, LIU, N,
     $                   NGIV
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            INDXQ( * ), IREC( * ), IWORK( * )
      DOUBLE PRECISION   BHI( * ), BLO( * ), D( * ), DREC( * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            DS, I, IDEL, IDLM, INEAR, IQ, IQ2, ITAU, IW,
     $                   IWP, IWSGN, IZ, J
      DOUBLE PRECISION   RHOW
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAED8_BR, DLAED9_BR, DLAMRG
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      K = 0
      NGIV = 0
      LIU = 0
      LDU = 0
      IF( N.LE.0 )
     $   RETURN
      IF( LIREC.LT.3*N .OR. LDREC.LT.7*N ) THEN
         INFO = -1
         RETURN
      END IF
*
*     Workspace: Z, DLAMBDA, W (N each), Q, Q2 (2*N each),
*     DELTA, WPROD, WSGN, TAU (N each), NEAR (2*N).
*
      IZ = 1
      IDLM = IZ + N
      IW = IDLM + N
      IQ = IW + N
      IQ2 = IQ + 2*N
      IDEL = IQ2 + 2*N
      IWP = IDEL + N
      IWSGN = IWP + N
      ITAU = IWSGN + N
      INEAR = ITAU + N
*
*     z = [ last rows of left child ; first rows of right child ].
*     Q(1,:) = [ first rows of left ; 0 ], Q(2,:) = [ 0 ; last rows of
*     right ]: the two parent boundary rows in the children's bases.
*
      DO 10 J = 1, CUTPNT
         WORK( IZ+J-1 ) = BHI( J )
         WORK( IQ+2*( J-1 ) ) = BLO( J )
         WORK( IQ+2*( J-1 )+1 ) = ZERO
   10 CONTINUE
      DO 20 J = CUTPNT + 1, N
         WORK( IZ+J-1 ) = BLO( J )
         WORK( IQ+2*( J-1 ) ) = ZERO
         WORK( IQ+2*( J-1 )+1 ) = BHI( J )
   20 CONTINUE
*
*     Sort and deflate.  PERM -> IREC(1:N), GIVCOL -> IREC(N+1:3N),
*     GIVNUM -> DREC(1:2N); these are kept as the record.
*
      RHOW = RHO
      CALL DLAED8_BR( 2, K, N, D, WORK( IQ ), 2, INDXQ, RHOW, CUTPNT,
     $                WORK( IZ ), WORK( IDLM ), WORK( IQ2 ), 2,
     $                WORK( IW ), IREC, NGIV, IREC( N+1 ), DREC,
     $                IWORK( 1 ), IWORK( N+1 ), INFO )
      IF( INFO.NE.0 )
     $   RETURN
      DS = 2*NGIV
      LIU = N + 2*NGIV
      LDU = DS
*
      IF( K.GT.0 ) THEN
         IF( K.LE.2 ) THEN
*
*           DLAED4/DLAED5 return the eigenvectors directly; keep the
*           inputs so that DLAES3 can regenerate them.
*
            CALL DCOPY( K, WORK( IDLM ), 1, DREC( DS+1 ), 1 )
            CALL DCOPY( K, WORK( IW ), 1, DREC( DS+K+1 ), 1 )
            DREC( DS+2*K+1 ) = RHOW
            LDU = DS + 2*K + 1
         END IF
*
*        Secular roots, corrected weights and the two parent rows.
*
         CALL DLAED9_BR( K, 2, D, WORK( IQ2 ), 2, WORK( IQ ), 2, RHOW,
     $                   WORK( IDLM ), WORK( IW ), WORK( IDEL ),
     $                   WORK( IWP ), WORK( IWSGN ), WORK( ITAU ),
     $                   WORK( INEAR ), INFO )
         IF( INFO.NE.0 )
     $      RETURN
*
         IF( K.GE.3 ) THEN
            CALL DCOPY( K, WORK( IDLM ), 1, DREC( DS+1 ), 1 )
            CALL DCOPY( K, WORK( IW ), 1, DREC( DS+K+1 ), 1 )
            CALL DCOPY( K, WORK( ITAU ), 1, DREC( DS+2*K+1 ), 1 )
            CALL DCOPY( 2*K, WORK( INEAR ), 1, DREC( DS+3*K+1 ), 1 )
            LDU = DS + 5*K
         END IF
*
         CALL DLAMRG( K, N-K, D, 1, -1, INDXQ )
      ELSE
         DO 30 I = 1, N
            INDXQ( I ) = I
   30    CONTINUE
      END IF
*
*     Boundary rows of the merged node (internal order).
*
      DO 40 J = 1, N
         BLO( J ) = WORK( IQ+2*( J-1 ) )
         BHI( J ) = WORK( IQ+2*( J-1 )+1 )
   40 CONTINUE
      RETURN
*
*     End of DLAES1
*
      END
