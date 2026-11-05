*> \brief \b DLAES0 boundary-row divide and conquer for one unreduced
*> block, recording every merge (first pass of DSTEDS).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAES0( N, D, E, SMLSIZ, TLVLS, NT, LDNT, INDXQ, BLO,
*                          BHI, DREC, LDREC, DPTR, IREC, LIREC, IPTR,
*                          WORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            DPTR, INFO, IPTR, LDNT, LDREC, LIREC, N,
*      $                   SMLSIZ, TLVLS
*       ..
*       .. Array Arguments ..
*       INTEGER            INDXQ( * ), IREC( * ), IWORK( * ),
*      $                   NT( LDNT, * )
*       DOUBLE PRECISION   BHI( * ), BLO( * ), D( * ), DREC( * ),
*      $                   E( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAES0 computes all eigenvalues of an unreduced symmetric
*> tridiagonal block by the boundary-row divide and conquer method of
*> DSTEBR (same splitting as DLAED0_BR), and records in O(N log N)
*> storage what is needed to apply each merge transform later:
*>   leaves: the dense eigenvector matrix U_p (DSTEQR);
*>   merges: the record of DLAES1.
*>
*> The tree is a perfect binary tree with 2**TLVLS leaves in heap
*> order: node 1 is the root, 2h and 2h+1 are the children of h, the
*> leaves are 2**TLVLS, ..., 2**(TLVLS+1)-1 from top to bottom.
*> NT(:,h) describes node h:
*>   NT(1,h) first row (relative to the block)
*>   NT(2,h) size
*>   NT(3,h) size of the left child (0 for a leaf)
*>   NT(4,h) active dimension K
*>   NT(5,h) number of deflation rotations
*>   NT(6,h) start of its integer record in IREC
*>   NT(7,h) start of its real record in DREC (leaf: U_p)
*>
*> On exit D holds the eigenvalues in the internal order of the root
*> and D(INDXQ(i)) is the i-th smallest one; column INDXQ(i) of the
*> (never formed) root eigenvector matrix is its eigenvector.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N        Order of the block, N >= 2.
*> \param[in,out] D    Diagonal on entry, eigenvalues on exit.
*> \param[in,out] E    Off-diagonal, destroyed.
*> \param[in] SMLSIZ   Maximum leaf size.
*> \param[out] TLVLS   Number of merge levels.
*> \param[out] NT      Node table, dimension (LDNT, 2**(TLVLS+1)-1).
*> \param[in] LDNT     LDNT >= 7.
*> \param[out] INDXQ   Sorting permutation, dimension (N).
*> \param[out] BLO     Workspace, dimension (N).
*> \param[out] BHI     Workspace, dimension (N).
*> \param[out] DREC    Real record storage, dimension (LDREC).
*> \param[in] LDREC
*> \param[in,out] DPTR Next free entry of DREC.
*> \param[out] IREC    Integer record storage, dimension (LIREC).
*> \param[in] LIREC
*> \param[in,out] IPTR Next free entry of IREC.
*> \param[out] WORK    Dimension (13*N).
*> \param[out] IWORK   Dimension (4*N+2).
*> \param[out] INFO    = 0: success; > 0: an eigenvalue computation
*>                     failed (encoded as DLAED0); 1 if the location cannot
*>                     be encoded in INTEGER; -1: record full.
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
      SUBROUTINE DLAES0( N, D, E, SMLSIZ, TLVLS, NT, LDNT, INDXQ, BLO,
     $                   BHI, DREC, LDREC, DPTR, IREC, LIREC, IPTR,
     $                   WORK, IWORK, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            DPTR, INFO, IPTR, LDNT, LDREC, LIREC, N,
     $                   SMLSIZ, TLVLS
*     ..
*     .. Array Arguments ..
      INTEGER            INDXQ( * ), IREC( * ), IWORK( * ),
     $                   NT( LDNT, * )
      DOUBLE PRECISION   BHI( * ), BLO( * ), D( * ), DREC( * ),
     $                   E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            CURLVL, CURPRB, H, I, J, K, LDU, LIU, MATSIZ,
     $                   MSD2, NGIV, SMM1, SPM1, SPM2, SUBMAT, SUBPBS
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAES1, DSTEQR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, HUGE
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Partition as DLAED0 / DLAED0_BR do.
*
      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE
*
      SPM1 = SUBPBS - 1
      DO 40 I = 1, SPM1
         SUBMAT = IWORK( I ) + 1
         SMM1 = SUBMAT - 1
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) )
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) )
   40 CONTINUE
*
*     Leaves: keep U_p and its first and last rows.
*
      DO 60 I = 0, SPM1
         IF( I.EQ.0 ) THEN
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         ELSE
            SUBMAT = IWORK( I ) + 1
            MATSIZ = IWORK( I+1 ) - IWORK( I )
         END IF
         H = 2**TLVLS + I
         IF( DPTR-1+MATSIZ*MATSIZ.GT.LDREC ) THEN
            INFO = -1
            RETURN
         END IF
         CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ),
     $                DREC( DPTR ), MATSIZ, WORK, INFO )
         IF( INFO.NE.0 ) THEN
            IF( SUBMAT.GT.( HUGE( 0 )-( SUBMAT+MATSIZ-1 ) ) /
     $          ( N+1 ) ) THEN
               INFO = 1
            ELSE
               INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
            END IF
            RETURN
         END IF
         NT( 1, H ) = SUBMAT
         NT( 2, H ) = MATSIZ
         NT( 3, H ) = 0
         NT( 4, H ) = MATSIZ
         NT( 5, H ) = 0
         NT( 6, H ) = 0
         NT( 7, H ) = DPTR
         DO 50 J = 1, MATSIZ
            BLO( SUBMAT+J-1 ) = DREC( DPTR+( J-1 )*MATSIZ )
            BHI( SUBMAT+J-1 ) = DREC( DPTR+( J-1 )*MATSIZ+MATSIZ-1 )
            INDXQ( SUBMAT+J-1 ) = J
   50    CONTINUE
         DPTR = DPTR + MATSIZ*MATSIZ
   60 CONTINUE
*
*     Merge level by level, left to right (as DLAED0_BR).
*
      CURLVL = 1
      CURPRB = 0
   70 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         SPM2 = SUBPBS - 2
         DO 80 I = 0, SPM2, 2
            IF( I.EQ.0 ) THEN
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0
            ELSE
               SUBMAT = IWORK( I ) + 1
               MATSIZ = IWORK( I+2 ) - IWORK( I )
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1
            END IF
            H = 2**( TLVLS-CURLVL ) + CURPRB
            CALL DLAES1( MATSIZ, MSD2, D( SUBMAT ), INDXQ( SUBMAT ),
     $                   E( SUBMAT+MSD2-1 ), BLO( SUBMAT ),
     $                   BHI( SUBMAT ), K, NGIV, IREC( IPTR ),
     $                   LIREC-IPTR+1, LIU, DREC( DPTR ), LDREC-DPTR+1,
     $                   LDU, WORK, IWORK( 2*N+3 ), INFO )
            IF( INFO.NE.0 ) THEN
               IF( INFO.GT.0 ) THEN
                  IF( SUBMAT.GT.( HUGE( 0 )-( SUBMAT+MATSIZ-1 ) ) /
     $                ( N+1 ) ) THEN
                     INFO = 1
                  ELSE
                     INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
                  END IF
               END IF
               RETURN
            END IF
            NT( 1, H ) = SUBMAT
            NT( 2, H ) = MATSIZ
            NT( 3, H ) = MSD2
            NT( 4, H ) = K
            NT( 5, H ) = NGIV
            NT( 6, H ) = IPTR
            NT( 7, H ) = DPTR
            IPTR = IPTR + LIU
            DPTR = DPTR + LDU
            IWORK( I / 2+1 ) = IWORK( I+2 )
   80    CONTINUE
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 70
      END IF
      RETURN
*
*     End of DLAES0
*
      END
