*> \brief \b DLAES2 top-down reconstruction of selected eigenvectors
*> of one unreduced block (second pass of DSTEDS).
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAES2( N, TLVLS, NT, LDNT, DREC, IREC, KB, ICOL,
*                          OUTCOL, Z, LDZ, C, TMP, PANEL, NB, INFO )
*
*> \par Purpose:
*> \verbatim
*> With X in root internal order and S = I(:,ICOL), evaluate
*>
*>    X*S = U * W_{L-1} * ... * W_1 * ( W_0 * S )
*>
*> from right to left.  At a merge node v, DLAES3 overwrites its
*> n_v-by-KB coefficient block by W_v*C_v.  The two row blocks pass
*> to its children.  At a leaf p, write U_p*C_p to the output rows.
*> The root selection is implicit; no full intermediate eigenvector
*> matrix is formed unless KB=N is itself requested.
*> \endverbatim
*
*> \param[in] N, TLVLS, NT, LDNT, DREC, IREC
*>          Records returned by DLAES0.
*> \param[in] KB
*>          Number of selected eigenvectors of this block.
*> \param[in] ICOL
*>          Root internal indices, dimension (KB).
*> \param[in] OUTCOL
*>          Destination columns of Z, dimension (KB).
*> \param[out] Z
*>          Output rows of this block, dimension (LDZ,*).
*> \param[in] LDZ
*>          Leading dimension of Z.
*> \param[out] C, TMP
*>          Workspace, dimension (N,KB) each.
*> \param[out] PANEL
*>          Workspace, dimension (N*NB).
*> \param[in] NB
*>          Secular-basis panel width.
*> \param[out] INFO
*>          Zero on success; positive if a small secular solve failed.
*>
*> \ingroup stebr
*
*> \par Contributors:
*  ==================
*>
*> Ruixin Jiang and Shaoshuai Zhang, University of Electronic Science
*> and Technology of China
*
      SUBROUTINE DLAES2( N, TLVLS, NT, LDNT, DREC, IREC, KB, ICOL,
     $                   OUTCOL, Z, LDZ, C, TMP, PANEL, NB, INFO )
      IMPLICIT NONE
      INTEGER            INFO, KB, LDNT, LDZ, N, NB, TLVLS
      INTEGER            ICOL( * ), IREC( * ), NT( LDNT, * ),
     $                   OUTCOL( * )
      DOUBLE PRECISION   C( N, * ), DREC( * ), PANEL( * ),
     $                   TMP( N, * ), Z( LDZ, * )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      LOGICAL            CONTIG
      INTEGER            H, I, NLEAF, NP, R0, T
      EXTERNAL           DGEMM, DLAES3
*
      INFO = 0
      IF( KB.LE.0 ) RETURN
      NLEAF = 2**TLVLS
*
*     A single leaf: copy selected columns of its checkpoint.
*
      IF( TLVLS.EQ.0 ) THEN
         NP = NT( 2, 1 )
         DO 20 T = 1, KB
            DO 10 I = 1, NP
               Z( I, OUTCOL( T ) ) = DREC( NT( 7, 1 )+
     $                               ( ICOL( T )-1 )*NP+I-1 )
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Root: generate only the columns selected by the implicit S.
*
      CALL DLAES3( 1, N, NT( 4, 1 ), NT( 5, 1 ), IREC( NT( 6, 1 ) ),
     $             DREC( NT( 7, 1 ) ), KB, ICOL, C, N, TMP, N,
     $             PANEL, NB, INFO )
      IF( INFO.NE.0 ) RETURN
*
*     Heap order visits every parent before its children. All KB
*     columns are propagated at each node, with a shared row stride N.
*
      DO 30 H = 2, NLEAF - 1
         R0 = NT( 1, H )
         CALL DLAES3( 0, NT( 2, H ), NT( 4, H ), NT( 5, H ),
     $                IREC( NT( 6, H ) ), DREC( NT( 7, H ) ), KB, ICOL,
     $                C( R0, 1 ), N, TMP, N, PANEL, NB, INFO )
         IF( INFO.NE.0 ) RETURN
   30 CONTINUE
*
*     Leaves. Consecutive destinations can be written by DGEMM
*     directly; interleaved destinations are scattered from TMP.
*
      CONTIG = .TRUE.
      DO 40 T = 2, KB
         IF( OUTCOL( T ).NE.OUTCOL( 1 )+T-1 ) CONTIG = .FALSE.
   40 CONTINUE
      DO 70 H = NLEAF, 2*NLEAF - 1
         R0 = NT( 1, H )
         NP = NT( 2, H )
         IF( CONTIG ) THEN
            CALL DGEMM( 'N', 'N', NP, KB, NP, ONE,
     $                  DREC( NT( 7, H ) ), NP, C( R0, 1 ), N,
     $                  ZERO, Z( R0, OUTCOL( 1 ) ), LDZ )
         ELSE
            CALL DGEMM( 'N', 'N', NP, KB, NP, ONE,
     $                  DREC( NT( 7, H ) ), NP, C( R0, 1 ), N,
     $                  ZERO, TMP, N )
            DO 60 T = 1, KB
               DO 50 I = 1, NP
                  Z( R0+I-1, OUTCOL( T ) ) = TMP( I, T )
   50          CONTINUE
   60       CONTINUE
         END IF
   70 CONTINUE
      RETURN
      END
