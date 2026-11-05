*> \brief \b DLAES3 applies one recorded merge transform to a block of
*> coefficient columns (used by DSTEDS).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAES3( MODE, N, K, NGIV, IREC, DREC, NC, ICOL, C,
*                          LDC, TMP, LDT, PANEL, NB, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDC, LDT, MODE, N, NB, NC, NGIV
*       ..
*       .. Array Arguments ..
*       INTEGER            ICOL( * ), IREC( * )
*       DOUBLE PRECISION   C( LDC, * ), DREC( * ), PANEL( * ),
*      $                   TMP( LDT, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAES3 overwrites the N-by-NC block C by W_v * C, where W_v is the
*> merge transform recorded by DLAES1 (see there), without forming
*> W_v:
*>   1. TMP(1:K,:) = Y * C(1:K,:)
*>   2. C(PERM(j),:) = TMP(j,:) for j <= K and the old C(j,:) for j > K
*>      (deflated rows pass through the identity block)
*>   3. apply the deflation rotations in reverse order.
*> Steps 2 and 3 are done one column at a time, with the deflated rows
*> of the column staged in PANEL, so that the column stays in cache.
*> On exit rows 1:CUTPNT of C are coordinates in the left child basis
*> and the remaining rows in the right child basis.
*>
*> MODE = 0: C is dense; Y is regenerated in panels of NB columns
*>           (DLAES4) and applied with DGEMM.
*> MODE = 1: column t of C is the unit vector e_{ICOL(t)} (root); only
*>           the selected columns of Y are generated and written
*>           directly into C.  C need not be set on entry.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] MODE   0 (dense block) or 1 (one-hot columns).
*> \param[in] N      Size of the merge node.
*> \param[in] K      Active dimension, 0 <= K <= N.
*> \param[in] NGIV   Number of deflation rotations.
*> \param[in] IREC   Integer record of DLAES1.
*> \param[in] DREC   Real record of DLAES1.
*> \param[in] NC     Number of columns.
*> \param[in] ICOL   (MODE = 1) the one-hot row of each column.
*> \param[in,out] C  Dimension (LDC,NC).
*> \param[in] LDC
*> \param[out] TMP   Workspace, dimension (LDT,NC), LDT >= MAX(K,1).
*> \param[in] LDT
*> \param[out] PANEL Workspace, dimension (MAX(K*NB,N)).
*> \param[in] NB     Panel width.
*> \param[out] INFO  = 0: success; > 0: DLAED4 failed (K <= 2).
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
      SUBROUTINE DLAES3( MODE, N, K, NGIV, IREC, DREC, NC, ICOL, C,
     $                   LDC, TMP, LDT, PANEL, NB, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDC, LDT, MODE, N, NB, NC, NGIV
*     ..
*     .. Array Arguments ..
      INTEGER            ICOL( * ), IREC( * )
      DOUBLE PRECISION   C( LDC, * ), DREC( * ), PANEL( * ),
     $                   TMP( LDT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CI, DS, I, IA, IB, JB, JJ, JW, T
      DOUBLE PRECISION   BETA, CS, DLAM, SN, XA, XB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLAED4, DLAES4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( N.LE.0 .OR. NC.LE.0 )
     $   RETURN
      DS = 2*NGIV
*
*     For K <= 2 form Y explicitly in PANEL with DLAED4, exactly as
*     DLAED9_BR does (DLAED4/DLAED5 return normalized eigenvectors).
*
      IF( K.GE.1 .AND. K.LE.2 ) THEN
         DO 10 JJ = 1, K
            CALL DLAED4( K, JJ, DREC( DS+1 ), DREC( DS+K+1 ),
     $                   PANEL( ( JJ-1 )*K+1 ), DREC( DS+2*K+1 ),
     $                   DLAM, INFO )
            IF( INFO.NE.0 )
     $         RETURN
   10    CONTINUE
      END IF
*
*     One-hot columns (root): column t of diag( Y, I ) * S is the unit
*     vector or one column of Y; scatter it directly into C and replay
*     the rotations on that column while it is in cache.
*
      IF( MODE.EQ.1 ) THEN
         DO 40 T = 1, NC
            DO 20 I = 1, N
               C( I, T ) = ZERO
   20       CONTINUE
            CI = ICOL( T )
            IF( CI.GT.K ) THEN
               C( IREC( CI ), T ) = ONE
            ELSE
               IF( K.LE.2 ) THEN
                  DO 25 I = 1, K
                     TMP( I, T ) = PANEL( ( CI-1 )*K+I )
   25             CONTINUE
               ELSE
                  CALL DLAES4( K, CI, DREC( DS+1 ), DREC( DS+K+1 ),
     $                         DREC( DS+2*K+1 ), DREC( DS+3*K+1 ),
     $                         TMP( 1, T ) )
               END IF
               DO 30 I = 1, K
                  C( IREC( I ), T ) = TMP( I, T )
   30          CONTINUE
            END IF
            DO 35 I = NGIV, 1, -1
               IA = IREC( N+2*I-1 )
               IB = IREC( N+2*I )
               CS = DREC( 2*I-1 )
               SN = -DREC( 2*I )
               XA = C( IA, T )
               XB = C( IB, T )
               C( IA, T ) = CS*XA + SN*XB
               C( IB, T ) = CS*XB - SN*XA
   35       CONTINUE
   40    CONTINUE
         RETURN
      END IF
*
*     Step 1: TMP(1:K,:) = Y * C(1:K,:).
*
      IF( K.GE.1 .AND. K.LE.2 ) THEN
         CALL DGEMM( 'N', 'N', K, NC, K, ONE, PANEL, K, C, LDC, ZERO,
     $               TMP, LDT )
      ELSE IF( K.GE.3 ) THEN
         DO 60 JB = 1, K, NB
            JW = MIN( NB, K-JB+1 )
            DO 50 JJ = 1, JW
               CALL DLAES4( K, JB+JJ-1, DREC( DS+1 ), DREC( DS+K+1 ),
     $                      DREC( DS+2*K+1 ), DREC( DS+3*K+1 ),
     $                      PANEL( ( JJ-1 )*K+1 ) )
   50       CONTINUE
            IF( JB.EQ.1 ) THEN
               BETA = ZERO
            ELSE
               BETA = ONE
            END IF
            CALL DGEMM( 'N', 'N', K, NC, JW, ONE, PANEL, K, C( JB, 1 ),
     $                  LDC, BETA, TMP, LDT )
   60    CONTINUE
      END IF
*
*     Steps 2 and 3, one column at a time so that the column stays in
*     cache: copy the deflated rows (identity block), scatter all rows
*     back to the children's concatenated order (PERM), then replay the
*     deflation rotations in reverse creation order.  DLAED8_BR applied
*     DROT( Q(:,a), Q(:,b), c, s ) to the columns of the children's
*     basis; on coordinates this is a rotation with (c, -s).
*
      DO 90 T = 1, NC
*
*        PANEL(1:N) is free after the products above and is small
*        enough to stay in cache: take a copy of the column there so
*        that the scatter reads the deflated rows from cache instead of
*        going through the N-by-NC array TMP.
*
         DO 65 I = K + 1, N
            PANEL( I ) = C( I, T )
   65    CONTINUE
         DO 70 I = 1, K
            C( IREC( I ), T ) = TMP( I, T )
   70    CONTINUE
         DO 75 I = K + 1, N
            C( IREC( I ), T ) = PANEL( I )
   75    CONTINUE
         DO 80 I = NGIV, 1, -1
            IA = IREC( N+2*I-1 )
            IB = IREC( N+2*I )
            CS = DREC( 2*I-1 )
            SN = -DREC( 2*I )
            XA = C( IA, T )
            XB = C( IB, T )
            C( IA, T ) = CS*XA + SN*XB
            C( IB, T ) = CS*XB - SN*XA
   80    CONTINUE
   90 CONTINUE
      RETURN
*
*     End of DLAES3
*
      END
