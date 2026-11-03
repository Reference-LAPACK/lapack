*> \brief \b DLAES4 regenerates one secular eigenvector from a DLAES1
*> merge record (used by DSTEDS).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAES4( K, J, DLAMBDA, W, TAU, NEAR, S )
*
*       .. Scalar Arguments ..
*       INTEGER            J, K
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   DLAMBDA( * ), NEAR( 2, * ), S( * ), TAU( * ),
*      $                   W( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAES4 forms column J of the K-by-K secular eigenvector matrix of a
*> merge recorded by DLAES1 (K >= 3), exactly as DLAED9_BR forms it to
*> update the boundary rows:
*>
*>    S(i) = W(i) / delta(i),   S = S / || S ||_2,
*>
*> with W the Lowner-corrected weights and delta(i) = DLAMBDA(i) -
*> lambda_J in the pole-relative representation of DLAED4_BR:
*> lambda_J = DLAMBDA(IORG) + TAU(J), where IORG = J if TAU(J) >= 0,
*> IORG = J+1 otherwise, IORG = K for J = K; the two pole-adjacent
*> differences are taken from NEAR(1:2,J).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] K       Active dimension, K >= 3.
*> \param[in] J       Column index, 1 <= J <= K.
*> \param[in] DLAMBDA Active poles, dimension (K).
*> \param[in] W       Corrected weights, dimension (K).
*> \param[in] TAU     Root offsets, dimension (K).
*> \param[in] NEAR    Pole-adjacent differences, dimension (2,K).
*> \param[out] S      Column J, dimension (K).
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
      SUBROUTINE DLAES4( K, J, DLAMBDA, W, TAU, NEAR, S )
*
*     .. Scalar Arguments ..
      INTEGER            J, K
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   DLAMBDA( * ), NEAR( 2, * ), S( * ), TAU( * ),
     $                   W( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IORG, N1, N2
      DOUBLE PRECISION   DLT, TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
      IF( J.EQ.K ) THEN
         IORG = K
         N1 = K - 1
         N2 = K
      ELSE IF( TAU( J ).GE.ZERO ) THEN
         IORG = J
         N1 = J
         N2 = J + 1
      ELSE
         IORG = J + 1
         N1 = J
         N2 = J + 1
      END IF
      TEMP = ZERO
      DO 10 I = 1, K
         IF( I.EQ.N1 ) THEN
            DLT = W( I ) / NEAR( 1, J )
         ELSE IF( I.EQ.N2 ) THEN
            DLT = W( I ) / NEAR( 2, J )
         ELSE
            DLT = W( I ) /
     $            ( ( DLAMBDA( I )-DLAMBDA( IORG ) ) - TAU( J ) )
         END IF
         TEMP = TEMP + DLT*DLT
         S( I ) = DLT
   10 CONTINUE
      TEMP = SQRT( TEMP )
      DO 20 I = 1, K
         S( I ) = S( I ) / TEMP
   20 CONTINUE
      RETURN
*
*     End of DLAES4
*
      END
