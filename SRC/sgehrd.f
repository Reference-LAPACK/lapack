*> \brief \b SGEHRD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SGEHRD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgehrd.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgehrd.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgehrd.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       REAL              A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGEHRD reduces a real general matrix A to upper Hessenberg form H by
*> an orthogonal similarity transformation:  Q**T * A * Q = H .
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>          It is assumed that A is already upper triangular in rows
*>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*>          set by a previous call to SGEBAL; otherwise they should be
*>          set to 1 and N respectively. See Further Details.
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the N-by-N general matrix to be reduced.
*>          On exit, the upper triangle and the first subdiagonal of A
*>          are overwritten with the upper Hessenberg matrix H, and the
*>          elements below the first subdiagonal, with the array TAU,
*>          represent the orthogonal matrix Q as a product of elementary
*>          reflectors. See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is REAL array, dimension (N-1)
*>          The scalar factors of the elementary reflectors (see Further
*>          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
*>          zero.
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
*>          The length of the array WORK.  LWORK >= max(1,N).
*>          For good performance, LWORK should generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
*> \ingroup gehrd
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of (ihi-ilo) elementary
*>  reflectors
*>
*>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
*>  exit in A(i+2:ihi,i), and tau in TAU(i).
*>
*>  The contents of A are illustrated by the following example, with
*>  n = 7, ilo = 2 and ihi = 6:
*>
*>  on entry,                        on exit,
*>
*>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
*>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
*>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
*>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
*>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
*>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
*>  (                         a )    (                          a )
*>
*>  where a denotes an element of the original matrix A, h denotes a
*>  modified element of the upper Hessenberg matrix H, and vi denotes an
*>  element of the vector defining H(i).
*>
*>  This file is a slight modification of LAPACK-3.0's SGEHRD
*>  subroutine incorporating improvements proposed by Quintana-Orti and
*>  Van de Geijn (2006). (See SLAHR2.)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE SGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1,
     $                     TSIZE = LDT*NBMAX )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0,
     $                     ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB,
     $                   NBMIN, NH, NX
      REAL               EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SGEHD2, SGEMM, SLAHR2, SLARFB,
     $                   STRMM,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      REAL               SROUNDUP_LWORK
      EXTERNAL           ILAENV, SROUNDUP_LWORK
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      NH = IHI - ILO + 1
      IF( INFO.EQ.0 ) THEN
*
*       Compute the workspace requirements
*
         IF( NH.LE.1 ) THEN
            LWKOPT = 1
         ELSE
            NB = MIN( NBMAX, ILAENV( 1, 'SGEHRD', ' ', N, ILO, IHI,
     $                              -1 ) )
            LWKOPT = N*NB + TSIZE
         ENDIF
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEHRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
*
      DO 10 I = 1, ILO - 1
         TAU( I ) = ZERO
   10 CONTINUE
      DO 20 I = MAX( 1, IHI ), N - 1
         TAU( I ) = ZERO
   20 CONTINUE
*
*     Quick return if possible
*
      IF( NH.LE.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Determine the block size
*
      NB = MIN( NBMAX, ILAENV( 1, 'SGEHRD', ' ', N, ILO, IHI, -1 ) )
      NBMIN = 2
      IF( NB.GT.1 .AND. NB.LT.NH ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code)
*
         NX = MAX( NB, ILAENV( 3, 'SGEHRD', ' ', N, ILO, IHI, -1 ) )
         IF( NX.LT.NH ) THEN
*
*           Determine if workspace is large enough for blocked code
*
            IF( LWORK.LT.LWKOPT ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code
*
               NBMIN = MAX( 2, ILAENV( 2, 'SGEHRD', ' ', N, ILO, IHI,
     $                 -1 ) )
               IF( LWORK.GE.(N*NBMIN + TSIZE) ) THEN
                  NB = (LWORK-TSIZE) / N
               ELSE
                  NB = 1
               END IF
            END IF
         END IF
      END IF
      LDWORK = N
*
      IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
*
*        Use unblocked code below
*
         I = ILO
*
      ELSE
*
*        Use blocked code
*
         IWT = 1 + N*NB
         DO 40 I = ILO, IHI - 1 - NX, NB
            IB = MIN( NB, IHI-I )
*
*           Reduce columns i:i+ib-1 to Hessenberg form, returning the
*           matrices V and T of the block reflector H = I - V*T*V**T
*           which performs the reduction, and also the matrix Y = A*V*T
*
            CALL SLAHR2( IHI, I, IB, A( 1, I ), LDA, TAU( I ),
     $                   WORK( IWT ), LDT, WORK, LDWORK )
*
*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
*           right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set
*           to 1
*
            EI = A( I+IB, I+IB-1 )
            A( I+IB, I+IB-1 ) = ONE
            CALL SGEMM( 'No transpose', 'Transpose',
     $                  IHI, IHI-I-IB+1,
     $                  IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE,
     $                  A( 1, I+IB ), LDA )
            A( I+IB, I+IB-1 ) = EI
*
*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
*           right
*
            CALL STRMM( 'Right', 'Lower', 'Transpose',
     $                  'Unit', I, IB-1,
     $                  ONE, A( I+1, I ), LDA, WORK, LDWORK )
            DO 30 J = 0, IB-2
               CALL SAXPY( I, -ONE, WORK( LDWORK*J+1 ), 1,
     $                     A( 1, I+J+1 ), 1 )
   30       CONTINUE
*
*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
*           left
*
            CALL SLARFB( 'Left', 'Transpose', 'Forward',
     $                   'Columnwise',
     $                   IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA,
     $                   WORK( IWT ), LDT, A( I+1, I+IB ), LDA,
     $                   WORK, LDWORK )
   40    CONTINUE
      END IF
*
*     Use unblocked code to reduce the rest of the matrix
*
      CALL SGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
*
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
*
      RETURN
*
*     End of SGEHRD
*
      END
