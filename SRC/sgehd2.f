*> \brief \b SGEHD2 reduces a general square matrix to upper Hessenberg form using an unblocked algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SGEHD2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgehd2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgehd2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgehd2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, ILO, INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGEHD2 reduces a real general matrix A to upper Hessenberg form H by
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
*>          1 <= ILO <= IHI <= max(1,N).
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the n by n general matrix to be reduced.
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
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
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
*> \ingroup gehd2
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
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE SGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARF1F, SLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEHD2', -INFO )
         RETURN
      END IF
*
      DO 10 I = ILO, IHI - 1
*
*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
*
         CALL SLARFG( IHI-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                TAU( I ) )
*
*        Apply H(i) to A(1:ihi,i+1:ihi) from the right
*
         CALL SLARF1F( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ),
     $                 A( 1, I+1 ), LDA, WORK )
*
*        Apply H(i) to A(i+1:ihi,i+1:n) from the left
*
         CALL SLARF1F( 'Left', IHI-I, N-I, A( I+1, I ), 1, TAU( I ),
     $                 A( I+1, I+1 ), LDA, WORK )
*
   10 CONTINUE
*
      RETURN
*
*     End of SGEHD2
*
      END
