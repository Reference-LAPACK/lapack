*> \brief \b CUNGKR computes the explicit Q factor from CGEQRF and CLARFT
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE CUNGKR(M, N, Q, LDQ)
*
*        .. Scalar Arguments ..
*        INTEGER           M, N, LDQ
*        ..
*        .. Array Arguments ..
*        COMPLEX           Q(LDQ,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CUNGKR generates an m by n complex matrix Q with orthonormal columns,
*> which is defined as the first n columns of the product of n
*> elementary reflectors
*>
*>       Q  =  I - V*T*V**H = H(1) H(2) . . . H(n)
*>
*> Where V is an m by n matrix whose columns are householder reflectors
*> as returned by CGEQRF and T is the n by n matrix returned by CLARFT
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix V. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix V, and the order of T.
*>          N >= 0.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>       Q is COMPLEX array, dimension (LDQ,N)
*>       On entry, the upper triangular part and diagonal contains
*>       The array T as returned from CLARFT. In addition, the
*>       strictly lower triangular portion of the i-th column contains
*>       the vector which defines the elementary reflector H(i),
*>       for i = 1,2,...,n, as returned by CGEQRF
*>       On exit, the m-by-n matrix Q.
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
*  =====================================================================
*     Cost: (2mn**2 + n**2 - n)/2
      SUBROUTINE CUNGKR(M, N, Q, LDQ)
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           M, N, LDQ
*     ..
*     .. Array Arguments ..
      COMPLEX           Q(LDQ,*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX           NEG_ONE, ONE
      PARAMETER(NEG_ONE=(-1.0E+0,0.0E+0), ONE=(1.0E+0,0.0E+0))
*     ..
*     .. Local Scalars ..
      INTEGER           I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL          CTRMM, CTRTRM, CLUMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC         MIN
*     ..
*     .. Executable Statements ..
*
*     Break Q apart as follows
*
*           |---|
*     Q =   | T |
*           | V |
*           |---|
*
*     Where T is an n-by-n upper triangular matrix, and V is an
*     m-by-n assumed unit lower trapezoidal matrix
*
*     In turn, break apart V as follows
*
*           |-----|
*     V =   | V_1 |
*           | V_2 |
*           |-----|
*
*     Where:
*
*     V_1 \in \C^{n\times n}   assumed unit lower triangular
*     V_2 \in \C^{m-n\times n}
*
*     Compute T = T*V_1**H
*
      CALL CTRTRM('Right', 'Upper', 'Conjugate Transpose', 
     $            'Non-unit', 'Unit', N, ONE, Q, LDQ, Q, LDQ)
*
*     Compute Q = -VT. This means that we need to break apart
*     Our computation in two parts
*
*           |--------|
*     Q =   | -V_1*T |
*           | -V_2*T |
*           |--------|
*
*     Q_2 = -V_2*T (TRMM) but only when necessary
*
      IF (M.GT.N) THEN
         CALL CTRMM('Right', 'Upper', 'No Transpose', 'Non-unit',
     $               M-N, N, NEG_ONE, Q, LDQ, Q(N+1,1), LDQ)
      END IF
*
*     Q_1 = -V_1*T (Lower-Upper Matrix-Matrix multiplication)
*
      CALL CLUMM('Left', 'Unit', 'Non-Unit', N, NEG_ONE, Q, LDQ)
*
*     Q = "I" + Q
*
      J = MIN(M,N)
      DO I = 1, J
         Q(I,I) = Q(I,I) + ONE
      END DO
      END SUBROUTINE
