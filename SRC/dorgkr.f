*> \brief \b DORGKR computes the explicit Q factor from DGEQRF and DLARFT
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE DORGKR(M, N, Q, LDQ)
*
*        .. Scalar Arguments ..
*        INTEGER           M, N, LDQ
*        ..
*        .. Array Arguments ..
*        DOUBLE PRECISION  Q(LDQ,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORGKR generates an m by n real matrix Q with orthonormal columns,
*> which is defined as the first n columns of the product of n
*> elementary reflectors
*>
*>       Q  =  I - V*T*V**T = H(1) H(2) . . . H(n)
*>
*> Where V is an m by n matrix whose columns are householder reflectors
*> as returned by DGEQRF and T is the n by n matrix returned by DLARFT
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
*>       Q is DOUBLE PRECISION array, dimension (LDQ,N)
*>       On entry, the upper triangular part and diagonal contains
*>       The array T as returned from DLARFT. In addition, the
*>       strictly lower triangular portion of the i-th column contains
*>       the vector which defines the elementary reflector H(i),
*>       for i = 1,2,...,n, as returned by DGEQRF
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
      SUBROUTINE DORGKR(M, N, Q, LDQ)
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           M, N, LDQ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION  NEG_ONE, ONE
      PARAMETER(NEG_ONE=-1.0D+0, ONE=1.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER           I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL          DTRMM, DTRTRM, DLUMM
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
*     V_1 \in \R^{n\times n}   assumed unit lower triangular
*     V_2 \in \R^{m-n\times n}
*
*     Compute T = T*V_1**T
*
      CALL DTRTRM('Right', 'Upper', 'Transpose', 'Non-unit', 'Unit',
     $            N, ONE, Q, LDQ, Q, LDQ)
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
         CALL DTRMM('Right', 'Upper', 'No Transpose', 'Non-unit',
     $               M-N, N, NEG_ONE, Q, LDQ, Q(N+1,1), LDQ)
      END IF
*
*     Q_1 = -V_1*T (Lower-Upper Matrix-Matrix multiplication)
*
      CALL DLUMM('Left', 'Unit', 'Non-Unit', N, NEG_ONE, Q, LDQ)
*
*     Q = "I" + Q
*
      J = MIN(M,N)
      DO I = 1, J
         Q(I,I) = Q(I,I) + ONE
      END DO
      END SUBROUTINE
