*> \brief \b SORGKL computes the explicit Q factor from SGEQLF and SLARFT
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE SORGKL(M, N, Q, LDQ)
*
*        .. Scalar Arguments ..
*        INTEGER           M, N, LDQ
*        ..
*        .. Array Arguments ..
*        REAL              Q(LDQ,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SORGKL generates an m by n real matrix Q with orthonormal columns,
*> which is defined as the last n columns of the product of n
*> elementary reflectors
*>
*>       Q  =  I - V*T*V**T = H(n) . . . H(2) H(1)
*>
*> Where V is an m by n matrix whose columns are householder reflectors
*> as returned by SGEQLF and T is the n by n matrix returned by SLARFT
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
*>       Q is REAL array, dimension (LDQ,N)
*>       On entry, Q(1:m-n+i-1,i) contains the vector which defines the
*>       elementary reflector H(i), for i=1,...,n as returned by SGEQLF.
*>       In addition, the lower triangular portion of the submatrix given
*>       by Q(m-n+1:m,1:n) will contain the arry T as returned by SLARFT.
*>       See further details for more information.
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
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*> The storage of the V and T components inside Q is best illustrated by
*> the following example with m = 5, n = 3.
*>
*> Q =   |----------|
*>       | V1 V2 V3 |
*>       | V1 V2 V3 |
*>       | T1 V2 V3 |
*>       | T1 T2 V3 |
*>       | T1 T2 T3 |
*>       |----------|
*>
*> \endverbatim
*>
*  =====================================================================

      SUBROUTINE SORGKL(M, N, Q, LDQ)
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           M, N, LDQ
*     ..
*     .. Array Arguments ..
      REAL              Q(LDQ,*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL              NEG_ONE, ONE
      PARAMETER(NEG_ONE=-1.0E+0, ONE=1.0E+0)
*     ..
*     .. Local Scalars ..
      INTEGER           I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL          STRMM, STRTRM, SLUMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC         MIN
*     ..
*     .. Executable Statements ..
*
*     Break Q apart as follows
*
*           |---|
*     Q =   | V |
*           | T |
*           |---|
*
*     Where T is an n-by-n lower triangular matrix, and V is as described
*     in the Further Details section
*
*     In turn, break apart V as follows
*
*           |-----|
*     V =   | V_2 |
*           | V_1 |
*           |-----|
*
*     Where:
*
*     V_1 \in \R^{n\times n}   assumed unit upper triangular
*     V_2 \in \R^{m-n\times n}
*
*     Compute T = T*V_1**T
*
      CALL STRTRM('Right', 'Lower', 'Transpose', 'Non-Unit', 'Unit',
     $         N, ONE, Q(M-N+1,1), LDQ, Q(M-N+1,1), LDQ)
*
*     Compute Q = -VT. This means that we need to break apart
*     Our computation in two parts
*
*           |--------|
*     Q =   | -V_2*T |
*           | -V_1*T |
*           |--------|
*
*     Q_2 = -V_2*T (TRMM) but only when necessary
*
      IF (M.GT.N) THEN
         CALL STRMM('Right', 'Lower', 'No Transpose', 'Non-Unit',
     $            M-N, N, NEG_ONE, Q(M-N+1,1), LDQ, Q, LDQ)
      END IF
*
*     Q_1 = -V_1*T (Lower-Upper Matrix-Matrix multiplication)
*
      CALL SLUMM('Right', 'Non-Unit', 'Unit', N, NEG_ONE,
     $         Q(M-N+1,1), LDQ)
*
*     Q = "I" + Q
*
      J = MIN(M,N)
      DO I = 1, J
         Q(M-N+I,I) = Q(M-N+I,I) + ONE
      END DO
      END SUBROUTINE
