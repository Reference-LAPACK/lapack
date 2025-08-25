*> \brief \b SORGRK computes the explicit Q factor from SGERQF and SLARFT
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE SORGRK(M, N, Q, LDQ)
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
*> SORGRK generates an m by n real matrix Q with orthonormal rows,
*> which is defined as the last m rows of the product of m
*> elementary reflectors
*>
*>       Q  =  I - V'*T*V = H(m) . . . H(2) H(1)
*>
*> Where V is an m by n matrix whose columns are householder reflectors
*> as returned by SGERQF and T is the n by n matrix returned by SLARFT
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix V, and the order of T.
*>          M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix V. N >= 0.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>       Q is REAL array, dimension (LDQ,N)
*>       On entry, Q(i,1:n-m-1+i) contains the vector which defines the
*>       elementary reflector H(i), for i=1,...,n as returned by SGERQF.
*>       In addition, the upper triangular portion of the submatrix given
*>       by Q(1:m,n-m:n) will contain the array T as returned by SLARFT.
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
*> the following example with m = 3, n = 5.
*>
*> Q =   |----------------|
*>       | V1 V1 T1 T1 T1 |
*>       | V2 V2 V2 T2 T2 |
*>       | V3 V3 V3 V3 T3 |
*>       |----------------|
*>
*> \endverbatim
*>
*  =====================================================================

      SUBROUTINE SORGRK(M, N, Q, LDQ)
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
*     .. Executable Statements ..
*
*     Break Q apart as follows
*
*           |-----|
*     Q =   | V T |
*           |-----|
*
*     Where T is an m-by-m upper triangular matrix, and V is as described
*     in the Further Details section
*
*     In turn, break apart V as follows
*
*           |---------|
*     V =   | V_2 V_1 |
*           |---------|
*
*     Where:
*
*     V_1 \in \R^{m\times m}   assumed unit lower triangular
*     V_2 \in \R^{m\times n-m}
*
*     Compute T = V_1'*T
*
      CALL STRTRM('Left', 'Upper', 'Transpose', 'Non-Unit', 'Unit',
     $         M, ONE, Q(1,N-M+1), LDQ, Q(1,N-M+1), LDQ)
*
*     Compute Q = -TV. This means that we need to break apart
*     Our computation in two parts
*
*           |---------------|
*     Q =   | -T*V_2 -T*V_1 |
*           |---------------|
*
*     Q_2 = -T*V_2 (TRMM) but only when necessary
*
      IF (N.GT.M) THEN
         CALL STRMM('Left', 'Upper', 'No Transpose', 'Non-Unit',
     $            M, N-M, NEG_ONE, Q(1,N-M+1), LDQ, Q, LDQ)
      END IF
*
*     Q_1 = -T*V_1 (Lower-Upper Matrix-Matrix multiplication)
*
      CALL SLUMM('Right', 'Unit', 'Non-Unit', M, NEG_ONE,
     $         Q(1,N-M+1), LDQ)
*
*     Q = "I" + Q
*
      J = MIN(M,N)
      DO I = 1, J
         Q(I,N-M+I) = Q(I,N-M+I) + ONE
      END DO
      END SUBROUTINE
