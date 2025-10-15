*> \brief \b DLARFT2 forms the triangular factor T of a block reflector H = I - vtvH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARFT2( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIRECT, STOREV
*       INTEGER            K, LDT, LDV, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARFT2 forms the triangular factor T of a real block reflector H
*> of order n, which is defined as a product of k elementary reflectors.
*>
*> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*>
*> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*>
*> If STOREV = 'C', the vector which defines the elementary reflector
*> H(i) is stored in the i-th column of the array V, and
*>
*>    H  =  I - V * T * V**T
*>
*> If STOREV = 'R', the vector which defines the elementary reflector
*> H(i) is stored in the i-th row of the array V, and
*>
*>    H  =  I - V**T * T * V
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIRECT
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          Specifies the order in which the elementary reflectors are
*>          multiplied to form the block reflector:
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*> \endverbatim
*>
*> \param[in] STOREV
*> \verbatim
*>          STOREV is CHARACTER*1
*>          Specifies how the vectors which define the elementary
*>          reflectors are stored (see also Further Details):
*>          = 'C': columnwise
*>          = 'R': rowwise
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the block reflector H. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The order of the triangular factor T (= the number of
*>          elementary reflectors). K >= 1.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension
*>                               (LDV,K) if STOREV = 'C'
*>                               (LDV,N) if STOREV = 'R'
*>          The matrix V. See further details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT,K)
*>          The k by k triangular factor T of the block reflector.
*>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*>          lower triangular. The rest of the array is not used.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= K.
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
*> \ingroup larft
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The shape of the matrix V and the storage of the vectors which define
*>  the H(i) is best illustrated by the following example with n = 5 and
*>  k = 3. The elements equal to 1 are not stored.
*>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*>                   ( v1  1    )                     (     1 v2 v2 v2 )
*>                   ( v1 v2  1 )                     (        1 v3 v3 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*>
*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*>                   (     1 v3 )
*>                   (        1 )
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLARFT2(DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT)
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*        .. Scalar Arguments
*
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
*
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*     .. Parameters ..
*
      DOUBLE PRECISION ONE, NEG_ONE, ZERO
      PARAMETER(ONE=1.0D+0, ZERO = 0.0D+0, NEG_ONE=-1.0D+0)
*
*     .. Local Scalars ..
*
      INTEGER           I,J,KMI,NMI,L,NX
      LOGICAL           QR,LQ,QL,RQ,LQT,RQT,DIRF,COLV,TDIRF,TCOLV
*
*     .. External Subroutines ..
*
      EXTERNAL          DTRMM,DGEMM,DLACPY
*
*     .. External Functions..
*
      LOGICAL           LSAME
      INTEGER           ILAENV
      EXTERNAL          LSAME, ILAENV
      ! temporary parameters
*     ..
*     .. Executable Statements ..
*
*     Determine what kind of Q we need to compute
*     We assume that if the user doesn't provide 'F' for DIRECT,
*     then they meant to provide 'B' and if they don't provide
*     'C' for STOREV, then they meant to provide 'R'
*
      DIRF = LSAME(DIRECT,'F')
      TDIRF = LSAME(DIRECT,'T')
      COLV = LSAME(STOREV,'C')
      TCOLV = LSAME(STOREV,'T')
*
*     QR happens when we have forward direction in column storage
*
      QR = DIRF.AND.COLV
*
*     LQT happens when we have forward direction in row storage and want to compute the transpose of
*     the T we would normally compute
*
      LQT = DIRF.AND.TCOLV
*
*     LQ happens when we have forward direction in row storage and want to compute the T we would
*     normally compute
*
      LQ = DIRF.AND.(.NOT.LQT)
*
*     QL happens when we have backward direction in column storage
*
      QL = (.NOT.DIRF).AND.COLV
*
*     RQT happens when we have backward direction in row storage and want to compute the transpose
*     of the T we would normally compute
*
      RQT = TDIRF.AND.(.NOT.COLV)
*
*     RQ happens when we have backward direction in row storage and want to compute the T that we
*     would normally compute
*
      RQ = (.NOT.RQT).AND.(.NOT.COLV)
*
*     If you want to see comments of how this subroutine works, we 
*     are essentially unrolling the recursion present in DLARFT, so
*     to see what we are doing in each step (for each switch) look at the
*     comments in dlarft.f by replacing K with I and L with I-1
*     
*
      IF (QR) THEN
*
*        Break V into 9 components
*        
*        V = |-----------------------|
*            |V_{1,1} 0       0      | i-1
*            |V_{2,1} V_{2,2} 0      | 1
*            |V_{3,1} V_{3,2} V_{3,3}| n-i
*            |-----------------------|
*             i-1     1       k-i
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit lower triangular
*
*        This is how we are going to view the matrix V at each step 
*        i=2,\dots,k, then we grow into V_{3,3} and repeat until we
*        reach the end. On each iteration V_{3,3} is not referenced
*
*        We will construct T one column at a time from left to right
*        after initializing T(1,1) = TAU(1)
*
*        T = |-------------------------|
*            | T_{1,1} T_{1,2} T_{1,3} | i-1
*            | 0       T_{2,2} T_{2,3} | 1
*            | 0       0       T_{3,3} | k-i
*            |-------------------------|
*              i-1     1       k-i
*
*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular 
*
*        Similarly as above, we will construct T_{1,2} and T_{2,2} at
*        each iteration i = 2, \dots k, and then grow into T_{1:3,3}. On
*        each iteration, T_{1:3,3} are not referenced. See dlarft.f
*        for details on how these formulae were constructed.
*
*        We now get
*
*        T_{1,2} = -T_{1,1}[V_{1,1}\\V_{2,1}\\V_{3,1}]'
*                       [0\\V_{2,2}\\V_{3,2}]T_{2,2}
*
*        T_{1,2} = -T_{1,1}(V_{2,1}' + V_{3,1}'V_{3,2})T_{2,2}
*
*        This means we will do the following
*
*        T_{1,2} = -V_{2,1}'T_{2,2} = -\tau_{i}V_{2,1}'
*        T_{1,2} = -\tau_{i}V_{3,2}' V_{3,1} + T_{1,2}
*        T_{1,2} = T_{1,1}T_{1,2}
*        T_{2,2} = \tau{i}
*
         T(1,1) = TAU(1)

         DO I = 2, K
*
*           T_{1,2} = -V_{2,1}'V_{2,2}T_{2,2} = -\tau_i V_{2,1}'
*           We must do this at copy time as otherwise gemv will do nothing
*           on the last column when n=k, but we neet to make sure we are
*           scaled by this value
*
            DO J = 1, I-1
               T(J,I) = -V(I,J)*TAU(I)
            END DO

*
*           T_{1,2} = -V_{3,1}'V_{3,2}T_{2,2} + T_{1,2}
*                   = -\tau{i} V_{3,2}'V_{3,1} + T_{1,2}
*
            CALL DGEMV('Transpose', N-I, I-1, -TAU(I), V(I+1,1),
     $            LDV, V(I+1,I), 1, ONE, T(1, I), 1)


*
*           T_{1,2} = T_{1,1}T_{1,2}
*
            CALL DTRMV('Upper', 'No Transpose', 'Non-unit', I-1, 
     $            T, LDT, T(1,I), 1)

*
*           T_{2,2} = \tau{i}
*
            T(I,I) = TAU(I)
         END DO
      ELSE IF (LQ) THEN
*
*        Break V into 9 components
*
*        V = |-------------------------|
*            | V_{1,1} V_{1,2} V_{1,3} | i-1
*            | 0       V_{2,2} V_{2,3} | 1
*            | 0       0       V_{3,3} | k-i
*            |-------------------------|
*              i-1     1       n-i
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit upper triangular
*
*        This is how we are going to view the matrix V at each step 
*        i=2,\dots,k, then we grow into V_{3,3} and repeat until we
*        reach the end. On each iteration V_{3,3} is not referenced
*
*        We will construct T one column at a time from left to right
*        after initializing T(1,1) = TAU(1)
*
*        T = |-------------------------|
*            | T_{1,1} T_{1,2} T_{1,3} | i-1
*            | 0       T_{2,2} T_{2,3} | 1
*            | 0       0       T_{3,3} | k-i
*            |-------------------------|
*              i-1     1       k-i
*
*        Similarly as above, we will construct T_{1,2} and T_{2,2} at
*        each iteration i = 2, \dots k, and then grow into T_{1:3,3}. On
*        each iteration, T_{1:3,3} are not referenced. See dlarft.f
*        for details on how these formulae were constructed.
*
*        We now get
*
*        T_{1,2} = -T_{1,1}[V_{1,1} V_{1,2} V_{1,3}][ 0 V_{2,2} V_{2,3} ]'T_{2,2}
*
*        T_{1,2} = -T_{1,1}(V_{1,2} + V_{1,3}V_{2,3}')T_{2,2}
*
*        This means we will do the following
*
*        T_{1,2} = -V_{1,2}T_{2,2} = -\tau_{i}V_{1,2}
*        T_{1,2} = -\tau_{i}V_{1,3}V_{2,3}' + T_{1,2}
*        T_{1,2} = T_{1,1}T_{1,2}
*        T_{2,2} = \tau{i}
*
         T(1,1) = TAU(1)

         DO I = 2, K

            ! T_{1,2} = -\tau_{i}V_{1,2}
            DO J = 1, I-1
               T(J, I) = -TAU(I)*V(J, I)
            END DO

            ! T_{1,2} = -\tau_{i}V_{1,3}V_{2,3}' + T_{1,2}
            CALL DGEMV('No Transpose', I-1, N-I, -TAU(I), V(1, I+1),
     $            LDV, V(I, I+1), LDV, ONE, T(1,I), 1)

            ! T_{1,2} = T_{1,1}T_{1,2}
            CALL DTRMV('Upper', 'No Transpose', 'Non-unit', I-1, 
     $            T, LDT, T(1,I), 1)

            ! T_{2,2} = \tau{i}
            T(I,I) = TAU(I)
         END DO
      ELSE IF (LQT) THEN
*
*        Break V into 9 components
*
*        V = |-------------------------|
*            | V_{1,1} V_{1,2} V_{1,3} | i-1
*            | 0       V_{2,2} V_{2,3} | 1
*            | 0       0       V_{3,3} | k-i
*            |-------------------------|
*              i-1     1       n-i
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit upper triangular
*
*        This is how we are going to view the matrix V at each step 
*        i=2,\dots,k, then we grow into V_{3,3} and repeat until we
*        reach the end. On each iteration V_{3,3} is not referenced
*
*        We will construct T one column at a time from left to right
*        after initializing T(1,1) = TAU(1)
*
*        T = |-------------------------|
*            | T_{1,1} 0       0       | i-1
*            | T_{2,1} T_{2,2} 0       | 1
*            | T_{3,1} T_{3,2} T_{3,3} | k-i
*            |-------------------------|
*              i-1     1       k-i
*
*        Similarly as above, we will construct T_{2,1} and T_{2,2} at
*        each iteration i = 2, \dots k, and then grow into T_{3,1:3}. On
*        each iteration, T_{3,1:3} are not referenced. See dlarft.f
*        for details on how these formulae were constructed.
*
*        We now get
*
*        T_{2,1} = -T_{2,2}[0 V_{2,2} V_{2,3}][V_{1,1} V_{1,2} V_{1,3}]'T_{1,1}
*
*        T_{2,1} = -T_{2,2}(V_{1,2}' + V_{2,3}V_{1,3}')T_{1,1}
*
*        This means we will do the following
*
*        T_{2,1} = -T_{2,2}V_{1,2}' = -\tau_{i}V_{1,2}'
*        T_{2,1} = -\tau_{i}V_{1,3}V_{2,3}' + T_{2,1}
*        T_{2,1} = T_{1,1}'T_{2,1}
*        T_{2,2} = \tau{i}
*
         T(1,1) = TAU(1)

         DO I = 2, K
*
*           T_{2,1} = -\tau_{i}V_{1,2}'
*
            DO J = 1, I-1
               T(I,J) = -TAU(I)*V(J,I)
            END DO
*
*           T_{2,1} = -\tau_{i}V_{1,3}V_{2,3}' + T_{2,1}
*
            CALL DGEMV('No transpose', I-1, N-I, -TAU(I), V(1, I+1),
     $            LDV, V(I, I+1), LDV, ONE, T(I, 1), LDT)

            ! T_{2,1}' = T_{1,1}'T_{2,1}'
            CALL DTRMV('Lower', 'Transpose', 'Non-unit', I-1, 
     $            T, LDT, T(I,1), LDT)

            T(I,I) = TAU(I)
         END DO
      ELSE IF (QL) THEN
*
*     Break V into 9 components
*     
*     V = |-------------------------|
*         | V_{1,1} V_{1,2} V_{1,3} | n-i
*         | 0       V_{2,2} V_{2,3} | 1
*         | 0       0       V_{3,3} | i-1
*         |-------------------------|
*           k-i     1       i-1
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit upper triangular
*
*        This is how we are going to view the matrix V at each step 
*        i=2,\dots,k, then we grow into V_{1,1} and repeat until we
*        reach the end. On each iteration V_{1,1} is not referenced
*
*        We will construct T one column at a time from right to left
*        after initializing T(K,K) = TAU(K)
* 
*     T = |-------------------------|
*         | T_{1,1} 0       0       | k-i
*         | T_{2,1} T_{2,2} 0       | 1
*         | T_{3,1} T_{3,2} T_{3,3} | i-1
*         |-------------------------|
*           k-i     1       i-1
*
*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular 
*
*        Similarly as above, we will construct T_{2,2} and T_{3,2} at
*        each iteration i = 2, \dots k, and then grow into T_{1:3,1}. On
*        each iteration, T_{1:3,1} are not referenced. See dlarft.f
*        for details on how these formulae were constructed.
*
*        We get that
*
*        T_{3,2} = -T_{3,3}[V_{1,3}\\V_{2,3}\\V_{3,3}]'
*        [V_{1,2}\\V_{2,2}\\0]T_{2,2}
*
*        T_{3,2} = -T_{3,3}(V_{1,3}'V_{1,2} + V_{2,3}')T_{2,2}
*
*        Thus, we will compute
*
*        T_{2,2} = \tau_{k-i+1}
*        T_{3,2} = -T_{3,3}V_{3,2}' = -\tau_{k-i+1}V_{3,2}'
*        T_{3,2} = -\tau_{k-i+1}V_{1,3}'V_{1,2} + T_{3,2}
*        T_{3,2} = T_{3,3}T_{3,2}
*
         T(K,K) = TAU(K)
         DO I = 2, K
            KMI = K-I+1
            NMI = N-I+1
*
*             T_{2,2} = \tau(k-i+1)
*
            T(KMI,KMI) = TAU(KMI)
*
*             T_{3,2} = -\tau(k-i+1)V_{2,3}'
*
            DO J = 1, I-1
               T(KMI + J, KMI) = -TAU(KMI)*V(NMI, KMI + J)
            END DO
*
*             T_{3,2} = -\tau(k-i+1)V_{1,3}'V_{1,2} + T_{3,2}
*
            CALL DGEMV('Transpose', N-I, I-1, -TAU(KMI), 
     $            V(1, KMI + 1), LDV, V(1, KMI), 1, ONE,
     $            T(KMI+1, KMI), 1)
*
*             T_{3,2} = T_{3,3}T_{3,2}
*
            CALL DTRMV('Lower', 'No Transpose', 'Non-unit', I-1,
     $            T(KMI + 1, KMI + 1), LDT, T(KMI + 1, KMI), 1)

         END DO
      ELSE IF (RQ) THEN
*
*     Break V into 9 components
*
*     V = |-------------------------|
*         | V_{1,1} 0       0       | k-i
*         | V_{2,1} V_{2,2} 0       | 1
*         | V_{3,1} V_{3,2} V_{3,3} | i-1
*         |-------------------------|
*           n-i     1       i-1
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit lower triangular
*
*        This is how we are going to view the matrix V at each step 
*        i=2,\dots,k, then we grow into V_{1,1} and repeat until we
*        reach the end. On each iteration V_{1,1} is not referenced
*
*        We will construct T one column at a time from right to left
*        after initializing T(K,K) = TAU(K)
*        
*     T = |-------------------------|
*         | T_{1,1} 0       0       | k-i
*         | T_{2,1} T_{2,2} 0       | 1
*         | T_{3,1} T_{3,2} T_{3,3} | i-1
*         |-------------------------|
*           k-i     1       i-1
*
*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular 
*
*        Similarly as above, we will construct T_{2,2} and T_{3,2} at
*        each iteration i = 2, \dots k, and then grow into T_{1:3,1}. On
*        each iteration, T_{1:3,1} are not referenced. See dlarft.f
*        for details on how these formulae were constructed.
*
*        We get that
*
*        T_{3,2} = -T_{3,3}[V_{3,1} V_{3,2} V_{3,3}][V_{2,1} V_{2,2} 0]'T_{2,2}
*
*        T_{3,2} = -T_{3,3}(V_{3,1}V_{2,1}' + V_{3,2})T_{2,2}
*
*        Thus, we will compute
*
*        T_{2,2} = \tau_{k-i+1}
*        T_{3,2} = -\tau_{k-i+1}V_{3,2}
*        T_{3,2} = -\tau_{k-i+1}V_{3,1}V_{2,1}' + T_{3,2}
*        T_{3,2} = T_{3,3}T_{3,2}
*
         T(K,K) = TAU(K)
         DO I = 2, K
            KMI = K-I+1
            NMI = N-I+1
*
*           T_{2,2} = \tau_{k-i+1}
*
            T(KMI,KMI) = TAU(KMI)
*
*           T_{3,2} = -\tau_{k-i+1}V_{3,2}
*
            DO J = 1, I-1
               T(KMI + J, KMI) = -TAU(KMI)*V(KMI + J, NMI)
            END DO
*
*           T_{3,2} = -\tau_{k-i+1}V_{3,1}V_{2,1}' + T_{3,2}
*
            CALL DGEMV('No Transpose', I-1, N-I, -TAU(KMI),
     $            V(KMI+1, 1), LDV, V(KMI, 1), LDV, ONE, T(KMI+1, KMI),
     $            1)
*
*           T_{3,2} = T_{3,3}T_{3,2}
*
            CALL DTRMV('Lower', 'No Transpose', 'Non-unit', I-1, 
     $            T(KMI+1, KMI+1), LDT, T(KMI+1, KMI), 1)
         END DO
      ELSE IF (RQT) THEN
*
*     Break V into 9 components
*
*     V = |-------------------------|
*         | V_{1,1} 0       0       | k-i
*         | V_{2,1} V_{2,2} 0       | 1
*         | V_{3,1} V_{3,2} V_{3,3} | i-1
*         |-------------------------|
*           n-i     1       i-1
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit lower triangular
*
*        This is how we are going to view the matrix V at each step 
*        i=2,\dots,k, then we grow into V_{1,1} and repeat until we
*        reach the end. On each iteration V_{1,1} is not referenced
*
*        We will construct T one column at a time from right to left
*        after initializing T(K,K) = TAU(K)
*        
*     T = |-------------------------|
*         | T_{1,1} T_{1,2} T_{1,3} | k-i
*         | 0       T_{2,2} T_{2,3} | 1
*         | 0       0       T_{3,3} | i-1
*         |-------------------------|
*           k-i     1       i-1
*
*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular 
*
*        Similarly as above, we will construct T_{2,2} and T_{2,3} at
*        each iteration i = 2, \dots k, and then grow into T_{1,1:3}. On
*        each iteration, T_{1,1:3} are not referenced. See dlarft.f
*        for details on how these formulae were constructed.
*
*        We get that
*
*        T_{2,3} = -T_{2,2}[V_{2,1} V_{2,2} 0][V_{3,1} V_{3,2} V_{3,3}]'T_{3,3}
*
*        After transposing when necessary to fit our blas routines' interface,
*           we get
*        T_{3,2} = -T_{2,2}(V_{3,1}V_{2,1} + V_{3,2})T_{3,3}
*
*        Thus, we will compute
*
*        T_{2,2} = \tau_{k-i+1}
*        T_{2,3} = -\tau_{k-i+1}V_{3,2}'
*        T_{2,3} = -\tau_{k-i+1}V_{3,1}V_{2,1} + T_{2,3}
*        T_{2,3} = T_{3,3}'T_{2,3}
*
         T(K,K) = TAU(K)
         DO I = 2, K
            KMI = K-I+1
            NMI = N-I+1
*
*           T_{2,2} = \tau_{k-i+1}
*
            T(KMI,KMI) = TAU(KMI)
*
*           T_{2,3} = -\tau_{k-i+1}V_{3,2}'
*
            DO J = 1, I-1
               T(KMI, KMI + J) = -TAU(KMI)*V(KMI + J, NMI)
            END DO
*
*           T_{2,3} = -\tau_{k-i+1}V_{3,1}V_{2,1} + T_{2,3}
*
            CALL DGEMV('No Transpose', I-1, N-I, -TAU(KMI),
     $            V(KMI+1, 1), LDV, V(KMI, 1), LDV, ONE, T(KMI, KMI+1),
     $            LDT)
*
*           T_{2,3} = T_{3,3}'T_{2,3}
*
            CALL DTRMV('Upper', 'Transpose', 'Non-unit', I-1, 
     $            T(KMI+1, KMI+1), LDT, T(KMI, KMI+1), LDT)
         END DO
      END IF
      END SUBROUTINE
