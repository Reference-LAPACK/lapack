c     Cost: n > k: 1/6 * (k^2-1)(2n+k)
c           n = k: 1/2 * (n^3-n)
      RECURSIVE SUBROUTINE MY_DLARFT_REC( DIRECT, STOREV, N, K, V, LDV,
     $                                    TAU, T, LDT)
         IMPLICIT NONE
         ! Arguemnts
         ! Scalars
         INTEGER           N, K, LDV, LDT
         CHARACTER         DIRECT, STOREV
         ! Matrix 
         DOUBLE PRECISION  V(LDV,*), T(LDT,*), TAU(N)

         ! Local variables
         INTEGER           I,J,L,MINNK
         LOGICAL           QR,LQ,QL,DIRF,COLV
         ! Parameters
         DOUBLE PRECISION ONE, NEG_ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO = 0.0, NEG_ONE=-1.0D+0)
         ! External functions
         LOGICAL           LSAME
         EXTERNAL          LSAME
         ! External subroutines
         EXTERNAL          DTRMM,DGEMM,DLACPY

         ! Break V apart into 6 components
         ! V = |---------------|
         !     |V_{1,1} 0      |
         !     |V_{2,1} V_{2,2}|
         !     |V_{3,1} V_{3,2}|
         !     |---------------|
         ! V_{1,1}\in\R^{k,k} unit lower triangular
         ! V_{2,1}\in\R^{n-k,k} rectangular
         ! V_{3,1}\in\R^{m-n,k} rectangular
         ! 
         ! V_{2,2}\in\R^{n-k,n-k} unit upper triangular
         ! V_{3,2}\in\R^{m-n,n-k} rectangular

         ! We will construct the T matrix 
         ! T = |---------------| =  |--------|
         !     |T_{1,1} T_{1,2}|    |T_1  T_3|
         !     |0       T_{2,2}|    |0    T_2|
         !     |---------------|    |--------|

         ! T is the triangular factor attained from block reflectors. 
         ! To motivate the structure, consider the product
         !
         ! (I - V_1T_1V_1^\top)(I - V_2T_2V_2^\top)
         ! = I - V_1T_1V_1^\top - V_2T_2V_2^\top + V_1T_1V_1^\topV_2T_2V_2^\top
         !
         ! Define T_3 = -T_1V_1^\topV_2T_2
         !
         ! Then, we can define the matrix V as 
         ! V = |-------|
         !     |V_1 V_2|
         !     |-------|
         !
         ! So, our product is equivalent to the matrix product
         ! I - VTV^\top
         ! So, we compute T_1, then T_2, then use these values to get T_3
         !
         ! The general scheme used is inspired by the approach inside DGEQRT3
         ! which was (at the time of writing this code):
         ! Based on the algorithm of Elmroth and Gustavson,
         ! IBM J. Res. Develop. Vol 44 No. 4 July 2000.

         IF(K.EQ.0.OR.N.EQ.0) THEN
            RETURN
         END IF
         ! Base case
         IF(K.EQ.1.OR.N.EQ.1) THEN
            T(1,1) = TAU(1)
            RETURN
         END IF

         ! Beginning of executable statements
!        MINNK = MIN(N,K)
!        L = MINNK / 2
         L = K / 2
         ! Determine what kind of Q we need to compute
         ! We assume that if the user doesn't provide 'F' for DIRECT,
         ! then they meant to provide 'B' and if they don't provide
         ! 'C' for STOREV, then they meant to provide 'R'
         DIRF = LSAME(DIRECT,'F')
         COLV = LSAME(STOREV,'C')
         ! QR happens when we have forward direction in column storage
         QR = DIRF.AND.COLV
         ! LQ happens when we have Forward direction in row storage
         LQ = DIRF.AND.(.NOT.COLV)
         ! QL happens when we have backward direction in column storage
         QL = (.NOT.DIRF).AND.COLV
         ! The last case is RQ. Due to how we strucutured this, if the
         ! above 3 are false, then RQ must be true, so we never store 
         ! this
         ! RQ happens when we have backward direction in row storage
         !RQ = (.NOT.DIRF).AND.(.NOT.COLV)


         ! Compute T3
         IF(QR) THEN
            ! If we are wide, then our 
            ! Compute T_1
            CALL MY_DLARFT_REC(DIRECT, STOREV, N, L, V, LDV, TAU, T, 
     $            LDT)
            ! Compute T_2
            CALL MY_DLARFT_REC(DIRECT, STOREV, N-L, K-L, V(L+1,L+1),
     $         LDV, TAU(L+1), T(L+1,L+1), LDT)
            ! Compute T_3 
            ! T_3 = V_{2,1}^\top
            DO J = 1, L
               DO I = 1, K-L
                  T(J,L+I) = V(L+I,J)
               END DO
            END DO
            ! T_3 = V_{2,1}^\top * V_{2,2}
            CALL DTRMM('Right', 'Lower', 'No transpose', 'Unit', 
     $            L, K - L, ONE, V(L+1, L+1), LDV, T(1, L+1), LDT)

            IF(N.GT.K) THEN
            ! T_3 = T_3 + V_{3,1}^\topV_{3,2}
               CALL DGEMM('Transpose', 'No transpose', L, K-L, N-K,
     $               ONE, V(K+1, 1), LDV, V(K+1,L+1), LDV, ONE, 
     $               T(1, L+1), LDT)
            END IF

            ! At this point, we have that T_3 = V_1^\top *V_2
            ! All that is left is to pre and post multiply by -T_1 and T_2
            ! respectively.

            ! T_3 = -T_1*T_3
            CALL DTRMM('Left', 'Upper', 'No transpose', 'Non-unit',
     $            L, K - L, NEG_ONE, T, LDT, T(1, L+1), LDT)
            ! T_3 = T_3*T_2
            CALL DTRMM('Right', 'Upper', 'No transpose', 'Non-unit',
     $            L, K - L, ONE, T(L+1,L+1), LDT, T(1, L+1), LDT)

         ELSE IF(LQ) THEN
            ! Compute T_1
            CALL MY_DLARFT_REC(DIRECT, STOREV, N, L, V, LDV, TAU, T,
     $         LDT)
            ! Compute T_2
            CALL MY_DLARFT_REC(DIRECT, STOREV, N-L, K-L, V(L+1,L+1),
     $         LDV, TAU(L+1), T(L+1,L+1), LDT)

            ! Begin computing T_3
            ! First, T_3 = V_1V_2^\top
            ! T_3 = V_{12}
            CALL DLACPY('All', L, K - L, V(1,L+1), LDV, T(1, L+1), LDT)

            ! T_3 = V_{12}V_{22}^\top = T_3V_{22}^\top
            CALL DTRMM('Right', 'Upper', 'Transpose', 'Unit', L, K-L,
     $         ONE, V(L+1, L+1), LDV, T(1, L+1), LDT)

            ! If needed, use the trailing components
            IF(N.GT.K) THEN
               CALL DGEMM('No transpose', 'Transpose', L, K-L, N-K, 
     $            ONE, V(1, K+1), LDV, V(L+1, K+1), LDV, ONE,
     $            T(1, L+1), LDT)
            END IF

            ! T_3 = -T_1T_3
            CALL DTRMM('Left', 'Upper', 'No transpose', 'Non-unit',
     $         L, K - L, NEG_ONE, T, LDT, T(1, L+1), LDT)

            ! T_3 = T_3T_1
            CALL DTRMM('Right', 'Upper', 'No transpose', 'Non-unit',
     $         L, K - L, ONE, T(L+1,L+1), LDT, T(1, L+1), LDT)
         ELSE IF(QL) THEN
            ! Compute T_1
            CALL MY_DLARFT_REC(DIRECT, STOREV, N-L, K-L, V, LDV, TAU,
     $         T, LDT)
            ! Compute T_2
            CALL MY_DLARFT_REC(DIRECT, STOREV, N, L, V(1, K-L+1), LDV,
     $         TAU(K-L+1), T(K-L+1,K-L+1), LDT)

            ! Begin computing T_3 = T_2V_2^\topV_1T_1

            ! T_3 = V_2^\top V_1

            ! T_3 = V_{2,2}^\top
            DO J = 1, K-L
               DO I = 1, L
                  T(K-L+I,J) = V(N-K+J, K-L+I)
               END DO
            END DO

            ! T_3 = V_{2,2}^\topV_{2,1} = T_3V_{2,1}
            CALL DTRMM('Right', 'Upper', 'No transpose', 'Unit',
     $         L, K - L, ONE, V(N-K+1,1), LDV, T(K-L+1,1), LDT)

            ! If needed, T_3 = V_{1,2}^\topV_{1,1} + T_3
            IF(N.GT.K) THEN
               CALL DGEMM('Transpose', 'No transpose', L, K-L, N-K,
     $            ONE, V(1,K-L+1), LDV, V, LDV, ONE, T(K-L+1,1), LDT)
            END IF

            ! T_3 = -T_2T_3
            CALL DTRMM('Left', 'Lower', 'No transpose', 'Non-unit',
     $         L, K-L, NEG_ONE, T(K-L+1,K-L+1), LDT, T(K-L+1,1), LDT)
            ! T_3 = T_3T_1
            CALL DTRMM('Right', 'Lower', 'No transpose', 'Non-unit',
     $         L, K-L, ONE, T, LDT, T(K-L+1,1), LDT)
         ELSE
            ! Else means RQ
            ! Compute T_1
            CALL MY_DLARFT_REC(DIRECT, STOREV, N-L, K-L, V, LDV, TAU,
     $         T, LDT)
            ! Compute T_2
            CALL MY_DLARFT_REC(DIRECT, STOREV, N, L, V(K-L+1,1), LDV,
     $         TAU(K-L+1), T(K-L+1,K-L+1), LDT)

            ! Begin computing T_3 = T_2V_2V_1^\topT_1

            ! T_3 = V_2V_1^\top

            ! T_3 = V_{2,2}
            CALL DLACPY('All', L, K-L, V(K-L+1,N-K+1), LDV, 
     $         T(K-L+1,1), LDT)

            ! T_3 = T_3V_{1,2}^\top
            CALL DTRMM('Right', 'Lower', 'Transpose', 'Unit',
     $         L, K-L, ONE, V(1, N-K+1), LDV, T(K-L+1,1), LDT)

            ! If needed, T_3 = V_{2,1}V_{1,1}^\top + T_3
            IF(N.GT.K) THEN
               CALL DGEMM('No transpose', 'Transpose', L, K-L, N-K,
     $            ONE, V(K-L+1,1), LDV, V, LDV, ONE, T(K-L+1,1), LDT)
            END IF

            ! T_3 = -T_2T_3
            CALL DTRMM('Left', 'Lower', 'No tranpose', 'Non-unit',
     $         L, K-L, NEG_ONE, T(K-L+1,K-L+1), LDT, T(K-L+1,1), LDT)

            ! T_3 = T_3T_1
            CALL DTRMM('Right', 'Lower', 'No tranpose', 'Non-unit',
     $         L, K-L, ONE, T, LDT, T(K-L+1,1), LDT)
         END IF
         
         ! Now, we have T in the correct form!
      END SUBROUTINE
