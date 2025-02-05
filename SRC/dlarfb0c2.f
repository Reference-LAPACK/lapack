      SUBROUTINE DLARFB0C2(SIDE, TRANS, DIRECT, STOREV, M, N, K, V, 
     $                     LDV, T, LDT, C, LDC)
         ! Scalar arguments
         INTEGER           M, N, K, LDV, LDC, LDT
         ! Note: We should probably get rid of SIDE and TRANS as these flags
         ! are not used as we are only using this routine inside
         ! x{OR,UN}G{QR,LQ,QR,LQ}, which have values that never change
         CHARACTER         SIDE, TRANS, DIRECT, STOREV 
         ! Array arguments
         DOUBLE PRECISION  V(LDV,*), C(LDC,*), T(LDT,*)
         ! Local scalars
         LOGICAL           QR, LQ, QL, DIRF, COLV
         ! External subroutines
         EXTERNAL          DGEMM, DTRMM
         ! External functions
         LOGICAL           LSAME
         EXTERNAL          LSAME
         ! Parameters
         DOUBLE PRECISION ONE, ZERO, NEG_ONE
         PARAMETER(ONE=1.0D+0, ZERO = 0.0D+0, NEG_ONE = -1.0D+0)

         ! Beginning of executable statements
         ! Convert our character flags to logical values
         DIRF = LSAME(DIRECT,'F')
         COLV = LSAME(STOREV,'C')

         ! Determine which of the 4 modes are using.
         ! QR is when we store the reflectors column by column and have the
         ! 'first' reflector stored in the first column
         QR = DIRF.AND.COLV

         ! LQ is when we store the reflectors row by row and have the
         ! 'first' reflector stored in the first row
         LQ = DIRF.AND.(.NOT.COLV)

         ! QL is when we store the reflectors column by column and have the
         ! 'first' reflector stored in the last column
         QL = (.NOT.DIRF).AND.COLV

         ! RQ is when we store the reflectors row by row and have the
         ! 'first' reflector stored in the last row
         ! RQ = (.NOT.DIRF).AND.(.NOT.COLV)
         ! Since we have exactly one of these 4 modes, we don't need to actually
         ! store the value of RQ, instead we assume this is the case if we fail
         ! the above 3 checks.

         IF (QR) THEN
            ! We are computing C = HC = (I - VTV')C
            ! Where: V = [ V1 ] and C = [ C1 ]
            !            [ V2 ]         [ C2 ]
            ! with the following dimensions:
            !     V1\in\R^{K\times K}
            !     V2\in\R^{M-K\times K}
            !     C1=0\in\R^{K\times N}
            !     C2\in\R^{M-K\times N}
            ! Since we are assuming that C1 is a zero matrix and it will be
            ! overwritten on exit, we can use this spot as a temporary workspace
            ! without having to allocate anything extra.
            ! This lets us simplify our above equation to get
            !
            ! C = HC = (I - [ V1 ]T [V1', V2'])[ 0  ]
            !               [ V2 ]             [ C2 ]
            !   = [ 0  ] - [ V1 ]T*V2'*C2
            !     [ C2 ]   [ V2 ]
            !
            !   = [ 0  ] - [ V1*T*V2'*C2 ]
            !     [ C2 ]   [ V2*T*V2'*C2 ]
            !
            !   = [      V1*T*V2'*C2 ]
            !     [ C2 - V2*T*V2'*C2 ]
            !
            ! So, we can order our computations as follows:
            !
            ! C1 = V2'*C2
            ! C1 = T*C1
            ! C2 = C2 - V2*C1
            ! C1 = -V1*C1
            !
            ! To achieve the same end result
            !
            ! C1 = V2'*C2
            !
            CALL DGEMM('Transpose', 'No Transpose', K, N, M - K,
     $                  ONE, V(K+1,1), LDV, C(K+1,1), LDC, ZERO, 
     $                  C, LDC)
            !
            ! C1 = T*C1
            !
            CALL DTRMM('Left', 'Upper', 'No Transpose', 'Non-unit', 
     $                  K, N, ONE, T, LDT, C, LDC)
            !
            ! C2 = C2 - V2*C1 = -V2*C1 + C2
            !
            CALL DGEMM('No Transpose', 'No Transpose', M-K, N, K,
     $                  NEG_ONE, V(K+1,1), LDV, C, LDC, ONE, 
     $                  C(K+1,1), LDC)
            !
            ! C1 = -V1*C1
            !
            CALL DTRMM('Left', 'Lower', 'No Transpose', 'Unit', 
     $                  K, N, NEG_ONE, V, LDV, C, LDC)
         ELSE IF (LQ) THEN
            ! We are computing C = CH' = C(I-V'T'V)
            ! Where: V = [ V1 V2 ] and C = [ C1 C2 ]
            ! with the following dimensions:
            !     V1\in\R^{K\times K}
            !     V2\in\R^{K\times N-K}
            !     C1=0\in\R^{M\times K}
            !     C2\in\R^{M\times N-K}
            ! Since we are assuming that C1 is a zero matrix and it will be
            ! overwritten on exit, we can use this spot as a temporary workspace
            ! without having to allocate anything extra.
            ! This lets us simplify our above equation to get
            !
            ! C = CH' = [ 0, C2 ](I - [ V1' ]T'[ V1, V2 ])
            !                         [ V2' ]
            !
            !   = [ 0, C2 ] - [ 0, C2 ][ V1' ]T'[ V1, V2 ]
            !                          [ V2' ]
            !
            !   = [ 0, C2 ] - C2*V2'*T'[ V1, V2 ]
            !
            !   = [ -C2*V2'*T'*V1, C2 - C2*V2'*T'*V2 ]
            !
            ! So, we can order our computations as follows:
            !
            ! C1 = C2*V2'
            ! C1 = C1*T'
            ! C2 = C2 - C1*V2
            ! C1 = -C1*V1
            !
            ! To achieve the same end result
            !
            ! C1 = C2*V2'
            !
            CALL DGEMM('No Transpose', 'Transpose', M, K, N-K, 
     $            ONE, C(1,K+1), LDC, V(1, K+1), LDV, ZERO, C, LDC)
            !
            ! C1 = C1*T'
            !
            CALL DTRMM('Right', 'Upper', 'Transpose', 'Non-unit',
     $            M, K, ONE, T, LDT, C, LDC)
            !
            ! C2 = C2 - C1*V2 = -C1*V2 + C2
            !
            CALL DGEMM('No Transpose', 'No Transpose', M, N-K, K,
     $            NEG_ONE, C, LDC, V(1,K+1), LDV, ONE, 
     $            C(1,K+1), LDC)
            !
            ! C1 = -C1*V1
            !
            CALL DTRMM('Right', 'Upper', 'No Transpose', 'Unit',
     $            M, K, NEG_ONE, V, LDV, C, LDC)
         ELSE IF (QL) THEN
            ! We are computing C = HC = (I - VTV')C
            ! Where: V = [ V2 ] and C = [ C2 ]
            !            [ V1 ]         [ C1 ]
            ! with the following dimensions:
            !     V1\in\R^{K\times K}
            !     V2\in\R^{M-K\times K}
            !     C1=0\in\R^{K\times N}
            !     C2\in\R^{M-K\times N}
            ! Since we are assuming that C1 is a zero matrix and it will be
            ! overwritten on exit, we can use this spot as a temporary workspace
            ! without having to allocate anything extra.
            ! This lets us simplify our above equation to get
            !
            ! C = HC = (I-[ V2 ]T[ V2' V1' ])[ C2 ]
            !             [ V1 ]             [ 0  ]
            !
            !   = [ C2 ] - [ V2 ]T*V2'*C2
            !     [ 0  ]   [ V1 ]
            !
            !   = [ C2 ] - [ V2*T*V2'*C2 ]
            !     [ 0  ]   [ V1*T*V2'*C2 ]
            !
            !   = [ C2 - V2*T*V2'*C2 ]
            !     [    - V1*T*V2'*C2 ]
            !
            ! So, we can order our computations as follows:
            !
            ! C1 = V2'*C2
            ! C1 = T*C1
            ! C2 = C2 - V2*C1
            ! C1 = -V1*C1
            !
            ! To achieve the same end result
            !
            ! C1 = V2'*C2
            !
            CALL DGEMM('Transpose', 'No Transpose', K, N, M-K, 
     $         ONE, V, LDV, C, LDC, ZERO, C(M-K+1, 1), LDC)
            !
            ! C1 = T*C1
            !
            CALL DTRMM('Left', 'Lower', 'No Transpose', 'Non-unit',
     $         K, N, ONE, T, LDT, C(M-K+1,1), LDC)
            !
            ! C2 = C2 - V2*C1 = -V2*C1 + C2
            !
            CALL DGEMM('No Transpose', 'No Transpose', M-K, N, K, 
     $         NEG_ONE, V, LDV, C(M-K+1,1), LDC, ONE, C, LDC)
            !
            ! C1 = -V1*C1
            !
            CALL DTRMM('Left', 'Upper', 'No Transpose', 'Unit', 
     $         K, N, NEG_ONE, V(M-K+1,1), LDV, C(M-K+1,1), LDC)
         ELSE ! IF (RQ) THEN
            ! We are computing C = CH' = C(I-V'T'V)
            ! Where: V = [ V2 V1] and C = [ C2 C1 ]
            ! with the following dimensions:
            !     V1\in\R^{K\times K}
            !     V2\in\R^{K\times N-K}
            !     C1=0\in\R^{M\times K}
            !     C2\in\R^{M\times N-K}
            ! Since we are assuming that C1 is a zero matrix and it will be
            ! overwritten on exit, we can use this spot as a temporary workspace
            ! without having to allocate anything extra.
            ! This lets us simplify our above equation to get
            !
            ! C = CH' = [ C2, 0 ] (I - [ V2' ]T'[ V2, V1 ]
            !                          [ V1' ]
            !
            !   = [ C2, 0 ] - [ C2, 0 ] [ V2' ]T'[ V2, V1 ]
            !                           [ V1' ]
            !
            !   = [ C2, 0 ] - C2*V2'*T'[ V2, V1 ]
            !
            !   = [ C2, 0 ] - [ C2*V2'*T'*V2, C2*V2'*T'*V1 ] 
            !
            !   = [ C2 - C2*V2'*T'*V2, -C2*V2'*T'*V1 ]
            !
            ! So, we can order our computations as follows:
            !
            ! C1 = C2*V2'
            ! C1 = C1*T'
            ! C2 = C2 - C1*V2
            ! C1 = -C1*V1
            ! 
            !
            ! To achieve the same end result
            !
            ! C1 = C2*V2'
            !
            CALL DGEMM('No Transpose', 'Transpose', M, K, N-K,
     $         ONE, C, LDC, V, LDV, ZERO, C(1, N-K+1), LDC)
            !
            ! C1 = C1*T'
            !
            CALL DTRMM('Right', 'Lower', 'Transpose', 'Non-unit', 
     $         M, K, ONE, T, LDT, C(1, N-K+1), LDC)
            !
            ! C2 = C2 - C1*V2 = -C1*V2 + C2
            !
            CALL DGEMM('No Transpose', 'No Transpose', M, N-K, K, 
     $         NEG_ONE, C(1, N-K+1), LDC, V, LDV, ONE, C, LDC)
            !
            ! C1 = -C1*V1
            !
            CALL DTRMM('Right', 'Lower', 'No Transpose', 'Unit', 
     $         M, K, NEG_ONE, V(1, N-K+1), LDV, C(1,N-K+1), LDC)
         END IF
      END SUBROUTINE
