      ! Essentially a recursive syrk but A is triangular
      RECURSIVE SUBROUTINE DST3RK(UPLOT, UPLOC, TRANS, DIAG, K,
     $            ALPHA, T, LDT, BETA, C, LDC)
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA,BETA
      INTEGER           K,LDT,LDT
      CHARACTER         UPLOT,UPLOC,TRANS,DIAG
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  T(LDT,*),C(LDC,*)
*     ..
*     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER(ZERO = 0.0D+0, ONE = 1.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER           L,NX
      LOGICAL           UPPERT,UPPERC,TRANSL,UNITT
*     ..
*     .. External Subroutines ..
      EXTERNAL          DSYRK,DTRMMOOP
*     ..
*     .. External Functions ..
      INTEGER           ILAENV
      LOGICAL           LSAME
      EXTERNAL          LSAME,ILAENV
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF(N.EQ.0) THEN
         RETURN
      END IF
*
*     Determine the crossover point into the unblocked variant
*
      NX = ILAENV(3, 'DST3RK', UPLOT // UPLOC // TRANS // DIAG,
     $      K, -1, -1, -1)
*
      IF(K.LT.NX) THEN
         CALL DSTRK(UPLOT, UPLOC, TRANS, UNIT, K, ALPHA, T, LDT,
     $         BETA, C, LDC)
         RETURN
      END IF
*
*     Convert our character inputs into logical variables
*
      UPPERT = LSAME(UPLOT,'U')
      UPPERC = LSAME(UPLOC,'U')
      TRANSL = LSAME(TRANS,'T')
      UNITT = LSAME(DIAG,'U')
*
*     Base case
*
      IF(N.EQ.1) THEN
         IF(BETA.EQ.ZERO) THEN
            C(1,1) = ZERO
         ELSE
            C(1,1) = BETA*C(1,1)
         END IF
         IF(UNITT) THEN
            C(1,1) = ALPHA + C(1,1)
         ELSE
            TMP = T(1,1)*T(1,1)
            C(1,1) = ALPHA*TMP + C(1,1)
         END IF
         RETURN
      END IF
*
*     Recursive case
*
      L = K/2
*
      IF(TRANSL) THEN
*
*        This means we are computing C = alpha*T**H*T + beta*C
*
         IF(UPPERT) THEN
*
*           This means T is upper triangular
*
*           Break C and T apart as follows
*               |-----------------|
*           T = | T_{1,1} T_{1,2} | l
*               | 0       T_{2,2} | k-l
*               |-----------------|
*                 l       k-l
*
*               |-----------------|
*           C = | C_{1,1} C_{1,2} | l
*               | C_{2,1} C_{2,2} | k-l
*               |-----------------|
*                 l       k-l
*
*           So, we are computing
*                      |-----------------------| |-----------------|         |-----------------|
*           C = \alpha | T_{1,1}**H 0          |*| T_{1,1} T_{1,2} | + \beta | C_{1,1} C_{1,2} |
*                      | T_{1,2}**H T_{2,2}**H | | 0       T_{2,2} |         | C_{2,1} C_{2,2} |
*                      |-----------------------| |-----------------|         |-----------------|
*
*           Which gives us the following componentwise representation of C
*
*           C_{1,1} = \alpha*T_{1,1}**H*T_{1,1} + \beta C_{1,1}
*           C_{1,2} = \alpha*T_{1,1}**H*T_{1,2} + \beta C_{1,2}
*           C_{2,1} = \alpha*T_{1,2}**H*T_{1,1} + \beta C_{2,1}
*           C_{2,2} = \alpha*T_{1,2}**H*T_{1,2} + \alpha*T_{2,2}**H*T_{2,2} + \beta C_{2,2}
*
*           Thus, we compute the following 
*
*           C_{1,1} = \alpha*T_{1,1}**H*T_{1,1} + \beta C_{1,1} (This routine)
*           C_{2,2} = \alpha*T_{2,2}**H*T_{2,2} + \beta C_{2,2} (This routine)
*           C_{2,2} = \alpha*T_{1,2}**H*T_{1,2} + C_{2,2}       (SYRK)
*
*           Compute C_{1,1} = \alpha*T_{1,1}**H*T_{1,1} + \beta C_{1,1}
*
            CALL DST3RK(UPLOT, UPLOC, TRANS, DIAG, L, ALPHA,
     $            T, LDT, BETA, C, LDC)
*
*           Compute C_{2,2}
*           C_{2,2} = \alpha*T_{2,2}**H*T_{2,2} + \beta C_{2,2}
*
            CALL DST3RK(UPLOT, UPLOC, TRANS, DIAG, K-L, ALPHA,
     $            T(L+1,L+1), LDT, BETA, C(L+1,L+1), LDC)
*
*           C_{2,2} = \alpha*T_{1,2}**H*T_{1,2} + C_{2,2}
*
            CALL DSYRK(UPLOC, TRANS, K-L, L, ALPHA, T(1,L+1), LDT,
     $            ONE, C(L+1,L+1), LDC)
            IF(UPPERC) THEN
*
*              Compute C_{1,2} = \alpha*T_{1,1}**H*T_{1,2} + \beta C_{1,2} (TRMMOOP)
*
               CALL DTRMMOOP('Left', UPLOT, 'Transpose',
     $               'No Transpose', DIAG, L, K-L, ALPHA, T, LDT,
     $               T(1, L+1), LDT, BETA, C(1,L+1), LDC)
            ELSE
*
*              Compute C_{2,1} = \alpha*T_{1,2}**H*T_{1,1} + \beta C_{2,1} (TRMMOOP)
*
               CALL DTRMMOOP('Right', UPLOT, 'No Transpose',
     $               'Transpose', DIAG, K-L, L, ALPHA, T, LDT,
     $               T(1, L+1), LDT, BETA, C(L+1,1), LDC)
            END IF
         ELSE
*
*           This means T is lower triangular
*
*           Break C and T apart as follows
*               |-----------------|
*           T = | T_{1,1} 0       | l
*               | T_{2,1} T_{2,2} | k-l
*               |-----------------|
*                 l       k-l
*
*               |-----------------|
*           C = | C_{1,1} C_{1,2} | l
*               | C_{2,1} C_{2,2} | k-l
*               |-----------------|
*                 l       k-l
*
*           So, we are computing
*                      |-----------------------| |-----------------|         |-----------------|
*           C = \alpha | T_{1,1}**H T_{2,1}**H |*| T_{1,1} 0       | + \beta | C_{1,1} C_{1,2} |
*                      | 0          T_{2,2}**H | | T_{2,1} T_{2,2} |         | C_{2,1} C_{2,2} |
*                      |-----------------------| |-----------------|         |-----------------|
*
*           Which gives us the following componentwise representation of C
*
*           C_{1,1} = \alpha*T_{1,1}**H*T_{1,1} + \alpha*T_{2,1}**H*T_{2,1} + \beta*C_{1,1}
*           C_{1,2} = \alpha*T_{2,1}**H*T_{2,2} + \beta*C_{1,2}
*           C_{2,1} = \alpha*T_{2,2}**H*T_{2,1} + \beta*C_{2,1}
*           C_{2,2} = \alpha*T_{2,2}**H*T_{2,2} + \beta*C_{2,2}
*
*           Thus, we compute the following 
*
*           C_{1,1} = \alpha*T_{1,1}**H*T_{1,1} + \beta*C_{1,1} (This routine)
*           C_{1,1} = \alpha*T_{2,1}**H*T_{2,1} + C_{1,1}       (SYRK)
*           C_{2,2} = \alpha*T_{2,2}**H*T_{2,2} + \beta*C_{2,2} (This routine)
*
*           Compute C_{1,1}
*           C_{1,1} = \alpha*T_{1,1}**H*T_{1,1} + \beta*C_{1,1}
*
            CALL DST3RK(UPLOT, UPLOC, TRANS, DIAG, L, ALPHA,
     $            T, LDT, BETA, C, LDC)
*
*           C_{1,1} = \alpha*T_{2,1}**H*T_{2,1} + C_{1,1}
*
            CALL DSYRK(UPLOC, TRANS, L, K-L, ALPHA, T(L+1,1), LDT,
     $            ONE, C, LDC)
            IF(UPPERC) THEN
*
*              Compute C_{1,2} = \alpha*T_{2,1}**H*T_{2,2} + \beta*C_{1,2}
*
               CALL DTRMMOOP('Right', UPLOT, 'No Transpose',
     $               'Transpose', DIAG, L, K-L, ALPHA, T(L+1,L+1), LDT,
     $                T(L+1,1), LDT, BETA, C(1,L+1), LDC)
            ELSE
*
*              Compute C_{2,1} = \alpha*T_{2,2}**H*T_{2,1} + \beta*C_{2,1}
*
               CALL DTRMMOOP('Left', UPLOT, 'Transpose',
     $               'No Transpose', DIAG, K-L, L, ALPHA,
     $               T(L+1,L+1), LDT, T(L+1, 1), LDT, BETA,
     $               C(L+1,1), LDC)
            END IF
         END IF
      ELSE
*
*        This means we are computing C = alpha*T*T**H + beta*C
*
         IF(UPPERT) THEN
*
*           This means T is upper triangular
*
*           Break C and T apart as follows
*               |-----------------|
*           T = | T_{1,1} T_{1,2} | l
*               | 0       T_{2,2} | k-l
*               |-----------------|
*                 l       k-l
*
*               |-----------------|
*           C = | C_{1,1} C_{1,2} | l
*               | C_{2,1} C_{2,2} | k-l
*               |-----------------|
*                 l       k-l
*
*           So, we are computing
*                      |-----------------| |-----------------------|         |-----------------|
*           C = \alpha | T_{1,1} T_{1,2} |*| T_{1,1}**H 0          | + \beta | C_{1,1} C_{1,2} |
*                      | 0       T_{2,2} | | T_{1,2}**H T_{2,2}**H |         | C_{2,1} C_{2,2} |
*                      |-----------------| |-----------------------|         |-----------------|
*
*           Which gives us the following componentwise representation of C
*
*           C_{1,1} = \alpha*T_{1,1}*T_{1,1}**H + \alpha*T_{1,2}*T_{1,2}**H + \beta*C_{1,1}
*           C_{1,2} = \alpha*T_{1,2}*T_{2,2}**H + \beta*C_{1,2}
*           C_{2,1} = \alpha*T_{2,2}*T_{1,2}**H + \beta*C_{2,1}
*           C_{2,2} = \alpha*T_{2,2}*T_{2,2}**H + \beta*C_{2,2}
*
*           Thus, we compute the following 
*
*           C_{1,1} = \alpha*T_{1,1}*T_{1,1}**H + \beta*C_{1,1} (This routine)
*           C_{1,1} = \alpha*T_{1,2}*T_{1,2}**H + C_{1,1}       (SYRK)
*           C_{2,2} = \alpha*T_{2,2}*T_{2,2}**H + \beta*C_{2,2} (This routine)
*
*           Compute C_{1,1}
*           C_{1,1} = \alpha*T_{1,1}*T_{1,1}**H + \beta*C_{1,1}
*
            CALL DST3RK(UPLOT, UPLOC, TRANS, DIAG, L, ALPHA,
     $            T, LDT, BETA, C, LDC)
*
*           C_{1,1} = \alpha*T_{1,2}*T_{1,2}**H + C_{1,1}
*
            CALL DSYRK(UPLOC, TRANS, L, K-L, ALPHA, T(1,L+1), LDT,
     $            ONE, C, LDC)
*
*           Compute C_{2,2} = \alpha*T_{2,2}*T_{2,2}**H + \beta*C_{2,2}
*
            CALL DST3RK(UPLOT, UPLOC, TRANS, DIAG, K-L, ALPHA,
     $            T(L+1,L+1), LDT, BETA, C(L+1,L+1), LDC)
            IF(UPPERC) THEN
*
*              Compute C_{1,2} = \alpha*T_{1,2}*T_{2,2}**H + \beta*C_{1,2}
*
               CALL DTRMMOOP('Right', UPLOT, 'Transpose',
     $               'No Transpose', DIAG, L, K-L, ALPHA,
     $               T(L+1,L+1), LDT, T(1, L+1), LDT, BETA,
     $               C(1,L+1), LDC)
            ELSE
*
*              Compute C_{2,1} = \alpha*T_{2,2}*T_{1,2}**H + \beta*C_{2,1}
*
               CALL DTRMMOOP('Left', UPLOT, 'No Transpose',
     $               'Transpose', DIAG, K-L, L, ALPHA,
     $               T(L+1,L+1), LDT, T(1, L+1), LDT, BETA,
     $               C(L+1,1), LDC)
            END IF
         ELSE
*
*           This means T is lower triangular
*
*           Break C and T apart as follows
*               |-----------------|
*           T = | T_{1,1} 0       | l
*               | T_{2,1} T_{2,2} | k-l
*               |-----------------|
*                 l       k-l
*
*               |-----------------|
*           C = | C_{1,1} C_{1,2} | l
*               | C_{2,1} C_{2,2} | k-l
*               |-----------------|
*                 l       k-l
*
*           So, we are computing
*                      |-----------------| |-----------------------|         |-----------------|
*           C = \alpha | T_{1,1} 0       |*| T_{1,1}**H T_{2,1}**H | + \beta | C_{1,1} C_{1,2} |
*                      | T_{2,1} T_{2,2} | | 0          T_{2,2}**H |         | C_{2,1} C_{2,2} |
*                      |-----------------| |-----------------------|         |-----------------|
*
*           Which gives us the following componentwise representation of C
*
*           C_{1,1} = \alpha*T_{1,1}*T_{1,1}**H + \beta*C_{1,1}
*           C_{1,2} = \alpha*T_{1,1}*T_{2,1}**H + \beta*C_{1,2}
*           C_{2,1} = \alpha*T_{2,1}*T_{1,1}**H + \beta*C_{2,1}
*           C_{2,2} = \alpha*T_{2,1}*T_{2,1}**H + \alpha*T_{2,2}*T_{2,2}**H + \beta*C_{2,2}
*
*           Thus, we compute the following 
*
*           C_{1,1} = \alpha*T_{1,1}*T_{1,1}**H + \beta*C_{1,1} (This routine)
*           C_{2,2} = \alpha*T_{2,2}*T_{2,2}**H + \beta*C_{2,2} (This routine)
*           C_{2,2} = \alpha*T_{2,1}*T_{2,1}**H + C_{2,2}       (SYRK)
*
*           Compute C_{1,1} = \alpha*T_{1,1}*T_{1,1}**H + \beta*C_{1,1}
*
            CALL DST3RK(UPLOT, UPLOC, TRANS, DIAG, L, ALPHA, T, LDT,
     $            BETA, C, LDC)
*
*           Compute C_{2,2}
*           C_{2,2} = \alpha*T_{2,2}*T_{2,2}**H + \beta*C_{2,2}
*
            CALL DST3RK(UPLOT, UPLOC, TRANS, DIAG, K-L, ALPHA,
     $            T(L+1,L+1), LDT, BETA, C(L+1,L+1), LDC)
*
*           C_{2,2} = \alpha*T_{2,1}*T_{2,1}**H + C_{2,2}
*
            CALL DSYRK(UPLOC, TRANS, K-L, L, ALPHA, T(L+1,1), LDT,
     $            ONE, C(L+1,L+1), LDC)
            IF(UPPERC) THEN
*
*              Compute C_{1,2} = \alpha*T_{1,1}*T_{2,1}**H + \beta*C_{1,2}
*
               CALL DTRMMOOP('Left', UPLOT, 'No Transpose',
     $               'Transpose', DIAG, L, K-L, ALPHA, T, LDT,
     $               T(L+1,1), LDT, BETA, C(1,L+1), LDC)
            ELSE
*
*              Compute C_{2,1} = \alpha*T_{2,1}*T_{1,1}**H + \beta*C_{2,1}
*
               CALL DTRMMOOP('Right', UPLOT, 'Transpose',
     $               'No Transpose', DIAG, K-L, L, ALPHA, T, LDT,
     $               T(L+1,1), LDT, BETA, C(L+1,1), LDC)
            END IF
         END IF
      END IF
      END SUBROUTINE
