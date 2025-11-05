! syrk but A is triangular. This is the base case we default to in dst3rk
      SUBROUTINE DSTRK(UPLOT, UPLOC, TRANS, DIAG, K, ALPHA, T, LDT,
     $            BETA, C, LDC)
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA,BETA
      INTEGER           K,LDT,LDC
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
      DOUBLE PRECISION  TMP
      INTEGER           I, L
      LOGICAL           UPPERT,UPPERC,TRANSL,UNITT
*     ..
*     .. External Subroutines ..
!     EXTERNAL          DSYRK,DTRMMOOP
*     ..
*     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DDOT
      EXTERNAL          LSAME,DDOT
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF(K.EQ.0) THEN
         RETURN
      END IF
*
*     Convert our character inputs into logical variables
*
      UPPERT = LSAME(UPLOT,'U')
      UPPERC = LSAME(UPLOC,'U')
      TRANSL = LSAME(TRANS,'T').OR.LSAME(TRANS,'C')
      UNITT = LSAME(DIAG,'U')
*
      IF(TRANSL) THEN
*
*        This means we are computing C = \alpha*T**H*T + \beta*C
*
         IF(UPPERC) THEN
*
*           This means we are only storing the upper triangular component of C
*
            IF (UPPERT) THEN
*
*              This means T is upper triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOT, 'Transpose',
     $                     'No Transpose', DIAG, I-1, 1, ALPHA, T, LDT, 
     $                     T(1,I), LDT, BETA, C(1,I), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(I-1, T(1,I), 1, T(1,I), 1) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOT, 'Transpose',
     $                     'No Transpose', DIAG, I, 1, ALPHA, T, LDT, 
     $                     T(1,I), LDT, BETA, C(1,I), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                 DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOT, 'No Transpose',
     $                     'Transpose', DIAG, 1, K-I, ALPHA,
     $                     T(I+1,I+1), LDT, T(I+1,I), LDT, BETA,
     $                     C(I,I+1), LDC)
*
                    C(I,I) = ALPHA * 
     $                  (DDOT(K-I, T(I+1,I), 1, T(I+1,I), 1) +
     $                  ONE) + BETA*C(I,I)
                 END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOT, 'No Transpose',
     $                     'Transpose', DIAG, 1, K-I+1, ALPHA,
     $                     T(I,I), LDT, T(I,I), LDT, BETA, C(I,I), LDC)
                  END DO
               END IF
            END IF
         ELSE
*
*           This means we are only storing the lower triangular component of C
*
            IF (UPPERT) THEN
*
*              This means T is upper triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOT, 'No Transpose', 
     $                     'Transpose', DIAG, 1, I-1, ALPHA, T, LDT, 
     $                     T(1,I), LDT, BETA, C(I,1), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(I-1, T(1,I), 1, T(1,I), 1) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOT, 'No Transpose', 
     $                     'Transpose', DIAG, 1, I, ALPHA, T, LDT, 
     $                     T(1,I), LDT, BETA, C(I,1), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOT, 'Transpose', 
     $                     'No Transpose', DIAG, K-I, 1, ALPHA,
     $                     T(I+1,I+1), LDT, T(I+1,I), LDT, BETA,
     $                     C(I+1,I), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(K-I, T(I+1,I), 1, T(I+1,I), 1) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOT, 'Transpose', 
     $                     'No Transpose', DIAG, K-I+1, 1, ALPHA,
     $                     T(I,I), LDT, T(I,I), LDT, BETA, C(I,I), LDC)
                  END DO
               END IF
            END IF
         END IF
      ELSE
*
*        This means we are computing C = \alpha*T*T**H + \beta*C
*
         IF(UPPERC) THEN
*
*           This means we are only storing the upper triangular component of C
*
            IF (UPPERT) THEN
*
*              This means T is upper triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOT, 'Transpose',
     $                     'No Transpose', DIAG, 1, K-I, ALPHA,
     $                     T(I+1,I+1), LDT, T(I,I+1), LDT, BETA,
     $                     C(I,I+1), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(K-I, T(I,I+1), LDT, T(I,I+1), LDT)
     $                  + ONE) + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOT, 'Transpose',
     $                     'No Transpose', DIAG, 1, K-I+1, ALPHA,
     $                     T(I,I), LDT, T(I,I), LDT, BETA, C(I,I), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOT, 'No Transpose',
     $                     'Transpose', DIAG, I-1, 1, ALPHA, T, LDT,
     $                     T(I,1), LDT, BETA, C(1,I), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(I-1, T(I,1), LDT, T(I,1), LDT) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOT, 'No Transpose',
     $                     'Transpose', DIAG, I, 1, ALPHA, T, LDT,
     $                     T(I,1), LDT, BETA, C(1,I), LDC)
                  END DO
               END IF
            END IF
         ELSE
*
*           This means we are only storing the lower triangular component of C
*
            IF (UPPERT) THEN
*
*              This means T is upper triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOT, 'No Transpose',
     $                     'Transpose', DIAG, K-I, 1, ALPHA,
     $                     T(I+1,I+1), LDT, T(I,I+1), LDT, BETA,
     $                     C(I+1,I), LDC)
*
                     C(I,I) = ALPHA *
     $                  (DDOT(K-I, T(I,I+1), LDT, T(I,I+1), LDT)
     $                  + ONE) + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOT, 'No Transpose',
     $                     'Transpose', DIAG, K-I+1, 1, ALPHA,
     $                     T(I,I), LDT, T(I,I), LDT, BETA, C(I,I), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOT, 'Transpose',
     $                     'No Transpose', DIAG, 1, I-1, ALPHA, T, LDT,
     $                     T(I,1), LDT, BETA, C(I,1), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(I-1, T(I,1), LDT, T(I,1), LDT) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOT, 'Transpose',
     $                     'No Transpose', DIAG, 1, I, ALPHA, T, LDT,
     $                     T(I,1), LDT, BETA, C(I,1), LDC)
                  END DO
               END IF
            END IF
         END IF
      END IF
      END SUBROUTINE
