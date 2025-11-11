*> \brief \b DSTRK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE DSTRK(UPLOA, UPLOC, TRANS, DIAG, K, ALPHA, A, LDA,
*    $            BETA, C, LDC)
*
*     .. Scalar Arguments ..
*     DOUBLE PRECISION  ALPHA,BETA
*     INTEGER           K,LDA,LDC
*     CHARACTER         UPLOA,UPLOC,TRANS,DIAG
*     ..
*     .. Array Arguments ..
*     DOUBLE PRECISION  A(LDA,*),C(LDC,*)
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSTRK  performs one of the symmetric rank k operations
*>
*>    C := alpha*A*A**T + beta*C,
*>
*> or
*>
*>    C := alpha*A**T*A + beta*C,
*>
*> where  alpha and beta  are scalars, C is a  k by k  symmetric matrix
*> and  A  is an  k by k  either upport or lower triangular matrix 
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLOA
*> \verbatim
*>          UPLOA is CHARACTER*1
*>           On  entry,   UPLOA specifies  whether  the  upper  or  lower
*>           triangular  part  of the  array  A  is to be  referenced  as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  A
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  A
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] UPLOC
*> \verbatim
*>          UPLOC is CHARACTER*1
*>           On  entry,   UPLOC specifies  whether  the  upper  or  lower
*>           triangular  part  of the  array  C  is to be  referenced  as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry,  TRANS  specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
*>
*>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
*>
*>              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry,  K specifies the rows and columns of the matrix C.
*>           K must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, k ).
*>          If UPLOA = 'U' or 'u', then the leading k by k upper triangular
*>          part of the array A must contain the upper triangular part of
*>          the triangular matrix, and the strictly lower triangular part of A
*>          is not referenced. If UPLOA = 'L' or 'l', then the leading k by k
*>          lower triangular part of the array A must contain the lower
*>          part of the triangular matrix, and the strictly upper triangular
*>          part of A is not referenced. If DIAG = 'U', then the diagonal
*>          component of A is not referenced and instead assumed to be unit
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in  the  calling  (sub)  program. LDA must be at least max( 1, k ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension ( LDC, N )
*>           Before entry  with  UPLOC = 'U' or 'u',  the leading  k by k
*>           upper triangular part of the array C must contain the upper
*>           triangular part  of the  symmetric matrix  and the strictly
*>           lower triangular part of C is not referenced.  On exit, the
*>           upper triangular part of the array  C is overwritten by the
*>           upper triangular part of the updated matrix.
*>           Before entry  with  UPLOC = 'L' or 'l',  the leading  n by n
*>           lower triangular part of the array C must contain the lower
*>           triangular part  of the  symmetric matrix  and the strictly
*>           upper triangular part of C is not referenced.  On exit, the
*>           lower triangular part of the array  C is overwritten by the
*>           lower triangular part of the updated matrix.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, k ).
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
*> \ingroup herk
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSTRK(UPLOA, UPLOC, TRANS, DIAG, K, ALPHA, A, LDA,
     $            BETA, C, LDC)
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA,BETA
      INTEGER           K,LDA,LDC
      CHARACTER         UPLOA,UPLOC,TRANS,DIAG
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*),C(LDC,*)
*     ..
*     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER(ZERO = 0.0D+0, ONE = 1.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION  TMP
      INTEGER           I, L
      LOGICAL           UPPERA,UPPERC,TRANSL,UNITT
*     ..
*     .. External Subroutines ..
      EXTERNAL          DTRMMOOP
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
      UPPERA = LSAME(UPLOA,'U')
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
            IF (UPPERA) THEN
*
*              This means T is upper triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOA, 'Transpose',
     $                     'No Transpose', DIAG, I-1, 1, ALPHA, A, LDA, 
     $                     A(1,I), LDA, BETA, C(1,I), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(I-1, A(1,I), 1, A(1,I), 1) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOA, 'Transpose',
     $                     'No Transpose', DIAG, I, 1, ALPHA, A, LDA, 
     $                     A(1,I), LDA, BETA, C(1,I), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                 DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOA, 'No Transpose',
     $                     'Transpose', DIAG, 1, K-I, ALPHA,
     $                     A(I+1,I+1), LDA, A(I+1,I), LDA, BETA,
     $                     C(I,I+1), LDC)
*
                    C(I,I) = ALPHA * 
     $                  (DDOT(K-I, A(I+1,I), 1, A(I+1,I), 1) +
     $                  ONE) + BETA*C(I,I)
                 END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOA, 'No Transpose',
     $                     'Transpose', DIAG, 1, K-I+1, ALPHA,
     $                     A(I,I), LDA, A(I,I), LDA, BETA, C(I,I), LDC)
                  END DO
               END IF
            END IF
         ELSE
*
*           This means we are only storing the lower triangular component of C
*
            IF (UPPERA) THEN
*
*              This means T is upper triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOA, 'No Transpose', 
     $                     'Transpose', DIAG, 1, I-1, ALPHA, A, LDA, 
     $                     A(1,I), LDA, BETA, C(I,1), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(I-1, A(1,I), 1, A(1,I), 1) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOA, 'No Transpose', 
     $                     'Transpose', DIAG, 1, I, ALPHA, A, LDA, 
     $                     A(1,I), LDA, BETA, C(I,1), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOA, 'Transpose', 
     $                     'No Transpose', DIAG, K-I, 1, ALPHA,
     $                     A(I+1,I+1), LDA, A(I+1,I), LDA, BETA,
     $                     C(I+1,I), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(K-I, A(I+1,I), 1, A(I+1,I), 1) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOA, 'Transpose', 
     $                     'No Transpose', DIAG, K-I+1, 1, ALPHA,
     $                     A(I,I), LDA, A(I,I), LDA, BETA, C(I,I), LDC)
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
            IF (UPPERA) THEN
*
*              This means T is upper triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOA, 'Transpose',
     $                     'No Transpose', DIAG, 1, K-I, ALPHA,
     $                     A(I+1,I+1), LDA, A(I,I+1), LDA, BETA,
     $                     C(I,I+1), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(K-I, A(I,I+1), LDA, A(I,I+1), LDA)
     $                  + ONE) + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOA, 'Transpose',
     $                     'No Transpose', DIAG, 1, K-I+1, ALPHA,
     $                     A(I,I), LDA, A(I,I), LDA, BETA, C(I,I), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOA, 'No Transpose',
     $                     'Transpose', DIAG, I-1, 1, ALPHA, A, LDA,
     $                     A(I,1), LDA, BETA, C(1,I), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(I-1, A(I,1), LDA, A(I,1), LDA) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOA, 'No Transpose',
     $                     'Transpose', DIAG, I, 1, ALPHA, A, LDA,
     $                     A(I,1), LDA, BETA, C(1,I), LDC)
                  END DO
               END IF
            END IF
         ELSE
*
*           This means we are only storing the lower triangular component of C
*
            IF (UPPERA) THEN
*
*              This means T is upper triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOA, 'No Transpose',
     $                     'Transpose', DIAG, K-I, 1, ALPHA,
     $                     A(I+1,I+1), LDA, A(I,I+1), LDA, BETA,
     $                     C(I+1,I), LDC)
*
                     C(I,I) = ALPHA *
     $                  (DDOT(K-I, A(I,I+1), LDA, A(I,I+1), LDA)
     $                  + ONE) + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Left', UPLOA, 'No Transpose',
     $                     'Transpose', DIAG, K-I+1, 1, ALPHA,
     $                     A(I,I), LDA, A(I,I), LDA, BETA, C(I,I), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOA, 'Transpose',
     $                     'No Transpose', DIAG, 1, I-1, ALPHA, A, LDA,
     $                     A(I,1), LDA, BETA, C(I,1), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DDOT(I-1, A(I,1), LDA, A(I,1), LDA) + ONE)
     $                  + BETA*C(I,I)
                  END DO
               ELSE
                  DO I = 1, K
                     CALL DTRMMOOP('Right', UPLOA, 'Transpose',
     $                     'No Transpose', DIAG, 1, I, ALPHA, A, LDA,
     $                     A(I,1), LDA, BETA, C(I,1), LDC)
                  END DO
               END IF
            END IF
         END IF
      END IF
      END SUBROUTINE
