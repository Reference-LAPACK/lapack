*> \brief \b ZHTRK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE ZHTRK(UPLOA, UPLOC, TRANS, DIAG, K, ALPHA, A, LDA,
*    $            BETA, C, LDC)
*
*     .. Scalar Arguments ..
*     DOUBLE PRECISION  ALPHA,BETA
*     INTEGER           K,LDA,LDC
*     CHARACTER         UPLOA,UPLOC,TRANS,DIAG
*     ..
*     .. Array Arguments ..
*     COMPLEX*16        A(LDA,*),C(LDC,*)
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZHTRK  performs one of the hermitian rank k operations
*>
*>    C := alpha*A*A**H + beta*C,
*>
*> or
*>
*>    C := alpha*A**H*A + beta*C,
*>
*> where  alpha and beta  are real scalars, C is a  k by k  hermitian matrix
*> and  A  is an  k by k  either upper or or lower triangular matrix 
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
*>              TRANS = 'N' or 'n'   C := alpha*A*A**H + beta*C.
*>
*>              TRANS = 'C' or 'c'   C := alpha*A**H*A + beta*C.
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
*>          ALPHA is REAL.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, k ).
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
*>          BETA is REAL.
*>           On entry, BETA specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension ( LDC, N )
*>           Before entry  with  UPLOC = 'U' or 'u',  the leading  k by k
*>           upper triangular part of the array C must contain the upper
*>           triangular part  of the  hermitian matrix  and the strictly
*>           lower triangular part of C is not referenced.  On exit, the
*>           upper triangular part of the array  C is overwritten by the
*>           upper triangular part of the updated matrix.
*>           Before entry  with  UPLOC = 'L' or 'l',  the leading  n by n
*>           lower triangular part of the array C must contain the lower
*>           triangular part  of the  hermitian matrix  and the strictly
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
      SUBROUTINE ZHTRK(UPLOA, UPLOC, TRANS, DIAG, K, ALPHA, A, LDA,
     $            BETA, C, LDC)
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA,BETA
      INTEGER           K,LDA,LDC
      CHARACTER         UPLOA,UPLOC,TRANS,DIAG
*     ..
*     .. Array Arguments ..
      COMPLEX*16        A(LDA,*),C(LDC,*)
*     ..
*     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER(ONE = 1.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER           I, L
      LOGICAL           UPPERA,UPPERC,TRANSL,UNITT
*     ..
*     .. External Subroutines ..
      EXTERNAL          ZTRMMOOP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC         DCMPLX,DBLE
*     ..
*     .. External Functions ..
      LOGICAL           LSAME
      COMPLEX*16        ZDOTC
      EXTERNAL          LSAME,ZDOTC
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
      TRANSL = LSAME(TRANS,'C')
      UNITT = LSAME(DIAG,'U')
      ! Consider erroring if (.NOT.TRANSL).AND.(.NOT.LSAME(TRANS,'N'))
      ! Also consider erroring if (.NOT.UPPERA).AND.(.NOT.LSAME(UPLOA,'L'))
      ! Also consider erroring if (.NOT.UPPERC).AND.(.NOT.LSAME(UPLOC,'L'))
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
                     CALL ZTRMMOOP('Left', UPLOA, 'Conjugate',
     $                     'No Transpose', DIAG, I-1, 1, DCMPLX(ALPHA), 
     $                     A, LDA, A(1,I), LDA, DCMPLX(BETA),
     $                     C(1,I), LDC)
*
                     C(I,I) = ALPHA * 
     $                  (DBLE(ZDOTC(I-1, A(1,I), 1, A(1,I), 1))
     $                  + ONE) + BETA*DBLE(C(I,I))
                  END DO
               ELSE
                  DO I = 1, K
                     CALL ZTRMMOOP('Left', UPLOA, 'Conjugate',
     $                     'No Transpose', DIAG, I, 1, DCMPLX(ALPHA), A,
     $                     LDA, A(1,I), LDA, DCMPLX(BETA), C(1,I), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                 DO I = 1, K
                     CALL ZTRMMOOP('Right', UPLOA, 'No Transpose',
     $                     'Conjugate', DIAG, 1, K-I, DCMPLX(ALPHA),
     $                     A(I+1,I+1), LDA, A(I+1,I), LDA, DCMPLX(BETA),
     $                     C(I,I+1), LDC)
 
                    C(I,I) = ALPHA * 
     $                  (DBLE(ZDOTC(K-I, A(I+1,I), 1, A(I+1,I), 1))
     $                  + ONE) + BETA*DBLE(C(I,I))
                 END DO
               ELSE
                  DO I = 1, K
                     CALL ZTRMMOOP('Right', UPLOA, 'No Transpose',
     $                     'Conjugate', DIAG, 1, K-I+1, DCMPLX(ALPHA),
     $                     A(I,I), LDA, A(I,I), LDA, DCMPLX(BETA),
     $                     C(I,I), LDC)
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
                     CALL ZTRMMOOP('Right', UPLOA, 'No Transpose', 
     $                     'Conjugate', DIAG, 1, I-1, DCMPLX(ALPHA), 
     $                     A, LDA, A(1,I), LDA, DCMPLX(BETA),
     $                     C(I,1), LDC)
 
                     C(I,I) = ALPHA * 
     $                  (DBLE(ZDOTC(I-1, A(1,I), 1, A(1,I), 1))
     $                  + ONE) + BETA*DBLE(C(I,I))
                  END DO
               ELSE
                  DO I = 1, K
                     CALL ZTRMMOOP('Right', UPLOA, 'No Transpose', 
     $                     'Conjugate', DIAG, 1, I, DCMPLX(ALPHA),
     $                     A, LDA, A(1,I), LDA, DCMPLX(BETA),
     $                     C(I,1), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL ZTRMMOOP('Left', UPLOA, 'Conjugate', 
     $                     'No Transpose', DIAG, K-I, 1, DCMPLX(ALPHA),
     $                     A(I+1,I+1), LDA, A(I+1,I), LDA, DCMPLX(BETA),
     $                     C(I+1,I), LDC)
 
                     C(I,I) = ALPHA * 
     $                  (DBLE(ZDOTC(K-I, A(I+1,I), 1, A(I+1,I), 1))
     $                  + ONE) + BETA*DBLE(C(I,I))
                  END DO
               ELSE
                  DO I = 1, K
                     CALL ZTRMMOOP('Left', UPLOA, 'Conjugate', 
     $                     'No Transpose', DIAG, K-I+1, 1,
     $                     DCMPLX(ALPHA), A(I,I), LDA, A(I,I), LDA,
     $                     DCMPLX(BETA), C(I,I), LDC)
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
                     CALL ZTRMMOOP('Right', UPLOA, 'Conjugate',
     $                     'No Transpose', DIAG, 1, K-I, DCMPLX(ALPHA),
     $                     A(I+1,I+1), LDA, A(I,I+1), LDA, DCMPLX(BETA),
     $                     C(I,I+1), LDC)
 
                     C(I,I) = ALPHA * 
     $                  (DBLE(ZDOTC(K-I, A(I,I+1), LDA,
     $                  A(I,I+1), LDA)) + ONE) + BETA*DBLE(C(I,I))
                  END DO
               ELSE
                  DO I = 1, K
                     CALL ZTRMMOOP('Right', UPLOA, 'Conjugate',
     $                     'No Transpose', DIAG, 1, K-I+1,
     $                     DCMPLX(ALPHA), A(I,I), LDA, A(I,I), LDA,
     $                     DCMPLX(BETA), C(I,I), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL ZTRMMOOP('Left', UPLOA, 'No Transpose',
     $                     'Conjugate', DIAG, I-1, 1, DCMPLX(ALPHA),
     $                     A, LDA, A(I,1), LDA, DCMPLX(BETA),
     $                     C(1,I), LDC)
 
                     C(I,I) = ALPHA * 
     $                  (DBLE(ZDOTC(I-1, A(I,1), LDA, A(I,1), LDA))
     $                  + ONE) + BETA*DBLE(C(I,I))
                  END DO
               ELSE
                  DO I = 1, K
                     CALL ZTRMMOOP('Left', UPLOA, 'No Transpose',
     $                     'Conjugate', DIAG, I, 1, DCMPLX(ALPHA),
     $                     A, LDA, A(I,1), LDA, DCMPLX(BETA),
     $                     C(1,I), LDC)
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
                     CALL ZTRMMOOP('Left', UPLOA, 'No Transpose',
     $                     'Conjugate', DIAG, K-I, 1, DCMPLX(ALPHA),
     $                     A(I+1,I+1), LDA, A(I,I+1), LDA, DCMPLX(BETA),
     $                     C(I+1,I), LDC)
 
                     C(I,I) = ALPHA *
     $                  (DBLE(ZDOTC(K-I, A(I,I+1), LDA,
     $                  A(I,I+1), LDA)) + ONE) + BETA*DBLE(C(I,I))
                  END DO
               ELSE
                  DO I = 1, K
                     CALL ZTRMMOOP('Left', UPLOA, 'No Transpose',
     $                     'Conjugate', DIAG, K-I+1, 1, DCMPLX(ALPHA),
     $                     A(I,I), LDA, A(I,I), LDA, DCMPLX(BETA),
     $                     C(I,I), LDC)
                  END DO
               END IF
            ELSE 
*
*              This means T is lower triangular
*
               IF(UNITT) THEN
                  DO I = 1, K
                     CALL ZTRMMOOP('Right', UPLOA, 'Conjugate',
     $                     'No Transpose', DIAG, 1, I-1, DCMPLX(ALPHA),
     $                     A, LDA, A(I,1), LDA, DCMPLX(BETA),
     $                     C(I,1), LDC)
 
                     C(I,I) = ALPHA * 
     $                  (DBLE(ZDOTC(I-1, A(I,1), LDA, A(I,1), LDA))
     $                  + ONE) + BETA*DBLE(C(I,I))
                  END DO
               ELSE
                  DO I = 1, K
                     CALL ZTRMMOOP('Right', UPLOA, 'Conjugate',
     $                     'No Transpose', DIAG, 1, I, DCMPLX(ALPHA),
     $                     A, LDA, A(I,1), LDA, DCMPLX(BETA),
     $                     C(I,1), LDC)
                  END DO
               END IF
            END IF
         END IF
      END IF
      END SUBROUTINE
