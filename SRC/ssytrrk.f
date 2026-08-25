*> \brief \b SSYTRRK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     RECURSIVE SUBROUTINE SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, K,
*    $            ALPHA, A, LDA, BETA, C, LDC)
*
*     .. Scalar Arguments ..
*     REAL  ALPHA,BETA
*     INTEGER           K,LDA,LDC
*     CHARACTER         UPLOA,UPLOC,TRANS,DIAG
*     ..
*     .. Array Arguments ..
*     REAL  A(LDA,*),C(LDC,*)
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SSYTRRK  performs one of the symmetric rank k operations
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
*>          ALPHA is REAL.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension ( LDA, k ).
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
*>          C is REAL array, dimension ( LDC, N )
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
*>  Level 3 Blas routine.
*>
*> \endverbatim
*>
*  =====================================================================
      RECURSIVE SUBROUTINE SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, K,
     $            ALPHA, A, LDA, BETA, C, LDC)
*
*     .. Scalar Arguments ..
      REAL  ALPHA,BETA
      INTEGER           K,LDA,LDC
      CHARACTER         UPLOA,UPLOC,TRANS,DIAG
*     ..
*     .. Array Arguments ..
      REAL  A(LDA,*),C(LDC,*)
*     ..
*     .. Parameters ..
      REAL  ZERO, ONE
      PARAMETER(ZERO = 0.0E+0, ONE = 1.0E+0)
*     ..
*     .. Local Scalars ..
      INTEGER           L,NX
      LOGICAL           UPPERA,UPPERC,TRANSL,UNITT
*     ..
*     .. External Subroutines ..
      EXTERNAL          SSYRK,STRMMOOP
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
      IF(K.EQ.0) THEN
         RETURN
      END IF
*
*     Determine the crossover point into the unblocked variant (and base case)
*
      NX = ILAENV(3, 'DSYTRRK', UPLOA // UPLOC // TRANS // DIAG,
     $      K, -1, -1, -1)
*
      IF(K.LT.NX.OR.K.EQ.1) THEN
         CALL SSYTRRK_LVL2(UPLOA, UPLOC, TRANS, DIAG, K,
     $         ALPHA, A, LDA, BETA, C, LDC)
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
*     Recursive case
*
      L = K/2
*
      IF(TRANSL) THEN
*
*        This means we are computing C = alpha*A**H*A + beta*C
*
         IF(UPPERA) THEN
*
*           This means A is upper triangular
*
*           Break C and A apart as follows
*               |-----------------|
*           A = | A_{1,1} A_{1,2} | l
*               | 0       A_{2,2} | k-l
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
*           C = \alpha | A_{1,1}**H 0          |*| A_{1,1} A_{1,2} | + \beta | C_{1,1} C_{1,2} |
*                      | A_{1,2}**H A_{2,2}**H | | 0       A_{2,2} |         | C_{2,1} C_{2,2} |
*                      |-----------------------| |-----------------|         |-----------------|
*
*           Which gives us the following componentwise representation of C
*
*           C_{1,1} = \alpha*A_{1,1}**H*A_{1,1} + \beta C_{1,1}
*           C_{1,2} = \alpha*A_{1,1}**H*A_{1,2} + \beta C_{1,2}
*           C_{2,1} = \alpha*A_{1,2}**H*A_{1,1} + \beta C_{2,1}
*           C_{2,2} = \alpha*A_{1,2}**H*A_{1,2} + \alpha*A_{2,2}**H*A_{2,2} + \beta C_{2,2}
*
*           Thus, we compute the following
*
*           C_{1,1} = \alpha*A_{1,1}**H*A_{1,1} + \beta C_{1,1} (This routine)
*           C_{2,2} = \alpha*A_{2,2}**H*A_{2,2} + \beta C_{2,2} (This routine)
*           C_{2,2} = \alpha*A_{1,2}**H*A_{1,2} + C_{2,2}       (SYRK)
*
*           Compute C_{1,1} = \alpha*A_{1,1}**H*A_{1,1} + \beta C_{1,1}
*
            CALL SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, L, ALPHA,
     $            A, LDA, BETA, C, LDC)
*
*           Compute C_{2,2}
*           C_{2,2} = \alpha*A_{2,2}**H*A_{2,2} + \beta C_{2,2}
*
            CALL SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, K-L, ALPHA,
     $            A(L+1,L+1), LDA, BETA, C(L+1,L+1), LDC)
*
*           C_{2,2} = \alpha*A_{1,2}**H*A_{1,2} + C_{2,2}
*
            CALL SSYRK(UPLOC, TRANS, K-L, L, ALPHA, A(1,L+1), LDA,
     $            ONE, C(L+1,L+1), LDC)
            IF(UPPERC) THEN
*
*              Compute C_{1,2} = \alpha*A_{1,1}**H*A_{1,2} + \beta C_{1,2} (TRMMOOP)
*
               CALL STRMMOOP('Left', UPLOA, 'Transpose',
     $               'No Transpose', DIAG, L, K-L, ALPHA, A, LDA,
     $               A(1, L+1), LDA, BETA, C(1,L+1), LDC)
            ELSE
*
*              Compute C_{2,1} = \alpha*A_{1,2}**H*A_{1,1} + \beta C_{2,1} (TRMMOOP)
*
               CALL STRMMOOP('Right', UPLOA, 'No Transpose',
     $               'Transpose', DIAG, K-L, L, ALPHA, A, LDA,
     $               A(1, L+1), LDA, BETA, C(L+1,1), LDC)
            END IF
         ELSE
*
*           This means A is lower triangular
*
*           Break C and A apart as follows
*               |-----------------|
*           A = | A_{1,1} 0       | l
*               | A_{2,1} A_{2,2} | k-l
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
*           C = \alpha | A_{1,1}**H A_{2,1}**H |*| A_{1,1} 0       | + \beta | C_{1,1} C_{1,2} |
*                      | 0          A_{2,2}**H | | A_{2,1} A_{2,2} |         | C_{2,1} C_{2,2} |
*                      |-----------------------| |-----------------|         |-----------------|
*
*           Which gives us the following componentwise representation of C
*
*           C_{1,1} = \alpha*A_{1,1}**H*A_{1,1} + \alpha*A_{2,1}**H*A_{2,1} + \beta*C_{1,1}
*           C_{1,2} = \alpha*A_{2,1}**H*A_{2,2} + \beta*C_{1,2}
*           C_{2,1} = \alpha*A_{2,2}**H*A_{2,1} + \beta*C_{2,1}
*           C_{2,2} = \alpha*A_{2,2}**H*A_{2,2} + \beta*C_{2,2}
*
*           Thus, we compute the following
*
*           C_{1,1} = \alpha*A_{1,1}**H*A_{1,1} + \beta*C_{1,1} (This routine)
*           C_{1,1} = \alpha*A_{2,1}**H*A_{2,1} + C_{1,1}       (SYRK)
*           C_{2,2} = \alpha*A_{2,2}**H*A_{2,2} + \beta*C_{2,2} (This routine)
*
*           Compute C_{1,1}
*           C_{1,1} = \alpha*A_{1,1}**H*A_{1,1} + \beta*C_{1,1}
*
            CALL SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, L, ALPHA,
     $            A, LDA, BETA, C, LDC)
*
*           C_{1,1} = \alpha*A_{2,1}**H*A_{2,1} + C_{1,1}
*
            CALL SSYRK(UPLOC, TRANS, L, K-L, ALPHA, A(L+1,1), LDA,
     $            ONE, C, LDC)
*
*           Compute C_{2,2} = \alpha*A_{2,2}**H*A_{2,2} + \beta*C_{2,2}
*
            CALL SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, K-L, ALPHA,
     $            A(L+1,L+1), LDA, BETA, C(L+1,L+1), LDC)
            IF(UPPERC) THEN
*
*              Compute C_{1,2} = \alpha*A_{2,1}**H*A_{2,2} + \beta*C_{1,2}
*
               CALL STRMMOOP('Right', UPLOA, 'No Transpose',
     $               'Transpose', DIAG, L, K-L, ALPHA, A(L+1,L+1), LDA,
     $                A(L+1,1), LDA, BETA, C(1,L+1), LDC)
            ELSE
*
*              Compute C_{2,1} = \alpha*A_{2,2}**H*A_{2,1} + \beta*C_{2,1}
*
               CALL STRMMOOP('Left', UPLOA, 'Transpose',
     $               'No Transpose', DIAG, K-L, L, ALPHA,
     $               A(L+1,L+1), LDA, A(L+1, 1), LDA, BETA,
     $               C(L+1,1), LDC)
            END IF
         END IF
      ELSE
*
*        This means we are computing C = alpha*T*T**H + beta*C
*
         IF(UPPERA) THEN
*
*           This means A is upper triangular
*
*           Break C and A apart as follows
*               |-----------------|
*           A = | A_{1,1} A_{1,2} | l
*               | 0       A_{2,2} | k-l
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
*           C = \alpha | A_{1,1} A_{1,2} |*| A_{1,1}**H 0          | + \beta | C_{1,1} C_{1,2} |
*                      | 0       A_{2,2} | | A_{1,2}**H A_{2,2}**H |         | C_{2,1} C_{2,2} |
*                      |-----------------| |-----------------------|         |-----------------|
*
*           Which gives us the following componentwise representation of C
*
*           C_{1,1} = \alpha*A_{1,1}*A_{1,1}**H + \alpha*A_{1,2}*A_{1,2}**H + \beta*C_{1,1}
*           C_{1,2} = \alpha*A_{1,2}*A_{2,2}**H + \beta*C_{1,2}
*           C_{2,1} = \alpha*A_{2,2}*A_{1,2}**H + \beta*C_{2,1}
*           C_{2,2} = \alpha*A_{2,2}*A_{2,2}**H + \beta*C_{2,2}
*
*           Thus, we compute the following
*
*           C_{1,1} = \alpha*A_{1,1}*A_{1,1}**H + \beta*C_{1,1} (This routine)
*           C_{1,1} = \alpha*A_{1,2}*A_{1,2}**H + C_{1,1}       (SYRK)
*           C_{2,2} = \alpha*A_{2,2}*A_{2,2}**H + \beta*C_{2,2} (This routine)
*
*           Compute C_{1,1}
*           C_{1,1} = \alpha*A_{1,1}*A_{1,1}**H + \beta*C_{1,1}
*
            CALL SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, L, ALPHA,
     $            A, LDA, BETA, C, LDC)
*
*           C_{1,1} = \alpha*A_{1,2}*A_{1,2}**H + C_{1,1}
*
            CALL SSYRK(UPLOC, TRANS, L, K-L, ALPHA, A(1,L+1), LDA,
     $            ONE, C, LDC)
*
*           Compute C_{2,2} = \alpha*A_{2,2}*A_{2,2}**H + \beta*C_{2,2}
*
            CALL SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, K-L, ALPHA,
     $            A(L+1,L+1), LDA, BETA, C(L+1,L+1), LDC)
            IF(UPPERC) THEN
*
*              Compute C_{1,2} = \alpha*A_{1,2}*A_{2,2}**H + \beta*C_{1,2}
*
               CALL STRMMOOP('Right', UPLOA, 'Transpose',
     $               'No Transpose', DIAG, L, K-L, ALPHA,
     $               A(L+1,L+1), LDA, A(1, L+1), LDA, BETA,
     $               C(1,L+1), LDC)
            ELSE
*
*              Compute C_{2,1} = \alpha*A_{2,2}*A_{1,2}**H + \beta*C_{2,1}
*
               CALL STRMMOOP('Left', UPLOA, 'No Transpose',
     $               'Transpose', DIAG, K-L, L, ALPHA,
     $               A(L+1,L+1), LDA, A(1, L+1), LDA, BETA,
     $               C(L+1,1), LDC)
            END IF
         ELSE
*
*           This means A is lower triangular
*
*           Break C and A apart as follows
*               |-----------------|
*           A = | A_{1,1} 0       | l
*               | A_{2,1} A_{2,2} | k-l
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
*           C = \alpha | A_{1,1} 0       |*| A_{1,1}**H A_{2,1}**H | + \beta | C_{1,1} C_{1,2} |
*                      | A_{2,1} A_{2,2} | | 0          A_{2,2}**H |         | C_{2,1} C_{2,2} |
*                      |-----------------| |-----------------------|         |-----------------|
*
*           Which gives us the following componentwise representation of C
*
*           C_{1,1} = \alpha*A_{1,1}*A_{1,1}**H + \beta*C_{1,1}
*           C_{1,2} = \alpha*A_{1,1}*A_{2,1}**H + \beta*C_{1,2}
*           C_{2,1} = \alpha*A_{2,1}*A_{1,1}**H + \beta*C_{2,1}
*           C_{2,2} = \alpha*A_{2,1}*A_{2,1}**H + \alpha*A_{2,2}*A_{2,2}**H + \beta*C_{2,2}
*
*           Thus, we compute the following
*
*           C_{1,1} = \alpha*A_{1,1}*A_{1,1}**H + \beta*C_{1,1} (This routine)
*           C_{2,2} = \alpha*A_{2,2}*A_{2,2}**H + \beta*C_{2,2} (This routine)
*           C_{2,2} = \alpha*A_{2,1}*A_{2,1}**H + C_{2,2}       (SYRK)
*
*           Compute C_{1,1} = \alpha*A_{1,1}*A_{1,1}**H + \beta*C_{1,1}
*
            CALL SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, L, ALPHA, A, LDA,
     $            BETA, C, LDC)
*
*           Compute C_{2,2}
*           C_{2,2} = \alpha*A_{2,2}*A_{2,2}**H + \beta*C_{2,2}
*
            CALL SSYTRRK(UPLOA, UPLOC, TRANS, DIAG, K-L, ALPHA,
     $            A(L+1,L+1), LDA, BETA, C(L+1,L+1), LDC)
*
*           C_{2,2} = \alpha*A_{2,1}*A_{2,1}**H + C_{2,2}
*
            CALL SSYRK(UPLOC, TRANS, K-L, L, ALPHA, A(L+1,1), LDA,
     $            ONE, C(L+1,L+1), LDC)
            IF(UPPERC) THEN
*
*              Compute C_{1,2} = \alpha*A_{1,1}*A_{2,1}**H + \beta*C_{1,2}
*
               CALL STRMMOOP('Left', UPLOA, 'No Transpose',
     $               'Transpose', DIAG, L, K-L, ALPHA, A, LDA,
     $               A(L+1,1), LDA, BETA, C(1,L+1), LDC)
            ELSE
*
*              Compute C_{2,1} = \alpha*A_{2,1}*A_{1,1}**H + \beta*C_{2,1}
*
               CALL STRMMOOP('Right', UPLOA, 'Transpose',
     $               'No Transpose', DIAG, K-L, L, ALPHA, A, LDA,
     $               A(L+1,1), LDA, BETA, C(L+1,1), LDC)
            END IF
         END IF
      END IF
      END SUBROUTINE
