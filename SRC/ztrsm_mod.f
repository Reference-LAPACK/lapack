*> \brief \b ZTRSM_MOD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     RECURSIVE SUBROUTINE ZTRSM_MOD(SIDE, UPLO, TRANSA, DIAG, M, N,
*    $      ALPHA, A, LDA, B, LDB)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA
*       INTEGER LDA,LDB,M,N
*       CHARACTER DIAG,SIDE,TRANSA,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),B(LDB,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRSM_MOD solves one of the matrix equations
*>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T.
*>
*> The matrix X is overwritten on B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether op( A ) appears on the left
*>           or right of X as follows:
*>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n'   op( A ) = A.
*>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c'   op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16.
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*>           zero then  A is not referenced and  B need not be set before
*>           entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, k ),
*>           where k is m when SIDE = 'L' or 'l'
*>             and k is n when SIDE = 'R' or 'r'.
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*>           upper triangular part of the array  A must contain the upper
*>           triangular matrix  and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*>           lower triangular part of the array  A must contain the lower
*>           triangular matrix  and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           A  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension ( LDB, N )
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain  the  right-hand  side  matrix  B,  and  on exit  is
*>           overwritten by the solution matrix  X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
*>           max( 1, m ).
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
*> \ingroup trsm
*  =====================================================================
      RECURSIVE SUBROUTINE ZTRSM_MOD(SIDE, UPLO, TRANSA, DIAG, M, N,
     $      ALPHA, A, LDA, B, LDB)
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA,NX
      LOGICAL LSIDE,NOUNIT,UPPER,NOTRANS
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE,ZERO,NEG_ONE
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0,NEG_ONE=-1.0D+0)
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
      NOTRANS = LSAME(TRANSA,'N')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
         INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
         INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
         INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND.
     +         (.NOT.LSAME(DIAG,'N'))) THEN
         INFO = 4
      ELSE IF (M.LT.0) THEN
         INFO = 5
      ELSE IF (N.LT.0) THEN
         INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
         INFO = 11
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DTRSM ',INFO)
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN

         CALL ZLASET('All', M, N, ZERO, ZERO, B, LDB)
         RETURN
      END IF
      ! Check if we are performing the base case or not
!     NX = ILAENV( 3, 'DTRSM_MOD', SIDE // UPLO // TRANSA // DIAG,
!    $   N, -1, -1, -1 )
      NX = 1 ! Add something else here later
      IF( MIN(M,N).LE.NX ) THEN
         CALL ZTRSM_LVL2_MOD(SIDE, UPLO, TRANSA, DIAG, M, N,
     $      ALPHA, A, LDA, B, LDB)
         RETURN
      END IF
*
*     Start the operations. Recursive Case
*
      K = MIN(M,N)/2
      IF( LSIDE ) THEN
         IF( UPPER ) THEN
            IF( NOTRANS ) THEN
*
*              Break A and B apart as follows
*                  |---------------|       |---------------|
*              A = | A_{11} A_{12} |   B = | B_{11} B_{12} |
*                  | 0      A_{22} |       | B_{21} B_{22} |
*                  |---------------|       |---------------|
*
*              Where
*
*              A_{11}\in\R^{k\times k} A_{12}\in\R^{  k\times m-k}
*                                      A_{22}\in\R^{m-k\times m-k}
*
*              B_{11}\in\R^{  k\times k} B_{12}\in\R^{  k\times n-k}
*              B_{21}\in\R^{m-k\times k} B_{22}\in\R^{m-k\times n-k}
*
*              We want to solve the system
*              |---------------|   |---------------|          |---------------|
*              | A_{11} A_{12} | * | X_{11} X_{12} | = \alpha | B_{11} B_{12} |
*              | 0      A_{22} |   | X_{21} X_{22} |          | B_{21} B_{22} |
*              |---------------|   |---------------|          |---------------|
*
*              Where X (and it's submatrix components) are the same shape as B
*              and its respective submatrices
*
*              This gives us the following system (overwriting X onto B)
*
*              A_{11}*X_{11} + A_{12}*X_{21} = \alpha B_{11}
*                              A_{22}*X_{21} = \alpha B_{21}
*
*              A_{11}*X_{12} + A_{12}*X_{22} = \alpha B_{12}
*                              A_{22}*X_{22} = \alpha B_{22}
*
*              We will compute this columnwise moving bottom to top
*              within each column. This leads to the following algorithm
*
*              solve A_{22}*X_{21} = \alpha B_{21}     (This routine)
*
*              B_{11} = -A_{12}*X_{21} + \alpha B_{11} (GEMM)
*              solve A_{11}*X_{11} = B_{11}            (This routine)
*
*              solve A_{22}*X_{22} = \alpha B_{22}     (This routine)
*
*              B_{12} = -A_{12}*X_{22} + \alpha B_{12} (GEMM)
*              solve A_{11}*X_{12} = B_{12}            (This routine)
*
               ! Compute X_{21} overwriting B_{21}
               CALL ZTRSM_MOD('Left', 'Upper', 'No Transpose', DIAG,
     $         M-K, K, ALPHA, A(K+1,K+1), LDA, B(K+1,1), LDB)

               ! Update B_{11} with X_{21}
               CALL ZGEMM('N', 'N', K, K, M-K, NEG_ONE,
     $         A(1,K+1), LDA, B(K+1,1), LDB, ALPHA, B, LDB)
               ! Compute X_{11} storing it in B_{11}
               CALL ZTRSM_MOD('Left', 'Upper', 'No Transpose', DIAG,
     $         K, K, ONE, A, LDA, B, LDB)

               ! Compute X_{22} overwriting B_{22}
               CALL ZTRSM_MOD('Left', 'Upper', 'No Transpose', DIAG,
     $         M-K, N-K, ALPHA, A(K+1,K+1), LDA, B(K+1,K+1), LDB)

               ! Update B_{12} with X_{22}
               CALL ZGEMM('N', 'N', K, N-K, M-K, NEG_ONE,
     $         A(1,K+1), LDA, B(K+1,K+1), LDB, ALPHA, B(1,K+1), LDB)
               !Compute X_{12} storing it in B_{12}
               CALL ZTRSM_MOD('Left', 'Upper', 'No Transpose', DIAG,
     $         K, N-K, ONE, A, LDA, B(1,K+1), LDB)
            ELSE ! A is transposed
*
*              Break A and B apart as follows
*                  |---------------|       |---------------|
*              A = | A_{11} A_{12} |   B = | B_{11} B_{12} |
*                  | 0      A_{22} |       | B_{21} B_{22} |
*                  |---------------|       |---------------|
*
*              Where
*
*              A_{11}\in\R^{k\times k} A_{12}\in\R^{  k\times m-k}
*                                      A_{22}\in\R^{m-k\times m-k}
*
*              B_{11}\in\R^{  k\times k} B_{12}\in\R^{  k\times n-k}
*              B_{21}\in\R^{m-k\times k} B_{22}\in\R^{m-k\times n-k}
*
*              We want to solve the system
*              |-------------------|   |---------------|          |---------------|
*              | A_{11}^T 0        | * | X_{11} X_{12} | = \alpha | B_{11} B_{12} |
*              | A_{12}^T A_{22}^T |   | X_{21} X_{22} |          | B_{21} B_{22} |
*              |-------------------|   |---------------|          |---------------|
*
*              Where X (and it's submatrix components) are the same shape as B
*              and its respective submatrices
*
*              This gives us the following system (overwriting X onto B)
*
*              A_{11}^T X_{11}                   = \alpha B_{11}
*              A_{12}^T X_{11} + A_{22}^T X_{21} = \alpha B_{21}
*
*              A_{11}^T X_{12}                   = \alpha B_{12}
*              A_{12}^T X_{12} + A_{22}^T X_{22} = \alpha B_{22}
*
*              Which we order as follows
*
*              Overwite B_{11} with solution to A_{11}^T X = \alpha B_{11} (This routine)
*
*              B_{21} = -A_{12}^T X_{11} + alpha B_{21} (GEMM)
*              Overwrite B_{21} with solution to A_{22}^T X = B_{21} (This routine)
*
*              Overwrite B_{12} with solution to A_{11}^T X = \alpha B_{12} (This routine)
*
*              B_{22} = -A_{12}^T X_{12} + alpha B_{22} (GEMM)
*              Overwrite B_{22} with solution to A_{22}^T X = B_{22} (This routine)
*

               ! Compute X_{11}
               CALL ZTRSM_MOD('left', 'upper', 'transpose', DIAG,
     $            K, K, ALPHA, A, LDA, B, LDB )
               ! Update B_{21}
               CALL ZGEMM( 'transpose', 'no transpose', M-K, K, K,
     $            NEG_ONE, A(1,K+1), LDA, B, LDB, ALPHA, B(K+1,1), LDB)
               ! Solve for X_{21}
               CALL ZTRSM_MOD( 'left', 'upper', 'transpose', DIAG,
     $            M-K, K, ONE, A(K+1,K+1), LDA, B(K+1,1), LDB )

               ! Compute X_{12}
               CALL ZTRSM_MOD( 'left', 'upper', 'transpose', DIAG,
     $            K, N-K, ALPHA, A, LDA, B(1,K+1), LDB )

               !Update B_{22}
               CALL ZGEMM( 'transpose', 'no transpose', M-K, N-K, K,
     $            NEG_ONE, A(1,K+1), LDA, B(1,K+1), LDB, ALPHA,
     $            B(K+1,K+1), LDB )
               ! Solve for B_{22}
               CALL ZTRSM_MOD( 'left', 'upper', 'transpose', DIAG,
     $            M-K, N-K, ONE, A(K+1,K+1), LDA, B(K+1,K+1), LDB )
            END IF
         ELSE ! A is lower triangular
            IF( NOTRANS ) THEN
*
*              Break A and B apart as follows
*                  |---------------|       |---------------|
*              A = | A_{11} 0      |   B = | B_{11} B_{12} |
*                  | A_{21} A_{22} |       | B_{21} B_{22} |
*                  |---------------|       |---------------|
*
*              Where
*
*              A_{11}\in\R^{  k\times k}
*              A_{21}\in\R^{m-k\times k} A_{22}\in\R^{m-k\times m-k}
*
*              B_{11}\in\R^{  k\times k} B_{12}\in\R^{  k\times n-k}
*              B_{21}\in\R^{m-k\times k} B_{22}\in\R^{m-k\times n-k}
*
*              We want to solve the system
*
*              |---------------|   |---------------|          |---------------|
*              | A_{11} 0      | * | X_{11} X_{12} | = \alpha | B_{11} B_{12} |
*              | A_{21} A_{22} |   | X_{21} X_{22} |          | B_{21} B_{22} |
*              |---------------|   |---------------|          |---------------|
*
*              Where X (and it's submatrix components) are the same shape as B
*              and its respective submatrices
*
*              This gives us the following system (overwriting X onto B)
*
*              A_{11}*X_{11}                 = \alpha B_{11}
*              A_{21}*X_{11} + A_{22}*X_{21} = \alpha B_{21}
*
*              A_{11}*X_{12}                 = \alpha B_{12}
*              A_{21}*X_{12} + A_{22}*X_{22} = \alpha B_{22}
*
*              Which we order as follows
*
*              Overwite B_{11} with solution to A_{11}^T X = \alpha B_{11} (This routine)
*
*              B_{21} = -A_{21}*X_{11} + \alpha B_{21} (GEMM)
*              Overwrite B_{21} with solution to A_{22} X = B_{21} (This routine)
*
*              Overwrite B_{12} with solution to A_{11} X = \alpha B_{12} (This routine)
*
*              B_{22} = -A_{21} X_{12} + \alpha B_{22} (GEMM)
*              Overwrite B_{22} with solution to A_{22} X = B_{22}
*
               ! Compute X_{11}
               CALL ZTRSM_MOD('left', 'lower', 'no transpose', DIAG,
     $            K, K, ALPHA, A, LDA, B, LDB)

               ! Update B_{21}
               CALL ZGEMM('no transpose', 'no transpose', M-K, K, K,
     $            NEG_ONE, A(K+1,1), LDA, B, LDB, ALPHA, B(K+1,1), LDB)
               ! Solve for X_{21}
               CALL ZTRSM_MOD('left', 'lower', 'no transpose',
     $            DIAG, M-K, K, ONE, A(K+1,K+1), LDA, B(K+1,1), LDB)

               ! Compute X_{12}
               CALL ZTRSM_MOD('left', 'lower', 'no transpose', DIAG,
     $            K, N-K, ALPHA, A, LDA, B(1,K+1), LDB )

               ! Update B_{22}
               CALL ZGEMM('no transpose', 'no transpose', M-K, N-K,
     $            K, NEG_ONE, A(K+1,1), LDA, B(1,K+1), LDB, ALPHA,
     $            B(K+1,K+1), LDB)
               ! Solve for X_{22}
               CALL ZTRSM_MOD('left', 'lower', 'no transpose', DIAG,
     $            M-K, N-K, ONE, A(K+1,K+1), LDA, B(K+1,K+1), LDB)
            ELSE ! A is transposed
*
*              Break A and B apart as follows
*                  |---------------|       |---------------|
*              A = | A_{11} 0      |   B = | B_{11} B_{12} |
*                  | A_{21} A_{22} |       | B_{21} B_{22} |
*                  |---------------|       |---------------|
*
*              Where
*
*              A_{11}\in\R^{  k\times k}
*              A_{21}\in\R^{m-k\times k} A_{22}\in\R^{m-k\times m-k}
*
*              B_{11}\in\R^{  k\times k} B_{12}\in\R^{  k\times n-k}
*              B_{21}\in\R^{m-k\times k} B_{22}\in\R^{m-k\times n-k}
*
*              We want to solve the system
*
*              |--------------------------|   |---------------|          |---------------|
*              | A_{11}^\top  A_{21}^\top | * | X_{11} X_{12} | = \alpha | B_{11} B_{12} |
*              | 0            A_{22}^\top |   | X_{21} X_{22} |          | B_{21} B_{22} |
*              |--------------------------|   |---------------|          |---------------|
*
*              Where X (and it's submatrix components) are the same shape as B
*              and its respective submatrices
*
*              This gives us the following system (overwriting X onto B)
*
*              A_{11}^\top*X_{11} + A_{21}^\top*X_{21} = \alpha B_{11}
*                                   A_{22}^\top*X_{21} = \alpha B_{21}
*
*              A_{11}^\top*X_{12} + A_{21}^\top*X_{22} = \alpha B_{12}
*                                   A_{22}^\top*X_{22} = \alpha B_{22}
*
               ! Compute X_{21}
               CALL ZTRSM_MOD('left', 'lower', 'transpose', DIAG,
     $            M-K, K, ALPHA, A(K+1,K+1), LDA, B(K+1,1), LDB)

               ! Update B_{11}
               CALL ZGEMM('transpose', 'no transpose', K, K, M-K,
     $            NEG_ONE, A(K+1,1), LDA, B(K+1,1), LDB, ALPHA,
     $            B, LDB)

               ! solve for X_{11}
               CALL ZTRSM_MOD('left', 'lower', 'transpose', DIAG,
     $            K, K, ONE, A, LDA, B, LDB)

               ! Solve for X_{22}
               CALL ZTRSM_MOD('left', 'lower', 'transpose', DIAG,
     $            M-K, N-K, ALPHA, A(K+1,K+1), LDA, B(K+1,K+1), LDB)

               ! Update B_{12}
               CALL ZGEMM('transpose', 'no transpose', K, N-K, M-K,
     $            NEG_ONE, A(K+1, 1), LDA, B(K+1,K+1), LDB, ALPHA,
     $            B(1,K+1), LDB)

               ! Solve for X_{12}
               CALL ZTRSM_MOD('left', 'lower', 'transpose', DIAG,
     $            K, N-K, ONE, A, LDA, B(1,K+1), LDB)
            END IF
         END IF
      ELSE ! A is on the right
         IF (UPPER) THEN
            IF( NOTRANS ) THEN
*
*              Break A and B apart as
*                    |-----------------|        |-----------------|
*              A =   | A_{11} A_{12}   |  B =   | B_{11} B_{12}   |
*                    | 0      A_{22}   |        | B_{21} B_{22}   |
*                    |-----------------|        |-----------------|
*
*              Where
*
*              A_{11}\in\R^{k\times k} A_{12}\in\R^{  k\times n-k}
*                                      A_{22}\in\R^{n-k\times n-k}
*
*              B_{11}\in\R^{  k\times k}  B_{12}\in\R^{  k\times n-k}
*              B_{21}\in\R^{m-k\times k}  B_{22}\in\R^{m-k\times n-k}
*
*              We want to solve the system below replacing B with X
*
*              |---------------| |---------------|          |---------------|
*              | X_{11} X_{12} | | A_{11} A_{12} | = \alpha | B_{11} B_{12} |
*              | X_{21} X_{22} | | 0      A_{22} |          | B_{21} B_{22} |
*              |---------------| |---------------|          |---------------|
*
*              Where X (and it's submatrix components) are the same shape as B
*              and its respective submatrices
*
*              This gives us the following system (overwriting X onto B)
*
*              X_{11} A_{11}                 = \alpha B_{11}
*              X_{11} A_{12} + X_{12} A_{22} = \alpha B_{12}
*
*              X_{21} A_{11}                 = \alpha B_{21}
*              X_{21} A_{12} + X_{22} A_{22} = \alpha B_{22}
*
               ! Solve for X_{11}
               CALL ZTRSM_MOD(SIDE, UPLO, TRANSA, DIAG, K, K, ALPHA,
     $            A, LDA, B, LDB)
               ! Update B_{12}
               CALL ZGEMM('No transpose', 'No transpose', K, N-K, K,
     $            NEG_ONE, B, LDB, A(1,K+1), LDA, ALPHA, B(1,K+1), LDB)
               ! Solve for X_{12}
               CALL ZTRSM_MOD(SIDE, UPLO, TRANSA, DIAG, K, N-K,
     $            ALPHA, A(K+1,K+1), LDA, B(1,K+1), LDB)
               ! Solve for X_{21}
               CALL ZTRSM_MOD(SIDE, UPLO, TRANSA, DIAG, M-K, K,
     $            ALPHA, A, LDA, B(K+1, 1), LDB)
               ! Update B_{22}
               CALL ZGEMM('No transpose', 'No transpose', M-K, N-K,
     $            K, NEG_ONE, B(K+1, 1), LDB, A(1,K+1), LDA, ALPHA,
     $            B(K+1, K+1), LDB)
               ! Solve for X_{22}
               CALL ZTRSM_MOD('Right', 'Upper', 'No Transpose', DIAG,
     $            M-K, N-K, ONE, A(K+1,K+1), LDA, B(K+1,K+1), LDB)
            ELSE ! A is transposed
*
*              Break A and B apart as
*                    |-----------------|        |-----------------|
*              A =   | A_{11} A_{12}   |  B =   | B_{11} B_{12}   |
*                    | 0      A_{22}   |        | B_{21} B_{22}   |
*                    |-----------------|        |-----------------|
*
*              Where
*
*              A_{11}\in\R^{k\times k} A_{12}\in\R^{  k\times n-k}
*                                      A_{22}\in\R^{n-k\times n-k}
*
*              B_{11}\in\R^{  k\times k}  B_{12}\in\R^{  k\times n-k}
*              B_{21}\in\R^{m-k\times k}  B_{22}\in\R^{m-k\times n-k}
*
*              We want to solve the system below replacing B with X
*
*              |---------------| |---------------------|          |---------------|
*              | X_{11} X_{12} | | A_{11}**T 0         | = \alpha | B_{11} B_{12} |
*              | X_{21} X_{22} | | A_{12}**T A_{22}**T |          | B_{21} B_{22} |
*              |---------------| |---------------------|          |---------------|
*
*              Where X (and it's submatrix components) are the same shape as B
*              and its respective submatrices
*
*              This gives us the following system (overwriting X onto B)
*
*              X_{11} A_{11}**T + X_{12} A_{12}**T = alpha B_{11}
*                                 X_{12} A_{22}**T = alpha B_{12}
*
*              X_{21} A_{11}**T + X_{22} A_{12}**T = alpha B_{21}
*                                 X_{22} A_{22}**T = alpha B_{22}
*
               ! Compute X_{12}
               CALL ZTRSM_MOD('Right', 'Upper', 'Transpose', DIAG,
     $            K, N-K, ALPHA, A(K+1,K+1), LDA, B(1, K+1), LDB)
               ! Update B_{11}
               CALL ZGEMM('No Transpose','Transpose', K, K, N-K,
     $            NEG_ONE, B(1,K+1), LDB, A(1,K+1), LDA, ALPHA, B, LDB)
               ! Solve for X_{11}
               CALL ZTRSM_MOD('Right', 'Upper', 'Transpose', DIAG,
     $            K, K, ONE, A, LDA, B, LDB)
               ! Compute X_{22}
               CALL ZTRSM_MOD('Right', 'Upper', 'Transpose', DIAG,
     $            M-K, N-K, ALPHA, A(K+1,K+1), LDA, B(K+1,K+1), LDB)
               ! Update B_{21}
               CALL ZGEMM('No Transpose', 'Transpose', M-K, K, K,
     $            NEG_ONE, B(K+1,K+1), LDB, A(1,K+1), LDA, ALPHA,
     $            B(K+1,1), LDB)
               ! Solve for X_{21}
               CALL ZTRSM_MOD('Right', 'Upper', 'Transpose', DIAG,
     $            M-K, K, ONE, A, LDA, B(K+1,1), LDB)
            END IF
         ELSE ! A is lower triangular
            IF( NOTRANS ) THEN
*
*              Break A and B apart as
*                    |-----------------|        |-----------------|
*              A =   | A_{11} 0        |  B =   | B_{11} B_{12}   |
*                    | A_{21} A_{22}   |        | B_{21} B_{22}   |
*                    |-----------------|        |-----------------|
*
*              Where
*
*              A_{11}\in\R^{  k\times k} 0
*              A_{21}\in\R^{n-k\times k} A_{22}\in\R^{n-k\times n-k}
*
*              B_{11}\in\R^{  k\times k}  B_{12}\in\R^{  k\times n-k}
*              B_{21}\in\R^{m-k\times k}  B_{22}\in\R^{m-k\times n-k}
*
*              We want to solve the system below replacing B with X
*
*              |---------------| |-----------------|          |---------------|
*              | X_{11} X_{12} | | A_{11} 0        | = \alpha | B_{11} B_{12} |
*              | X_{21} X_{22} | | A_{21} A_{22}   |          | B_{21} B_{22} |
*              |---------------| |-----------------|          |---------------|
*
*              Where X (and it's submatrix components) are the same shape as B
*              and its respective submatrices
*
*              This gives us the following system (overwriting X onto B)
*
*              X_{11} A_{11} + X_{12} A_{21} = \alpha B_{11}
*                              X_{12} A_{22} = \alpha B_{12}
*
*              X_{21} A_{11} + X_{22} A_{21} = \alpha B_{21}
*                              X_{22} A_{22} = \alpha B_{22}
*
               ! Solve for X_{12}
               CALL ZTRSM_MOD('Right', 'Lower', 'No Transpose', DIAG,
     $            K, N-K, ALPHA, A(K+1,K+1), LDA, B(1,K+1), LDB)
               ! Update B_{11}
               CALL ZGEMM('No Transpose', 'No Transpose', K, K, N-K,
     $            NEG_ONE, B(1,K+1), LDB, A(K+1,1), LDA, ALPHA, B, LDB)
               ! Solve for X_{11}
               CALL ZTRSM_MOD('Right', 'Lower', 'No Transpose', DIAG,
     $            K, K, ONE, A, LDA, B, LDB)
               ! Solve for X_{22}
               CALL ZTRSM_MOD('Right', 'Lower', 'No Transpose', DIAG,
     $            M-K, N-K, ALPHA, A(K+1,K+1), LDA, B(K+1,K+1), LDB)
               ! Update B_{21}
               CALL ZGEMM('No Transpose', 'Transpose', M-K, K, N-K,
     $            NEG_ONE, B(K+1,K+1), LDB, A(K+1,1), LDA,
     $            ALPHA, B(K+1,1), LDB)
               ! Solve for X_{21}
               CALL ZTRSM_MOD('Right', 'Lower', 'No Transpose', DIAG,
     $            M-K, K, ONE, A, LDA, B(K+1,1), LDB)
            ELSE  ! A is transposed
*
*              Break A and B apart as
*                    |-----------------|        |-----------------|
*              A =   | A_{11} 0        |  B =   | B_{11} B_{12}   |
*                    | A_{21} A_{22}   |        | B_{21} B_{22}   |
*                    |-----------------|        |-----------------|
*
*              Where
*
*              A_{11}\in\R^{  k\times k} 0
*              A_{21}\in\R^{n-k\times k} A_{22}\in\R^{n-k\times n-k}
*
*              B_{11}\in\R^{  k\times k}  B_{12}\in\R^{  k\times n-k}
*              B_{21}\in\R^{m-k\times k}  B_{22}\in\R^{m-k\times n-k}
*
*              We want to solve the system below replacing B with X
*
*              |---------------| |-----------------------|          |---------------|
*              | X_{11} X_{12} | | A_{11}**T A_{21}**T   | = \alpha | B_{11} B_{12} |
*              | X_{21} X_{22} | | 0         A_{22}**T   |          | B_{21} B_{22} |
*              |---------------| |-----------------------|          |---------------|
*
*              Where X (and it's submatrix components) are the same shape as B
*              and its respective submatrices
*
*              This gives us the following system (overwriting X onto B)
*
*              X_{11} A_{11}**T                    = \alpha B_{11}
*              X_{11} A_{21}**T + X_{12} A_{22}**T = \alpha B_{12}
*
*              X_{21} A_{11}**T                    = \alpha B_{21}
*              X_{21} A_{21}**T + X_{22} A_{22}**T = \alpha B_{22}
*
               ! Solve for X_{11}
               CALL ZTRSM_MOD('Right', 'Lower', 'Transpose', DIAG,
     $            K, K, ALPHA, A, LDA, B, LDB)
               ! Update B_{12}
               CALL ZGEMM('No Transpose', 'Transpose', K, N-K, K,
     $            NEG_ONE, B, LDB, A(K+1,1), LDA,
     $            ALPHA, B(1,K+1), LDB)
               ! Solve for X_{12}
               CALL ZTRSM_MOD('Right', 'Lower', 'Transpose', DIAG,
     $            K, N-K, ONE, A(K+1,K+1), LDA, B(1,K+1), LDB)
               ! Solve for X_{21}
               CALL ZTRSM_MOD('Right', 'Lower', 'Transpose', DIAG,
     $            M-K, K, ALPHA, A, LDA, B(K+1,1), LDB)
               ! Update B_{22}
               CALL ZGEMM('No Transpose', 'Transpose', M-K, N-K, K,
     $            NEG_ONE, B(K+1,1), LDB, A(K+1,1), LDA,
     $            ALPHA, B(K+1,K+1), LDB)
            END IF
         END IF
      END IF
      RETURN
*
*     End of ZTRSM_MOD
*
      END
