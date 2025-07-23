*> \brief \b DKYTRS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DKYTRS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dkytrs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dkytrs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dkytrs.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DKYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DKYTRS solves a system of linear equations A*X = B with a real
*> skew-symmetric matrix A using the factorization A = U*D*U**T or
*> A = L*D*L**T computed by SKYTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by SKYTRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by SKYTRF.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
*> \ingroup hetrs
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Derived from subroutine dsytrs.
*>
*> \endverbatim
*
*> \par Contributors:
*  ==================
*>
*>    Shuo Zheng, China, Jul 2025 \n
*>
*  =====================================================================
      SUBROUTINE DKYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, K, KP
      DOUBLE PRECISION   AK, AKM1, AKM1K, BK, BKM1, DENOM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( (N.LT.0) .OR. (MOD(N,2).NE.0) ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DKYTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B, where A = U*D*U**T.
*
*        First solve U*D*X = B, overwriting B with X.
*
*        K is the main loop index, decreasing from N to 1
*        in steps of 2.
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 30
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           2 x 2 diagonal block
*
*           Interchange rows K-1 and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K-1 ) THEN
               CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
*
*           Multiply by inv(U(K)), where U(K) is the transformation
*           stored in columns K-1 and K of A.
*
            CALL DGER( K-2, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ),
     $                 LDB, B( 1, 1 ), LDB )
            CALL DGER( K-2, NRHS, -ONE, A( 1, K-1 ), 1, B( K-1, 1 ),
     $                 LDB, B( 1, 1 ), LDB )
*
*           Multiply by the inverse of the diagonal block.
*
            CALL DSCAL( NRHS, -ONE / A( K-1, K ), B( K, 1 ), LDB )
            CALL DSCAL( NRHS, ONE / A( K-1, K ), B( K-1, 1 ), LDB )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K-1, 1 ), LDB )
*
            K = K - 2
         ELSEIF( IPIV( K ).LT.0 ) THEN
*
*           2 x 2 diagonal block
*
*           Interchange rows K and K-1, then K-1 and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K-1 ) THEN
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K-1, 1 ), LDB )
               CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
*
*           Multiply by inv(U(K)), where U(K) is the transformation
*           stored in columns K-1 and K of A.
*
            CALL DGER( K-2, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ),
     $                 LDB, B( 1, 1 ), LDB )
            CALL DGER( K-2, NRHS, -ONE, A( 1, K-1 ), 1, B( K-1, 1 ),
     $                 LDB, B( 1, 1 ), LDB )
*
*           Multiply by the inverse of the diagonal block.
*
            CALL DSCAL( NRHS, -ONE / A( K-1, K ), B( K, 1 ), LDB )
            CALL DSCAL( NRHS, ONE / A( K-1, K ), B( K-1, 1 ), LDB )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K-1, 1 ), LDB )
*
            K = K - 2
         ELSE
*
*           2 x 2 diagonal block
*
*           Multiply by inv(U(K)), where U(K) is the transformation
*           stored in columns K-1 and K of A.
*
            CALL DGER( K-2, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ),
     $                 LDB, B( 1, 1 ), LDB )
            CALL DGER( K-2, NRHS, -ONE, A( 1, K-1 ), 1, B( K-1, 1 ),
     $                 LDB, B( 1, 1 ), LDB )
*
*           Multiply by the inverse of the diagonal block.
*
            CALL DSCAL( NRHS, -ONE / A( K-1, K ), B( K, 1 ), LDB )
            CALL DSCAL( NRHS, ONE / A( K-1, K ), B( K-1, 1 ), LDB )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K-1, 1 ), LDB )
*
            K = K - 2
         END IF
*
         GO TO 10
   30    CONTINUE
*
*        Next solve U**T *X = B, overwriting B with X.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
   40    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 50
*
         IF( IPIV( K+1 ).GT.0 ) THEN
*
*           2 x 2 diagonal block
*
*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
*           stored in columns K and K+1 of A.
*
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1,
     $                  K ),
     $                  1, ONE, B( K, 1 ), LDB )
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB,
     $                  A( 1, K+1 ), 1, ONE, B( K+1, 1 ), LDB )
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K+1 )
            IF( KP.NE.K ) THEN
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
            K = K + 2
         ELSEIF( IPIV( K+1 ).LT.0 ) THEN
*
*           2 x 2 diagonal block
*
*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
*           stored in columns K and K+1 of A.
*
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1,
     $                  K ),
     $                  1, ONE, B( K, 1 ), LDB )
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB,
     $                  A( 1, K+1 ), 1, ONE, B( K+1, 1 ), LDB )
*
*           Interchange rows K and -IPIV(K), then K and K+1.
*
            KP = -IPIV( K+1 )
            IF( KP.NE.K ) THEN
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( K, 1 ), LDB )
            END IF
            K = K + 2
         ELSE
*
*           2 x 2 diagonal block
*
*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
*           stored in columns K and K+1 of A.
*
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB,
     $                  A( 1, K ), 1, ONE, B( K, 1 ), LDB )
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB,
     $                  A( 1, K+1 ), 1, ONE, B( K+1, 1 ), LDB )
            K = K + 2
         END IF
*
         GO TO 40
   50    CONTINUE
*
      ELSE
*
*        Solve A*X = B, where A = L*D*L**T.
*
*        First solve L*D*X = B, overwriting B with X.
*
*        K is the main loop index, increasing from 1 to N
*        in steps of 2.
*
         K = 1
   60    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 80
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           2 x 2 diagonal block
*
*           Interchange rows K+1 and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K+1 ) THEN
               CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
*
*           Multiply by inv(L(K)), where L(K) is the transformation
*           stored in columns K and K+1 of A.
*
            IF( K.LT.N-1 ) THEN
               CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K ), 1,
     $                    B( K, 1 ), LDB, B( K+2, 1 ), LDB )
               CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K+1 ), 1,
     $                    B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            END IF
*
*           Multiply by the inverse of the diagonal block.
*
            CALL DSCAL( NRHS, -ONE / A( K+1, K ), B( K, 1 ), LDB )
            CALL DSCAL( NRHS, ONE / A( K+1, K ), B( K+1, 1 ), LDB )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
*
            K = K + 2
         ELSEIF( IPIV( K ).LT.0 ) THEN
*
*           2 x 2 diagonal block
*
*           Interchange rows K and K+1, then K+1 and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K+1 ) THEN
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
               CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
*
*           Multiply by inv(L(K)), where L(K) is the transformation
*           stored in columns K and K+1 of A.
*
            IF( K.LT.N-1 ) THEN
               CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K ), 1,
     $                    B( K, 1 ), LDB, B( K+2, 1 ), LDB )
               CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K+1 ), 1,
     $                    B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            END IF
*
*           Multiply by the inverse of the diagonal block.
*
            CALL DSCAL( NRHS, -ONE / A( K+1, K ), B( K, 1 ), LDB )
            CALL DSCAL( NRHS, ONE / A( K+1, K ), B( K+1, 1 ), LDB )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
*
            K = K + 2
         ELSE
*
*           2 x 2 diagonal block
*
*           Multiply by inv(L(K)), where L(K) is the transformation
*           stored in columns K and K+1 of A.
*
            IF( K.LT.N-1 ) THEN
               CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K ), 1,
     $                    B( K, 1 ), LDB, B( K+2, 1 ), LDB )
               CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K+1 ), 1,
     $                    B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            END IF
*
*           Multiply by the inverse of the diagonal block.
            CALL DSCAL( NRHS, -ONE / A( K+1, K ), B( K, 1 ), LDB )
            CALL DSCAL( NRHS, ONE / A( K+1, K ), B( K+1, 1 ), LDB )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
*
            K = K + 2
         END IF
*
         GO TO 60
   80    CONTINUE
*
*        Next solve L**T *X = B, overwriting B with X.
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = N
   90    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 100
*
         IF( IPIV( K-1 ).GT.0 ) THEN
*
*           2 x 2 diagonal block
*
*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
*           stored in columns K-1 and K of A.
*
            IF( K.LT.N ) THEN
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE,
     $                     B( K+1, 1 ), LDB, A( K+1, K ),
     $                     1, ONE, B( K, 1 ), LDB )
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE,
     $                     B( K+1, 1 ), LDB, A( K+1, K-1 ), 1,
     $                     ONE, B( K-1, 1 ), LDB )
            END IF
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K-1 )
            IF( KP.NE.K ) THEN
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END IF
            K = K - 2
         ELSEIF( IPIV( K-1 ).LT.0 ) THEN
*
*           2 x 2 diagonal block
*
*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
*           stored in columns K-1 and K of A.
*
            IF( K.LT.N ) THEN
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE,
     $                     B( K+1, 1 ), LDB, A( K+1, K ),
     $                     1, ONE, B( K, 1 ), LDB )
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE,
     $                     B( K+1, 1 ), LDB, A( K+1, K-1 ), 1,
     $                     ONE, B( K-1, 1 ), LDB )
            END IF
*
*           Interchange rows K and -IPIV(K), then K and K-1.
*
            KP = -IPIV( K-1 )
            IF( KP.NE.K ) THEN
               CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( K, 1 ), LDB )
            END IF
            K = K - 2
         ELSE
*
*           2 x 2 diagonal block
*
*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
*           stored in columns K-1 and K of A.
*
            IF( K.LT.N ) THEN
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE,
     $                     B( K+1, 1 ), LDB, A( K+1, K ),
     $                     1, ONE, B( K, 1 ), LDB )
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE,
     $                     B( K+1, 1 ), LDB, A( K+1, K-1 ), 1,
     $                     ONE, B( K-1, 1 ), LDB )
            END IF
            K = K - 2
         END IF
*
         GO TO 90
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of DKYTRS
*
      END
