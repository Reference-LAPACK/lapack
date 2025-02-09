*> \brief \b CHPTRS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CHPTRS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chptrs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chptrs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chptrs.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX            AP( * ), B( LDB, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CHPTRS solves a system of linear equations A*X = B with a complex
*> Hermitian matrix A stored in packed format using the factorization
*> A = U*D*U**H or A = L*D*L**H computed by CHPTRF.
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;
*>          = 'L':  Lower triangular, form is A = L*D*L**H.
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
*> \param[in] AP
*> \verbatim
*>          AP is COMPLEX array, dimension (N*(N+1)/2)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by CHPTRF, stored as a
*>          packed triangular matrix.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by CHPTRF.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,NRHS)
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
*>          < 0: if INFO = -i, the i-th argument had an illegal value
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
*> \ingroup hptrs
*
*  =====================================================================
      SUBROUTINE CHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            AP( * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, K, KC, KP
      REAL               S
      COMPLEX            AK, AKM1, AKM1K, BK, BKM1, DENOM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV, CGERU, CLACGV, CSSCAL, CSWAP,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, REAL
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHPTRS', -INFO )
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
*        Solve A*X = B, where A = U*D*U**H.
*
*        First solve U*D*X = B, overwriting B with X.
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = N
         KC = N*( N+1 ) / 2 + 1
   10    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 30
*
         KC = KC - K
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K )
     $         CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*           Multiply by inv(U(K)), where U(K) is the transformation
*           stored in column K of A.
*
            CALL CGERU( K-1, NRHS, -ONE, AP( KC ), 1, B( K, 1 ), LDB,
     $                  B( 1, 1 ), LDB )
*
*           Multiply by the inverse of the diagonal block.
*
            S = REAL( ONE ) / REAL( AP( KC+K-1 ) )
            CALL CSSCAL( NRHS, S, B( K, 1 ), LDB )
            K = K - 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Interchange rows K-1 and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K-1 )
     $         CALL CSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
*
*           Multiply by inv(U(K)), where U(K) is the transformation
*           stored in columns K-1 and K of A.
*
            CALL CGERU( K-2, NRHS, -ONE, AP( KC ), 1, B( K, 1 ), LDB,
     $                  B( 1, 1 ), LDB )
            CALL CGERU( K-2, NRHS, -ONE, AP( KC-( K-1 ) ), 1,
     $                  B( K-1, 1 ), LDB, B( 1, 1 ), LDB )
*
*           Multiply by the inverse of the diagonal block.
*
            AKM1K = AP( KC+K-2 )
            AKM1 = AP( KC-1 ) / AKM1K
            AK = AP( KC+K-1 ) / CONJG( AKM1K )
            DENOM = AKM1*AK - ONE
            DO 20 J = 1, NRHS
               BKM1 = B( K-1, J ) / AKM1K
               BK = B( K, J ) / CONJG( AKM1K )
               B( K-1, J ) = ( AK*BKM1-BK ) / DENOM
               B( K, J ) = ( AKM1*BK-BKM1 ) / DENOM
   20       CONTINUE
            KC = KC - K + 1
            K = K - 2
         END IF
*
         GO TO 10
   30    CONTINUE
*
*        Next solve U**H *X = B, overwriting B with X.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
         KC = 1
   40    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 50
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Multiply by inv(U**H(K)), where U(K) is the transformation
*           stored in column K of A.
*
            IF( K.GT.1 ) THEN
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', K-1, NRHS, -ONE, B,
     $                     LDB, AP( KC ), 1, ONE, B( K, 1 ), LDB )
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
            END IF
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K )
     $         CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            KC = KC + K
            K = K + 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Multiply by inv(U**H(K+1)), where U(K+1) is the transformation
*           stored in columns K and K+1 of A.
*
            IF( K.GT.1 ) THEN
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', K-1, NRHS, -ONE, B,
     $                     LDB, AP( KC ), 1, ONE, B( K, 1 ), LDB )
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
*
               CALL CLACGV( NRHS, B( K+1, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', K-1, NRHS, -ONE, B,
     $                     LDB, AP( KC+K ), 1, ONE, B( K+1, 1 ), LDB )
               CALL CLACGV( NRHS, B( K+1, 1 ), LDB )
            END IF
*
*           Interchange rows K and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K )
     $         CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            KC = KC + 2*K + 1
            K = K + 2
         END IF
*
         GO TO 40
   50    CONTINUE
*
      ELSE
*
*        Solve A*X = B, where A = L*D*L**H.
*
*        First solve L*D*X = B, overwriting B with X.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
         KC = 1
   60    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 80
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K )
     $         CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*           Multiply by inv(L(K)), where L(K) is the transformation
*           stored in column K of A.
*
            IF( K.LT.N )
     $         CALL CGERU( N-K, NRHS, -ONE, AP( KC+1 ), 1, B( K, 1 ),
     $                     LDB, B( K+1, 1 ), LDB )
*
*           Multiply by the inverse of the diagonal block.
*
            S = REAL( ONE ) / REAL( AP( KC ) )
            CALL CSSCAL( NRHS, S, B( K, 1 ), LDB )
            KC = KC + N - K + 1
            K = K + 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Interchange rows K+1 and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K+1 )
     $         CALL CSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
*
*           Multiply by inv(L(K)), where L(K) is the transformation
*           stored in columns K and K+1 of A.
*
            IF( K.LT.N-1 ) THEN
               CALL CGERU( N-K-1, NRHS, -ONE, AP( KC+2 ), 1, B( K,
     $                     1 ),
     $                     LDB, B( K+2, 1 ), LDB )
               CALL CGERU( N-K-1, NRHS, -ONE, AP( KC+N-K+2 ), 1,
     $                     B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            END IF
*
*           Multiply by the inverse of the diagonal block.
*
            AKM1K = AP( KC+1 )
            AKM1 = AP( KC ) / CONJG( AKM1K )
            AK = AP( KC+N-K+1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 70 J = 1, NRHS
               BKM1 = B( K, J ) / CONJG( AKM1K )
               BK = B( K+1, J ) / AKM1K
               B( K, J ) = ( AK*BKM1-BK ) / DENOM
               B( K+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
   70       CONTINUE
            KC = KC + 2*( N-K ) + 1
            K = K + 2
         END IF
*
         GO TO 60
   80    CONTINUE
*
*        Next solve L**H *X = B, overwriting B with X.
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = N
         KC = N*( N+1 ) / 2 + 1
   90    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 100
*
         KC = KC - ( N-K+1 )
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Multiply by inv(L**H(K)), where L(K) is the transformation
*           stored in column K of A.
*
            IF( K.LT.N ) THEN
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', N-K, NRHS, -ONE,
     $                     B( K+1, 1 ), LDB, AP( KC+1 ), 1, ONE,
     $                     B( K, 1 ), LDB )
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
            END IF
*
*           Interchange rows K and IPIV(K).
*
            KP = IPIV( K )
            IF( KP.NE.K )
     $         CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Multiply by inv(L**H(K-1)), where L(K-1) is the transformation
*           stored in columns K-1 and K of A.
*
            IF( K.LT.N ) THEN
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', N-K, NRHS, -ONE,
     $                     B( K+1, 1 ), LDB, AP( KC+1 ), 1, ONE,
     $                     B( K, 1 ), LDB )
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
*
               CALL CLACGV( NRHS, B( K-1, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', N-K, NRHS, -ONE,
     $                     B( K+1, 1 ), LDB, AP( KC-( N-K ) ), 1, ONE,
     $                     B( K-1, 1 ), LDB )
               CALL CLACGV( NRHS, B( K-1, 1 ), LDB )
            END IF
*
*           Interchange rows K and -IPIV(K).
*
            KP = -IPIV( K )
            IF( KP.NE.K )
     $         CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            KC = KC - ( N-K+2 )
            K = K - 2
         END IF
*
         GO TO 90
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of CHPTRS
*
      END
