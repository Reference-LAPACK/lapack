      SUBROUTINE SSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, 
     $                    WORK, INFO )
*
*  -- LAPACK PROTOTYPE routine (version 3.2) --
*
*  -- Written by Julie Langou of the Univ. of TN    --
*     May 2010
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SSYTRS2 solves a system of linear equations A*X = B with a real
*  symmetric matrix A using the factorization A = U*D*U**T or
*  A = L*D*L**T computed by SSYTRF and converted by SSYCONV.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the details of the factorization are stored
*          as an upper or lower triangular matrix.
*          = 'U':  Upper triangular, form is A = U*D*U**T;
*          = 'L':  Lower triangular, form is A = L*D*L**T.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The block diagonal matrix D and the multipliers used to
*          obtain the factor U or L as computed by SSYTRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D
*          as determined by SSYTRF.
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  WORK    (workspace) REAL array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, J, K, KP
      REAL               AK, AKM1, AKM1K, BK, BKM1, DENOM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SSWAP, STRSM, XERBLA
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
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTRS2', -INFO )
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
*        Solve A*X = B, where A = U*D*U'.
*
*       P' * B  
        K=N
        DO WHILE ( K .GE. 1 )
         IF( IPIV( K ).GT.0 ) THEN
*           1 x 1 diagonal block
*           Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K )
     $         CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-1
         ELSE
*           2 x 2 diagonal block
*           Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( KP.EQ.-IPIV( K-1 ) )
     $         CALL SSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-2
         END IF
        END DO
*
*  Compute (U \P' * B) -> B    [ (U \P' * B) ]
*
        CALL STRSM('L','U','N','U',N,NRHS,ONE,A,N,B,N)
*
*  Compute D \ B -> B   [ D \ (U \P' * B) ]
*       
         I=N
         DO WHILE ( I .GE. 1 )
            IF( IPIV(I) .GT. 0 ) THEN
              CALL SSCAL( NRHS, ONE / A( I, I ), B( I, 1 ), N )
            ELSEIF ( I .GT. 1) THEN
               IF ( IPIV(I-1) .EQ. IPIV(I) ) THEN
                  AKM1K = WORK(I)
                  AKM1 = A( I-1, I-1 ) / AKM1K
                  AK = A( I, I ) / AKM1K
                  DENOM = AKM1*AK - ONE
                  DO 15 J = 1, NRHS
                     BKM1 = B( I-1, J ) / AKM1K
                     BK = B( I, J ) / AKM1K
                     B( I-1, J ) = ( AK*BKM1-BK ) / DENOM
                     B( I, J ) = ( AKM1*BK-BKM1 ) / DENOM
 15              CONTINUE
               I = I - 1
               ENDIF
            ENDIF
            I = I - 1
         END DO
*
*      Compute (U' \ B) -> B   [ U' \ (D \ (U \P' * B) ) ]
*
         CALL STRSM('L','U','T','U',N,NRHS,ONE,A,N,B,N)
*
*       P * B  [ P * (U' \ (D \ (U \P' * B) )) ]
*
        K=1
        DO WHILE ( K .LE. N )
         IF( IPIV( K ).GT.0 ) THEN
*           1 x 1 diagonal block
*           Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K )
     $         CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+1
         ELSE
*           2 x 2 diagonal block
*           Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( K .LT. N .AND. KP.EQ.-IPIV( K+1 ) )
     $         CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+2
         ENDIF
        END DO
*
      ELSE
*
*        Solve A*X = B, where A = L*D*L'.
*
*       P' * B  
        K=1
        DO WHILE ( K .LE. N )
         IF( IPIV( K ).GT.0 ) THEN
*           1 x 1 diagonal block
*           Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K )
     $         CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+1
         ELSE
*           2 x 2 diagonal block
*           Interchange rows K and -IPIV(K+1).
            KP = -IPIV( K+1 )
            IF( KP.EQ.-IPIV( K ) )
     $         CALL SSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+2
         ENDIF
        END DO
*
*  Compute (L \P' * B) -> B    [ (L \P' * B) ]
*
        CALL STRSM('L','L','N','U',N,NRHS,ONE,A,N,B,N)
*
*  Compute D \ B -> B   [ D \ (L \P' * B) ]
*       
         I=1
         DO WHILE ( I .LE. N )
            IF( IPIV(I) .GT. 0 ) THEN
              CALL SSCAL( NRHS, ONE / A( I, I ), B( I, 1 ), N )
            ELSE
                  AKM1K = WORK(I)
                  AKM1 = A( I, I ) / AKM1K
                  AK = A( I+1, I+1 ) / AKM1K
                  DENOM = AKM1*AK - ONE
                  DO 25 J = 1, NRHS
                     BKM1 = B( I, J ) / AKM1K
                     BK = B( I+1, J ) / AKM1K
                     B( I, J ) = ( AK*BKM1-BK ) / DENOM
                     B( I+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
 25              CONTINUE
                  I = I + 1
            ENDIF
            I = I + 1
         END DO
*
*  Compute (L' \ B) -> B   [ L' \ (D \ (L \P' * B) ) ]
* 
        CALL STRSM('L','L','T','U',N,NRHS,ONE,A,N,B,N)
*
*       P * B  [ P * (L' \ (D \ (L \P' * B) )) ]
*
        K=N
        DO WHILE ( K .GE. 1 )
         IF( IPIV( K ).GT.0 ) THEN
*           1 x 1 diagonal block
*           Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K )
     $         CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-1
         ELSE
*           2 x 2 diagonal block
*           Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( K.GT.1 .AND. KP.EQ.-IPIV( K-1 ) )
     $         CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-2
         ENDIF
        END DO
*
      END IF
*
      RETURN
*
*     End of SSYTRS2
*
      END
