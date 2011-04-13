      SUBROUTINE ZSYSWAPR( UPLO, N, A, LDA, I1, I2)
*
*  -- LAPACK auxiliary routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      CHARACTER        UPLO
      INTEGER          I1, I2, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX   A( LDA, N )
*
*  Purpose
*  =======
*
*  ZSYSWAPR applies an elementary permutation on the rows and the columns of
*  a symmetric matrix.
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
*  A       (input/output) DOUBLE COMPLEX array, dimension (LDA,N)
*          On entry, the NB diagonal matrix D and the multipliers
*          used to obtain the factor U or L as computed by ZSYTRF.
*
*          On exit, if INFO = 0, the (symmetric) inverse of the original
*          matrix.  If UPLO = 'U', the upper triangular part of the
*          inverse is formed and the part of A below the diagonal is not
*          referenced; if UPLO = 'L' the lower triangular part of the
*          inverse is formed and the part of A above the diagonal is
*          not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  I1      (input) INTEGER
*          Index of the first row to swap
*
*  I2      (input) INTEGER
*          Index of the second row to swap
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE COMPLEX     TMP
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZSWAP
*     ..
*     .. Executable Statements ..
*
      UPPER = LSAME( UPLO, 'U' )
      IF (UPPER) THEN
*
*         UPPER
*         first swap
*          - swap column I1 and I2 from I1 to I1-1 
         CALL ZSWAP( I1-1, A(1,I1), 1, A(1,I2), 1 )
*
*          second swap :
*          - swap A(I1,I1) and A(I2,I2)
*          - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1     
         TMP=A(I1,I1)
         A(I1,I1)=A(I2,I2)
         A(I2,I2)=TMP
*
         DO I=1,I2-I1-1
            TMP=A(I1,I1+I)
            A(I1,I1+I)=A(I1+I,I2)
            A(I1+I,I2)=TMP
         END DO
*
*          third swap
*          - swap row I1 and I2 from I2+1 to N
         DO I=I2+1,N
            TMP=A(I1,I)
            A(I1,I)=A(I2,I)
            A(I2,I)=TMP
         END DO
*
        ELSE
*
*         LOWER
*         first swap
*          - swap row I1 and I2 from I1 to I1-1 
         CALL ZSWAP( I1-1, A(I1,1), LDA, A(I2,1), LDA )
*
*         second swap :
*          - swap A(I1,I1) and A(I2,I2)
*          - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1     
          TMP=A(I1,I1)
          A(I1,I1)=A(I2,I2)
          A(I2,I2)=TMP
*
          DO I=1,I2-I1-1
             TMP=A(I1+I,I1)
             A(I1+I,I1)=A(I2,I1+I)
             A(I2,I1+I)=TMP
          END DO
*
*         third swap
*          - swap col I1 and I2 from I2+1 to N
          DO I=I2+1,N
             TMP=A(I,I1)
             A(I,I1)=A(I,I2)
             A(I,I2)=TMP
          END DO
*
      ENDIF
      END SUBROUTINE ZSYSWAPR

