      SUBROUTINE ZLAPMR( FORWRD, M, N, X, LDX, K )
      IMPLICIT NONE
*
*     Originally ZLAPMT
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     Adapted to ZLAPMR
*     July 2010
*
*     .. Scalar Arguments ..
      LOGICAL            FORWRD
      INTEGER            LDX, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            K( * )
      COMPLEX*16         X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLAPMR rearranges the rows of the M by N matrix X as specified
*  by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.
*  If FORWRD = .TRUE.,  forward permutation:
*
*       X(K(I),*) is moved X(I,*) for I = 1,2,...,M.
*
*  If FORWRD = .FALSE., backward permutation:
*
*       X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.
*
*  Arguments
*  =========
*
*  FORWRD  (input) LOGICAL
*          = .TRUE., forward permutation
*          = .FALSE., backward permutation
*
*  M       (input) INTEGER
*          The number of rows of the matrix X. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix X. N >= 0.
*
*  X       (input/output) COMPLEX*16 array, dimension (LDX,N)
*          On entry, the M by N matrix X.
*          On exit, X contains the permuted matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X, LDX >= MAX(1,M).
*
*  K       (input/output) INTEGER array, dimension (M)
*          On entry, K contains the permutation vector. K is used as
*          internal workspace, but reset to its original value on
*          output.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IN, J, JJ
      COMPLEX*16         TEMP
*     ..
*     .. Executable Statements ..
*
      IF( M.LE.1 )
     $   RETURN
*
      DO 10 I = 1, M
         K( I ) = -K( I )
   10 CONTINUE
*
      IF( FORWRD ) THEN
*
*        Forward permutation
*
         DO 50 I = 1, M
*
            IF( K( I ).GT.0 )
     $         GO TO 40
*
            J = I
            K( J ) = -K( J )
            IN = K( J )
*
   20       CONTINUE
            IF( K( IN ).GT.0 )
     $         GO TO 40
*
            DO 30 JJ = 1, N
               TEMP = X( J, JJ )
               X( J, JJ ) = X( IN, JJ )
               X( IN, JJ ) = TEMP
   30       CONTINUE
*
            K( IN ) = -K( IN )
            J = IN
            IN = K( IN )
            GO TO 20
*
   40       CONTINUE
*
   50    CONTINUE
*
      ELSE
*
*        Backward permutation
*
         DO 90 I = 1, M
*
            IF( K( I ).GT.0 )
     $         GO TO 80
*
            K( I ) = -K( I )
            J = K( I )
   60       CONTINUE
            IF( J.EQ.I )
     $         GO TO 80
*
            DO 70 JJ = 1, N
               TEMP = X( I, JJ )
               X( I, JJ ) = X( J, JJ )
               X( J, JJ ) = TEMP
   70       CONTINUE
*
            K( J ) = -K( J )
            J = K( J )
            GO TO 60
*
   80       CONTINUE
*
   90    CONTINUE
*
      END IF
*
      RETURN
*
*     End of ZLAPMT
*
      END

