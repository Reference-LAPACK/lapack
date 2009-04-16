      INTEGER FUNCTION ILAUPLO( UPLO )
*
*  -- LAPACK routine (version 3.2) --
*
*  -- April 2009                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
*     ..
*
*  Purpose
*  =======
*
*  This subroutine translated from a character string specifying a
*  upper- or lower-triangular matrix to the relevant BLAST-specified
*  integer constant.
*
*  ILAUPLO returns an INTEGER.  If ILAUPLO < 0, then the input is not
*  a character indicating an upper- or lower-triangular matrix.
*  Otherwise ILAUPLO returns the constant value corresponding to UPLO.
*
*  Arguments
*  =========
*  UPLO    (input) CHARACTER
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*  =====================================================================
*
*     .. Parameters ..
      INTEGER BLAS_UPPER, BLAS_LOWER
      PARAMETER ( BLAS_UPPER = 121, BLAS_LOWER = 122 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
      IF( LSAME( UPLO, 'U' ) ) THEN
         ILAUPLO = BLAS_UPPER
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         ILAUPLO = BLAS_LOWER
      ELSE
         ILAUPLO = -1
      END IF
      RETURN
*
*     End of ILAUPLO
*
      END
