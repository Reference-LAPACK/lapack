      INTEGER FUNCTION ILADIAG( DIAG )
*
*  -- LAPACK routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2008
*     .. Scalar Arguments ..
      CHARACTER          DIAG
*     ..
*
*  Purpose
*  =======
*
*  This subroutine translated from a character string specifying if a
*  matrix has unit diagonal or not to the relevant BLAST-specified
*  integer constant.
*
*  ILADIAG returns an INTEGER.  If ILADIAG < 0, then the input is not a
*  character indicating a unit or non-unit diagonal.  Otherwise ILADIAG
*  returns the constant value corresponding to DIAG.
*
*  Arguments
*  =========
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*  =====================================================================
*
*     .. Parameters ..
      INTEGER BLAS_NON_UNIT_DIAG, BLAS_UNIT_DIAG
      PARAMETER ( BLAS_NON_UNIT_DIAG = 131, BLAS_UNIT_DIAG = 132 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
      IF( LSAME( DIAG, 'N' ) ) THEN
         ILADIAG = BLAS_NON_UNIT_DIAG
      ELSE IF( LSAME( DIAG, 'U' ) ) THEN
         ILADIAG = BLAS_UNIT_DIAG
      ELSE
         ILADIAG = -1
      END IF
      RETURN
*
*     End of ILADIAG
*
      END
