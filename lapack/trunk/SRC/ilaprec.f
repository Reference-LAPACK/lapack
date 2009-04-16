      INTEGER FUNCTION ILAPREC( PREC )
*
*  -- LAPACK routine (version 3.2.1)                                    --
*
*  -- April 2009                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          PREC
*     ..
*
*  Purpose
*  =======
*
*  This subroutine translated from a character string specifying an
*  intermediate precision to the relevant BLAST-specified integer
*  constant.
*
*  ILAPREC returns an INTEGER.  If ILAPREC < 0, then the input is not a
*  character indicating a supported intermediate precision.  Otherwise
*  ILAPREC returns the constant value corresponding to PREC.
*
*  Arguments
*  =========
*  PREC   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'S':  Single
*          = 'D':  Double
*          = 'I':  Indigenous
*          = 'X', 'E':  Extra
*  =====================================================================
*
*     .. Parameters ..
      INTEGER BLAS_PREC_SINGLE, BLAS_PREC_DOUBLE, BLAS_PREC_INDIGENOUS,
     $           BLAS_PREC_EXTRA
      PARAMETER ( BLAS_PREC_SINGLE = 211, BLAS_PREC_DOUBLE = 212,
     $     BLAS_PREC_INDIGENOUS = 213, BLAS_PREC_EXTRA = 214 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
      IF( LSAME( PREC, 'S' ) ) THEN
         ILAPREC = BLAS_PREC_SINGLE
      ELSE IF( LSAME( PREC, 'D' ) ) THEN
         ILAPREC = BLAS_PREC_DOUBLE
      ELSE IF( LSAME( PREC, 'I' ) ) THEN
         ILAPREC = BLAS_PREC_INDIGENOUS
      ELSE IF( LSAME( PREC, 'X' ) .OR. LSAME( PREC, 'E' ) ) THEN
         ILAPREC = BLAS_PREC_EXTRA
      ELSE
         ILAPREC = -1
      END IF
      RETURN
*
*     End of ILAPREC
*
      END
