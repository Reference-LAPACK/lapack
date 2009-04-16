      INTEGER FUNCTION ILATRANS( TRANS )
*
*  -- LAPACK routine (version 3.2) --
*
*  -- April 2009                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
*     ..
*
*  Purpose
*  =======
*
*  This subroutine translates from a character string specifying a
*  transposition operation to the relevant BLAST-specified integer
*  constant.
*
*  ILATRANS returns an INTEGER.  If ILATRANS < 0, then the input is not
*  a character indicating a transposition operator.  Otherwise ILATRANS
*  returns the constant value corresponding to TRANS.
*
*  Arguments
*  =========
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  No transpose
*          = 'T':  Transpose
*          = 'C':  Conjugate transpose
*  =====================================================================
*
*     .. Parameters ..
      INTEGER BLAS_NO_TRANS, BLAS_TRANS, BLAS_CONJ_TRANS
      PARAMETER ( BLAS_NO_TRANS = 111, BLAS_TRANS = 112,
     $     BLAS_CONJ_TRANS = 113 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
      IF( LSAME( TRANS, 'N' ) ) THEN
         ILATRANS = BLAS_NO_TRANS
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
         ILATRANS = BLAS_TRANS
      ELSE IF( LSAME( TRANS, 'C' ) ) THEN
         ILATRANS = BLAS_CONJ_TRANS
      ELSE
         ILATRANS = -1
      END IF
      RETURN
*
*     End of ILATRANS
*
      END
