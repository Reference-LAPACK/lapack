      CHARACTER*1 FUNCTION CHLA_TRANSTYPE( TRANS )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     October 2008
*     .. Scalar Arguments ..
      INTEGER            TRANS
*     ..
*
*  Purpose
*  =======
*
*  This subroutine translates from a BLAST-specified integer constant to
*  the character string specifying a transposition operation.
*
*  CHLA_TRANSTYPE returns an CHARACTER*1.  If CHLA_TRANSTYPE is 'X',
*  then input is not an integer indicating a transposition operator.
*  Otherwise CHLA_TRANSTYPE returns the constant value corresponding to
*  TRANS.
*
*  Arguments
*  =========
*  TRANS   (input) INTEGER
*          Specifies the form of the system of equations:
*          = BLAS_NO_TRANS   = 111 :  No Transpose
*          = BLAS_TRANS      = 112 :  Transpose
*          = BLAS_CONJ_TRANS = 113 :  Conjugate Transpose
*  =====================================================================
*
*     .. Parameters ..
      INTEGER BLAS_NO_TRANS, BLAS_TRANS, BLAS_CONJ_TRANS
      PARAMETER ( BLAS_NO_TRANS = 111, BLAS_TRANS = 112,
     $     BLAS_CONJ_TRANS = 113 )
*     ..
*     .. Executable Statements ..
      IF( TRANS.EQ.BLAS_NO_TRANS ) THEN
         CHLA_TRANSTYPE = 'N'
      ELSE IF( TRANS.EQ.BLAS_TRANS ) THEN
         CHLA_TRANSTYPE = 'T'
      ELSE IF( TRANS.EQ.BLAS_CONJ_TRANS ) THEN
         CHLA_TRANSTYPE = 'C'
      ELSE
         CHLA_TRANSTYPE = 'X'
      END IF
      RETURN
*
*     End of CHLA_TRANSTYPE
*
      END
