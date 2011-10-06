*> \brief \b ILAPREC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition
*  ==========
*
*       INTEGER FUNCTION ILAPREC( PREC )
* 
*       .. Scalar Arguments ..
*       CHARACTER          PREC
*       ..
*  
*  Purpose
*  =======
*
*>\details \b Purpose:
*>\verbatim
*>
*> This subroutine translated from a character string specifying an
*> intermediate precision to the relevant BLAST-specified integer
*> constant.
*>
*> ILAPREC returns an INTEGER.  If ILAPREC < 0, then the input is not a
*> character indicating a supported intermediate precision.  Otherwise
*> ILAPREC returns the constant value corresponding to PREC.
*>
*>\endverbatim
*
*  Arguments
*  =========
*
*
*  Authors
*  =======
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup auxOTHERcomputational
*
*  =====================================================================
      INTEGER FUNCTION ILAPREC( PREC )
*
*  -- LAPACK computational routine (version 3.2.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          PREC
*     ..
*
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
