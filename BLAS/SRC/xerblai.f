*> \brief \b XERBLAI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE XERBLAI( SRNAME, INFO, INDX )
*
*       .. Scalar Arguments ..
*       CHARACTER*(*)      SRNAME
*       INTEGER            INFO, INDX
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> XERBLAI is an error handler for the BLAS/LAPACK routines.
*> It is called by an BLAS/LAPACK routine if an input parameter has an
*> invalid value.  A message is printed and execution stops.
*>
*> Installers may consider modifying the STOP statement in order to
*> call system-specific exception-handling facilities.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SRNAME
*> \verbatim
*>          SRNAME is CHARACTER*(*)
*>          The name of the routine which called XERBLAI.
*> \endverbatim
*>
*> \param[in] INFO
*> \verbatim
*>          INFO is INTEGER
*>          The position of the invalid parameter in the parameter list
*>          of the calling routine.
*> \endverbatim
*>
*> \param[in] INDX
*> \verbatim
*>          INDX is INTEGER
*>          The position at the invalid parameter in the parameter list
*>          of the calling routine.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*> \author Igor S. Gerasimov
*
*> \ingroup xerblai
*
*  =====================================================================
      SUBROUTINE XERBLAI( SRNAME, INFO, INDX )
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO, INDX
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
*     ..
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO, INDX
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ',
     $      'an illegal value at position ', I0, '.' )
*
*     End of XERBLAI
*
      END
