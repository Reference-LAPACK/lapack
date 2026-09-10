*> \brief \b XERBLA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE XERBLA( SRNAME, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER*(*)      SRNAME
*       INTEGER            INFO
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> XERBLA  is an error handler for the LAPACK routines.
*> It is called by an LAPACK routine if an input parameter has an
*> invalid value.  A message is printed and execution stops.
*>
*> Users can replace the system XERBLA by calling SET_XERBLA
*> with a replacement handler that takes the same arguments, which
*> will then be called on error instead of the above behaviour. This
*> can then be negated by calling SET_XERBLA with NULL().
*> The current handler value can be retried by calling GET_XERBLA
*> with an output argument:
*>
*>   PROGRAM HELLO
*>     EXTERNAL GET_XERBLA
*>     PROCEDURE(), POINTER :: ALREADY_CB
*>     INTERFACE
*>       SUBROUTINE GET_XERBLA(CB_RET)
*>         PROCEDURE(), POINTER :: CB_RET
*>       END SUBROUTINE
*>     END INTERFACE
*>     CALL GET_XERBLA(ALREADY_CB)
*>   END PROGRAM HELLO
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
*>          The name of the routine which called XERBLA.
*> \endverbatim
*>
*> \param[in] INFO
*> \verbatim
*>          INFO is INTEGER
*>          The position of the invalid parameter in the parameter list
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
*
*> \ingroup xerbla
*
*  =====================================================================
      SUBROUTINE XERBLA( SRNAME, INFO )
#ifdef LAPACK_ILP64
#define CB_MODULE XERBLA_CALLBACKS_64
#define CB_CALLBACK ACTIVE_CALLBACK_64
#else
#define CB_MODULE XERBLA_CALLBACKS
#define CB_CALLBACK ACTIVE_CALLBACK
#endif
      USE CB_MODULE, ONLY: CB_CALLBACK
      IMPLICIT NONE
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
*     ..
*     .. Executable Statements ..
*
      IF (ASSOCIATED(CB_CALLBACK)) THEN
          CALL CB_CALLBACK(SRNAME, INFO)
          RETURN
      END IF
*
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
