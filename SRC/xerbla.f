*> \brief \b XERBLA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download XERBLA + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/xerbla.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/xerbla.f">
*> [TXT]</a>
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
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
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
      INTRINSIC          LEN_TRIM, NULL
*     ..
      PROCEDURE(XERBLA_INTERFACE), POINTER :: ACTIVE_CALLBACK => NULL()
      PROCEDURE(XERBLA_INTERFACE), POINTER :: CB_RET
      PROCEDURE(XERBLA_INTERFACE) :: CB
      ABSTRACT INTERFACE
        SUBROUTINE XERBLA_INTERFACE(SRNAME, INFO)
          CHARACTER*(*), INTENT(IN) :: SRNAME
          INTEGER, INTENT(IN) :: INFO
        END SUBROUTINE
      END INTERFACE
*     ..
*     .. Executable Statements ..
*
      IF (ASSOCIATED(ACTIVE_CALLBACK)) THEN
        CALL ACTIVE_CALLBACK(SRNAME, INFO)
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
      ENTRY SET_XERBLA(CB)
        ACTIVE_CALLBACK => CB
      RETURN
*
      ENTRY GET_XERBLA(CB_RET)
        CB_RET => ACTIVE_CALLBACK
      RETURN
*
      END
