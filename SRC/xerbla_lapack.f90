!> \brief \b SET_LAPACK_XERBLA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SET_LAPACK_XERBLA(CB)
!       PROCEDURE(XERBLA_INTERFACE) :: CB
!       ABSTRACT INTERFACE
!         SUBROUTINE XERBLA_INTERFACE(SRNAME, INFO)
!           CHARACTER*(*), INTENT(IN) :: SRNAME
!           INTEGER, INTENT(IN) :: INFO
!         END SUBROUTINE
!       END INTERFACE
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SET_LAPACK_XERBLA overrides the LAPACK XERBLA with a replacement subroutine.
!> This can then be negated by calling SET_LAPACK_XERBLA with NULL().
!> The current handler value can be retrieved by calling GET_LAPACK_XERBLA:
!>
!>   PROGRAM HELLO
!>     PROCEDURE(XERBLA_INTERFACE), POINTER :: ALREADY_CB
!>     INTERFACE
!>       SUBROUTINE XERBLA_INTERFACE(SRNAME, INFO)
!>         CHARACTER*(*), INTENT(IN) :: SRNAME
!>         INTEGER, INTENT(IN) :: INFO
!>       END SUBROUTINE
!>       FUNCTION GET_LAPACK_XERBLA() RESULT(CB_RET)
!>         IMPLICIT NONE
!>         PROCEDURE(XERBLA_INTERFACE), POINTER :: CB_RET
!>       END FUNCTION
!>     END INTERFACE
!>     ALREADY_CB => GET_LAPACK_XERBLA()
!>   END PROGRAM HELLO
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CB
!> \verbatim
!>          CB is a pointer to a PROCEDURE that takes the same
!>          arguments as XERBLA.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Ed J
!
!> \date September 2026
!
!> \ingroup xerbla
!
!  =====================================================================
module lapack_xerbla
  private
  public :: active_callback, xerbla_interface
  intrinsic          null
  procedure(xerbla_interface), pointer :: active_callback => null()
  abstract interface
    subroutine xerbla_interface(srname, info)
      character*(*), intent(in) :: srname
      integer, intent(in) :: info
    end subroutine
  end interface
end module lapack_xerbla

subroutine xerbla_lapack(srname, info)
  use lapack_xerbla
  character*(*)      srname
  integer            info
  if (.not. associated(active_callback)) then
    print *, 'Error: LAPACK XERBLA called but no callback registered'
    stop
  end if
  call active_callback(srname, info)
end

subroutine set_lapack_xerbla(cb)
  use lapack_xerbla
  implicit none
  procedure(xerbla_interface) :: cb
  active_callback => cb
end

function get_lapack_xerbla() result(cb_ret)
  use lapack_xerbla
  implicit none
  procedure(xerbla_interface), pointer :: cb_ret
  cb_ret => active_callback
end
