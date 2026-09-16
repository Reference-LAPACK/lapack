!> \brief \b SET_BLAS_XERBLA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SET_BLAS_XERBLA(CB)
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
!> SET_BLAS_XERBLA overrides the BLAS XERBLA with a replacement subroutine.
!> This can then be negated by calling SET_BLAS_XERBLA with NULL().
!> The current handler value can be retrieved by calling GET_BLAS_XERBLA:
!>
!>   PROGRAM HELLO
!>     PROCEDURE(XERBLA_INTERFACE), POINTER :: ALREADY_CB
!>     INTERFACE
!>       SUBROUTINE XERBLA_INTERFACE(SRNAME, INFO)
!>         CHARACTER*(*), INTENT(IN) :: SRNAME
!>         INTEGER, INTENT(IN) :: INFO
!>       END SUBROUTINE
!>       FUNCTION GET_BLAS_XERBLA() RESULT(CB_RET)
!>         IMPLICIT NONE
!>         PROCEDURE(XERBLA_INTERFACE), POINTER :: CB_RET
!>       END FUNCTION
!>     END INTERFACE
!>     ALREADY_CB => GET_BLAS_XERBLA()
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
module blas_xerbla
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
end module blas_xerbla

subroutine xerbla_blas(srname, info)
  use blas_xerbla
  character*(*)      srname
  integer            info
  if (.not. associated(active_callback)) then
    print *, 'Error: BLAS XERBLA called but no callback registered'
    stop
  end if
  call active_callback(srname, info)
end

subroutine set_blas_xerbla(cb)
  use blas_xerbla
  implicit none
  procedure(xerbla_interface) :: cb
  active_callback => cb
end

function get_blas_xerbla() result(cb_ret)
  use blas_xerbla
  implicit none
  procedure(xerbla_interface), pointer :: cb_ret
  cb_ret => active_callback
end
