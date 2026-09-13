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
!       ABSTRACT INTERFACE
!         SUBROUTINE XERBLA_INTERFACE(SRNAME, INFO)
!           CHARACTER*(*), INTENT(IN) :: SRNAME
!           INTEGER, INTENT(IN) :: INFO
!         END SUBROUTINE
!       END INTERFACE
!       PROCEDURE(XERBLA_INTERFACE) :: CB
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
!>     INTERFACE
!>       SUBROUTINE XERBLA_INTERFACE(SRNAME, INFO)
!>         CHARACTER*(*), INTENT(IN) :: SRNAME
!>         INTEGER, INTENT(IN) :: INFO
!>       END SUBROUTINE
!>       SUBROUTINE GET_LAPACK_XERBLA(CB)
!>         IMPORT :: XERBLA_INTERFACE
!>         IMPLICIT NONE
!>         PROCEDURE(XERBLA_INTERFACE), POINTER, INTENT(OUT) :: CB
!>       END SUBROUTINE
!>     END INTERFACE
!>     PROCEDURE(XERBLA_INTERFACE), POINTER :: ALREADY_CB
!>     CALL GET_LAPACK_XERBLA(ALREADY_CB)
!>   END PROGRAM HELLO
!> \endverbatim
!>
!> Finally, there is a SET_XERBLA in LAPACK which will override
!> both handlers.
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
module xerbla_lapack
  private
  public :: active_callback, xerbla_interface
  intrinsic          null
  abstract interface
    subroutine xerbla_interface(srname, info)
      character*(*), intent(in) :: srname
      integer, intent(in) :: info
    end subroutine
  end interface
  procedure(xerbla_interface), pointer :: active_callback => null()
end module xerbla_lapack

subroutine xerbla_lapack_sub(srname, info)
  use xerbla_lapack
  character*(*)      srname
  integer            info
  if (.not. associated(active_callback)) then
    print *, 'Error: LAPACK XERBLA called but no callback registered'
    stop
  end if
  call active_callback(srname, info)
end

subroutine set_lapack_xerbla(cb)
  use xerbla_lapack
  implicit none
  procedure(xerbla_interface) :: cb
  active_callback => cb
end

! A subroutine, not a function returning a procedure pointer: ifx
! miscompiles references to such a function when the caller declares
! it through an interface body (ICE in 2025.1, a call through a null
! pointer in 2026.1).
subroutine get_lapack_xerbla(cb)
  use xerbla_lapack
  implicit none
  procedure(xerbla_interface), pointer, intent(out) :: cb
  cb => active_callback
end

subroutine set_xerbla(cb)
  implicit none
  procedure() :: cb
  external set_lapack_xerbla
  external set_blas_xerbla
  call set_blas_xerbla(cb)
  call set_lapack_xerbla(cb)
end
