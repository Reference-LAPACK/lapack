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
!> SET_BLAS_XERBLA overrides the BLAS XERBLA with a replacement subroutine.
!> This can then be negated by calling SET_BLAS_XERBLA with a NULL value,
!> although this is only possible via a C interface as it is forbidden
!> in the Fortran standard.
!>
!> Users should usually also use SET_LAPACK_XERBLA since the error
!> handlers may be separate due to the vagaries of linkage.
!> GET_BLAS_XERBLA returns the current BLAS XERBLA. It always returns an
!> associated value, so can always be passed to SET_BLAS_XERBLA.
!> The current handler value can be retrieved by calling GET_BLAS_XERBLA:
!>
!>   PROGRAM HELLO
!>     INTERFACE
!>       SUBROUTINE XERBLA_INTERFACE(SRNAME, INFO)
!>         CHARACTER*(*), INTENT(IN) :: SRNAME
!>         INTEGER, INTENT(IN) :: INFO
!>       END SUBROUTINE
!>       SUBROUTINE GET_BLAS_XERBLA(CB)
!>         IMPORT :: XERBLA_INTERFACE
!>         IMPLICIT NONE
!>         PROCEDURE(XERBLA_INTERFACE), POINTER, INTENT(OUT) :: CB
!>       END SUBROUTINE
!>     END INTERFACE
!>     PROCEDURE(XERBLA_INTERFACE), POINTER :: ALREADY_CB
!>     CALL GET_BLAS_XERBLA(ALREADY_CB)
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
module xerbla_blas
  private
  public :: active_callback, default_blas_xerbla, xerbla_interface
  intrinsic          null
  abstract interface
    subroutine xerbla_interface(srname, info)
      character*(*), intent(in) :: srname
      integer, intent(in) :: info
    end subroutine
  end interface
  procedure(xerbla_interface), pointer :: active_callback => null()
contains
  ! The default handler, the historical XERBLA. It is a module procedure,
  ! so it has no fixed external name; GET_BLAS_XERBLA returns it while no
  ! other handler is installed.
  subroutine default_blas_xerbla(srname, info)
    character*(*), intent(in) :: srname
    integer, intent(in) :: info
    intrinsic          len_trim
    write (*, fmt = 9999) srname(1:len_trim(srname)), info
    stop
9999 format (' ** On entry to ', a, ' parameter number ', i2, ' had ', &
             'an illegal value')
  end subroutine
end module xerbla_blas

! Called by XERBLA: runs the installed handler, or the default.
subroutine xerbla_blas_sub(srname, info)
  use xerbla_blas
  character*(*)      srname
  integer            info
  if (associated(active_callback)) then
    call active_callback(srname, info)
  else
    call default_blas_xerbla(srname, info)
  end if
end

subroutine set_blas_xerbla(cb)
  use xerbla_blas
  implicit none
  procedure(xerbla_interface) :: cb
  active_callback => cb
end

! A subroutine, not a function returning a procedure pointer: ifx
! miscompiles references to such a function when the caller declares
! it through an interface body (ICE in 2025.1, a call through a null
! pointer in 2026.1).
subroutine get_blas_xerbla(cb) ! never returns disassociated pointer
  use xerbla_blas
  implicit none
  procedure(xerbla_interface), pointer :: cb
  if (associated(active_callback)) then
    cb => active_callback
  else
    cb => default_blas_xerbla
  end if
end
