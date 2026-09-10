#ifdef LAPACK_ILP64
#define CB_MODULE xerbla_callbacks_64
#define CB_CALLBACK active_callback_64
#else
#define CB_MODULE xerbla_callbacks
#define CB_CALLBACK active_callback
#endif

module CB_MODULE
  implicit none
  private
  public :: set_xerbla, get_xerbla, CB_CALLBACK
  procedure(xerbla_interface), pointer :: CB_CALLBACK => null()
  abstract interface
    subroutine xerbla_interface(srname, info)
      character*(*), intent(in) :: srname
      integer, intent(in) :: info
    end subroutine
  end interface
contains
  subroutine set_xerbla(cb)
    implicit none
    procedure(xerbla_interface) :: cb
    CB_CALLBACK => cb
  end subroutine set_xerbla
  subroutine get_xerbla(cb_ret)
    implicit none
    procedure(xerbla_interface), pointer :: cb_ret
    cb_ret => CB_CALLBACK
  end subroutine get_xerbla
end module CB_MODULE
