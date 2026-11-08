!> \brief \b GET_LAPACK_LIBNAME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE GET_LAPACK_LIBNAME(BUFFER, OUT_LEN)
!       CHARACTER(KIND=C_CHAR, LEN=*), INTENT(OUT) :: BUFFER
!       INTEGER,                       INTENT(OUT) :: OUT_LEN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> GET_LAPACK_LIBNAME writes a free-form, human-orientated string
!> identifying the active LAPACK library implementation into the
!> passed-in character buffer. The content is at the provider's
!> discretion. It has no mandated format, and should reflect the active
!> execution layer. For example, it may identify a dispatcher or
!> trampoline itself, or the underlying implementation providing
!> routines such as DGETRF, and may vary across configurations or
!> hardware targets as appropriate.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] BUFFER
!> \verbatim
!>          BUFFER will be filled with the string, including a
!>          terminating NUL byte for C compatibility. It will
!>          not go past the end of the provided buffer, and will
!>          always be NUL-terminated. However if a length-zero
!>          buffer is passed in, nothing will be written. If a
!>          longer-than-needed buffer is passed in, the rest will
!>          be filled with NULs.
!> \endverbatim
!
!> \param[out] OUT_LEN
!> \verbatim
!>          OUT_LEN will be set to the number of bytes actually
!>          written, including the terminating NUL byte but not
!>          any extra NULs. If the buffer had length zero, OUT_LEN
!>          will be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Ed J
!
!> \date November 2026
!
!> \ingroup ilaver
!
!  =====================================================================
subroutine get_lapack_libname(buffer, out_len)
  use, intrinsic :: iso_c_binding, only: c_char, c_null_char
  implicit none
  character(kind=c_char, len=*), intent(out) :: buffer
  integer,                       intent(out) :: out_len
  character(kind=c_char, len=*), parameter :: local_data = c_char_"Reference LAPACK"
  out_len = min(len(local_data) + 1, len(buffer))
  if (out_len > 1) buffer(1:out_len-1) = local_data(1:out_len-1)
  if (out_len > 0) buffer(out_len:) = c_null_char
end
