subroutine lapack_c_slasrt(id, n, d, info) bind(c, name='lapack_c_slasrt')
    use iso_c_binding, only: c_char, c_int, c_float
    implicit none
    character(kind=c_char), intent(in), value :: id
    integer(kind=c_int), intent(in), value :: n
    real(kind=c_float), intent(inout) :: d(*)
    integer(kind=c_int), intent(out) :: info

    character :: id_fortran

    external slasrt

    ! Convert the character to fortran
    id_fortran = id
    ! Forward call to fortran implementation
    call slasrt( id_fortran, n, d, info)

end subroutine lapack_c_slasrt

subroutine lapack_c_dlasrt(id, n, d, info) bind(c, name='lapack_c_dlasrt')
    use iso_c_binding, only: c_char, c_int, c_double
    implicit none
    character(kind=c_char), intent(in), value :: id
    integer(kind=c_int), intent(in), value :: n
    real(kind=c_double), intent(inout) :: d(*)
    integer(kind=c_int), intent(out) :: info

    character :: id_fortran

    external dlasrt

    ! Convert the character to fortran
    id_fortran = id
    ! Forward call to fortran implementation
    call dlasrt( id_fortran, n, d, info)

end subroutine lapack_c_dlasrt