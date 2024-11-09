subroutine lapack_c_slae2(a, b, c, rt1, rt2) bind(c, name='lapack_c_slae2')
    use iso_c_binding, only: c_float
    implicit none
    real(kind=c_float), intent(in), value :: a, b, c
    real(kind=c_float), intent(out) :: rt1, rt2

    external slae2

    ! Forward call to fortran implementation
    call slae2( a, b, c, rt1, rt2)

end subroutine lapack_c_slae2

subroutine lapack_c_dlae2(a, b, c, rt1, rt2) bind(c, name='lapack_c_dlae2')
    use iso_c_binding, only: c_double
    implicit none
    real(kind=c_double), intent(in), value :: a, b, c
    real(kind=c_double), intent(out) :: rt1, rt2

    external dlae2

    ! Forward call to fortran implementation
    call dlae2( a, b, c, rt1, rt2)

end subroutine lapack_c_dlae2