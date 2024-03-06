subroutine lapack_c_slaev2(a, b, c, rt1, rt2, cs1, sn1) bind(c, name='lapack_c_slaev2')
    use iso_c_binding, only: c_float
    implicit none
    real(kind=c_float), intent(in), value :: a, b, c
    real(kind=c_float), intent(out) :: rt1, rt2, cs1, sn1

    external slaev2

    ! Forward call to fortran implementation
    call slaev2( a, b, c, rt1, rt2, cs1, sn1)

end subroutine lapack_c_slaev2

subroutine lapack_c_dlaev2(a, b, c, rt1, rt2, cs1, sn1) bind(c, name='lapack_c_dlaev2')
    use iso_c_binding, only: c_double
    implicit none
    real(kind=c_double), intent(in), value :: a, b, c
    real(kind=c_double), intent(out) :: rt1, rt2, cs1, sn1

    external dlaev2

    ! Forward call to fortran implementation
    call dlaev2( a, b, c, rt1, rt2, cs1, sn1)

end subroutine lapack_c_dlaev2