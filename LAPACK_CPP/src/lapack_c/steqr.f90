subroutine lapack_c_ssteqr(compz, n, d, e, Z, ldz, work, info) bind(c, name='lapack_c_ssteqr')
    use iso_c_binding
    implicit none
    character(kind=c_char), intent(in), value :: compz
    integer(kind=c_int), intent(in), value :: n, ldz
    real(kind=c_float), intent(inout) :: d(*), e(*), Z(ldz,*)
    real(kind=c_float), intent(out) :: work(*)
    integer(kind=c_int), intent(out) :: info

    external ssteqr

    ! Forward call to fortran implementation
    call ssteqr( compz, n, d, e, Z, ldz, work, info )

end subroutine lapack_c_ssteqr

subroutine lapack_c_dsteqr(compz, n, d, e, Z, ldz, work, info) bind(c, name='lapack_c_dsteqr')
    use iso_c_binding
    implicit none
    character(kind=c_char), intent(in), value :: compz
    integer(kind=c_int), intent(in), value :: n, ldz
    real(kind=c_double), intent(inout) :: d(*), e(*), Z(ldz,*)
    real(kind=c_double), intent(out) :: work(*)
    integer(kind=c_int), intent(out) :: info

    external dsteqr

    ! Forward call to fortran implementation
    call dsteqr( compz, n, d, e, Z, ldz, work, info )

end subroutine lapack_c_dsteqr

subroutine lapack_c_csteqr(compz, n, d, e, Z, ldz, work, info) bind(c, name='lapack_c_csteqr')
    use iso_c_binding
    implicit none
    character(kind=c_char), intent(in), value :: compz
    integer(kind=c_int), intent(in), value :: n, ldz
    real(kind=c_float), intent(inout) :: d(*), e(*)
    complex(kind=c_float_complex), intent(inout) :: Z(ldz,*)
    real(kind=c_float), intent(out) :: work(*)
    integer(kind=c_int), intent(out) :: info

    external csteqr

    ! Forward call to fortran implementation
    call csteqr( compz, n, d, e, Z, ldz, work, info )

end subroutine lapack_c_csteqr

subroutine lapack_c_zsteqr(compz, n, d, e, Z, ldz, work, info) bind(c, name='lapack_c_zsteqr')
    use iso_c_binding
    implicit none
    character(kind=c_char), intent(in), value :: compz
    integer(kind=c_int), intent(in), value :: n, ldz
    real(kind=c_double), intent(inout) :: d(*), e(*)
    complex(kind=c_double_complex), intent(inout) :: Z(ldz,*)
    real(kind=c_double), intent(out) :: work(*)
    integer(kind=c_int), intent(out) :: info

    external zsteqr

    ! Forward call to fortran implementation
    call zsteqr( compz, n, d, e, Z, ldz, work, info )

end subroutine lapack_c_zsteqr