#ifdef USE_FORTRAN_BLAS
subroutine lapack_c_drot(n, dx ,incx, dy, incy, c, s) bind(c, name='lapack_c_drot')
    use iso_c_binding
    implicit none
    integer(kind=c_int), value, intent(in) :: n, incx, incy
    real(kind=c_double), value, intent(in) :: c, s
    real(kind=c_double), dimension(*), intent(inout) :: dx, dy

    external drot

    ! Forward call to fortran implementation
    call drot( n, dx ,incx, dy, incy, c, s )

end subroutine lapack_c_drot

subroutine lapack_c_srot(n, dx ,incx, dy, incy, c, s) bind(c, name='lapack_c_srot')
    use iso_c_binding
    implicit none
    integer(kind=c_int), value, intent(in) :: n, incx, incy
    real(kind=c_float), value, intent(in) :: c, s
    real(kind=c_float), dimension(*), intent(inout) :: dx, dy

    external srot

    ! Forward call to fortran implementation
    call srot( n, dx ,incx, dy, incy, c, s )

end subroutine lapack_c_srot

subroutine lapack_c_crot(n, dx ,incx, dy, incy, c, s) bind(c, name='lapack_c_crot')
    use iso_c_binding
    implicit none
    integer(kind=c_int), value, intent(in) :: n, incx, incy
    real(kind=c_float), value, intent(in) :: c
    complex(kind=c_float_complex), value, intent(in) :: s
    complex(kind=c_float_complex), dimension(*), intent(inout) :: dx, dy

    external crot

    ! Forward call to fortran implementation
    call crot( n, dx ,incx, dy, incy, c, s )

end subroutine lapack_c_crot

subroutine lapack_c_zrot(n, dx ,incx, dy, incy, c, s) bind(c, name='lapack_c_zrot')
    use iso_c_binding
    implicit none
    integer(kind=c_int), value, intent(in) :: n, incx, incy
    real(kind=c_double), value, intent(in) :: c
    complex(kind=c_double_complex), value, intent(in) :: s
    complex(kind=c_double_complex), dimension(*), intent(inout) :: dx, dy

    external zrot

    ! Forward call to fortran implementation
    call zrot( n, dx ,incx, dy, incy, c, s )

end subroutine lapack_c_zrot
#endif