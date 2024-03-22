subroutine dsteqr3 ( compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, info )
    use iso_c_binding, only: c_char
    implicit none
    ! Interface to C
    interface
        subroutine lapack_c_dsteqr3(compz_c, n_c, d_c, e_c, z_c, ldz_c, work_c,&
             lwork_c, rwork_c, lrwork_c, info_c) bind(c, name="lapack_c_dsteqr3")
            use iso_c_binding, only: c_char, c_double, c_double_complex, c_int
            implicit none
            integer(c_int), intent(in) :: n_c, ldz_c, lwork_c, lrwork_c
            integer(c_int), intent(out) :: info_c
            character(c_char), intent(in) :: compz_c
            real(c_double), intent(inout) :: d_c(*), e_c(*), rwork_c(*), z_c(ldz_c,*), work_c(*)
        end subroutine lapack_c_dsteqr3
    end interface

    ! Arguments
    integer, intent(in) :: n, ldz, lwork, lrwork
    integer, intent(out) :: info
    character, intent(in) :: compz
    double precision, intent(inout) :: d(*), e(*), rwork(*), z(ldz,*), work(*)
    character(c_char) :: compz_to_c

    ! Convert characters to C characters
    compz_to_c = compz

    ! Call C function
    call lapack_c_dsteqr3(compz_to_c, n, d, e, z, ldz, work, lwork, rwork, lrwork, info)
end subroutine dsteqr3
