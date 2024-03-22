subroutine lapack_c_slartg(f, g, c, s, r) bind(c, name='lapack_c_slartg')
    use iso_c_binding
    implicit none
    real(kind=c_float), intent(in), value :: f, g
    real(kind=c_float), intent(out) :: c, s, r

    external slartg

    ! Forward call to fortran implementation
    call slartg( f, g, c, s, r )

end subroutine lapack_c_slartg

subroutine lapack_c_dlartg(f, g, c, s, r) bind(c, name='lapack_c_dlartg')
    use iso_c_binding
    implicit none
    real(kind=c_double), intent(in), value :: f, g
    real(kind=c_double), intent(out) :: c, s, r

    external dlartg

    ! Forward call to fortran implementation
    call dlartg( f, g, c, s, r )

end subroutine lapack_c_dlartg

subroutine lapack_c_clartg(f, g, c, s, r) bind(c, name='lapack_c_clartg')
    use iso_c_binding
    implicit none
    complex(kind=c_float_complex), intent(in), value :: f, g
    complex(kind=c_float_complex), intent(out) :: s, r
    real(kind=c_float), intent(out) :: c

    external clartg

    ! Forward call to fortran implementation
    call clartg( f, g, c, s, r )

end subroutine lapack_c_clartg

subroutine lapack_c_zlartg(f, g, c, s, r) bind(c, name='lapack_c_zlartg')
    use iso_c_binding
    implicit none
    complex(kind=c_double_complex), intent(in), value :: f, g
    complex(kind=c_double_complex), intent(out) :: s, r
    real(kind=c_double), intent(out) :: c

    external zlartg

    ! Forward call to fortran implementation
    call zlartg( f, g, c, s, r )

end subroutine lapack_c_zlartg