subroutine clasr3 ( side, direct, m, n, k, C, ldc, S, lds, A, lda, work, lwork )
    use iso_c_binding, only: c_char
    implicit none
    ! Interface to C
    interface
        subroutine lapack_c_clasr3(side_c, direct_c, m_c, n_c, k_c, C_c,&
             ldc_c, S_c, lds_c, A_c, lda_c, work_c, lwork_c) bind(c, name="lapack_c_clasr3")
            use iso_c_binding, only: c_char, c_float, c_float_complex, c_int
            implicit none
            integer(c_int), intent(in) :: m_c,n_c,k_c,ldc_c, lds_c, lda_c, lwork_c
            character(c_char), intent(in) :: side_c, direct_c
            real(c_float), intent(in) :: C_c(ldc_c,*)
            complex(c_float_complex), intent(in) :: S_c(lds_c,*)
            complex(c_float_complex), intent(inout) :: A_c(lda_c,*), work_c(*)
        end subroutine lapack_c_clasr3
    end interface

    ! Arguments
    integer, intent(in) :: m,n,k,ldc, lds, lda, lwork
    character, intent(in) :: side, direct
    real, intent(in) :: C(ldc,*)
    complex, intent(in) :: S(lds,*)
    complex, intent(inout) :: A(lda,*), work(*)
    character(c_char) :: side_to_c, direct_to_c

    ! Convert characters to C characters
    side_to_c = side
    direct_to_c = direct

    ! Call C function
    call lapack_c_clasr3(side_to_c, direct_to_c, m, n, k, C, ldc, S, lds, A, lda, work, lwork)
end subroutine clasr3
