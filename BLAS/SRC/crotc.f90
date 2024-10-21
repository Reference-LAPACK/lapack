!> \brief \b CROTC applies a chain of rotation sequences to a matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!   subroutine crotc(side, dir, startup, shutdown, m, n, k,&
!                    A, lda, C, ldc, S, lds)
!    .. Scalar Arguments ..
!    integer, intent(in) :: m, n, k
!    ...
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CROTC applies a chain of k rotation sequences of length n to a matrix A.
!>
!> Each rotation is specified by a cosine and a sine, stored in the
!> matrices C and S respectively. Rotation G(i,j) is formed by
!> C(i,j) and S(i,j).
!>
!> If side = 'L', rotation G(i,j) is applied to rows i and i+1 of A.
!> [ A(i,j)   ] = [  C(i,j)        S(i,j) ] [ A(i,j)   ]
!> [ A(i+1,j) ]   [ -conjg(S(i,j))  C(i,j) ] [ A(i+1,j) ]
!> If side = 'R', rotation G(i,j) is applied to columns j and j+1 of A.
!> [ A(i,j)   A(i,j+1)   ] = [ A(i,j)   A(i,j+1)   ] [  C(i,j) -conjg(S(i,j)) ]
!> [ A(i+1,j) A(i+1,j+1) ]   [ A(i+1,j) A(i+1,j+1) ] [  S(i,j)  C(i,j)       ]
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] side
!> \verbatim
!>          side is CHARACTER*1
!>          If side = 'L', the rotations are applied to A from the left.
!>          If side = 'R', the rotations are applied to A from the right.
!> \endverbatim
!>
!> \param[in] dir
!> \verbatim
!>          dir is CHARACTER*1
!>          If dir = 'F', the rotations are applied in sequence from the
!>          first column/row to the last column/row.
!>          If dir = 'B', the rotations are applied in sequence from the
!>          last column/row to the first column/row.
!> \endverbatim
!>
!> \param[in] startup
!> \verbatim
!>          startup is LOGICAL
!>          If startup = .FALSE., the first (k-1) x (k-1) triangle
!>          of rotations is not applied.
!> \endverbatim
!>
!> \param[in] shutdown
!> \verbatim
!>          shutdown is LOGICAL
!>          If shutdown = .FALSE., the last (k-1) x (k-1) triangle
!>          of rotations is not applied.
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is INTEGER
!>          If side = 'L', m is the number of columns of A.
!>          If side = 'R', m is the number of rows of A.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The number of rotations in one sequence.
!> \endverbatim
!>
!> \param[in] k
!> \verbatim
!>          k is INTEGER
!>          The number of sequences of rotations.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array
!>          If side = 'L', A has dimension (n+1,m).
!>          If side = 'R', A has dimension (m,n+1).
!>          The matrix to which the rotations are applied.
!> \endverbatim
!>
!> \param[in] lda
!> \verbatim
!>          lda is INTEGER
!>          The leading dimension of A.
!>          If side = 'L', lda >= n+1.
!>          If side = 'R', lda >= m.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is REAL array, dimension (ldc,k)
!>          The matrix containing the cosines of the rotations.
!> \endverbatim
!>
!> \param[in] ldc
!> \verbatim
!>          ldc is INTEGER
!>          The leading dimension of C.
!>          ldc >= n.
!> \endverbatim
!>
!> \param[in,out] S
!> \verbatim
!>          S is COMPLEX array, dimension (lds,k)
!>          The matrix containing the sines of the rotations.
!> \endverbatim
!>
!> \param[in] lds
!> \verbatim
!>          lds is INTEGER
!>          The leading dimension of S.
!>          lds >= n.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Thijs Steel, KU Leuven, Belgium
!
!> \date October 2024
!
!> \ingroup rotc
!
subroutine crotc(side, dir, startup, shutdown, m, n, k,&
    A, lda, C, ldc, S, lds)
!    .. Scalar Arguments ..
    integer, intent(in) :: m, n, k, lda, ldc, lds
    character, intent(in) :: dir, side
    logical, intent(in) :: startup, shutdown
!    .. Array Arguments ..
    complex, intent(inout) :: A(lda,*)
    complex, intent(in) :: S(lds,*)
    real, intent(in) :: C(ldc,*)
!   .. Local Scalars ..
    integer i, j, l, j1, j2, incj, incj1, incj2, info
    complex temp, sn
    real cs
!    .. Executable Statements ..

!   Test the input parameters
    info = 0
    if(.not. (side .eq. 'L' .or. side .eq. 'R')) then
        info = 1
    end if
    if(.not. (dir .eq. 'F' .or. dir .eq. 'B')) then
        info = 2
    end if
    if(m .lt. 0) then
        info = 5
    end if
    if(n .lt. 0) then
        info = 6
    end if
    if(k .lt. 0) then
        info = 7
    end if
    if(side .eq. 'L') then
        if(lda .lt. n+1) then
            info = 9
        end if
    else
        if(lda .lt. m) then
            info = 9
        end if
    end if
    if(ldc .lt. n) then
        info = 11
    end if
    if(lds .lt. n) then
        info = 13
    end if

    if(info .ne. 0) then
        call xerbla('CROTC ', info)
        return
    end if

!   Determine ranges for loops around C and S
!   The range for sequence l is:
!   j1+(l-1)*incj1:incj:j2+(l-1)*incj2
    if( dir .eq. 'F') then
        incj = 1
        if(startup) then
            j1 = 1
            incj1 = 0
        else
            j1 = k
            incj1 = -1
        end if
        j2 = n
        if(shutdown) then
            incj2 = 0
        else
            incj2 = -1
        end if
    else
        incj = -1
        j1 = 1
        if(startup) then
            incj1 = 1
        else
            incj1 = 0
        end if
        if(shutdown) then
            j2 = 0
            incj2 = 0
        else
            j2 = n-k+1
            incj2 = 1
        end if
    end if

!   Apply the rotations
    if(side .eq. 'L') then
        do l = 1, k
            do j = j1+(l-1)*incj1, j2+(l-1)*incj2, incj
                cs = C(j,l)
                sn = S(j,l)
                do i = 1, m
                    temp = cs*A(i,j) + sn*A(i,j+1)
                    A(i,j+1) = -conjg(sn*A(i,j)) + cs*A(i,j+1)
                    A(i,j) = temp
                end do
            end do
        end do
    else
        do l = 1, k
            do j = j1+(l-1)*incj1, j2+(l-1)*incj2, incj
                cs = C(l,j)
                sn = S(l,j)
                do i = 1, m
                    temp = cs*A(j,i) + sn*A(j+1,i)
                    A(j+1,i) = -conjg(sn*A(j,i)) + cs*A(j+1,i)
                    A(j,i) = temp
                end do
            end do
        end do
    end if

end subroutine crotc
