*> \brief \b DLAQZ3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAQZ3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqz3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqz3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqz3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAQZ3 normalizes a generalized Schur decomposition and calculates the eigenvalues of the pencil.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ILSCHUR
*> \verbatim
*>          ILSCHUR is LOGICAL
*>              Determines whether or not to update the full Schur form
*> \endverbatim
*> \param[in] ILQ
*> \verbatim
*>          ILQ is LOGICAL
*>              Determines whether or not to update the matrix Q
*> \endverbatim
*>
*> \param[in] ILZ
*> \verbatim
*>          ILZ is LOGICAL
*>              Determines whether or not to update the matrix Z
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, Q, and Z.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          ILO and IHI mark the rows and columns of (A,B) which
*>          are to be normalized
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max( 1, N ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB, N)
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max( 1, N ).
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*> \endverbatim
*>
*> \param[out] ALPHAR
*> \verbatim
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)
*>          The real parts of each scalar alpha defining an eigenvalue
*>          of GNEP.
*> \endverbatim
*>
*> \param[out] ALPHAI
*> \verbatim
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)
*>          The imaginary parts of each scalar alpha defining an
*>          eigenvalue of GNEP.
*>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
*>          positive, then the j-th and (j+1)-st eigenvalues are a
*>          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j).
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>          The scalars beta that define the eigenvalues of GNEP.
*>          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
*>          beta = BETA(j) represent the j-th eigenvalue of the matrix
*>          pair (A,B), in one of the forms lambda = alpha/beta or
*>          mu = beta/alpha.  Since either lambda or mu may overflow,
*>          they should not, in general, be computed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          = 1: Normalization of a 2x2 block with real eigenvalues failed
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Thijs Steel
*
*> \date May 2020
*
*> \ingroup doubleGEcomputational
*>
*  =====================================================================
      subroutine dlaqz3(ilschur,ilq,ilz,n,ilo,ihi,A,ldA,B,ldB,Q,ldQ,Z,
     $   ldZ,alphar,alphai,beta,info)
      implicit none

*     Arguments
      logical,intent(in) :: ilschur,ilq,ilz
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ
      integer,intent(inout) :: info

      double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $   Z(ldZ,*),alphar(*),alphai(*),beta(*)

*     Parameters
      double precision :: zero,one,half
      parameter(zero=0.0d0,one=1.0d0,half=0.5d0)

*     Local scalars
      double precision :: smlnum,ulp,safmin,safmax,temp,c1,s1,c2,s2,H(2,
     $   3),b11,b22
      integer :: istartm,istopm,k,k2,i
      logical :: bulge

*     External Functions
      double precision,external :: dlamch

*     External Subroutines
      external            :: dlartg

      info = 0

*     Get machine constants
      safmin = dlamch('SAFE MINIMUM')
      safmax = one/safmin
      call dlabad(safmin,safmax)
      ulp = dlamch('precision')
      smlnum = safmin*(dble(n)/ulp)


*     Store eigenvalues and standardize the blocks
      k = ilo
      do while (k .le. ihi)

*        Determine range
         if (ilschur) then
            istartm = 1
            istopm = n
         else
            istartm = k
            istopm = k+1
         end if

         bulge = .false.
         if (k .lt. ihi) then
            if (A(k+1,k) .ne. zero) then
               bulge = .true.
            end if
         end if
         if (.not. bulge) then
*           1x1 eigenvalue block
            if( B(k,k) .lt. zero ) then
               do k2 = istartm,k
                  A(k2,k)=-A(k2,k)
                  B(k2,k)=-B(k2,k)
               end do
               if (ilz) then
                  do k2 = 1,n
                     Z(k2,k)=-Z(k2,k)
                  end do
               end if
            end if

            alphar(k) = A(k,k)
            alphai(k) = zero
            beta(k) = B(k,k)
            k = k+1
         else
*           2x2 eigenvalue block
            call dlag2(A(k,k),ldA,B(k,k),ldB,safmin,beta(k),beta(k+1),
     $         alphar(k),alphar(k+1),alphai(k))
            alphai(k+1) =-alphai(k)

            if (alphai(k) .eq. zero) then
*              2x2 block has real eigenvalues

               do i = 1,10

                  if(i .eq. 10) then
*                    Maximum number of iterations reached,
*                    semi-normalize the block and indicate an error
                     if(B(k+1,k).ne.zero)then
*                       Make B upper triangular
                        call dlartg(B(k,k),B(k+1,k),c1,s1,temp)
                        B(k,k) = temp
                        B(k+1,k) = zero

                        call drot(istopm-k+1,A(k,k),ldA,A(k+1,k),ldA,c1,
     $                     s1)
                        call drot(istopm-(k+1)+1,B(k,k+1),ldB,B(k+1,
     $                     k+1),ldB,c1,s1)
                        if (ilq) then
                           call drot(n,Q(1,k),1,Q(1,k+1),1,c1,s1)
                        end if
                     end if
                     info = 1
                     k = k+2
                     exit

                  end if

*                 Recalculate eigenvalues
                  call dlag2(A(k,k),ldA,B(k,k),ldB,safmin,beta(k),
     $               beta(k+1),alphar(k),alphar(k+1),alphai(k))
                  alphai(k+1) =-alphai(k)

                  H(1:2,1:2) = beta(k)*A(k:k+1,k:k+1)-alphar(k)*B(k:k+1,
     $               k:k+1)
                  if (H(2,2) .ge. H(1,2)) then
                     call dlartg(H(2,2),H(2,1),c1,s1,temp)
                  else
                     call dlartg(H(1,2),H(1,1),c1,s1,temp)
                  end if
                  call drot(2,H(1,2),1,H(1,1),1,c1,s1)

                  call drot(k+1-istartm+1,A(istartm,k+1),1,A(istartm,k),
     $               1,c1,s1)
                  call drot(k+1-istartm+1,B(istartm,k+1),1,B(istartm,k),
     $               1,c1,s1)
                  if (ilz) then
                     call drot(n,Z(1,k+1),1,Z(1,k),1,c1,s1)
                  end if

                  if (A(k+1,k) .ge. B(k+1,k)) then
                     call dlartg(A(k,k),A(k+1,k),c1,s1,temp)
                  else
                     call dlartg(B(k,k),B(k+1,k),c1,s1,temp)
                  end if

                  call drot(istopm-k+1,A(k,k),ldA,A(k+1,k),ldA,c1,s1)
                  call drot(istopm-k+1,B(k,k),ldB,B(k+1,k),ldB,c1,s1)
                  if (ilq) then
                     call drot(n,Q(1,k),1,Q(1,k+1),1,c1,s1)
                  end if

                  if (abs(A(k+1,k)) .lt. ulp*(abs(A(k,k))+abs(A(k+1,
     $               k+1))) .and. abs(B(k+1,k)) .lt. ulp*(abs(B(k,
     $               k))+abs(B(k+1,k+1)))) then
                     A(k+1,k) = zero
                     B(k+1,k) = zero
                     exit
                  end if
               end do
            else

*              Standardize, that is, rotate so that
*
*                       ( B11  0  )
*                   B = (         )  with B11 non-negative.
*                       (  0  B22 )
*
               if(B(k+1,k).ne.zero)then
*                 Make B upper triangular 
*                 (should already be the case, this is just to be sure)
                  call dlartg(B(k,k),B(k+1,k),c1,s1,temp)
                  B(k,k) = temp
                  B(k+1,k) = zero

                  call drot(istopm-k+1,A(k,k),ldA,A(k+1,k),ldA,c1,s1)
                  call drot(istopm-(k+1)+1,B(k,k+1),ldB,B(k+1,k+1),ldB,
     $               c1,s1)
                     if (ilq) then
                        call drot(n,Q(1,k),1,Q(1,k+1),1,c1,s1)
                     end if
               end if

               call dlasv2(B(k,k),B(k,k+1),B(k+1,k+1),b22,b11,s2,c2,s1,
     $            c1)

               if( b11.lt.zero ) then
                  c2 =-c2
                  s2 =-s2
                  b11 =-b11
                  b22 =-b22
               end if
               call drot(istopm-(k+2)+1,B(k,k+2),ldB,B(k+1,k+2),ldB,c1,
     $            s1)
               call drot(istopm-k+1,A(k,k),ldA,A(k+1,k),ldA,c1,s1)
               if(ilq) then
                  call drot(n,Q(1,k),1,Q(1,k+1),1,c1,s1)
               end if

               call drot(k-1-istartm+1,B(istartm,k),1,B(istartm,k+1),1,
     $            c2,s2)
               call drot(k+1-istartm+1,A(istartm,k),1,A(istartm,k+1),1,
     $            c2,s2)
               if(ilz) then
                  call drot(n,Z(1,k),1,Z(1,k+1),1,c2,s2)
               end if

               B(k,k) = b11
               B(k,k+1) = zero
               B(k+1,k) = zero
               B(k+1,k+1) = b22

               if( b22.lt.zero ) then
                  do k2 = istartm,k+1
                     A( k2,k+1 ) =-A( k2,k+1 )
                     B( k2,k+1 ) =-B( k2,k+1 )
                  end do
                  if (ilz) then
                     do k2 = 1,n
                        Z(k2,k+1)=-Z(k2,k+1)
                     end do
                  end if
                  b22 =-b22
               end if
               k = k+2

            end if
         end if
      end do

      end subroutine
