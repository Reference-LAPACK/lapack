*> \brief \b DLAQZ4
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAQZ4 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqz4.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqz4.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqz4.f">
*> [TXT]</a>
*> \endhtmlonly
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAQZ4 performs AED
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
*> \param[in] NW
*> \verbatim
*>          NW is INTEGER
*>          The desired size of the deflation window.
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
*> \param[out] NS
*> \verbatim
*>          NS is INTEGER
*>          The number of unconverged eigenvalues available to
*>          use as shifts.
*> \endverbatim
*>
*> \param[out] ND
*> \verbatim
*>          ND is INTEGER
*>          The number of converged eigenvalues found.
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
*> \param[in,out] QC
*> \verbatim
*>          QC is DOUBLE PRECISION array, dimension (LDQC, NW)
*> \endverbatim
*>
*> \param[in] LDQC
*> \verbatim
*>          LDQC is INTEGER
*> \endverbatim
*>
*> \param[in,out] ZC
*> \verbatim
*>          ZC is DOUBLE PRECISION array, dimension (LDZC, NW)
*> \endverbatim
*>
*> \param[in] LDZC
*> \verbatim
*>          LDZ is INTEGER
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
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
      subroutine dlaqz4(ilschur,ilq,ilz,n,ilo,ihi,nw,A,ldA,B,ldB,Q,ldQ,
     $   Z,ldZ,ns,nd,alphar,alphai,beta,Qc,ldQc,Zc,ldZc,work,lwork,
     $   info)
      implicit none

*     Arguments
      logical,intent(in) :: ilschur,ilq,ilz
      integer,intent(in) :: n,ilo,ihi,nw,ldA,ldB,ldQ,ldZ,ldQc,ldZc,
     $   lwork

      double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $   Z(ldZ,*),alphar(*),alphai(*),beta(*)
      integer,intent(out) :: ns,nd,info
      double precision :: Qc(ldQc,*),Zc(ldZc,*),work(*)

*     Parameters
      double precision :: zero,one,half
      parameter(zero=0.0d0,one=1.0d0,half=0.5d0)

*     Local Scalars
      logical :: bulge
      integer :: jw,kwtop,kwbot,istopm,istartm,k,k2,dtgexc_info,ifst,
     $   ilst,lworkreq,n_shifts,qz_small_info
      double precision :: s,smlnum,ulp,safmin,safmax,c1,s1,temp

*     External Functions
      double precision,external :: dlamch

      info = 0

*     Set up deflation window
      jw = min(nw,ihi-ilo+1)
      kwtop = ihi-jw+1
      if (kwtop .eq. ilo) then
         s = zero
      else
         s = A(kwtop,kwtop-1)
      end if

*     Determine required workspace
      ifst = 1
      ilst = jw
      call dtgexc(.true.,.true.,jw,A,ldA,B,ldB,Qc,ldQc,Zc,ldZc,ifst,
     $   ilst,work,-1,dtgexc_info)
      lworkreq = int(work(1))
      call dlaqz6('S','V','V',jw,1,jw,A(kwtop,kwtop),ldA,B(kwtop,kwtop),
     $   ldB,alphar,alphai,beta,Qc,ldQc,Zc,ldZc,work,-1,qz_small_info)
      lworkreq = max(lworkreq,int(work(1))+2*jw**2)
      lworkreq = max(lworkreq,n*nw,2*nw**2+n)
      if (lwork .eq.-1) then
*        workspace query, quick return
         work(1) = lworkreq
         return
      else if (lwork .lt. lworkreq) then
         info =-26
      end if

      if( info.NE.0 ) then
         CALL xerbla( 'DLAQZ4',-info )
         return
      end if

*     Get machine constants
      safmin = dlamch('SAFE MINIMUM')
      safmax = one/safmin
      call dlabad(safmin,safmax)
      ulp = dlamch('precision')
      smlnum = safmin*(dble(n)/ulp)

      if (ihi .eq. kwtop) then
*        1 by 1 deflation window, just try a regular deflation
         alphar(kwtop) = A(kwtop,kwtop)
         alphai(kwtop) = zero
         beta(kwtop) = B(kwtop,kwtop)
         ns = 1
         nd = 0
         if (abs(s) .le. max(smlnum,ulp*abs(A(kwtop,kwtop)))) then
            ns = 0
            nd = 1
            if (kwtop .gt. ilo) then
               A(kwtop,kwtop-1) = zero
            end if
         end if
      end if


*     Store window in case of convergence failure
      call dlacpy('ALL',jw,jw,A(kwtop,kwtop),ldA,work,jw)
      call dlacpy('ALL',jw,jw,B(kwtop,kwtop),ldB,work(jw**2+1),jw)

*     Transform window to real schur form
      call dlaset('Full',jw,jw,zero,one,Qc,ldQc)
      call dlaset('Full',jw,jw,zero,one,Zc,ldZc)
      call dlaqz6('S','V','V',jw,1,jw,A(kwtop,kwtop),ldA,B(kwtop,kwtop),
     $   ldB,alphar,alphai,beta,Qc,ldQc,Zc,ldZc,work(2*jw**2+1),
     $   lwork-2*jw**2,qz_small_info)

      if(qz_small_info .ne. 0) then
*        Convergence failure, restore the window and exit
         nd = 0
         ns = jw-qz_small_info
         call dlacpy('ALL',jw,jw,work,jw,A(kwtop,kwtop),ldA)
         call dlacpy('ALL',jw,jw,work(jw**2+1),jw,B(kwtop,kwtop),ldB)
         return
      end if 

*     Deflation detection loop
      if (kwtop .eq. ilo .or. s .eq. zero) then
         kwbot = kwtop-1
      else
         kwbot = ihi
         k = 1
         k2 = 1
         do while (k .le. jw)
            bulge = .false.
            if (kwbot-kwtop+1 .ge. 2) then
               bulge = A(kwbot,kwbot-1) .ne. zero
            end if
            if (bulge) then

*              Try to deflate complex conjugate eigenvalue pair
               temp = abs(A(kwbot,kwbot))+sqrt( abs(A(kwbot,kwbot-1)) )*
     $            sqrt( abs(A(kwbot-1,kwbot)) )
               if(temp .eq. zero)then
                  temp = abs(s)
               end if
               if (max(abs(s*Qc(1,kwbot-kwtop)),abs(s*Qc(1,kwbot-kwtop+
     $            1))) .le. max(smlnum,ulp*temp) ) then
*                 Deflatable
                  kwbot = kwbot-2
               else
*                 Not deflatable, move out of the way
                  ifst = kwbot-kwtop+1
                  ilst = k2
                  call dtgexc(.true.,.true.,jw,A(kwtop,kwtop),ldA,
     $               B(kwtop,kwtop),ldB,Qc,ldQc,Zc,ldZc,ifst,ilst,work,
     $               lwork,dtgexc_info)
                  k2 = k2+2
               end if
               k = k+2
            else

*              Try to deflate real eigenvalue
               temp = abs(A(kwbot,kwbot))
               if(temp .eq. zero) then
                  temp = abs(s)
               end if
               if ((abs(s*Qc(1,kwbot-kwtop+1))) .le. max(ulp*temp,
     $            smlnum)) then
*                 Deflatable
                  kwbot = kwbot-1
               else
*                 Not deflatable, move out of the way
                  ifst = kwbot-kwtop+1
                  ilst = k2
                  call dtgexc(.true.,.true.,jw,A(kwtop,kwtop),ldA,
     $              B(kwtop,kwtop),ldB,Qc,ldQc,Zc,ldZc,ifst,ilst,work,
     $              lwork,dtgexc_info)
                  k2 = k2+1
               end if

               k = k+1

            end if
         end do
      end if

*     Store eigenvalues
      nd = ihi-kwbot
      ns = jw-nd
      k = kwtop
      do while (k .le. ihi)
         bulge = .false.
         if (k .lt. ihi) then
            if (A(k+1,k) .ne. zero) then
               bulge = .true.
            end if
         end if
         if (bulge) then
*           2x2 eigenvalue block
            call dlag2(A(k,k),ldA,B(k,k),ldB,safmin,beta(k),beta(k+1),
     $         alphar(k),alphar(k+1),alphai(k))
            alphai(k+1) =-alphai(k)
            k = k+2
         else
*           1x1 eigenvalue block
            alphar(k) = A(k,k)
            alphai(k) = zero
            beta(k) = B(k,k)
            k = k+1
         end if
      end do

      if (kwtop .ne. ilo .and. s .ne. zero) then
*        Reflect spike back, this will create optimally packed bulges
         A(kwtop:kwbot,kwtop-1) = A(kwtop,kwtop-1)*Qc(1,1:jw-nd)
         do k = kwbot-1,kwtop,-1
            call dlartg(A(k,kwtop-1),A(k+1,kwtop-1),c1,s1,temp)
            A(k,kwtop-1) = temp
            A(k+1,kwtop-1) = zero
            k2 = max(kwtop,k-1)
            call drot(ihi-k2+1,A(k,k2),ldA,A(k+1,k2),ldA,c1,s1)
            call drot(ihi-(k-1)+1,B(k,k-1),ldB,B(k+1,k-1),ldB,c1,s1)
            call drot(jw,Qc(1,k-kwtop+1),1,Qc(1,k+1-kwtop+1),1,c1,s1)
         end do

*        Chase bulges down
         istartm = kwtop
         istopm = ihi
         k = kwbot-1
         do while (k .ge. kwtop)
            if ((k .ge. kwtop+1) .and. A(k+1,k-1) .ne. zero) then

*              Move double pole block down and remove it
               do k2 = k-1,kwbot-2
                  call dlaqz2(.true.,.true.,k2,kwtop,kwtop+jw-1,kwbot,A,
     $               ldA,B,ldB,jw,kwtop,Qc,ldQc,jw,kwtop,Zc,ldZc)
               end do

               k = k-2
            else

*              k points to single shift
               do k2 = k,kwbot-2

*                 Move shift down
                  call dlartg(B(k2+1,k2+1),B(k2+1,k2),c1,s1,temp)
                  B(k2+1,k2+1) = temp
                  B(k2+1,k2) = zero
                  call drot(k2+2-istartm+1,A(istartm,k2+1),1,A(istartm,
     $               k2),1,c1,s1)
                  call drot(k2-istartm+1,B(istartm,k2+1),1,B(istartm,
     $               k2),1,c1,s1)
                  call drot(jw,Zc(1,k2+1-kwtop+1),1,Zc(1,k2-kwtop+1),1,
     $               c1,s1)
            
                  call dlartg(A(k2+1,k2),A(k2+2,k2),c1,s1,temp)
                  A(k2+1,k2) = temp
                  A(k2+2,k2) = zero
                  call drot(istopm-k2,A(k2+1,k2+1),ldA,A(k2+2,k2+1),ldA,
     $               c1,s1)
                  call drot(istopm-k2,B(k2+1,k2+1),ldB,B(k2+2,k2+1),ldB,
     $               c1,s1)
                  call drot(jw,Qc(1,k2+1-kwtop+1),1,Qc(1,k2+2-kwtop+1),
     $               1,c1,s1)

               end do

*              Remove the shift
               call dlartg(B(kwbot,kwbot),B(kwbot,kwbot-1),c1,s1,temp)
               B(kwbot,kwbot) = temp
               B(kwbot,kwbot-1) = zero
               call drot(kwbot-istartm,B(istartm,kwbot),1,B(istartm,
     $            kwbot-1),1,c1,s1)
               call drot(kwbot-istartm+1,A(istartm,kwbot),1,A(istartm,
     $            kwbot-1),1,c1,s1)
               call drot(jw,Zc(1,kwbot-kwtop+1),1,Zc(1,kwbot-1-kwtop+1),
     $            1,c1,s1)

               k = k-1
            end if
         end do

      end if

*     Apply Qc and Zc to rest of the matrix
      if (ilschur) then
         istartm = 1
         istopm = n
      else
         istartm = ilo
         istopm = ihi
      end if

      if (istopm-ihi > 0) then
         call dgemm('T','N',jw,istopm-ihi,jw,one,Qc,ldQc,A(kwtop,ihi+1),
     $      ldA,zero,work,jw)
         call dlacpy('ALL',jw,istopm-ihi,work,jw,A(kwtop,ihi+1),ldA)
         call dgemm('T','N',jw,istopm-ihi,jw,one,Qc,ldQc,B(kwtop,ihi+1),
     $      ldB,zero,work,jw)
         call dlacpy('ALL',jw,istopm-ihi,work,jw,B(kwtop,ihi+1),ldB)
      end if
      if (ilq) then
         call dgemm('N','N',n,jw,jw,one,Q(1,kwtop),ldQ,Qc,ldQc,zero,
     $      work,n)
         call dlacpy('ALL',n,jw,work,n,Q(1,kwtop),ldQ)
      end if

      if (kwtop-1-istartm+1 > 0) then
         call dgemm('N','N',kwtop-istartm,jw,jw,one,A(istartm,kwtop),
     $      ldA,Zc,ldZc,zero,work,kwtop-istartm)
        call dlacpy('ALL',kwtop-istartm,jw,work,kwtop-istartm,A(istartm,
     $     kwtop),ldA)
         call dgemm('N','N',kwtop-istartm,jw,jw,one,B(istartm,kwtop),
     $      ldB,Zc,ldZc,zero,work,kwtop-istartm)
        call dlacpy('ALL',kwtop-istartm,jw,work,kwtop-istartm,B(istartm,
     $     kwtop),ldB)
      end if
      if (ilz) then
         call dgemm('N','N',n,jw,jw,one,Z(1,kwtop),ldZ,Zc,ldZc,zero,
     $      work,n)
         call dlacpy('ALL',n,jw,work,n,Z(1,kwtop),ldZ)
      end if

      end subroutine
