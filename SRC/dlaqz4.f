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
*> DLAQZ4 Executes a single multishift QZ sweep
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
*> \endverbatim
*>
*> \param[in] NSHIFTS
*> \verbatim
*>          NSHIFTS is INTEGER
*>          The desired number of shifts to use
*> \endverbatim
*>
*> \param[in] NBLOCK_DESIRED
*> \verbatim
*>          NBLOCK_DESIRED is INTEGER
*>          The desired size of the computational windows
*> \endverbatim
*>
*> \param[in] SR
*> \verbatim
*>          SR is DOUBLE PRECISION array. SR contains
*>          the real parts of the shifts to use.
*> \endverbatim
*>
*> \param[in] SI
*> \verbatim
*>          SI is DOUBLE PRECISION array. SI contains
*>          the imaginary parts of the shifts to use.
*> \endverbatim
*>
*> \param[in] SS
*> \verbatim
*>          SS is DOUBLE PRECISION array. SS contains
*>          the scale of the shifts to use.
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
*> \param[in,out] QC
*> \verbatim
*>          QC is DOUBLE PRECISION array, dimension (LDQC, NBLOCK_DESIRED)
*> \endverbatim
*>
*> \param[in] LDQC
*> \verbatim
*>          LDQC is INTEGER
*> \endverbatim
*>
*> \param[in,out] ZC
*> \verbatim
*>          ZC is DOUBLE PRECISION array, dimension (LDZC, NBLOCK_DESIRED)
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
      subroutine dlaqz4(ilschur,ilq,ilz,n,ilo,ihi,nshifts,
     $   nblock_desired,sr,si,ss,A,ldA,B,ldB,Q,ldQ,Z,ldZ,Qc,ldQc,Zc,
     $   ldZc,work,lwork,info)
      implicit none

*     Function arguments
      logical,intent(in) :: ilschur,ilq,ilz
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ,lwork,nshifts,
     $   nblock_desired,ldQc,ldZc

      double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $   Z(ldZ,*),Qc(ldQc,*),Zc(ldZc,*),work(*),sr(*),si(*),ss(*)

      integer,intent(out) :: info

*     Parameters
      double precision :: zero,one,half
      parameter(zero=0.0d0,one=1.0d0,half=0.5d0)

*     Local scalars
      integer :: i,j,ns,istartm,istopm,sheight,swidth,k,np,istartb,
     $   istopb,ishift,nblock,npos
      double precision :: temp,v(3),c1,s1,c2,s2,H(2,3),swap

      info = 0
      if (nblock_desired .lt. nshifts+1) then
         info =-8
      end if
      if (lwork .eq.-1) then
*        workspace query, quick return
         work(1) = n*nblock_desired
         return
      else if (lwork .lt. n*nblock_desired) then
         info =-25
      end if

      if( info.NE.0 ) then
         call xerbla( 'DLAQZ4',-info )
         return
      end if

*     Executable statements

      if (nshifts .lt. 2) then
         return
      end if

      if (ilo .ge. ihi) then
         return
      end if

      if (ilschur) then
         istartm = 1
         istopm = n
      else
         istartm = ilo
         istopm = ihi
      end if

*     Shuffle shifts into pairs of real shifts and pairs
*     of complex conjugate shifts assuming complex
*     conjugate shifts are already adjacent to one
*     another

      do i = 1,nshifts-2,2
         if( si( i ).NE.-si( i+1 ) ) then
*
            swap = sr( i )
            sr( i ) = sr( i+1 )
            sr( i+1 ) = sr( i+2 )
            sr( i+2 ) = swap

            swap = si( i )
            si( i ) = si( i+1 )
            si( i+1 ) = si( i+2 )
            si( i+2 ) = swap
            
            swap = ss( i )
            ss( i ) = ss( i+1 )
            ss( i+1 ) = ss( i+2 )
            ss( i+2 ) = swap
         end if
      end do

*     NSHFTS is supposed to be even, but if it is odd,
*     then simply reduce it by one.  The shuffle above
*     ensures that the dropped shift is real and that
*     the remaining shifts are paired.

      ns = nshifts-mod(nshifts,2)
      npos = max(nblock_desired-ns,1)

*     The following block introduces the shifts and chases
*     them down one by one just enough to make space for
*     the other shifts. The near-the-diagonal block is
*     of size (ns+1) x ns.

      call dlaset('Full',ns+1,ns+1,zero,one,Qc,ldQc)
      call dlaset('Full',ns,ns,zero,one,Zc,ldZc)

      do i = 1,ns,2
*        Introduce the shift
         call dlaqz1(A(ilo,ilo),ldA,B(ilo,ilo),ldB,sr(i),sr(i+1),si(i),
     $      ss(i),ss(i+1),v)

         temp = v(2)
         call dlartg(temp,v(3),c1,s1,v(2))
         call dlartg(v(1),v(2),c2,s2,temp)

         call drot(ns,A(ilo+1,ilo),ldA,A(ilo+2,ilo),ldA,c1,s1)
         call drot(ns,A(ilo,ilo),ldA,A(ilo+1,ilo),ldA,c2,s2)
         call drot(ns,B(ilo+1,ilo),ldB,B(ilo+2,ilo),ldB,c1,s1)
         call drot(ns,B(ilo,ilo),ldB,B(ilo+1,ilo),ldB,c2,s2)
         call drot(ns+1,Qc(1,2),1,Qc(1,3),1,c1,s1)
         call drot(ns+1,Qc(1,1),1,Qc(1,2),1,c2,s2)

*        Chase the shift down
         do j = 1,ns-1-i

            call dlaqz2(.true.,.true.,j,1,ns,ihi-ilo+1,A(ilo,ilo),ldA,
     $         B(ilo,ilo),ldB,ns+1,1,Qc,ldQc,ns,1,Zc,ldZc)

         end do

      end do

*     Update the rest of the pencil

*     Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm)
*     from the left with Qc(1:ns+1,1:ns+1)'
      sheight = ns+1
      swidth = istopm-(ilo+ns)+1
      if (swidth > 0) then
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,A(ilo,
     $      ilo+ns),ldA,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,A(ilo,ilo+ns),
     $      ldA)
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,B(ilo,
     $      ilo+ns),ldB,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,B(ilo,ilo+ns),
     $      ldB)
      end if
      if (ilq) then
        call dgemm('N','N',n,sheight,sheight,one,Q(1,ilo),ldQ,Qc,ldQc,
     $     zero,work,n)
         call dlacpy('ALL',n,sheight,work,n,Q(1,ilo),ldQ)
      end if

*     Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1)
*     from the right with Zc(1:ns,1:ns)
      sheight = ilo-1-istartm+1
      swidth = ns
      if (sheight > 0) then
         call dgemm('N','N',sheight,swidth,swidth,one,A(istartm,ilo),
     $      ldA,Zc,ldZc,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,A(istartm,ilo),
     $      ldA)
         call dgemm('N','N',sheight,swidth,swidth,one,B(istartm,ilo),
     $      ldB,Zc,ldZc,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,B(istartm,ilo),
     $      ldB)
      end if
      if (ilz) then
         call dgemm('N','N',n,swidth,swidth,one,Z(1,ilo),ldZ,Zc,ldZc,
     $      zero,work,n)
         call dlacpy('ALL',n,swidth,work,n,Z(1,ilo),ldZ)
      end if

*     The following block chases the shifts down to the bottom
*     right block. If possible, a shift is moved down npos
*     positions at a time

      k = ilo
      do while (k < ihi-ns)
         np = min(ihi-ns-k,npos)
*        Size of the near-the-diagonal block
         nblock = ns+np
*        istartb points to the first row we will be updating
         istartb = k+1
*        istopb points to the last column we will be updating
         istopb = k+nblock-1

         call dlaset('Full',ns+np,ns+np,zero,one,Qc,ldQc)
         call dlaset('Full',ns+np,ns+np,zero,one,Zc,ldZc)

*        Near the diagonal shift chase
         do i = ns-1,1,-2
            do j = 0,np-1
*              Move down the block with index k+i+j-1, updating
*              the (ns+np x ns+np) block:
*              (k:k+ns+np,k:k+ns+np-1)
               call dlaqz2(.true.,.true.,k+i+j-1,istartb,istopb,ihi,A,
     $            ldA,B,ldB,nblock,k+1,Qc,ldQc,nblock,k,Zc,ldZc)
            end do
         end do

*        Update rest of the pencil

*        Update A(k+1:k+ns+np, k+ns+np:istopm) and
*        B(k+1:k+ns+np, k+ns+np:istopm)
*        from the left with Qc(1:ns+np,1:ns+np)'
         sheight = ns+np
         swidth = istopm-(k+ns+np)+1
         if (swidth > 0) then
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,A(k+1,
     $      k+ns+np),ldA,zero,work,sheight)
            call dlacpy('ALL',sheight,swidth,work,sheight,A(k+1,
     $         k+ns+np),ldA)
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,B(k+1,
     $      k+ns+np),ldB,zero,work,sheight)
            call dlacpy('ALL',sheight,swidth,work,sheight,B(k+1,
     $         k+ns+np),ldB)
         end if
         if (ilq) then
        call dgemm('N','N',n,nblock,nblock,one,Q(1,k+1),ldQ,Qc,ldQc,
     $     zero,work,n)
            call dlacpy('ALL',n,nblock,work,n,Q(1,k+1),ldQ)
         end if

*        Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1)
*        from the right with Zc(1:ns+np,1:ns+np)
         sheight = k-istartm+1
         swidth = nblock
         if (sheight > 0) then
            call dgemm('N','N',sheight,swidth,swidth,one,A(istartm,k),
     $         ldA,Zc,ldZc,zero,work,sheight)
            call dlacpy('ALL',sheight,swidth,work,sheight,A(istartm,k),
     $         ldA)
            call dgemm('N','N',sheight,swidth,swidth,one,B(istartm,k),
     $         ldB,Zc,ldZc,zero,work,sheight)
            call dlacpy('ALL',sheight,swidth,work,sheight,B(istartm,k),
     $         ldB)
         end if
         if (ilz) then
            call dgemm('N','N',n,nblock,nblock,one,Z(1,k),ldZ,Zc,ldZc,
     $         zero,work,n)
            call dlacpy('ALL',n,nblock,work,n,Z(1,k),ldZ)
         end if

         k = k+np

      end do

*     The following block removes the shifts from the bottom right corner
*     one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi).

      call dlaset('Full',ns,ns,zero,one,Qc,ldQc)
      call dlaset('Full',ns+1,ns+1,zero,one,Zc,ldZc)

*     istartb points to the first row we will be updating
      istartb = ihi-ns+1
*     istopb points to the last column we will be updating
      istopb = ihi

      do i = 1,ns,2
*        Chase the shift down to the bottom right corner
         do ishift = ihi-i-1,ihi-2
            call dlaqz2(.true.,.true.,ishift,istartb,istopb,ihi,A,ldA,B,
     $         ldB,ns,ihi-ns+1,Qc,ldQc,ns+1,ihi-ns,Zc,ldZc)
         end do
         
      end do

*     Update rest of the pencil

*     Update A(ihi-ns+1:ihi, ihi+1:istopm)
*     from the left with Qc(1:ns,1:ns)'
      sheight = ns
      swidth = istopm-(ihi+1)+1
      if (swidth > 0) then
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,
     $      A(ihi-ns+1,ihi+1),ldA,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,A(ihi-ns+1,
     $      ihi+1),ldA)
         call dgemm('T','N',sheight,swidth,sheight,one,Qc,ldQc,
     $      B(ihi-ns+1,ihi+1),ldB,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,B(ihi-ns+1,
     $      ihi+1),ldB)
      end if
      if (ilq) then
         call dgemm('N','N',n,ns,ns,one,Q(1,ihi-ns+1),ldQ,Qc,ldQc,zero,
     $      work,n)
         call dlacpy('ALL',n,ns,work,n,Q(1,ihi-ns+1),ldQ)
      end if

*     Update A(istartm:ihi-ns,ihi-ns:ihi)
*     from the right with Zc(1:ns+1,1:ns+1)
      sheight = ihi-ns-istartm+1
      swidth = ns+1
      if (sheight > 0) then
         call dgemm('N','N',sheight,swidth,swidth,one,A(istartm,ihi-ns),
     $      ldA,Zc,ldZc,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,A(istartm,
     $      ihi-ns),ldA)
         call dgemm('N','N',sheight,swidth,swidth,one,B(istartm,ihi-ns),
     $      ldB,Zc,ldZc,zero,work,sheight)
         call dlacpy('ALL',sheight,swidth,work,sheight,B(istartm,
     $      ihi-ns),ldB)
      end if
      if (ilz) then
      call dgemm('N','N',n,ns+1,ns+1,one,Z(1,ihi-ns),ldZ,Zc,ldZc,zero,
     $   work,n)
         call dlacpy('ALL',n,ns+1,work,n,Z(1,ihi-ns),ldZ)
      end if

      end subroutine
