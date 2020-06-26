*> \brief \b SLAQZ2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SLAQZ2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqz2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqz2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqz2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>      SLAQZ2 chases a 2x2 shift bulge in a matrix pencil down a single position
*> \endverbatim
*
*
*  Arguments:
*  ==========
*
*>
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
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>              Index indicating the position of the bulge.
*>              On entry, the bulge is located in
*>              (A(k+1:k+2,k:k+1),B(k+1:k+2,k:k+1)).
*>              On exit, the bulge is located in
*>              (A(k+2:k+3,k+1:k+2),B(k+2:k+3,k+1:k+2)).
*> \endverbatim
*>
*> \param[in] ISTARTM
*> \verbatim
*>          ISTARTM is INTEGER
*> \endverbatim
*>
*> \param[in] ISTOPM
*> \verbatim
*>          ISTOPM is INTEGER
*>              Updates to (A,B) are restricted to
*>              (istartm:k+3,k:istopm). It is assumed
*>              without checking that istartm <= k+1 and
*>              k+2 <= istopm
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*> \endverbatim
*>
*> \param[inout] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>              The leading dimension of A as declared in
*>              the calling procedure.
*> \endverbatim
*
*> \param[inout] B
*> \verbatim
*>          B is REAL array, dimension (LDB,N)
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>              The leading dimension of B as declared in
*>              the calling procedure.
*> \endverbatim
*>
*> \param[in] NQ
*> \verbatim
*>          NQ is INTEGER
*>              The order of the matrix Q
*>
*> \param[in] QSTART
*> \verbatim
*>          QSTART is INTEGER
*>              Start index of the matrix Q. Rotations are applied
*>              To columns k+2-qStart:k+4-qStart of Q.
*> \endverbatim
*
*> \param[inout] Q
*> \verbatim
*>          Q is REAL array, dimension (LDQ,NQ)
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>              The leading dimension of Q as declared in
*>              the calling procedure.
*> \endverbatim
*>
*> \param[in] NZ
*> \verbatim
*>          NZ is INTEGER
*>              The order of the matrix Z
*>
*> \param[in] ZSTART
*> \verbatim
*>          ZSTART is INTEGER
*>              Start index of the matrix Z. Rotations are applied
*>              To columns k+1-qStart:k+3-qStart of Z.
*> \endverbatim
*
*> \param[inout] Z
*> \verbatim
*>          Z is REAL array, dimension (LDZ,NZ)
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>              The leading dimension of Q as declared in
*>              the calling procedure.
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
      subroutine slaqz2(ilq,ilz,k,istartm,istopm,ihi,A,ldA,B,ldB,nQ,
     $   qStart,Q,ldQ,nZ,zStart,Z,ldZ)
      implicit none
*
*     Arguments
      logical,intent(in) :: ilq,ilz
      integer,intent(in) :: k,ldA,ldB,ldQ,ldZ,istartm,istopm,nQ,nZ,
     $   qStart,zStart,ihi
      real :: A(ldA,*),B(ldB,*),Q(ldQ,*),Z(ldZ,*)
*
*     Parameters
      real :: zero,one,half
      parameter(zero=0.0,one=1.0,half=0.5)
*
*     Local variables
      real :: H(2,3),c1,s1,c2,s2,temp
*
      if(k+2 .eq. ihi) then
*        Shift is located on the edge of the matrix, remove it
         H = B(ihi-1:ihi,ihi-2:ihi)
*        Make H upper triangular
         call slartg(H(1,1),H(2,1),c1,s1,temp)
         H(2,1) = zero
         H(1,1) = temp
         call srot(2,H(1,2),2,H(2,2),2,c1,s1)
*
         call slartg(H(2,3),H(2,2),c1,s1,temp)
         call srot(1,H(1,3),1,H(1,2),1,c1,s1)
         call slartg(H(1,2),H(1,1),c2,s2,temp)
*
         call srot(ihi-istartm+1,B(istartm,ihi),1,B(istartm,ihi-1),1,c1,
     $      s1)
         call srot(ihi-istartm+1,B(istartm,ihi-1),1,B(istartm,ihi-2),1,
     $      c2,s2)
         B(ihi-1,ihi-2) = zero
         B(ihi,ihi-2) = zero
         call srot(ihi-istartm+1,A(istartm,ihi),1,A(istartm,ihi-1),1,c1,
     $      s1)
         call srot(ihi-istartm+1,A(istartm,ihi-1),1,A(istartm,ihi-2),1,
     $      c2,s2)
         if (ilz) then
            call srot(nz,Z(1,ihi-zStart+1),1,Z(1,ihi-1-zStart+1),1,c1,
     $         s1)
            call srot(nz,Z(1,ihi-1-zStart+1),1,Z(1,ihi-2-zStart+1),1,c2,
     $         s2)
         end if
*
         call slartg(A(ihi-1,ihi-2),A(ihi,ihi-2),c1,s1,temp)
         A(ihi-1,ihi-2) = temp
         A(ihi,ihi-2) = zero
         call srot(istopm-ihi+2,A(ihi-1,ihi-1),ldA,A(ihi,ihi-1),ldA,c1,
     $      s1)
         call srot(istopm-ihi+2,B(ihi-1,ihi-1),ldB,B(ihi,ihi-1),ldB,c1,
     $      s1)
         if (ilq) then
            call srot(nq,Q(1,ihi-1-qStart+1),1,Q(1,ihi-qStart+1),1,c1,
     $         s1)
         end if
*
         call slartg(B(ihi,ihi),B(ihi,ihi-1),c1,s1,temp)
         B(ihi,ihi) = temp
         B(ihi,ihi-1) = zero
         call srot(ihi-istartm,B(istartm,ihi),1,B(istartm,ihi-1),1,c1,
     $      s1)
         call srot(ihi-istartm+1,A(istartm,ihi),1,A(istartm,ihi-1),1,c1,
     $      s1)
         if (ilz) then
            call srot(nz,Z(1,ihi-zStart+1),1,Z(1,ihi-1-zStart+1),1,c1,
     $         s1)
         end if
*
      else
*
*        Normal operation, move bulge down
*
         H = B(k+1:k+2,k:k+2)
*
*        Make H upper triangular
*
         call slartg(H(1,1),H(2,1),c1,s1,temp)
         H(2,1) = zero
         H(1,1) = temp
         call srot(2,H(1,2),2,H(2,2),2,c1,s1)
*
*        Calculate Z1 and Z2
*
         call slartg(H(2,3),H(2,2),c1,s1,temp)
         call srot(1,H(1,3),1,H(1,2),1,c1,s1)
         call slartg(H(1,2),H(1,1),c2,s2,temp)
*
*        Apply transformations from the right
*
         call srot(k+3-istartm+1,A(istartm,k+2),1,A(istartm,k+1),1,c1,
     $      s1)
         call srot(k+3-istartm+1,A(istartm,k+1),1,A(istartm,k),1,c2,s2)
         call srot(k+2-istartm+1,B(istartm,k+2),1,B(istartm,k+1),1,c1,
     $      s1)
         call srot(k+2-istartm+1,B(istartm,k+1),1,B(istartm,k),1,c2,s2)
         if (ilz) then
            call srot(nZ,Z(1,k+2-zStart+1),1,Z(1,k+1-zStart+1),1,c1,s1)
            call srot(nZ,Z(1,k+1-zStart+1),1,Z(1,k-zStart+1),1,c2,s2)
         end if
         B(k+1,k) = zero
         B(k+2,k) = zero
*
*        Calculate Q1 and Q2
*
         call slartg(A(k+2,k),A(k+3,k),c1,s1,temp)
         A(k+2,k) = temp
         A(k+3,k) = zero
         call slartg(A(k+1,k),A(k+2,k),c2,s2,temp)
         A(k+1,k) = temp
         A(k+2,k) = zero
*
*     Apply transformations from the left
*
         call srot(istopm-k,A(k+2,k+1),ldA,A(k+3,k+1),ldA,c1,s1)
         call srot(istopm-k,A(k+1,k+1),ldA,A(k+2,k+1),ldA,c2,s2)
*
         call srot(istopm-k,B(k+2,k+1),ldB,B(k+3,k+1),ldB,c1,s1)
         call srot(istopm-k,B(k+1,k+1),ldB,B(k+2,k+1),ldB,c2,s2)
         if (ilq) then
            call srot(nQ,Q(1,k+2-qStart+1),1,Q(1,k+3-qStart+1),1,c1,s1)
            call srot(nQ,Q(1,k+1-qStart+1),1,Q(1,k+2-qStart+1),1,c2,s2)
         end if
*
      end if
*
*     End of SLAQZ2
*
      end subroutine