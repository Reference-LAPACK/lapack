*> \brief \b SLAQZ1
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SLAQZ1 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqz1.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqz1.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqz1.f">
*> [TXT]</a>
*> \endhtmlonly
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>      Given a 3-by-3 matrix pencil (A,B), SLAQZ1 sets v to a
*>      scalar multiple of the first column of the product
*>
*>      (*)  K = (A - (beta2*sr2 - i*si)*B)*B^(-1)*(beta1*A - (sr2 + i*si2)*B)*B^(-1).
*>
*>      It is assumed that either
*>
*>              1) sr1 = sr2
*>          or
*>              2) si = 0.
*>
*>      This is useful for starting double implicit shift bulges
*>      in the QZ algorithm.
*> \endverbatim
*
*
*  Arguments:
*  ==========
*
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>              The 3-by-3 matrix A in (*).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>              The leading dimension of A as declared in
*>              the calling procedure.
*> \endverbatim
*
*> \param[in] B
*> \verbatim
*>          B is REAL array, dimension (LDB,N)
*>              The 3-by-3 matrix B in (*).
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>              The leading dimension of B as declared in
*>              the calling procedure.
*> \endverbatim
*>
*> \param[in] SR1
*> \verbatim
*>          SR1 is REAL
*> \endverbatim
*>
*> \param[in] SR2
*> \verbatim
*>          SR2 is REAL
*> \endverbatim
*>
*> \param[in] SI
*> \verbatim
*>          SI is REAL
*> \endverbatim
*>
*> \param[in] BETA1
*> \verbatim
*>          BETA1 is REAL
*> \endverbatim
*>
*> \param[in] BETA2
*> \verbatim
*>          BETA2 is REAL
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is REAL array, dimension (N)
*>              A scalar multiple of the first column of the
*>              matrix K in (*).
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
      subroutine slaqz1(A,ldA,B,ldB,sr1,sr2,si,beta1,beta2,v)
      implicit none
*
*     Arguments
      integer,intent(in) :: ldA,ldB
      real,intent(in) :: A(ldA,*),B(ldB,*)
      real,intent(out) :: sr1,sr2,si,beta1,beta2,v(*)
*
*     Parameters
      real :: zero,one,half
      parameter(zero=0.0,one=1.0,half=0.5)
*
*     Local scalars
      real :: w(2),safmin,safmax,scale
*
*     External Functions
      real,external :: slamch
*
      safmin = slamch('SAFE MINIMUM')
      safmax = one/safmin
*
*     Calculate first shifted vector
*
      w(1) = beta1*A(1,1)-sr1*B(1,1)
      w(2) = beta1*A(2,1)-sr1*B(2,1)
      scale = sqrt( abs(w(1)) ) * sqrt( abs(w(2)) )
      if(scale .ge. safmin .and. scale .le. safmax) then
         w(1) = w(1)/scale
         w(2) = w(2)/scale
      end if
*
*     Solve linear system
*
      w(2) = w(2)/B(2,2)
      w(1) = (w(1)-B(1,2)*w(2))/B(1,1)
      scale = sqrt( abs(w(1)) ) * sqrt( abs(w(2)) )
      if(scale .ge. safmin .and. scale .le. safmax) then
         w(1) = w(1)/scale
         w(2) = w(2)/scale
      end if
*
*     Apply second shift
*
      v(1) = beta2*(A(1,1)*w(1)+A(1,2)*w(2))-sr2*(B(1,1)*w(1)+B(1,
     $   2)*w(2))
      v(2) = beta2*(A(2,1)*w(1)+A(2,2)*w(2))-sr2*(B(2,1)*w(1)+B(2,
     $   2)*w(2))
      v(3) = beta2*(A(3,1)*w(1)+A(3,2)*w(2))-sr2*(B(3,1)*w(1)+B(3,
     $   2)*w(2))
*
*     Account for imaginary part
*
      v(1) = v(1)+si*si*B(1,1)
*
*     Check for overflow
*
      if( abs(v(1)).gt.safmax .or. abs(v(2)) .gt. safmax .or. abs(v(3)).
     $   gt.safmax .or. v(1).ne.v(1) .or. v(2).ne.v(2) .or. v(3).ne.v(3)
     $   ) then
         v(1) = zero
         v(2) = zero
         v(3) = zero
      end if
*
*     End of SLAQZ1
*
      end subroutine
