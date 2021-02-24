module LA_CONSTANTS
!
!  -- BLAS/LAPACK module --
!     May 06, 2016
!
!  Standard constants
!
   double precision, parameter :: zero = 0.0
   double precision, parameter :: half = 0.5
   double precision, parameter :: one = 1.0
   double precision, parameter :: two = 2.0
   double precision, parameter :: three = 3.0
   double precision, parameter :: four = 4.0
   double precision, parameter :: eight = 8.0
   double precision, parameter :: ten = 10.0
   complex*16, parameter :: czero = ( 0.0, 0.0 )
   complex*16, parameter :: chalf = ( 0.5, 0.0 )
   complex*16, parameter :: cone = ( 1.0, 0.0 )
   character*1, parameter :: sprefix = 'D'
   character*1, parameter :: cprefix = 'Z'
!
!  Model parameters
!
   double precision, parameter :: eps =  0.11102230246251565404D-015
   double precision, parameter :: ulp =  0.22204460492503130808D-015
   double precision, parameter :: safmin =  0.22250738585072013831D-307
   double precision, parameter :: safmax =  0.44942328371557897693D+308
   double precision, parameter :: smlnum =  0.10020841800044863890D-291
   double precision, parameter :: bignum =  0.99792015476735990583D+292
   double precision, parameter :: rtmin =  0.10010415475915504622D-145
   double precision, parameter :: rtmax =  0.99895953610111751404D+146
!
!  Blue's scaling constants
!
   double precision, parameter :: tsml =  0.14916681462400413487D-153
   double precision, parameter :: tbig =  0.19979190722022350281D+147
   double precision, parameter :: ssml =  0.44989137945431963828D+162
   double precision, parameter :: sbig =  0.11113793747425387417D-161
end module LA_CONSTANTS
