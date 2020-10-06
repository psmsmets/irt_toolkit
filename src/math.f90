!*************************************************************************************************************************
!
!                                                       M A T H
!
!  MODULE:       math.f90
!
!  Programmer:   Pieter S. M. Smets
!                Seismology Devision - Koninklijk Nederlands Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         August 29, 2014
!
!  Language:     Fortran-90
!
!  Description:  . 
!
!*************************************************************************************************************************

Module math

   use syms
   use io_types

!
!..define dms types
!
   type type_dms
      integer(int16) :: deg
      integer(int8)  :: min
      real(single)   :: sec
   end type type_dms
!

Contains

! =========================================================================================
!
!..Subroutine sph2llh
!
!  convert ph (azimuth), th (colatitude), and r (radius) to lon (longitude), lat
!  (latitude), and h (altitude) above the ellipsoid WGS84
!
!****************************
   Subroutine sph2llh ( ph, th, r, lat, lon, h, lonshift )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in  )  :: ph, th, r  
      double precision, intent( out )  :: lat, lon, h
      double precision, intent(in), optional  :: lonshift
!
! ----
!          
!.....Spherical convertion: radial distance r, polar angle θ (theta), and azimuthal angle φ (phi)
!
      if (present(lonshift)) then
        lon = lonrangef(ph * r2d) - lonshift ! longitude (deg)
      else
        lon = lonrangef(ph * r2d) ! longitude (deg)
      end if
      lat = 90.d0 - th * r2d    ! latitude (deg)
      h   = r - RE              ! altitude (km)
!
      return
!
   End subroutine sph2llh
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine llh2sph
!
!  convert lon (longitude), lat (latitude), and h (altitude) to ph (azimuth), 
!  th (colatitude), and r (radius) to 
!
!****************************
   Subroutine llh2sph ( lat, lon, h, ph, th, r, lonshift )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in  )  :: lat, lon, h
      double precision, intent( out )  :: ph, th, r
      double precision, intent(in), optional  :: lonshift
!
! ----
!          
!.....Spherical convertion: radial distance r, polar angle θ (theta), and azimuthal angle φ (phi)
!
      if (present(lonshift)) then
        ph = lonrangef(lon+lonshift) * d2r   ! longitude (deg)
      else
        ph = lonrangef(lon) * d2r   ! longitude (deg)
      end if
      th = ( 90.d0 - lat ) * d2r  ! latitude (deg)
      r   = re + h                ! altitude (km)
!
      return
!
   End subroutine llh2sph 
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine distance
!
!  get distance between to coordinates
!
!****************************
   Subroutine distance ( lat1, lon1, lat2, lon2, d, earthRadius )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )     :: lat1, lon1, lat2, lon2
      double precision, intent( out )    :: d
      integer, intent(in), optional      :: earthRadius
!
!.....local variables
!
      double precision                   :: latA, lonA, latB, lonB
!
! ----
!        
!.....Degrees to radians transform
!
      latA = lat1 * d2r
      lonA = lon1 * d2r
      latB = lat2 * d2r
      lonB = lon2 * d2r
!
!.....Great circle distance, less subjective to rounding errors for small distances
!
      d = 2 * dasin(   dsqrt(  ( dsin(  ( latA - latB ) / 2  ) )**2 + & 
         & dcos( latA ) * dcos( latB ) * ( dsin(  ( lonA - lonB ) / 2  ) )**2  )   )
!
!.....Great circle distance to km
!
      if (present(earthRadius)) then
          d = d * earthRadius
      else
          d = d * RE
      end if
!
      return
!
   End subroutine distance
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine euclidean_distance 
!
!  get euclidean distance between two cartesian coordinates
!
!****************************
   Subroutine euclidean_distance ( x1, y1, z1, x2, y2, z2, d )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )   :: x1, y1, z1, x2, y2, z2 
      double precision, intent( out )  :: d
!
! ----
!
      d = dsqrt(  ( x1 - x2 ) ** 2 + ( y1 - y2 ) ** 2 + ( z1 - z2 ) ** 2  )
!
      return
!
   End subroutine euclidean_distance 
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine euclidean_distance_llh 
!
!  get euclidean distance between two geodetic coordinates
!
!****************************
   Subroutine euclidean_distance_llh ( lat1, lon1, h1, lat2, lon2, h2, d )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )   :: lat1, lon1, h1, lat2, lon2, h2 
      double precision, intent( out )  :: d
!          
!.....local variables
!
      double precision   :: x1, y1, z1, x2, y2, z2
!
! ----
!
      call llh2cart( lat1, lon1, h1, x1, y1, z1 )
      call llh2cart( lat2, lon2, h2, x2, y2, z2 )
!
      call euclidean_distance( x1, y1, z1, x2, y2, z2, d )
!
      return
!
   End subroutine euclidean_distance_llh 
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine euclidean_distance_sph 
!
!  get euclidean distance between two spherical coordinates
!
!****************************
   Subroutine euclidean_distance_sph ( ph1, th1, r1, ph2, th2, r2, d )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )   :: ph1, th1, r1, ph2, th2, r2 
      double precision, intent( out )  :: d
!          
!.....local variables
!
      double precision   :: x1, y1, z1, x2, y2, z2
!
! ----
!
!.....Spherical convertion: radial distance r, polar angle θ (theta), and azimuthal angle φ (phi)
!
      call sph2cart( ph1, th1, r1, x1, y1, z1 )
      call sph2cart( ph2, th2, r2, x2, y2, z2 )
!
      call euclidean_distance( x1, y1, z1, x2, y2, z2, d )
!
      return
!
   End subroutine euclidean_distance_sph 
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine sph2cart
!
!  Convert shperical to cartesian coordinates
!
!****************************
   Subroutine sph2cart ( ph, th, r, x, y, z )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )   :: ph, th, r
      double precision, intent( out )  :: x, y, z
!
! ----
!
!.....Spherical convertion: radial distance r, polar angle θ (theta), and azimuthal angle φ (phi)
!
      x = r * dsin( th ) * dcos( ph )
      y = r * dsin( th ) * dsin( ph )
      z = r * dcos( th )
!
      return
!
   End subroutine sph2cart 
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine cart2sph
!
!  Convert cartesian to spherical coordinates 
!
!****************************
   Subroutine cart2sph ( x, y, z, ph, th, r )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )   :: x, y, z
      double precision, intent( out )  :: ph, th, r
!
! ----
!
!.....Spherical convertion: radial distance r, polar angle θ (theta), and azimuthal angle φ (phi)
!
      r  = dsqrt( x ** 2 + y ** 2 + z ** 2 )
      th = dacos( z / r )
      ph = datan2( y, x )
!
      return
!
   End subroutine cart2sph 
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine llh2cart 
!
!  Convert lat lon to cartesian coordinates
!
!****************************
   Subroutine llh2cart ( lat, lon, h, x, y, z )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )   :: lat, lon, h
      double precision, intent( out )  :: x, y, z
!          
!.....local variables
!
      double precision  :: ph, th, r
!
! ----
!
      call llh2sph( lat, lon, h, ph, th, r )
      call sph2cart( ph, th, r, x, y, z )
!
      return
!
   End subroutine llh2cart 
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine cart2llh
!
!  Convert cartesian to lat lon coordinates 
!
!****************************
   Subroutine cart2llh ( x, y, z, lat, lon, h )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )   :: x, y, z
      double precision, intent( out )  :: lat, lon, h
!          
!.....local variables
!
      double precision  :: ph, th, r
!
! ----
!
      call cart2sph( x, y, z, ph, th, r )
      call sph2llh( ph, th, r, lat, lon, h )
!
      return
!
   End subroutine cart2llh 
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine angle
!
!  get angle between two coordinates
!
!****************************
   Subroutine angle ( lat1, lon1, lat2, lon2, a )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )     :: lat1, lon1, lat2, lon2
      double precision, intent( out )    :: a
!
!.....local variables
!
      double precision                 :: latA, lonA, latB, lonB
!
! ----
!        
!.....Degrees to radians transform
!
      latA = lat1 * d2r
      lonA = lon1 * d2r
      latB = lat2 * d2r
      lonB = lon2 * d2r
!
!.....Bearing angle
!
      a = datan2(  &
         &   dsin( lonB - lonA ) * dcos( latB ) , &
         &   dcos( latA ) * dsin( latB ) - dsin( latA ) * dcos( latB ) * dcos( lonB - lonA ) &
         & )
!
!.....Course to degrees
!
      a = a * r2d
!
      return
!
   End subroutine angle
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine course
!
!  get the great circle distance and bearing angle between two coordinates
!
!****************************
   Subroutine course ( lat1, lon1, lat2, lon2, d, a )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )     :: lat1, lon1, lat2, lon2
      double precision, intent( out )    :: d, a
!
!.....local variables
!
      double precision                 :: latA, lonA, latB, lonB
!
! ----
!        
!.....Degrees to radians transform
!
      latA = lat1 * d2r
      lonA = lon1 * d2r
      latB = lat2 * d2r
      lonB = lon2 * d2r
!
!.....Great circle distance, less subjective to rounding errors for small distances
!
      d = 2 * dasin( dsqrt(   &
         &   (  dsin( ( latA - latB ) / 2 )  )**2 + & 
         &   dcos( latA ) * dcos( latB ) * (  dsin( ( lonA - lonB ) / 2 )  )**2 &
         & ) )
!
!.....Bearing angle
!
      a = datan2(  &
         &   dsin( lonB - lonA ) * dcos( latB ) , &
         &   dcos( latA ) * dsin( latB ) - dsin( latA ) * dcos( latB ) * dcos( lonB - lonA )  &
         & )
!
!.....Great circle distance to km
!
      d = d * RE
!
!.....Course to degrees
!
      a = a * r2d
!
      return
!
   End subroutine course
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine radial_distance 
!
!  get the end coordinate given a start coordinate and a distance and bearing
!
!****************************
   Subroutine radial_distance ( lat1, lon1, d, a, lat2, lon2, a2 )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )            :: lat1, lon1, a, d
      double precision, intent( out )           :: lat2, lon2
      double precision, intent( out ), optional :: a2
!
!.....local variables
!
      double precision                 :: az, di, lat, lon, az2
      double precision, parameter      :: EPS = 1.0d0 / 10**6
!
! ----
!
!.....Degrees to radians transform
!
      lat = lat1 * d2r
      lon = lon1 * d2r
      az  = a * d2r
      di  = d * d2r
!
!.....A point (lat1,lon1) is a distance d out on the tc radial from point (lat0,lon0)
!
      lat2 = dasin(  dsin( lat ) * dcos( di ) + dcos( lat ) * dsin( di ) * dcos( az )  )
!
      if ( dcos( lat2 ) .LE. EPS ) then
         lon2 = lon ! endpoint is a pole
      else
         lon2 = lon + datan2(  &
            &   dsin( az ) * dsin( di ) * dcos( lat ), &
            &   dcos( di ) - dsin( lat ) * dsin( lat2 ) &
            & )
      end if
!
!.....Radians to degrees
!
      lat2 = lat2 * r2d
      lon2 = lon2 * r2d
!
!.....bearing angle of final point
!
      if ( present( a2 ) ) then
         CALL ANGLE ( lat2, lon2, lat1, lon1, az2 )
         a2 = lonrangeF ( az2 + 180.d0 )
      end if
!
      return
!
   End subroutine radial_distance
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine sphdistance
!
!  get distance between to coordinates
!
!****************************
   Subroutine sphdistance ( ph1, th1, ph2, th2, d )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )     :: ph1, th1, ph2, th2
      double precision, intent( out )    :: d
!
! ----
!
!.....Great circle distance, less subjective to rounding errors for small distances
!
      d = 2 * dasin(   dsqrt(  ( dsin(  ( th1 - th2 ) / 2  ) )**2 + & 
         & dsin( th1 ) * dsin( th2 ) * ( dsin(  ( ph1 - ph2 ) / 2  ) )**2  )   )
!
!.....Great circle distance to km
!
      d = d * RE
!
      return
!
   End subroutine sphdistance
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine twopirange
!
!  check and correct longitude to be within [0,twopi) rad interval.
!
!****************************
   Subroutine twopirange ( ph, offset )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( inout )        :: ph
      double precision, intent( in ), optional :: offset
!
!.....local variables
!
      double precision :: o = 0.d0
!
! ----
!
!.....Set offset
!
      if ( present( offset ) ) o = offset
!
!.....Degrees to radians transform
!
      ph = ph - floor(  ( ph - o ) / twopi  ) * twopi
!
      return
!
   End subroutine twopirange
!
! =========================================================================================


! =========================================================================================
!
!..Function twopirangef
!
!  check and correct longitude to be within [0,twopi) rad interval.
!
!****************************
   Function twopirangef ( ph, offset ) RESULT ( ph2 )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )           :: ph
      double precision                         :: ph2
      double precision, intent( in ), optional :: offset
!
!.....local variables
!
      double precision :: o = 0.d0
!
! ----
!
!.....Set offset
!
      if ( present( offset ) ) o = offset
!
!.....Degrees to radians transform
!
      ph2 = ph - floor(  ( ph - o ) / twopi  ) * twopi
!
      return
!
   End function twopirangef
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine npole_range
!
!  check and correct latitude at north pole.
!
!****************************
   Subroutine npole_range ( ph, th, lpole )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      logical, intent( out )            :: lpole
      double precision, intent( inout ) :: ph, th
!
!.....local variables
!
      double precision :: ph2
!
! ----
!
!.....Correct th
!
      if ( th .LT. 0.d0 ) then
         th  = dabs( th )
         ph2 = ph + pi
         ph = twopirangef( ph2 )
         lpole = .true.
      else
         lpole = .false.
      end if
!
      return
!
   End subroutine npole_range
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine spole_range
!
!  check and correct latitude at south pole.
!
!****************************
   Subroutine spole_range ( ph, th, lpole )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      logical, intent( out )            :: lpole
      double precision, intent( inout ) :: ph, th
!
!.....local variables
!
      double precision :: ph2
!
! ----
!
!.....Correct th
!
      if ( th .gt. pi ) then
         th = twopi - th
         ph2 = ph + pi
         ph = twopirangef( ph2 )
         lpole = .true.
      else
         lpole = .false.
      end if
!
      return
!
   End subroutine spole_range
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine lonrange
!
!  check and correct longitude to be within [0,360) degrees interval.
!
!****************************
   Subroutine lonrange ( lon, offset )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( inout )        :: lon
      double precision, intent( in ), optional :: offset
!
!.....local variables
!
      double precision :: o
!
! ----
!
!.....Set offset
!
      if ( present(offset) ) then
        o = offset
      else
        o = 0.d0
      end if
!
!.....Degrees to radians transform
!
      lon = lon - floor(  ( lon - o ) / 360.d0  ) * 360.d0
!
      return
!
   End subroutine lonrange
!
! =========================================================================================


! =========================================================================================
!
!..Function lonrangeF
!
!  check and correct longitude to be within [0,360) degrees interval.
!
!****************************
   Function lonrangeF ( lon, offset ) RESULT ( lon2 )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in )           :: lon
      double precision                         :: lon2
      double precision, intent( in ), optional :: offset
!
!.....local variables
!
      double precision :: o
!
! ----
!
      lon2=lon
      call lonrange( lon2, offset )
!
      return
!
   End function lonrangeF
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine latrange 
!
!  check and correct latitdue to be within [-90,90] degrees interval.
!
!****************************
   Subroutine latrange ( lat, lon )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( inout ) :: lat, lon
!
! ----
!
!.....Degrees to radians transform
!
      if (  dabs( lat ) .gt. 90.d0 ) then
         lat = sign( 1.d0, lat ) * 180.d0 - lat
         lon = lonrangeF( lon + 180.d0 )
      end if
!
      return
!
   End subroutine latrange
!
! =========================================================================================



! =========================================================================================
!
!..Function degrees2dms
!
!  Convert decimal degrees to degrees-minutus-seconds
!
!****************************
   Subroutine degrees2dms ( deg, dms, round_sec )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      real(double)  , intent(in)             :: deg
      type(type_dms), intent(out)            :: dms
      logical       , intent(in) , optional  :: round_sec
!
!.....Local variables
!
      logical       :: rnd
      real(single)  :: dummy
!
! ----
!
      rnd=.false.
      if (present(round_sec)) rnd=round_sec
!
      if (rnd) then
         dummy=real(anint(deg*3600._double,kind=double)/3600._double,kind=single)
      else
         dummy=real(deg,kind=single)
      end if
!
      dms%deg=int(dummy,kind=int16)
      dummy=(dummy-dms%deg)*60
      dms%min=int(dummy,kind=int8)
      dms%sec=real((dummy-dms%min)*60,kind=single)
!
      return
!
   End subroutine degrees2dms
!
! =========================================================================================


! =========================================================================================
!
!..Function deg2dms
!
!  Convert decimal degrees to degrees-minutus-seconds
!
!****************************
   Function deg2dms ( deg ) result( dms )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      type(type_dms)                        :: dms
      real(double)  , intent(in)            :: deg
!
! ----
!
      call degrees2dms( deg, dms, .false. )
!
      return
!
   End function deg2dms
!
! =========================================================================================


! =========================================================================================
!
!..Function deg2dmsr
!
!  Convert decimal degrees to degrees-minutus-seconds
!
!****************************
   Function deg2dmsr ( deg ) result( dms )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      type(type_dms)                        :: dms
      real(double)  , intent(in)            :: deg
!
! ----
!
      call degrees2dms( deg, dms, .true. )
!
      return
!
   End function deg2dmsr
!
! =========================================================================================


! =========================================================================================
!
!..Function dms2deg
!
!  Convert degrees-minutus-seconds to decimal degrees
!
!****************************
   Function dms2deg( dms ) result( deg )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      real(double)                :: deg
      type(type_dms), intent(in)  :: dms
!
! ----
!
      deg = real(dms%deg,kind=double)+real(dms%min,double)/60._double+dms%sec/3600._double
!
      return
!
   End function dms2deg
!
! =========================================================================================


! =========================================================================================
!
!..Function dms2str
!
!   Print degrees-minutus-seconds to string
!
!****************************
   Function dms2str( dms ) result( str )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      character(len=20)           :: str
      type(type_dms), intent(in)  :: dms
!
!.....Local variables
!
      character  :: d*4, m*2, s*6
!
! ----
!
      write(d,"(i4)") dms%deg
      write(m,"(i2)") dms%min
      write(s,"(f6.3)") dms%sec
!
      str=d//'d '//m//'m '//s//'s'
!
      return
!
   End function dms2str
!
! =========================================================================================

! =========================================================================================
!
!..Function sec 
!
!  Secans
!
!****************************
   double precision Function sec ( a )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in ) :: a
!
! ----
!
      sec = 1.d0 / dsin( a )
!
      return
!
   End function sec
!
! =========================================================================================


! =========================================================================================
!
!..Function csc 
!
!  Cosecans
!
!****************************
   double precision Function csc ( a )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in ) :: a
!
! ----
!
      csc = 1.d0 / dcos( a )
!
      return
!
   End function csc
!
! =========================================================================================


! =========================================================================================
!
!..Function cot 
!
!  Cotangens
!
!****************************
   double precision Function cot ( a )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, intent( in ) :: a
!
! ----
!
      cot = 1.d0 / dtan( a )
!
      return
!
   End function cot
!
! =========================================================================================


! =========================================================================================
!
!..Subroutine smooth 
!
!****************************
   Subroutine smooth ( y, span )
!****************************
!
      implicit none
!          
!.....dummy variables
!
      double precision, dimension ( : ), intent( inout )            :: y
      integer                          , intent( in    ), optional  :: span
!
!.....local variables
!
      integer  :: nofy, i, s, nofb
      double precision, dimension ( size( y, 1 ) ) :: yy
!
! ----
!
!.....set span
!
      if ( present( span ) ) then
         if ( mod( span, 2 ) .EQ. 0 ) then
            s = span + 1
         else
            s = 5
         end if
      else
         s = 5
      end if
!
!.....get length
!
      nofy = size( y, 1 )
!
      if ( nofy .LE. span .or. span .LT. 3 ) return
!
!.....boundary points
!
      nofb = ( s - 1 ) / 2
!
      do i = 0, nofb - 1
        yy(    1 + i ) = sum(  y( 1 : 2 * i + 1 )  ) / ( 2 * i + 1 ) 
        yy( nofy - i ) = sum(  y( nofy - 2 * i : nofy )  ) / ( 2 * i + 1 )
      end do
!
!.....regular values
!
      do i = nofb + 1, nofy - nofb
         yy( i ) = sum( y ( i - 2 : i + 2 ) ) / 5
      end do
!
!.....copy
!
      y = yy
!
      return
!
   End subroutine smooth 
!
! =========================================================================================

End module math
