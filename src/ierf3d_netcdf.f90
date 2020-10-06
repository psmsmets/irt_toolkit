!*************************************************************************************************************************
!
!                                              I E R F 3 D _ N E T C D F
!
!  Module:       IERF3D_NETCDF 
!
!  Programmer:   Pieter S. M. Smets
!                Seismology Devision - Koninklijk Nederlands Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         October 21, 2014
!
!  Language:     Fortran-90
!
!  Description:  This module contains writing modules for writing netcdf files for the eigen ray
!                finder IERF3D.
!
!                   isnetcd         check if file is of type netcdf
!                   write_ierf3d    write IERF3D results to netcdf file
!                   check           netcdf error check
!
!*************************************************************************************************************************

Module ierf3d_netcdf

   use netcdf

!
!..variable unit description
!
   character (len = *), parameter, private :: units       = "units"
   character (len = *), parameter, private :: att_missing = "missing_value"
   real(kind=4)       , parameter, private :: val_missing = -9999.
   character (len = *), parameter, private :: time_units = "-"
   character (len = *), parameter, private :: lat__units = "deg N"
   character (len = *), parameter, private :: lon__units = "deg E"
   character (len = *), parameter, private :: alt__units = "km"
   character (len = *), parameter, private :: azi__units = "deg"
   character (len = *), parameter, private :: elev_units = "deg"
   character (len = *), parameter, private :: capp_units = "m/s"
   character (len = *), parameter, private :: incl_units = "deg"
   character (len = *), parameter, private :: d____units = "km"
   character (len = *), parameter, private :: loss_units = "dB RE 1"
   character (len = *), parameter, private :: refl_units = "-"
   character (len = *), parameter, private :: tt___units = "s"
!
!..variable short name
!
   character (len = *), parameter, private :: time_short_name = "time"
   character (len = *), parameter, private :: azi__short_name = "azi"
   character (len = *), parameter, private :: elev_short_name = "elev"
   character (len = *), parameter, private :: capp_short_name = "capp"
   character (len = *), parameter, private :: bazi_short_name = "bazi"
   character (len = *), parameter, private :: bdev_short_name = "bazidev"
   character (len = *), parameter, private :: apoi_short_name = "azi_poi"
   character (len = *), parameter, private :: bpoi_short_name = "bazi_poi"
   character (len = *), parameter, private :: incl_short_name = "incl"
   character (len = *), parameter, private :: lat__short_name = "lat"
   character (len = *), parameter, private :: lon__short_name = "lon"
   character (len = *), parameter, private :: alt__short_name = "h"
   character (len = *), parameter, private :: lat0_short_name = "lat0"
   character (len = *), parameter, private :: lon0_short_name = "lon0"
   character (len = *), parameter, private :: alt0_short_name = "h0"
   character (len = *), parameter, private :: lat1_short_name = "lat1"
   character (len = *), parameter, private :: lon1_short_name = "lon1"
   character (len = *), parameter, private :: alt1_short_name = "h1"
   character (len = *), parameter, private :: d____short_name = "d"
   character (len = *), parameter, private :: loss_short_name = "tloss"
   character (len = *), parameter, private :: hmin_short_name = "hmin"
   character (len = *), parameter, private :: hmax_short_name = "hmax"
   character (len = *), parameter, private :: hmn__short_name = "hmean"
   character (len = *), parameter, private :: refl_short_name = "nofrefl"
   character (len = *), parameter, private :: tt___short_name = "tt"
!
!..variable long name
!
   character (len = *), parameter, private :: long_name      = "long_name"
   character (len = *), parameter, private :: time_long_name = "Time"
   character (len = *), parameter, private :: azi__long_name = "Azimuth"
   character (len = *), parameter, private :: elev_long_name = "Elevation"
   character (len = *), parameter, private :: capp_long_name = "Trace velocity"
   character (len = *), parameter, private :: bazi_long_name = "Back-azimuth"
   character (len = *), parameter, private :: bdev_long_name = "Back-azimuth deviation"
   character (len = *), parameter, private :: apoi_long_name = "Source-receiver azimuth"
   character (len = *), parameter, private :: bpoi_long_name = "Receiver-source azimuth"
   character (len = *), parameter, private :: incl_long_name = "Inclination"
   character (len = *), parameter, private :: lat__long_name = "Ray latitude"
   character (len = *), parameter, private :: lon__long_name = "Ray longitude"
   character (len = *), parameter, private :: alt__long_name = "Ray altitude"
   character (len = *), parameter, private :: lat0_long_name = "Source latitude"
   character (len = *), parameter, private :: lon0_long_name = "Source longitude"
   character (len = *), parameter, private :: alt0_long_name = "Source altitude"
   character (len = *), parameter, private :: lat1_long_name = "Receiver latitude"
   character (len = *), parameter, private :: lon1_long_name = "Receiver longitude"
   character (len = *), parameter, private :: alt1_long_name = "Receiver altitude"
   character (len = *), parameter, private :: d____long_name = "Euclidean distance to receiver"
   character (len = *), parameter, private :: loss_long_name = "Transmission loss"
   character (len = *), parameter, private :: hmin_long_name = "Minimum ray altitude"
   character (len = *), parameter, private :: hmax_long_name = "Maximum ray altitude"
   character (len = *), parameter, private :: hmn__long_name = "Mean ray altitude"
   character (len = *), parameter, private :: refl_long_name = "Number of reflections"
   character (len = *), parameter, private :: tt___long_name = "Travel time"

Contains

! =========================================================================================
!
!..Function isnetcdf
!
!  Check if file is of netcdf type
!
!****************************
   Logical Function isnetcdf ( file1 )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character ( LEN = * ), intent ( in )  :: file1
!
!.....Local variables
!
      integer(kind=4)  :: ncid, ierr
!
!  ---
!
!.....Open the netCDF file.
!
      ierr = nf90_open ( file1, 0, ncid )
!
!.....If no error, file is netcdf
!
      isnetcdf = ierr .EQ. nf90_noerr
!
!.....If open succesfull, close file again
!
      if ( isnetcdf ) call check ( nf90_close ( ncid ) )
!
      return
!
   End function isnetcdf
!
! =========================================================================================


! =========================================================================================
!
!..subroutine write_ierf3d 
!
!  Write IERF3D results to netcdf grid file.
!
!****************************
   Subroutine write_ierf3d ( file1, source, receiver, azi_poi, bazi_poi, azi, elev, lat, lon, alt, &
      & capp, bazi, bdev, incl, d, tt, tloss, re_unit, hmin, hmax, hmean, nofrefl, tstr  )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character (len=*)                , intent (in)  :: file1
      double precision, dimension (3)  , intent (in)  :: source, receiver
      double precision                 , intent (in)  :: azi_poi, bazi_poi
      double precision, dimension (:)  , intent (in)  :: azi, elev
      double precision, dimension (:,:), intent (in)  :: lat, lon, alt, capp, bazi, bdev
      double precision, dimension (:,:), intent (in)  :: incl, d, tt, tloss, hmin, hmax, hmean
      integer(kind=4)       , dimension (:,:), intent (in)  :: nofrefl
      character (len=23)               , intent (in)  :: tstr
      character (len=2)               , intent (in)  :: re_unit
!
!.....Local variables
!
      integer(kind=4)  :: ncid 
      integer(kind=4)  :: lon0_var, lat0_var, h0_var, lon1_var, lat1_var, h1_var
      integer(kind=4)  :: lon_var, lat_var, h_var, t_var, apoi_var, bpoi_var
      integer(kind=4)  :: azi_dim, elev_dim, t_dim
      integer(kind=4)  :: azi_var, elev_var, capp_var, bazi_var, incl_var, d_var
      integer(kind=4)  :: tt_var, tloss_var, hmax_var, hmin_var, hmean_var, refl_var, bdev_var
!
  integer(kind=4), parameter  :: shuffle = 1 ! If non-zero, turn on the shuffle filter.
  integer(kind=4), parameter  :: deflate = 1 ! If non-zero, turn on deflate specified by level.
  integer(kind=4), parameter  :: level   = 9 ! Set deflate level [1,9]
!
      integer(kind=4), dimension ( 2 )  :: dimids
      character(len=40) :: tloss_unit
!
!  ---
!
!.....Create the netCDF file.
!
      call check ( nf90_create ( path=file1, cmode=NF90_NETCDF4, ncid=ncid ) )
!
!.....Define time
!
      call check ( nf90_def_dim ( ncid, time_short_name, len=23, dimid = t_dim ) )
      call check ( nf90_def_var ( ncid, time_short_name, nf90_char, t_dim, t_var ) )
      call check ( nf90_put_att ( ncid, t_var, long_name, time_long_name ) )
      call check ( nf90_put_att ( ncid, t_var, units, time_units ) )
!
!.....Define source coordinates
!
      call check ( nf90_def_var ( ncid, lat0_short_name, nf90_float, lat0_var ) )
      call check ( nf90_put_att ( ncid, lat0_var, long_name, lat0_long_name ) )
      call check ( nf90_put_att ( ncid, lat0_var, units, lat__units ) )
!
      call check ( nf90_def_var ( ncid, lon0_short_name, nf90_float, lon0_var ) )
      call check ( nf90_put_att ( ncid, lon0_var, long_name, lon0_long_name ) )
      call check ( nf90_put_att ( ncid, lon0_var, units, lon__units ) )
!
      call check ( nf90_def_var ( ncid, alt0_short_name, nf90_float, h0_var ) )
      call check ( nf90_put_att ( ncid, h0_var, long_name, alt0_long_name ) )
      call check ( nf90_put_att ( ncid, h0_var, units, alt__units ) )
!
!.....Define receiver coordinates
!
      call check ( nf90_def_var ( ncid, lat1_short_name, nf90_float, lat1_var ) )
      call check ( nf90_put_att ( ncid, lat1_var, long_name, lat1_long_name ) )
      call check ( nf90_put_att ( ncid, lat1_var, units, lat__units ) )
!
      call check ( nf90_def_var ( ncid, lon1_short_name, nf90_float, lon1_var ) )
      call check ( nf90_put_att ( ncid, lon1_var, long_name, lon1_long_name ) )
      call check ( nf90_put_att ( ncid, lon1_var, units, lon__units ) )
!
      call check ( nf90_def_var ( ncid, alt1_short_name, nf90_float, h1_var ) )
      call check ( nf90_put_att ( ncid, h1_var, long_name, alt1_long_name ) )
      call check ( nf90_put_att ( ncid, h1_var, units, alt__units ) )
!
!.....Define bearing angles
!
      call check ( nf90_def_var ( ncid, apoi_short_name, nf90_float, apoi_var ) )
      call check ( nf90_put_att ( ncid, apoi_var, long_name, apoi_long_name ) )
      call check ( nf90_put_att ( ncid, apoi_var, units, azi__units ) )
!
      call check ( nf90_def_var ( ncid, bpoi_short_name, nf90_float, bpoi_var ) )
      call check ( nf90_put_att ( ncid, bpoi_var, long_name, bpoi_long_name ) )
      call check ( nf90_put_att ( ncid, bpoi_var, units, azi__units ) )
!
!.....Define azimuth and elevation dimension
!
      call check ( nf90_def_dim ( ncid, azi__short_name, size( azi  ), azi_dim  ) )
      call check ( nf90_def_dim ( ncid, elev_short_name, size( elev ), elev_dim ) )
!
!.....Define azimuth and elevation variables
!
      call check ( nf90_def_var ( ncid, azi__short_name, nf90_float, azi_dim, azi_var ) )
      call check ( nf90_put_att ( ncid, azi_var, long_name, azi__long_name ) )
      call check ( nf90_put_att ( ncid, azi_var, units    , azi__units ) )
!
      call check ( nf90_def_var ( ncid, elev_short_name, nf90_float, elev_dim, elev_var ) )
      call check ( nf90_put_att ( ncid, elev_var, long_name, elev_long_name ) )
      call check ( nf90_put_att ( ncid, elev_var, units    , elev_units ) )!
!
      dimids = (/ azi_dim, elev_dim /)
!
!.....Define variables
!
      call check ( nf90_def_var ( ncid, capp_short_name, nf90_float, dimids, capp_var ) )
      call check ( nf90_put_att ( ncid, capp_var, long_name, capp_long_name ) )
      call check ( nf90_put_att ( ncid, capp_var, units, capp_units ) )
      call check ( nf90_put_att ( ncid, capp_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, bazi_short_name, nf90_float, dimids, bazi_var ) )
      call check ( nf90_put_att ( ncid, bazi_var, long_name, bazi_long_name ) )
      call check ( nf90_put_att ( ncid, bazi_var, units, azi__units ) )
      call check ( nf90_put_att ( ncid, bazi_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, bdev_short_name, nf90_float, dimids, bdev_var ) )
      call check ( nf90_put_att ( ncid, bdev_var, long_name, bdev_long_name ) )
      call check ( nf90_put_att ( ncid, bdev_var, units, azi__units ) )
      call check ( nf90_put_att ( ncid, bdev_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, incl_short_name, nf90_float, dimids, incl_var ) )
      call check ( nf90_put_att ( ncid, incl_var, long_name, incl_long_name ) )
      call check ( nf90_put_att ( ncid, incl_var, units, incl_units ) )
      call check ( nf90_put_att ( ncid, incl_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, lon__short_name, nf90_float, dimids, lon_var ) )
      call check ( nf90_put_att ( ncid, lon_var, long_name, lon__long_name ) )
      call check ( nf90_put_att ( ncid, lon_var, units, lon__units ) )
      call check ( nf90_put_att ( ncid, lon_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, lat__short_name, nf90_float, dimids, lat_var ) )
      call check ( nf90_put_att ( ncid, lat_var, long_name, lat__long_name ) )
      call check ( nf90_put_att ( ncid, lat_var, units, lat__units ) )
      call check ( nf90_put_att ( ncid, lon_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, alt__short_name, nf90_float, dimids, h_var ) )
      call check ( nf90_put_att ( ncid, h_var, long_name, alt__long_name ) )
      call check ( nf90_put_att ( ncid, h_var, units, alt__units ) )
      call check ( nf90_put_att ( ncid, h_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, tt___short_name, nf90_float, dimids, tt_var ) )
      call check ( nf90_put_att ( ncid, tt_var, long_name, tt___long_name ) )
      call check ( nf90_put_att ( ncid, tt_var, units, tt___units ) )
      call check ( nf90_put_att ( ncid, tt_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, d____short_name, nf90_float, dimids, d_var ) )
      call check ( nf90_put_att ( ncid, d_var, long_name, d____long_name ) )
      call check ( nf90_put_att ( ncid, d_var, units, d____units ) )
      call check ( nf90_put_att ( ncid, d_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, loss_short_name, nf90_float, dimids, tloss_var ) )
      call check ( nf90_put_att ( ncid, tloss_var, long_name, loss_long_name ) )
      tloss_unit=trim(loss_units)//trim(re_unit)
      call check ( nf90_put_att ( ncid, tloss_var, units, trim(tloss_unit) ) )
      call check ( nf90_put_att ( ncid, tloss_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, hmin_short_name, nf90_float, dimids, hmin_var ) )
      call check ( nf90_put_att ( ncid, hmin_var, long_name, hmin_long_name ) )
      call check ( nf90_put_att ( ncid, hmin_var, units, alt__units ) )
      call check ( nf90_put_att ( ncid, hmin_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, hmax_short_name, nf90_float, dimids, hmax_var ) )
      call check ( nf90_put_att ( ncid, hmax_var, long_name, hmax_long_name ) )
      call check ( nf90_put_att ( ncid, hmax_var, units, alt__units ) )
      call check ( nf90_put_att ( ncid, hmax_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, hmn__short_name, nf90_float, dimids, hmean_var ) )
      call check ( nf90_put_att ( ncid, hmean_var, long_name, hmn__long_name ) )
      call check ( nf90_put_att ( ncid, hmean_var, units, alt__units ) )
      call check ( nf90_put_att ( ncid, hmean_var, att_missing, val_missing ) )
!
      call check ( nf90_def_var ( ncid, refl_short_name, nf90_int, dimids, refl_var ) )
      call check ( nf90_put_att ( ncid, refl_var, long_name, refl_long_name ) )
      call check ( nf90_put_att ( ncid, refl_var, units, refl_units ) )
      call check ( nf90_put_att ( ncid, refl_var, att_missing, val_missing ) )
!
!.....Set compression
!
      call check( nf90_def_var_deflate ( ncid, capp_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, bazi_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, bdev_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, lat_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, lon_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, incl_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, h_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, tt_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, d_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, tloss_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, hmin_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, hmax_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, hmean_var, shuffle, deflate, level ) )
      call check( nf90_def_var_deflate ( ncid, refl_var, shuffle, deflate, level ) )
!
!.....End definitions
!
      call check ( nf90_enddef ( ncid ) )
!
!.....Write data
!
      call check ( nf90_put_var ( ncid, t_var, tstr ) )
!
      call check ( nf90_put_var ( ncid, lat0_var, real(source(1) )) )
      call check ( nf90_put_var ( ncid, lon0_var, real(source(2) )) )
      call check ( nf90_put_var ( ncid, h0_var  , real(source(3) )) )
!
      call check ( nf90_put_var ( ncid, lat1_var, real(receiver(1) )) )
      call check ( nf90_put_var ( ncid, lon1_var, real(receiver(2) )) )
      call check ( nf90_put_var ( ncid, h1_var  , real(receiver(3) )) )
!
      call check ( nf90_put_var ( ncid, apoi_var , real( azi_poi)) )
      call check ( nf90_put_var ( ncid, bpoi_var , real(bazi_poi)) )
!
      call check ( nf90_put_var ( ncid, azi_var  , real(azi  )) )
      call check ( nf90_put_var ( ncid, elev_var , real(elev )) )
      call check ( nf90_put_var ( ncid, lat_var  , real(lat  )) )
      call check ( nf90_put_var ( ncid, lon_var  , real(lon  )) )
      call check ( nf90_put_var ( ncid, h_var    , real(alt  )) )
      call check ( nf90_put_var ( ncid, capp_var , real(capp )) )
      call check ( nf90_put_var ( ncid, bazi_var , real(bazi )) )
      call check ( nf90_put_var ( ncid, bdev_var , real(bdev )) )
      call check ( nf90_put_var ( ncid, incl_var , real(incl )) )
      call check ( nf90_put_var ( ncid, d_var    , real(d    )) )
      call check ( nf90_put_var ( ncid, tt_var   , real(tt   )) )
      call check ( nf90_put_var ( ncid, tloss_var, real(tloss)) )
      call check ( nf90_put_var ( ncid, hmin_var , real(hmin )) )
      call check ( nf90_put_var ( ncid, hmax_var , real(hmax )) )
      call check ( nf90_put_var ( ncid, hmean_var, real(hmean)) )
      call check ( nf90_put_var ( ncid, refl_var , nofrefl    ) )
!
!.....Close netcdf file
!
      call check ( nf90_close ( ncid ) )
!
      return
!
   End subroutine write_ierf3d
!
! =========================================================================================


! =========================================================================================
!
!..subroutine check 
!
!  Check ierr of NETCDF routines and print error message and terminate program.
!
!****************************
   Subroutine check ( istatus )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      integer, intent ( in )  :: istatus
!
!  ---
!
      if ( istatus .ne. nf90_noerr ) then
         print *, trim(  adjustl( nf90_strerror ( istatus ) )  )
         stop
      end if
!
      return
!
   End subroutine check
!
! =========================================================================================
!
End module ierf3d_netcdf

