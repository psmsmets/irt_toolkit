! ***************************************************************************
!
!                        Infrasound Eigenray Ray Finder 3D
!                                     IERF3D
!
! ***************************************************************************
!
! Main program to simulate infrasound rays in the atmosphere, using 3D ECMWF
! atmospheric profiles of temperature and wind (91 level model).
!
!
! Run the ray tracer by executing the run.csh file (version 1.2, March 2012)
!
! ---
! Pieter Smets (smets@knmi.nl)
!
! Seismology Division,
! Royal Netherlands MeteoroLOGICAL Institute (KNMI)
! De Bilt, The Netherlands

! ==================================================================================================
! //////////////////////////////////////////////////////////////////////////////////////////////////
! ==================================================================================================
!
!****************************
Program ierf3d 
!****************************
!
  use syms
  use string_utility
  use io_types
  use tools_atmo
  use tools_topo
  use tools_netcdf, only: isnetcdf
  use ierf3d_netcdf, only : write_ierf3d
  use math
  use eom
  use cubicweight
  use time
  use omp_lib
!!  use, intrinsic :: iso_fortran_env
!
   implicit none
! 
! Main parameters
!
  integer(int32), parameter :: max_len=256_int32
  integer(int32), parameter :: missing_int32=-9999_int32
  real(double)  , parameter :: missing_double=-9999._double
! 
! Main variables
!
  character(len=max_len)  :: inf, atmo_path, topo_path, prefix, suffix, filename
  character(len=max_len)  :: arg_file, arg_type, arg_value, arg_line, dummy, errpref, warpref
  character(len=2) :: tloss_re_unit
!
  logical, target  :: verb
  logical   :: arg_parsed, exists, debug, overwrite, clip, parallel, eof
  logical   :: atmo_grib, atmo_nc, topo_from_dem, topo_from_atmo, rk_adaptive
  logical   :: tloss_attenuation, tloss_geometrical, tloss_absorption, tloss_coherent, tloss_density
  logical   :: jacobian_det_is_zero
!
  integer(int8)   :: inu, verblvl
  integer(int32)  :: error, i, j, ia, ie, cnt, nofrefl, nofs, arg_cnt 
  integer(int32)  :: hours, minutes, progress, dprogress, threads, max_threads
  integer(int32)  :: azimuth_count, elevation_count, rays_count, nofCaustics
!
  real(double), target  :: bearing_sr, bearing_rs, bearing_d
  real(double)  :: lat0, lon0, h0, lat1, lon1, h1, tloss_frequency, clip_max_offset
  real(double)  :: azimuth, azimuth_first, azimuth_last, azimuth_step, azimuth_range
  real(double)  :: azimuth_offset, azimuth_offset_first, azimuth_offset_last
  real(double)  :: elevation, elevation_first, elevation_last, elevation_step, elevation_range
  real(double)  :: ph1, th1, r1, dd, lat, lon, px, py, pz, inclPrev
  real(double)  :: aziMin, aziMax, aziStep, rk_stepsize, rk_stepvar, re_unit_conversion
  real(double)  :: phMin, phMax, thMin, thMax, rMin, rMax, refl_hmin, refl_hmax, altmin, altmax, altmean
  real(double)  :: azi1, elev1, h, hdLon, hdLat, s, factor, A0, d, ud, tloss, frq, A
  real(double)  :: runtime, seconds, jacobian_det, rho, rho0, c, c0, absorption, alpha 
!
! arrays
!
  real(double), allocatable, dimension(:)  :: azimuth_vec, elevation_vec
  integer(int32) , allocatable, dimension ( :, : )    :: store_refl
  real(double), allocatable, dimension ( :, : )    :: store_lat, store_lon, store_alt
  real(double), allocatable, dimension ( :, : )    :: store_tt, store_d, store_tloss
  real(double), allocatable, dimension ( :, : )    :: store_bazi, store_bazidev, store_capp, store_incl
  real(double), allocatable, dimension ( :, : )    :: store_hmax, store_hmin, store_hmean
!
  real(double), dimension ( 6 )       :: atmo_specs
  real(double), dimension ( 19 )      :: f0, f, df, k1, k2, k3, k4
!
! pointers
!
  logical     , pointer :: verbose
  real(double), pointer :: azi_poi, bazi_poi
!
! ----
!
! ----------------------------------------------------------------------------80
!
! Initialize 
!
! ----------------------------------------------------------------------------80
!
! Get start time
!
  runtime = omp_get_wtime()
!
! Set defaults parameters
!
  include 'ierf3d_defaults.f90'
!
! Set error strings
!
  errpref = 'Error :'
  warpref = 'Warning :'
!
! ----------------------------------------------------------------------------80
!
! Get parameters from command line and input file
!
! ----------------------------------------------------------------------------80
!
! Parse command line arguments
!
  include 'ierf3d_parse_command_args.f90'
!
! Check input parameters
!
  include 'ierf3d_check_args.f90'
!
! ----------------------------------------------------------------------------80
!
! Parallelize
!
! ----------------------------------------------------------------------------80
!
  if (parallel) call omp_set_num_threads( threads )
!
! ----------------------------------------------------------------------------80
!
! Print header
!
! ----------------------------------------------------------------------------80
!
  verbose=>verb
  if (verb) then
!
    print "(a)", ''
    print "(a)", 'IERF3D - Infrasound Eigenray Ray Finder 3D'
    print "(a)", ''
!
    if (debug) print "(a)", '// Debug mode //'
!
    print "('> ',a)" , 'Atmosphere'
    print "(4x,2a)"  , 'File name              = ', trim(atmo_path)
    print "(4x,2a)"  , 'Centre                 = ', trim(atmo%ranges%centre)
    print "(4x,2a)"  , 'Edition                = ', trim(atmo%ranges%editionnumber)
    print "(4x,2a)"  , 'Date time              = ', trim(atmo%ranges%tstr)
    print "(4x,2a)"  , 'Stream                 = ', trim(atmo%ranges%stream)
    print "(4x,2a)"  , 'Datatype               = ', trim(atmo%ranges%datatype)
    print "(4x,2a)"  , 'Gridtype               = ', trim(atmo%ranges%gridtype)
    print "(4x,2a)"  , 'TypeOfLevel            = ', trim(atmo%ranges%typeOfLevel)
    print "(4x,a,i4)", 'Levels             (-) = ', atmo%ranges%noflvl
    print "(4x,a,i4)", 'Steps              (-) = ', atmo%ranges%nofstp
    if (atmo%ranges%nofens.gt.0) &
      print "(4x,a,i4)", 'Ensembles          (-) = ', atmo%ranges%nofens
    print "(4x,2(a,f6.2),a,f5.2,a,i4,a)", "Latitude       (deg N) = ", &
      & atmo%latMin, " to ", atmo%latMax, ", ", atmo%dlat, " (#", atmo%noflat, ")"
    if ( atmo%l360 ) then
      print "(4x,2(a,f6.2),a,f5.2,a,i4,a)", "Longitude      (deg E) = ", &
        0.0d0, " to ", 360.0d0, ", ", atmo%ranges%dlon, " (#", atmo%ranges%noflon, ")"
    else
      print "(4x,2(a,f6.2),a,f5.2,a,i4,a)", "Longitude      (deg E) = ", &
        atmo%lonMin, " to ", atmo%lonMax, ", ", atmo%dlon, " (#", atmo%noflon, ")"
    end if
    print "(4x,2(a,f6.2),a,f5.2,a,i4,a)", "Altitude          (km) = ", &
      atmo%hMin, " to ", atmo%hmax, ", ", atmo%dH, " (#", atmo%nofh, ")"
    print "(4x,a,1l)"  , 'Reverse atmosphere     = ', atmo%lreverse
!
    print "('> ',a)", 'Topography'
    print "(4x,a,1l)", 'Topography             = ', ltopo
    if (ltopo) then
      if (topo_from_atmo) print "(4x,a,l1)",   'Topography_from_atmo   = ', topo_from_atmo
      if (topo_from_dem ) print "(4x,a,l1)",   'Topography_from_dem    = ', topo_from_dem
      if (topo_from_dem ) print "(4x,a,2a)", 'Topography             = ', trim(topo_path)
    end if
!
    print "('> ',a)", 'Source location'
    print "(4x,a,f8.4)", 'Latitude       (deg N) = ', lat0
    print "(4x,a,f8.4)", 'Longitude      (deg E) = ', lon0
    print "(4x,a,f8.4)", 'Altitude          (km) = ', h0
    print "(4x,a,f8.4)", 'Start bearing    (deg) = ', bearing_rs
!
    print "('> ',a)", 'Receiver location'
    print "(4x,a,f8.4)", 'Latitude       (deg N) = ', lat1
    print "(4x,a,f8.4)", 'Longitude      (deg E) = ', lon1
    print "(4x,a,f8.4)", 'Altitude          (km) = ', h1
    print "(4x,a,f8.4)", 'Final bearing    (deg) = ', bearing_sr
    print "(4x,a,f8.4)", 'Bearing distance  (km) = ', bearing_d
!
    print "('> ',a)", 'Ray domain'
    print "(4x,2(a,f8.3),a,f8.3,a,i4,a)", 'Elevation        (deg) = ', &
      elevation_first, " to ", elevation_last, ", ", &
      elevation_step, " (#", elevation_count, ")"
    print "(4x,2(a,f8.3),a,f8.3,a,i4,a)", 'Azimuth          (deg) = ', &
      azimuth_first, " to ", azimuth_last, ", ", &
      azimuth_step, " (#", azimuth_count, ")"
!
    print "('> ',a)", 'Transmission loss'
    print "(4x,a,1l)", 'Attenuation            = ', tloss_attenuation
    if (tloss_attenuation) then
      print "(4x,a,1l)",   'Geometrical loss       = ', tloss_geometrical
      print "(4x,a,1l)",   'Absorption             = ', tloss_absorption
      print "(4x,a,f5.2)", 'Frequency         (Hz) = ', tloss_frequency
      print "(4x,2a)",     'Reference unit         = dB re 1', trim(tloss_re_unit)
    end if
!
    print "('> ',a)", 'Solver options'
    print "(4x,a,f5.2)", "RK4 stepsize       (s) = ", rk_stepsize
    print "(4x,a,1l)"  , "Adaptive               = ", rk_adaptive
!
    print "('> ',a)", 'Output'
    print "(4x,2a)"    , 'filename               = ', trim(prefix)//trim(filename)//trim(suffix)//'.nc' 
    print "(4x,a,l1)"  , 'Overwrite              = ', overwrite
    print "(4x,a,1l)"  , "Clip                   = ", clip
    if (clip) print "(4x,a,f5.2)", "Clip max offset   (km) = ", clip_max_offset
!
    print "('> ',a)", 'Miscellaneous parameters'
    print "(4x,a,i1)"        , 'Verbose level          = ', verblvl
    print "(4x,a,i3)"        , 'OpenMP threads         = ', threads
!
  end if
!
! --------------------------------------------------------------------------------------------------
!
! Get atmospheric profile and topography
!
! --------------------------------------------------------------------------------------------------
!
    include 'ierf3d_grib.f90'
!
! --------------------------------------------------------------------------------------------------
!
! Set ranges
!
! --------------------------------------------------------------------------------------------------
!
    include 'ierf_ranges.f90'
!
! --------------------------------------------------------------------------------------------------
!
! RAYTRACE (rk4 integration)
!
! --------------------------------------------------------------------------------------------------
!
! Get back-azimuth deviation receiver-source
!
  bazi_poi=>bearing_rs
   azi_poi=>bearing_sr
!
! Allocate
!
    allocate (  &
      & azimuth_vec( azimuth_count ), &
      & elevation_vec( elevation_count ), &
      & store_lat( azimuth_count, elevation_count ), &
      & store_lon( azimuth_count, elevation_count ), &
      & store_alt( azimuth_count, elevation_count ), &
      & store_capp( azimuth_count, elevation_count ), &
      & store_bazi( azimuth_count, elevation_count ), &
      & store_bazidev( azimuth_count, elevation_count ), &
      & store_incl( azimuth_count, elevation_count ), &
      & store_tt( azimuth_count, elevation_count ), &
      & store_tloss( azimuth_count, elevation_count ), &
      & store_d( azimuth_count, elevation_count ), &
      & store_hmin( azimuth_count, elevation_count ), &
      & store_hmax( azimuth_count, elevation_count ), &
      & store_hmean( azimuth_count, elevation_count ), &
      & store_refl( azimuth_count, elevation_count ) &
      &  )
!
! ..Init
!
      store_lat=missing_double
      store_lon=missing_double
      store_alt=missing_double
      store_capp=missing_double
      store_bazi=missing_double
      store_bazidev=missing_double
      store_incl=missing_double
      store_tt=missing_double
      store_tloss=missing_double
      store_d=missing_double
      store_hmin=missing_double
      store_hmax=missing_double
      store_hmean=missing_double
      store_refl=missing_double
!
! Parallelize
!
    call omp_set_num_threads ( threads )
!
! Ray trace code
!
    include 'ierf3d_raytrace.f90'
!
! Calculate back-azimuth deviation for rays
!
    where(store_bazi.gt.missing_double) store_bazidev = store_bazi - lonrangef(bazi_poi)
!
!   Calculate altitude above sea level for source and receiver
!
    if ( ltopo ) then
      call TOPO_GET ( lon0, lat0, h, hdLon, hdLat )
      if ( h .gt. h0 ) h0 = h
      call TOPO_GET ( lon1, lat1, h, hdLon, hdLat )
      if ( h .gt. h1 ) h1 = h
    end if
!     
! Write netcdf output file
!
    call write_ierf3d ( trim(prefix)//trim(filename)//trim(suffix)//'.nc', &
      & (/ lat0, lon0, h0 /), (/ lat1, lon1, h1 /), &
      & lonrangef(azi_poi), lonrangef(bazi_poi), azimuth_vec, elevation_vec, & 
      & store_lat, store_lon, store_alt, &
      & store_capp, store_bazi, store_bazidev, &
      & store_incl, store_d, store_tt, store_tloss, tloss_re_unit, & 
      & store_hmin, store_hmax, store_hmean, store_refl, atmo%ranges%tstr &
      & )
!
! -----------------------------------------------------------------------------------------
!
! CLEAN UP
!
! -----------------------------------------------------------------------------------------
!
!   Deallocate raytrace vectors and arrays
!
    deallocate( &
      azimuth_vec, elevation_vec, store_lon, store_lat, store_alt, store_capp, store_bazi, &
      store_bazidev, store_incl, store_tt, store_tloss, store_d, store_hmin, store_hmax,   &
      store_hmean, store_refl &
      & )
! 
!   Deallocate all atmo arrays
!
    call atmo_deallocate ( 'all' )
    call topo_deallocate ()
!
! -----------------------------------------------------------------------------------------
!
! GET END TIME
!
! -----------------------------------------------------------------------------------------
!
    runtime = omp_get_wtime() - runtime
!
    hours   = floor( runtime / 3600.d0 )
    minutes = floor( runtime / 60.d0 - hours * 60.d0 )
    seconds = runtime - hours * 3600.d0 - minutes * 60.d0
!
! -----------------------------------------------------------------------------------------
!
! DISPLAY TERMINATE MESSAGE
!
! -----------------------------------------------------------------------------------------
!   
    if ( verbose ) then
      print "(a,i2,a,i2,a,f6.3,a)", 'IERF3D succesfully terminated. Elapsed time is ', &
       & hours, ' hours ', minutes, ' minutes ', seconds, ' seconds.'
    end if
!
Contains

  include 'ierf3d_help.f90'
  include 'irt_toolkit_version.f90'
  include 'irt_filename.f90'

End program ierf3d
