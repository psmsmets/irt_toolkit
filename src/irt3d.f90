Program irt3d
!*****************************************************************************80
!
!  Infrasound Ray Tracer 3D - IRT3D 
!
!  Description:
!
!    Main program to simulate infrasound rays in the atmosphere, using 3D ECMWF
!    atmospheric profiles of temperature and wind.
!
!    Run the ray tracer by executing the run.sh file
!
!  Modified:
!
!    2 August 2017
!
!  Author:
!
!    Pieter Smets
!
!    Seismology Division,
!    Royal Netherlands Meteorological Institute (KNMI)
!    De Bilt, The Netherlands
!
!*****************************************************************************80

  use syms
  use string_utility
  use tools_atmo
  use tools_topo
  use tools_netcdf, only : isnetcdf
  use math
  use eom
  use cubicweight
  use omp_lib
!
   IMPLICIT NONE
! 
! GLOBAL VARIABLES
!
!.....input file variables
!
      CHARACTER :: in_inputfile*256, in_atmo*256, in_topo*256, in_reverse*10
      CHARACTER :: out_prefix*256, out_refl*256, out_rays*256, out_angle*256, dummy*256
      CHARACTER :: verb*128, cross_ver*128, re_unit*2
!
      CHARACTER, DIMENSION ( 999 ) :: marker_sym*10
!
      LOGICAL     :: lexist, verbose, lreflection, lsaveRays, lsaveRefl 
      LOGICAL     :: jacobian, jacobian_density, jacobian_det_is_zero 
!
      INTEGER(kind=4) :: i, j, ierr, cnt, out_type, io, marker_cnt 
      INTEGER(kind=4) :: nofElev, nofAzi, nofRays, nofRefl, nofCaustics
      INTEGER(kind=4) :: hours, minutes, progress, dprogress
!
      DOUBLE PRECISION :: lat0, lon0, h0, elevMin, elevMax, elevStep, re_unit_conversion
      DOUBLE PRECISION :: aziMin, aziMax, aziStep, rk_stepsize
      DOUBLE PRECISION :: phMin, phMax, thMin, thMax, rMin, rMax, altRefl, altMax, refl_hmin
      DOUBLE PRECISION :: azi1, elev1, h, hdLon, hdLat, s, factor, A0, d, tloss, frq, A
      DOUBLE PRECISION :: t1, t2, time, seconds, jacobian_det, c, c0, rho, rho0
      DOUBLE PRECISION :: px, py, pz, inclPrev, rk_stepvar, absorption, alpha
!
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION ( :, : ) :: angles
!
      DOUBLE PRECISION, DIMENSION ( 6 )         :: atmo_specs
      DOUBLE PRECISION, DIMENSION ( 19 )        :: f0, f, df, k1, k2, k3, k4
      DOUBLE PRECISION, DIMENSION (   999, 6 )  :: marker
      DOUBLE PRECISION, DIMENSION (  9999, 7 )  :: store_refl
      DOUBLE PRECISION, DIMENSION ( 99999, 7 )  :: store_rays
!
! ----
!
! --------------------------------------------------------------------------------------------------
!
! GET START TIME
!
! --------------------------------------------------------------------------------------------------
!
     call cpu_time( t1 )
!
! --------------------------------------------------------------------------------------------------
!
! GET PARAMETERS FROM INPUT FILE
!
! --------------------------------------------------------------------------------------------------
!
      INCLUDE 'irt_argc.f90'
!
! --------------------------------------------------------------------------------------------------
!
! DEFINE TRUE OR FALSE
!
! --------------------------------------------------------------------------------------------------
!
!.....verbose true or false
!
      verbose = ( verb( 1 : 1 ) .EQ. 'v' .OR. verb .EQ. 'y' .OR. verb .EQ. 'yes' )
!
!.....reverse horizontal wind
!
      atmo%lreverse = ( in_reverse .EQ. 'y' .OR. in_reverse .EQ. 'yes' )
!
!.....topography
!
      dummy = StrLowCase ( TRIM ( in_topo ) )
      ltopo = ( dummy .NE. 'n' .AND. dummy .NE. 'no' .AND. dummy .NE. '-' )
!
!.....prepare output
!
      lsaveRefl = ( out_type .EQ. 1 .OR. out_type .EQ. 3 )
      lsaveRays = ( out_type .EQ. 2 .OR. out_type .EQ. 3 )
!
!.....enable transmission loss
!
      jacobian = .true.
!
!.....enable transmission loss by density variations
!
      jacobian_density = .true.
!
!.....enable absorption
!
!     absorption = .false.
!
! --------------------------------------------------------------------------------------------------
!
! PRINT HEADER
!
! --------------------------------------------------------------------------------------------------
!
      IF ( verbose ) THEN
         PRINT "(A)", ' ______________________________________________________________________________ '
         PRINT "(A)", '|                                                                              |'
         PRINT "(A)", '|                                                                              |'
         PRINT "(A)", '|                       IRT3D - Infrasound Ray Tracer 3D                       |'
         PRINT "(A)", '|                                                                              |'
         PRINT "(A)", '|______________________________________________________________________________|'
         PRINT "(A)", ' '
      ELSE
         PRINT "(A)", 'IRT3D started.'
      ENDIF
!     
! --------------------------------------------------------------------------------------------------
!
! CHECK INPUT
!
! --------------------------------------------------------------------------------------------------
!
      INCLUDE 'irt_check.f90'
!
! --------------------------------------------------------------------------------------------------
!
! PREPARE OUTPUT FILENAMES
!
! --------------------------------------------------------------------------------------------------
!
      out_angle = TRIM ( out_prefix ) // '_angle.out'
      out_rays  = TRIM ( out_prefix ) // '_rays.out'
      out_refl  = TRIM ( out_prefix ) // '_refl.out'
!     
! --------------------------------------------------------------------------------------------------
!
! SET NUMBER OFF ...
!
! --------------------------------------------------------------------------------------------------
!
      INCLUDE 'irt_nof.f90'
!
! --------------------------------------------------------------------------------------------------
!
! PRINT PARAMETERS
!
! --------------------------------------------------------------------------------------------------
!
      INCLUDE 'irt_params.f90'
!
! --------------------------------------------------------------------------------------------------
!
! GET ATMOSPHERIC PROFILE AND TOPOGRAPHY
!
! --------------------------------------------------------------------------------------------------
!
      INCLUDE 'irt_grib.f90'
!
! --------------------------------------------------------------------------------------------------
!
! SET RANGES
!
! --------------------------------------------------------------------------------------------------
!
      INCLUDE 'irt_ranges.f90'
!
! --------------------------------------------------------------------------------------------------
!
! RAYTRACE (rk4 integration)
!
! --------------------------------------------------------------------------------------------------
!
      INCLUDE 'irt3d_raytrace.f90'
!
! -----------------------------------------------------------------------------------------
!
! CLEAN UP
!
! -----------------------------------------------------------------------------------------
!
! 
!.....Deallocate all atmo arrays
!
      CALL ATMO_DEALLOCATE ( 'all' )
      CALL TOPO_DEALLOCATE ()
!
! -----------------------------------------------------------------------------------------
!
! GET END TIME
!
! -----------------------------------------------------------------------------------------
!
      call cpu_time( t2 )
      time = t2 - t1
!
      hours   = FLOOR ( time / 3600.d0 )
      minutes = FLOOR ( time / 60.d0 - hours * 60.d0 )
      seconds = time - hours * 3600.d0 - minutes * 60.d0
!
! -----------------------------------------------------------------------------------------
!
! DISPLAY TERMINATE MESSAGE
!
! -----------------------------------------------------------------------------------------
!   
      IF ( verbose ) THEN
         PRINT "(A,I2,A,I2,A,F6.3,A)", 'IRT3D succesfully terminated. Elapsed time is ', &
            & hours, ' hours ', minutes, ' minutes ', seconds, ' seconds.'
         PRINT "(A)", ''
      ELSE
!!         PRINT "(A,I2,A,I2,A,F6.3,A)", 'IRT3D succesfully terminated. Elapsed time is ', &
!!            & hours, ' hours ', minutes, ' minutes ', seconds, ' seconds.'
      ENDIF
!

Contains

!  include 'irt3d_help.f90'
  include 'irt_toolkit_version.f90'

End program irt3d
