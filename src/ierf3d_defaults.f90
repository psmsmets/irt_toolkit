!*****************************************************************************80
!
!                                 I E R F _ D E F A U L T S
!
!  INCLUDE:      ierf_defaults.f90
!
!  Programmer:   Pieter S. M. Smets(1,2)
!                (1) Seismology Devision - Koninklijk Nederlands
!                    Meteorologisch Instituut (KNMI), De Bilt, The Netherlands
!                (2) Geoscience and Engineering - Delft University of Technology
!                    Delft, The Netherlands
!
!  Modified:     August 2, 2017
!
!  Language:     Fortran-90
!
!  Description:  ierf3d default values 
!
!*****************************************************************************80

!
! general
!
  inf                  = ''
  inu                  = 1_int8
!
! atmosphere
!
  atmo_path            = ''
  atmo_grib            = .true.
  atmo_nc              = .false.
  atmo%hmin            = 0._double
  atmo%hmax            = 60._double
  atmo%dh              = 0.5_double
  atmo%ens             = 0_int16
  atmo%step            = 0_int16
  atmo%lreverse        = .false.
!
! topography
!
  ltopo                = .true.
  topo_path            = ''
  topo_from_atmo       = .false.
  topo_from_dem        = .false.
!
! source
!
  lat0                 = missing_double
  lon0                 = missing_double
  h0                   = missing_double
!
! receiver
!
  lat1                 = missing_double
  lon1                 = missing_double
  h1                   = missing_double
!
! ray parameters
!
  azimuth              = missing_double
  azimuth_first        = missing_double
  azimuth_last         = missing_double
  azimuth_step         =   0.5_double
  azimuth_count        = missing_int32
  azimuth_range        = missing_double
  azimuth_offset       = missing_double
  azimuth_offset_first = missing_double
  azimuth_offset_last  = missing_double
!
  elevation            = missing_double
  elevation_first      =  0.0_double
  elevation_last       = 90.0_double
  elevation_step       =  0.5_double
  elevation_count      = missing_int32
!
! transmission loss
!
  tloss_attenuation    = .true.
  tloss_geometrical    = .true.
  tloss_density        = .true. ! density depending impedance
  tloss_absorption     = .false. ! Sutherland and Bass
  tloss_frequency      = 1.0_double
  tloss_re_unit        = 'm'
!
! solver options
!
  rk_adaptive          = .true.
  rk_stepsize          = 2.0_double
!
! output
!
  prefix               = ''
  suffix               = ''
  filename             = ''
  overwrite            = .false.
  clip                 = .true.
  clip_max_offset      = 25.0_double
!
! miscellaneous parameters
!
  verb                 = .true.
  verblvl              = 1_int8
  debug                = .false.
  parallel             = .false.
  threads              = -1_int32
  max_threads          = omp_get_max_threads()
!
