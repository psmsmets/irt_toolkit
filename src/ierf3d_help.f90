!*****************************************************************************80
!
!                               I E R F 3 D _ H E L P
!
!  Module:       /
!
!  Programmer:   Pieter S. M. Smets
!                R&D depart. of Seismology and Acoustics - Koninklijk Nederlands
!                Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         May 20, 2016
!
!  Language:     Fortran-90
!
!  Description:  This module contains a subroutine to print the ierf3d help text
!                to the stdout screen.
!
!
!*****************************************************************************80


!*****************************************************************************80
Subroutine ierf3d_show_help ( version )
!*****************************************************************************80
!
! Subroutine ierf3d_show_help
!
! Description
!
!   Show ierf3d help text
!
!*****************************************************************************80
!
  character(len=*), intent(in)  :: version

  print "(a)", ''
  print "(3a)", 'ierf3d ', trim(version), ' - Infrasound Eigen Ray Finder 3D'
  print "(a)", ''
  print "(a)", 'Usage: ierf3d parameter_file [OPTIONS]'
  print "(a)", ''    
  print "(a)", 'Options and arguments (except files and paths) are not case sensitive.'    
  print "(a)", 'Mandatory arguments to long options are mandatory for short options too.'
  print "(a)", 'General:'
  print "(2x,2a)",'-?, --help              ','Display this help text.'
! print "(2x,2a)",'-q, --quick             ','Quick and dirty, avoid all checks.'
  print "(2x,2a)",'-v, --verbose           ','Verbose information and warnings'
  print "(2x,2a)",'-n, --silent            ','No verbose.'
  print "(2x,2a)",'-p, --prefix ..         ','Add a prefix to the output file name.'
  print "(2x,2a)",'-u, --suffix ..         ','Add a suffix to the output file name.'
  print "(2x,2a)",'-o, --overwrite         ','Overwrite existing output files, by default not'
  print "(28x,a)"                           ,'allowed.'
  print "(2x,2a)",'-t, --threads [..]      ','Max number of openmp processing threads. Default = 1.'
  print "(2x,2a)",'    --debug             ','Verbose everything.'                  
  print "(2x,2a)",'    --version           ','Print version number.'                  
  print "(2x,2a)",'-i, --input             ','Print input argument file options.'                  
  print "(2x,2a)",'-a, --atmosphere        ','Set atmospheric conditions file.'
! print "(2x,2a)",'-r, --reverse           ','Reverse source-receiver and flip atmospheric wind.'
  print "(2x,2a)",'-f, --filename          ','Set filename. Default is "eigenrays_%T_%Y%M%D_%H".'
  print "(28x,a)"                           ,'Filename variables are: '
  print "(28x,a)"                           ,' %C = centre'
  print "(28x,a)"                           ,' %T = datatype'
  print "(28x,a)"                           ,' %Y = year'
  print "(28x,a)"                           ,' %M = month'
  print "(28x,a)"                           ,' %D = day of month'
  print "(28x,a)"                           ,' %O = day of year'
  print "(28x,a)"                           ,' %H = hour'
  print "(28x,a)"                           ,' %S = forecast step'
  print "(28x,a)"                           ,' %N = ensemble number'
!                                             ***************************************************80
  print "(2x,a)", ''
  stop
!
!*****************************************************************************80
End subroutine ierf3d_show_help
!*****************************************************************************80


!*****************************************************************************80
Subroutine ierf3d_show_arguments ( version )
!*****************************************************************************80
!
! Subroutine ierf3d_show_arguments
!
! Description
!
!   Show ierf3d arguments
!
!*****************************************************************************80
!
  character(len=*), intent(in)  :: version

  print "(a)", '#'
  print "(3a)",'# ierf3d v', trim(version), ' - Infrasound Eigen Ray Finder 3D'
  print "(a)", '#'
  print "(a)", '#            Input argument file'
  print "(a)", '#'
  print "(a)", '#--------------- Atmosphere -----------------'
  print "(a)", 'atmo_file		= test.grib'
  print "(a)", 'atmo_alt_min		= 0 #km'
  print "(a)", 'atmo_alt_max		= 0 #km'
  print "(a)", 'atmo_alt_step		= 0 #km'
  print "(a)", 'atmo_fc_step		= 0'
  print "(a)", 'atmo_ens_num		= 0'
  print "(a)", 'atmo_reverse		= F'
!  print "(a)", 'atmo_taper		= T # taper wind from ground to taper_height'
!  print "(a)", 'atmo_taper_height	= 3.0 # km'
  print "(a)", '#--------------- Topography -----------------'
  print "(a)", 'topo			= T'
  print "(a)", 'topo_file		= atmosphere # atmosphere | path_to_dem'
  print "(a)", '#------------- Source location --------------'
  print "(a)", 'source_lat		= #deg N'
  print "(a)", 'source_lon		= #deg E'
  print "(a)", 'source_alt		= #km'
  print "(a)", '#------------ Receiver location -------------'
  print "(a)", 'receiver_lat		= #deg N'
  print "(a)", 'receiver_lon		= #deg E'
  print "(a)", 'receiver_alt		= #km'
  print "(a)", '#---------------- Ray domain ----------------'
  print "(a)", 'azimuth			= #deg'
  print "(a)", 'azimuth_first		= #deg'
  print "(a)", 'azimuth_last		= #deg'
  print "(a)", 'azimuth_step		= #deg'
  print "(a)", 'azimuth_count		= #deg'
  print "(a)", 'azimuth_range		= #deg; set azimuth range centered on bearing'
  print "(a)", 'azimuth_offset		= #deg; set azimuth range centered on bearing'
  print "(a)", 'azimuth_offset_first	= #deg; set azimuth range left from bearing'
  print "(a)", 'azimuth_offset_last	= #deg; set azimuth range right from bearing'
  print "(a)", 'elevation		= #deg'
  print "(a)", 'elevation_first		= #deg'
  print "(a)", 'elevation_last		= #deg'
  print "(a)", 'elevation_step		= #deg'
  print "(a)", 'elevation_count		= #deg'
  print "(a)", '#------- Transmission loss / Jacobian -------'
  print "(a)", 'tloss_geometrical	= T'
  print "(a)", 'tloss_absorption	= F'
!  print "(a)", 'tloss_density		= T'
!  print "(a)", 'tloss_coherent		= T'
  print "(a)", 'tloss_frequency		= 1.5 #Hz'
  print "(a)", 'tloss_re_unit		= m'
  print "(a)", '#------------- Solver options ---------------'
  print "(a)", 'rk_stepsize		= 2.0 #s'
  print "(a)", 'rk_adaptive		= T'
!  print "(a)", 'rk_stepvar		= 10.0 #-'
  print "(a)", '#------------------ Output ------------------'
  print "(a)", 'prefix			= '
  print "(a)", 'suffix			= '
  print "(a)", 'filename		= '
  print "(a)", 'overwrite               = F'
  print "(a)", 'clip 			= T'
  print "(a)", 'clip_max_offset		= 25.0 #km'
  print "(a)", '#-------- Miscellaneous Parameters ----------'
  print "(a)", 'verbose                 = F'
  print "(a)", 'verbose_level           = 1'
  print "(a)", 'debug			= F'
  print "(a)", 'parallel_threads        = 2'
  stop
!
!*****************************************************************************80
End subroutine ierf3d_show_arguments
!*****************************************************************************80
