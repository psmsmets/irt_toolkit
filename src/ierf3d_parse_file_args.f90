!*****************************************************************************80
!
!                   I E R F 3 D _ P A R S E _ F I L E _ A R G S
!
!  Module:       / (include)
!
!  Programmer:   Pieter S. M. Smets
!                R&D depart. of Seismology and Acoustics - Koninklijk Nederlands
!                Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         August 2, 2017
!
!  Language:     Fortran-90
!
!  Description:  Parse ierf3d input file arguments
!
!
!*****************************************************************************80


!
! get file_input from command-line argument (first argument)
!
  call get_command_argument( arg_cnt, arg_file )
!
! check length inputfile
!
  if ( len_trim(arg_file).ge.max_len ) then
    print "(a,1x,a,i3,a)", trim(errpref), 'input file string too long (len>= ', max_len, &
       '). Use symbolic links to avoid too long strings.'
    stop
  end if
!
! help?
!
  if (strlowcase(trim(arg_file)).eq.'--help') call ierf3d_show_help( trim(irt_toolkit_version()) )
!
! open tmp input file
!
  inquire ( file=arg_file, exist=exists ) 
!
  if (exists) then
    open( unit=inu, file=arg_file, status='old', action='read', iostat=error )
    if ( error.ne.0_int32 ) then
      print "(a,1x,3a)", trim(errpref), 'something wrong with input file "',trim(arg_file),'.'
      stop
    end if
    eof=.false.
  else
    print "(a,1x,3a)", trim(errpref), 'input file "', trim(arg_file),' not found.'
    stop
    eof=.true.
  end if
!
! loop over file
!
  do while (.not.eof.and.exists)
!
!   read line
!
    read(inu,"(a)",iostat=error) arg_line
!
!   eof?
!
    if (error.lt.0) then
      eof=.true.
      exit 
    else if (error.gt.1) then   
      print "(a,1x,a)", trim(errpref), '@ ierf_parse_arg while reading in input parameter file.'
      stop
    end if
!
!   skip #comments
!
    if (arg_line(1:1).eq.'#') cycle
!
!   split argument / value
!
    i=scan(arg_line,'=',.false.)
    j=scan(arg_line,'#',.false.)
    arg_type  = strcompress(arg_line(1:i-1))
    if (j.eq.0_int32) then
      arg_value = strcompressedges(arg_line(i+1:len(arg_line)))
    else
      arg_value = strcompressedges(arg_line(i+1:j-1))
    end if
!
!   no argument value present? Skip!
!
    if ( len_trim(arg_value).eq.0 ) cycle
!
!   overwrite default values
!
    select case ( strlowcase(arg_type) )
!
!--------------- Atmosphere -----------------
!
      case('atmo_file')
        atmo_path = arg_value
      case('atmo_alt_min')
        read(arg_value,*,iostat=error) atmo%hmin
      case('atmo_alt_max')
        read(arg_value,*,iostat=error) atmo%hmax
      case('atmo_alt_step')
        read(arg_value,*,iostat=error) atmo%dh
      case('atmo_fc_step')
        read(arg_value,*,iostat=error) atmo%step
      case('atmo_ens_num')
        read(arg_value,*,iostat=error) atmo%ens
      case('atmo_rev','atmo_reverse')
        read(arg_value,*,iostat=error) atmo%lreverse
!
!--------------- Topography -----------------
!
      case('topo')
        read(arg_value,*,iostat=error) ltopo
      case('topo_file')
        read(arg_value,*,iostat=error) topo_path
!
!------------- Source location --------------
!
      case('source_lat','source_latitude')
        read(arg_value,*,iostat=error) lat0
      case('source_lon','source_longitude')
        read(arg_value,*,iostat=error) lon0
      case('source_alt','source_altitude')
        read(arg_value,*,iostat=error) h0
!
!------------ Receiver location -------------
!
      case('receiver_lat','receiver_latitude')
        read(arg_value,*,iostat=error) lat1
      case('receiver_lon','receiver_longitude')
        read(arg_value,*,iostat=error) lon1
      case('receiver_alt','receiver_altitude')
        read(arg_value,*,iostat=error) h1
!
!--------------- Ray options ----------------
!
!     azimuth
!
      case('azimuth')
        read(arg_value,*,iostat=error) azimuth
      case('azimuth_first')
        read(arg_value,*,iostat=error) azimuth_first
      case('azimuth_last')
        read(arg_value,*,iostat=error) azimuth_last
      case('azimuth_step')
        read(arg_value,*,iostat=error) azimuth_step
      case('azimuth_count')
        read(arg_value,*,iostat=error) azimuth_count
      case('azimuth_range')
        read(arg_value,*,iostat=error) azimuth_range
      case('azimuth_offset')
        read(arg_value,*,iostat=error) azimuth_offset
      case('azimuth_offset_first')
        read(arg_value,*,iostat=error) azimuth_offset_first
      case('azimuth_offset_last')
        read(arg_value,*,iostat=error) azimuth_offset
!
!     elevation
!
      case('elevation')
        read(arg_value,*,iostat=error) elevation
      case('elevation_first')
        read(arg_value,*,iostat=error) elevation_first
      case('elevation_last')
        read(arg_value,*,iostat=error) elevation_last
      case('elevation_step')
        read(arg_value,*,iostat=error) elevation_step
      case('elevation_count')
        read(arg_value,*,iostat=error) elevation_count
!
!------- Transmission loss / Jacobian -------
!
      case('tloss_jacobian','tloss_geometrical')
        read(arg_value,*,iostat=error) tloss_geometrical
      case('tloss_absorption')
        read(arg_value,*,iostat=error) tloss_absorption
      case('tloss_density')
        read(arg_value,*,iostat=error) tloss_density
      case('tloss_coherent')
        read(arg_value,*,iostat=error) tloss_coherent
      case('tloss_frequency','tloss_freq')
        read(arg_value,*,iostat=error) tloss_frequency
      case('tloss_re_unit')
        read(arg_value,*,iostat=error) dummy
        tloss_re_unit = strlowcase(dummy(1:2))
!
!------------- Solver options ---------------
!
      case('rk_stepsize','rk_step')
        read(arg_value,*,iostat=error) rk_stepsize
      case('rk_adaptive')
        read(arg_value,*,iostat=error) rk_adaptive
      case('rk_stepvar')
        read(arg_value,*,iostat=error) rk_stepvar
!
!------------------ Output -----------------
!
      case('prefix')
        read(arg_value,*,iostat=error) prefix
      case('suffix')
        read(arg_value,*,iostat=error) suffix
      case('filename')
        read(arg_value,*,iostat=error) filename 
      case('overwrite')
        read(arg_value,*,iostat=error) overwrite
      case('clip')
        read(arg_value,*,iostat=error) clip
      case('clip_max_offset')
        read(arg_value,*,iostat=error) clip_max_offset
!
!-------- Miscellaneous Parameters ----------
!
      case('verbose')
        read(arg_value,*,iostat=error) verb
      case('verbose_level')
        read(arg_value,*,iostat=error) verblvl
      case('debug')
        read(arg_value,*,iostat=error) debug
      case('parallel_threads')
        read(arg_value,*,iostat=error) threads
!
!     otherwise ignore
!
      case default
!
        cycle
!
    end select
!
!   Error?
!
    if (error.ne.0) then
      print "(a,1x,7a)", trim(errpref), 'reading parameter "', trim(arg_type), '" with value "', &
        trim(arg_value), '" from file ', trim(arg_file), '.'
      stop
   end if
!
  end do
!
! close argument file
!
  if (exists) then
    close(inu,iostat=error)
    if (error.ne.0) then
       print "(a,1x,3a)", trim(errpref), 'cannot close input file "',trim(arg_file),'.'
       stop
    else
      arg_parsed=.true.
    end if
  end if
!
