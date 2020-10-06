!*****************************************************************************80
!
!                        I E R F 3 D _ C H E C K _ A R G S 
!
!  INCLUDE:      ierf3d_check_args.f90
!
!  Programmer:   Pieter S. M. Smets(1,2)
!                (1) Seismology Devision - Koninklijk Nederlands
!                    Meteorologisch Instituut (KNMI), De Bilt, The Netherlands
!                (2) Geoscience and Engineering - Delft University of Technology
!                    Delft, The Netherlands
!
!  Modified:     August 3, 2017
!
!  Language:     Fortran-90
!
!  Description:  check ierf3d input parameter arguments 
!
!*****************************************************************************80
 
!
! check length input strings
!
  if ( len_trim(atmo_path).ge.max_len ) then
    print "(a,1x,2a,i3,1a)", trim(errpref), 'atmosphere file string contains too many ', &
      & 'characters (>= ', max_len, '). Use symbolic links to avoid too long strings.'
    stop
  end if
  if ( len_trim(topo_path).ge.max_len ) then
    print "(a,1x,2a,i3,1a)", trim(errpref), 'topography file string contains too many ', &
      & 'characters (>= ', max_len, '). Use symbolic links to avoid too long strings.'
    stop
  end if
  if ( len_trim(prefix).ge.max_len ) then
    print "(a,1x,2a,i3,1a)", trim(errpref), 'output prefix contains too many ', &
      & 'characters (>= ', max_len, ').'
    stop
  end if
  if ( len_trim(suffix).ge.max_len ) then
    print "(a,1x,2a,i3,1a)", trim(errpref), 'output suffix contains too many ', &
      & 'characters (>= ', max_len, ').'
    stop
  end if
!     
! atmosphere
!
  inquire(file=trim(atmo_path),exist=exists)
  if (.not.exists) then
    print "(a,1x,3a)", trim(errpref), 'atmospheric profile "', trim(atmo_path), &
      & '" does not exist.'
    stop
  end if
  call atmo_set_ranges( trim(atmo_path) )
!
  if ( atmo%hmin.lt.0._double ) then
    atmo%hmin = 0.d0
    print "(a,1x,2a)", trim(warpref), 'minimal altitude cannot be smaller ', &
      & 'than zero! --> atmo_alt_min = 0.0km.'
  end if
  if ( atmo%hmax.lt.atmo%hmin ) then
    atmo%hmax = 60.d0
    print "(a,1x,2a)", trim(warpref), 'maximal altitude cannot be smaller ', &
      & 'than minimal altitude! --> atmo_alt_max = 70km.'
  end if
  if ( atmo%hmax.gt.75._double ) then
    atmo%hmax = 60.d0
    print "(a,1x,2a)", trim(warpref), 'maximal altitude cannot be larger ', &
      & 'than 75km! --> atmo_alt_max = 75km.'
  end if
  if ( atmo%dh.lt.0.1_double ) then
    atmo%dh = .1d0
    print "(a,1x,2a)", trim(warpref), 'altitude stepsize cannot be smaller ', &
      & 'than 100m --> atmo_alt_step = 0.1km.'
  end if
  if ( atmo%hmax.lt.atmo%hmin ) then
    atmo%hmax = 60.d0
    print "(a,1x,2a)", trim(warpref), 'maximal altitude cannot be smaller ', &
      & 'than minimal altitude! --> atmo_alt_max = 70km.'
  end if
  if ( atmo%ens.lt.0) then
    atmo%ens  = 0
    print "(a,1x,2a)", trim(warpref), 'ensemble number cannot be smaller than zero! ', &
      & ' --> atmo_ens_num = 0'
  end if
  if ( atmo%step.lt.0 ) then
    atmo%step  = 0
    print "(a,1x,2a)", trim(warpref), 'forcast step cannot be smaller than zero! ', &
     & ' --> atmo_fc_step = 0'
  end if
!
! topography ( orography or sea-level based atmospheric profile)
!
  if (ltopo) then
!
    dummy = strlowcase(topo_path)
    select case (dummy)
      case ('ecmwf','atmosphere','','-') ! default due to ''
        topo_from_atmo=.true.
      case default
        topo_from_dem =.true.
    end select
!
    if (topo_from_dem) then
      inquire(file=trim(topo_path),exist=exists)
      if (.not.exists) then
        print "(a,1x,3a)", trim(errpref), 'topography file "', trim(topo_path), &
          & '" does not exist.'
        stop
      end if
    end if
!
  end if
!
! source
!
  if ( h0.lt.atmo%hmin ) then
    h0 = atmo%hmin
    print "(a,1x,2a)", trim(warpref), 'source altitude cannot be smaller ', &
     & 'than minimal atmosphere altitude! --> source_alt = atmo_alt_min.'
  end if
  if ( h0.gt.atmo%hmax ) then
    h0 = atmo%hmax
    print "(a,1x,2a)", trim(errpref), 'source altitude cannot be larger ', &
     & 'than maximal atmosphere altitude! --> source_alt = atmo_alt_max.'
  end if 
  if ( lonrangef(lon0,atmo%lonmin).lt.atmo%lonmin &
    .or.lonrangef(lonrangef(lon0,atmo%lonmin)-atmo%lonmin).gt.atmo%lonrange ) then
    print "(a,1x,1a)", trim(errpref), 'source longitude is not part of the atmospheric range.'
    stop
  end if
  if ( lat0.lt.atmo%latmin.or.lat0.gt.atmo%latmax ) then
    print "(a,1x,1a)", trim(errpref), 'source latitude is not part of the atmospheric range.'
    stop
  end if
  call lonrange(lon0)
!
! receiver
!
  if ( h1.lt.atmo%hmin ) then
    h1 = atmo%hmin
    print "(a,1x,2a)", trim(errpref), 'receiver altitude cannot be smaller ', &
     & 'than minimal atmosphere altitude! --> receiver_alt = atmo_alt_min.'
  end if
  if ( h1.gt.atmo%hmax ) then
    h1 = atmo%hmax
    print "(a,1x,2a)", trim(errpref), 'receiver altitude cannot be larger ', &
     & 'than maximal atmosphere altitude! --> receiver_alt = atmo_alt_max.'
  end if
  if ( lonrangef(lon1,atmo%lonmin).lt.atmo%lonmin &
    .or.lonrangef(lonrangef(lon1,atmo%lonmin)-atmo%lonmin).gt.atmo%lonrange ) then
    print "(a,1x,1a)", trim(errpref), 'receiver longitude is not part of the atmospheric range.'
    stop
  end if
  if ( lat1.lt.atmo%latmin.or.lat1.gt.atmo%latmax ) then
    print "(a,1x,1a)", trim(warpref), 'receiver latitude is not part of the atmospheric range.'
    stop
  end if
  call lonrange(lon1)
!
! source - receiver bearing and distance
!
  call course(lat0,lon0,lat1,lon1,bearing_d,bearing_sr)
  bearing_sr=lonrangef(bearing_sr)
  call angle(lat1,lon1,lat0,lon0,bearing_rs)
  bearing_rs=lonrangef(bearing_rs+180._double)
!
! elevation
!
  if (elevation_count.ne.missing_int32.and.elevation_count.lt.0_int32) then
    elevation_count = missing_int32
    print "(a,1x,2a)", trim(warpref), 'elevation count cannot be smaller than 0! ', &
      & '--> elevation_count = missing_value.'
  end if
!
  if ( abs(elevation).gt.90._double.and.abs(elevation-missing_double).gt.1.d-8 ) then
    elevation = missing_double
    print "(a,1x,2a)", trim(warpref), 'elevation angle should be within +/-90 deg!', &
      & ' --> elevation = missing_value.'
  end if
  if ( abs(elevation_first).gt.90._double ) then
    elevation_first = -90._double
    print "(a,1x,2a)", trim(warpref), 'elevation angle should be within +/-90 deg!', &
      & ' --> elevation_first = 90.0 deg.'
  end if
  if ( abs(elevation_last).gt.90._double ) then
    elevation_last = 90._double
    print "(a,1x,2a)", trim(warpref), 'elevation angle should be within +/-90 deg!', &
      & ' --> elevation_last = 90.0 deg.'
  end if
  if (abs(elevation_step).lt.1.d-6) then
    elevation_step = 1.d-6
    print "(a,1x,2a)", trim(warpref), 'elevation step cannot be smaller than 10^-6! ', &
      & '--> elevation_step = 10^-6 deg.'
  end if
! elevation (single)
  if (abs(elevation-missing_double).gt.1.d-8) then
    elevation_range = 0._double
    elevation_first = elevation
    elevation_last  = elevation
    elevation_step  = missing_double
    elevation_count = 1_int32
  else
    elevation_range = elevation_last-elevation_first
    elevation_step  = sign(elevation_step,elevation_range)
    if (abs(elevation_step).gt.elevation_range) elevation_step=sign(elevation_range,elevation_step)
    elevation_count = elevation_range/elevation_step+1_int32
  end if
!
! azimuth
!
  if (azimuth_count.ne.missing_int32.and.azimuth_count.lt.0_int32) then
    azimuth_count = missing_int32
    print "(a,1x,2a)", trim(warpref), 'azimuth count cannot be smaller than 0! ', &
      & '--> azimuth_count = missing_value.'
  end if
  if (abs(azimuth_step).lt.1.d-6) then
    azimuth_step = 1.d-6
    print "(a,1x,2a)", trim(warpref), 'azimuth step cannot be smaller than 10^-6! ', &
      & '--> azimuth_step = 10^-6 deg.'
  end if
!
  if (abs(azimuth_first-missing_double).gt.1.d-8 &
    & .and.abs(azimuth_first-missing_double).gt.1.d-8) then
    azimuth_range = lonrangef(azimuth_last-azimuth_first)
    azimuth_first = azimuth
    azimuth_last  = azimuth
  elseif (abs(azimuth_first-missing_double).lt.1.d-8 &
    & .and.abs(azimuth_first-missing_double).lt.1.d-8) then
    ! bearing set?
    if (abs(azimuth-missing_double).lt.1.d-8) then
      azimuth=bearing_sr
    else
      call lonrange(azimuth)
    end if
    ! offset_first & offset_last
    if ( abs(azimuth_offset_first-missing_double).gt.1.d-8 &
       .and. abs(azimuth_offset_last-missing_double).gt.1.d-8 ) then
      azimuth_range = azimuth_offset_last - azimuth_offset_first
      azimuth_first = azimuth + azimuth_offset_first
      azimuth_last  = azimuth + azimuth_offset_last
    ! offset
    else if ( abs(azimuth_offset-missing_double).gt.1.d-8 ) then
      azimuth_range = 2 * azimuth_offset
      azimuth_first = azimuth - azimuth_offset
      azimuth_last  = azimuth + azimuth_offset
    ! range
    else if ( abs(azimuth_range-missing_double).gt.1.d-8 ) then
      azimuth_first = azimuth - azimuth_range/2
      azimuth_last  = azimuth + azimuth_range/2
    ! single
    else
      azimuth_range = 0._double
      azimuth_first = azimuth
      azimuth_last  = azimuth
      azimuth_step  = missing_double
      azimuth_count = 1_int32
    end if
  else
    print "(a,1x,a)", trim(errpref), 'azimuth parameters not set properly!'
  end if
!
  if (abs(azimuth_range).gt.1.d-8) then
    if (azimuth_count.eq.missing_int32) then
      azimuth_step  = sign(azimuth_step,azimuth_range)
      azimuth_count = azimuth_range/azimuth_step+1_int32
    else
      azimuth_step = azimuth_range/(azimuth_count-1_int32)
    end if
  end if
!
  rays_count=azimuth_count*elevation_count
!
! transmission loss
!
  if ( tloss_frequency.lt.0.01_double ) then
    tloss_frequency = 0.01_double
    print "(a,1x,2a)", trim(warpref), 'transmission loss frequency too low!', &
      & '--> tloss_frequency = 0.01Hz.'
  end if
!
  if (tloss_re_unit.eq.'km') then
    re_unit_conversion=1._double
  else
    tloss_re_unit='m'
    re_unit_conversion = 1000._double
  end if
!
! solver options - rk stepsize
!
  if ( rk_stepsize.lt.0.1_double ) then
    rk_stepsize = 1.0_double
    print "(a,1x,2a)", trim(warpref), 'solver stepsize cannot be smaller than 0.1s! ', &
      & '--> rk_stepsize = 1s'
  end if
!
! output
!
  if ( clip_max_offset.lt.1.0_double ) then
    clip_max_offset = 1.0_double
    print "(a,1x,2a)", trim(warpref), 'clip maximum offset cannot be smaller than 1km! ', &
      & '--> clip_max_offset = 1km.'
  end if
!
  if (len_trim(filename).eq.0) then
    filename = 'eigenrays_%T_%Y%M%D_%H'
  else
    prefix=''
    suffix=''
  end if
  filename = irt_substitute_filename( string=trim(filename), e=atmo%ranges%epoch,&
    c=atmo%ranges%centre, t=atmo%ranges%dataType, s=atmo%step, n=atmo%ens )
!
! File already exists?
!
  if (.not.overwrite) then
    inquire(file=trim(prefix)//trim(filename)//trim(suffix)//'.nc',exist=exists)
    if (exists) then
      print "(a,1x,3a)", trim(errpref), 'output file "', &
        trim(prefix)//trim(filename)//trim(suffix)//'.nc', '"already exists and overwriting is off.'
      stop
    end if
  end if
! 
