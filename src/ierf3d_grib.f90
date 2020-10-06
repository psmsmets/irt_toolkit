!*****************************************************************************80
!
!                                    I R T _ G R I B 
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
!  Description:  Include file with atmo and topo grib calls.
!
!
!*****************************************************************************80
!
! Read atmospheric profile (grib)
!
  if (verb) print "('> ',a)", 'Get atmospheric conditions'
! Allocate atmo arrays
  if (verblvl.gt.1) print "(4x,a)" , 'Allocate raw atmospheric grid'
  call atmo_allocate ( 'raw' )
! Set atmospheric profile
  if (verblvl.gt.1) print "(4x,a)" , 'Read raw atmosphere grid'
  call atmo_read_raw ( trim(atmo_path), ltopo )
! Allocate atmo grd arrays
  if (verblvl.gt.1) print "(4x,a)" , 'Allocate regular atmospheric grid'
  call atmo_allocate ( 'grd' )
!
! Get topography
!
! check topography file, otherwise use ECMWF height of level 1
!
  if ( ltopo ) then
!
    if (verb) print "('> ',a)", 'Get topography'
!
    if (topo_from_atmo) then
      if (verblvl.gt.1) print "(4x,a)" , &
        'Get topography from ecmwf grib field'
      if ( atmo%ranges%ltopo ) then
        allocate ( topo%h( 0:atmo%ranges%noflon+1, 0:atmo%ranges%noflat+1 ) )
        topo%h=0.d0
        topo%ranges%dlat   = atmo%ranges%dlat
        topo%ranges%latmin = atmo%ranges%latmin - topo%ranges%dlat
        topo%ranges%latmax = atmo%ranges%latmax + topo%ranges%dlat
        topo%ranges%noflat = atmo%ranges%noflat 
        topo%ranges%dlon   = atmo%ranges%dlon
        topo%ranges%lonmin = atmo%ranges%lonmin - topo%ranges%dlon
        topo%ranges%lonmax = atmo%ranges%lonmax + topo%ranges%dlon
        topo%ranges%noflon = atmo%ranges%noflon
        topo%ranges%nofpoints = topo%ranges%noflat * topo%ranges%noflon
        topo%circ = (/ atmo%l360, .FALSE. /)
        topo%h( 1:topo%ranges%noflon, 1:topo%ranges%noflat ) = atmo%raw%topo
        call CW_SET_BOUND2 ( topo%ranges%noflon, topo%ranges%noflat, topo%h, topo%circ )
      else
        ltopo = .false.
        print *, 'Warning: no topography present in grib file!'
      end if
!            
    else if (topo_from_dem) then
!      call topo_set_new (  atmo%ranges%latmin, atmo%ranges%latmax, &
!        & atmo%ranges%lonmin, atmo%ranges%lonmax, frq, atmo%l360  )
      call topo_read (  trim(topo_path), atmo%ranges%latmin, atmo%ranges%latmax, &
        & atmo%ranges%lonmin, atmo%ranges%lonmax, atmo%l360  )
    end if
!
!   Update grib missing topography?
!
    if ( ltopo .and. .not. atmo%ranges%ltopo ) then
      do i = 1, atmo%ranges%noflon
        do j = 1, atmo%ranges%noflat
          call topo_get ( &
            & atmo%ranges%lonmin + ( i - 1 ) * atmo%ranges%dlon, &
            & atmo%ranges%latmin + ( j - 1 ) * atmo%ranges%dlat, h )
          atmo%raw%h( i, j, : ) = atmo%raw%h( i, j, : ) + h
          atmo%raw%topo( i, j ) = h
        end do
      end do
    end if
!
  end if
!
! Grid atmospheric profile
!
  if (verb) print "('> ',a)", 'Regrid atmospheric conditions'
!
! Regularize raw atmo data to grd
!
  if (verblvl.gt.1) print "(4x,a)" , 'Regularize raw atmospheric grid'
  call atmo_raw2grd ()
! 
! Deallocate atmo raw arrays
!
  if (verblvl.gt.1) print "(4x,a)" , 'Deallocate raw atmospheric grid'
  call atmo_deallocate ( 'raw' )
!
