!*************************************************************************************************************************
!
!                                                    I R T _ G R I B
!
!  INCLUDE:      irt_grib.f90
!
!  Programmer:   Pieter S. M. Smets
!                Seismology Devision - Koninklijk Nederlands Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         October 30, 2014
!
!  Language:     Fortran-90
!
!  Description:  Include file with atmo and topo grib calls. 
!
!*************************************************************************************************************************

! --------------------------------------------------------------------------------------------------
!
! READ ATMOSPHERIC PROFILE
!
! --------------------------------------------------------------------------------------------------

      if ( verbose ) print "(A)", 'READ ATMOSPHERIC PROFILE'
! 
!.....Allocate atmo arrays
!
      call atmo_allocate ( 'raw' )
!
!.....Set atmospheric profile
!
      call atmo_read_raw ( trim( in_atmo ), ltopo )
! 
!.....Allocate atmo grd arrays
!
      call atmo_allocate ( 'grd' )

! --------------------------------------------------------------------------------------------------
!
! GET TOPOGRAPHY
!
! --------------------------------------------------------------------------------------------------

!
!.....check topography file, otherwise use ECMWF height of level 1
!
      if ( ltopo ) then
!
         if ( verbose ) print "(A)", 'SET TOPOGRAPHY'
!
         dummy = StrLowCase (  trim( in_topo )  )
!
         if ( dummy .eq. 'a' .or. dummy .eq. 'e' .or. &
            & dummy .eq. 'atmo' .or. dummy .eq. 'ecmwf'  ) then
!
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
         else if ( dummy .eq. 'y' .or. dummy .eq. 'yes' .AND. &
            & .NOT. ISNETCDF( trim( in_topo ) ) ) then
!
            call topo_set_new (  atmo%ranges%latmin, atmo%ranges%latmax, &
               & atmo%ranges%lonmin, atmo%ranges%lonmax, frq, atmo%l360  )
!
         else
!
            call topo_read (  trim( in_topo ), atmo%ranges%latmin, atmo%ranges%latmax, &
               & atmo%ranges%lonmin, atmo%ranges%lonmax, atmo%l360  )
!
         end if
!
!........update grib missing topography?
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

! --------------------------------------------------------------------------------------------------
!
! GRID ATMOSPHERIC PROFILE
!
! --------------------------------------------------------------------------------------------------

      if ( verbose ) print "(A)", 'GRID ATMOSPHERIC PROFILE'

!
!.....Regularize raw atmo data to grd
!
      call atmo_raw2grd ()
! 
!.....Deallocate atmo raw arrays
!
      call atmo_deallocate ( 'raw' )
!
      if ( verbose ) print "(A)", ''
!
