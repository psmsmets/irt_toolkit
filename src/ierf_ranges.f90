!*************************************************************************************************************************
!
!                                                I E R F _ R A N G E S
!
!  INCLUDE:      ierf_ranges.f90
!
!  Programmer:   Pieter S. M. Smets
!                Seismology Devision - Koninklijk Nederlands Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         August 29, 2014
!
!  Language:     Fortran-90
!
!  Description:  Include file with set range procedures. 
!
!*************************************************************************************************************************

!
!.....print header
!
      if ( verbose ) print "(A)", 'SET RANGES'
!
!.....Check stop conditions with atmospheric profile
!
      sc%hmin = atmo%hmin
      sc%hmax = atmo%hmax
!
!.....Copy atmo ranges to stop conditions
!
      sc%latmin = atmo%latmin
      sc%latmax = atmo%latmax
      sc%lonmin = atmo%lonmin
      sc%lonmax = atmo%lonmax
!
!.....Set starting altitude
!
      if ( ltopo ) then
        call topo_get ( lon0, lat0, h, hdLon, hdLat )
        if ( h.gt.h0 ) h0=h
        call topo_get ( lon1, lat1, h )
        if ( h.gt.h1 ) h1=h
      end if
      call llh2sph( lat1, lon1, h1, ph1, th1, r1, atmo%lonshift )
!
!.....Calculate Spherical and Euclidean distance between source and poi 
!
      call distance( lat0, lon0, lat1, lon1, d )
      call euclidean_distance_llh( lat0, lon0, h0, lat1, lon1, h1, ud )
!
      sc%d = ud
!
!.....Get longitude and latitude range
!
      sc%lonrange = atmo%dlon * ( atmo%noflon - 1 )
      sc%latrange = atmo%dlat * ( atmo%noflat - 1 )
!
!.....Print stop conditons
!
      if ( verbose ) THEN
         print "(3x,2(A,F6.2))", &
            & "  o Latitude range     (deg) : ", sc%latmin, " to ", sc%latmax
         if ( atmo%l360 ) THEN
            print "(3x,2(A,F6.2))", &
               & "  o Longitude range    (deg) : ", 0.d0, " to ", 360.d0
         else
            print "(3x,2(A,F6.2))", &
               & "  o Longitude range    (deg) : ", sc%lonmin, " to ", sc%lonmax
         end if
         print "(3x,1(A,F9.2))", &
            & "  o Euclidean distance  (km) : ", sc%d 
         print "(A)", ''
      end if
!
