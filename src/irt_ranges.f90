!*************************************************************************************
!
!                                   I R T _ R A N G E S
!
!  INCLUDE:      irt_ranges.f90
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
!***************************************************************************************

!
!.....print header
!
      if ( verbose ) print "(A)", 'SET RANGES'
!
!.....Set starting altitude
!
      if ( ltopo ) then
        call topo_get ( lon0, lat0, h )
        if ( h.gt.h0 ) h0=h
      end if
!
!.....Check stop conditions with atmospheric profile
!
      sc%hmin = atmo%hmin
      sc%hmax = atmo%hmax
!
      if ( sc%latmin .lt. atmo%latmin ) sc%latmin = atmo%latmin
      if ( sc%latmax .gt. atmo%latmax ) sc%latmax = atmo%latmax
      if ( sc%lonmin .lt. atmo%lonmin ) sc%lonmin = atmo%lonmin
      if ( sc%lonmax .gt. atmo%lonmax ) sc%lonmax = atmo%lonmax
!
!.....Correct logicals based on type
!
      select case ( sc%stype )
         case ( 'll' )
            sc%d        = 40000.d0
         case ( 'ld' )
         case ( 'dl' )
            if ( atmo%l360 ) THEN
               sc%lonmin =    0.d0
               sc%lonmax =  360.d0
            end if
         case ( 'd ' )
            if ( atmo%l360 ) THEN
               sc%lonmin =    0.d0
               sc%lonmax =  360.d0
            end if
      end select
!
!.....Get longitude and latitude range
!
      sc%lonrange = atmo%dlon * ( atmo%noflon - 1 )
      sc%latrange = atmo%dlat * ( atmo%noflat - 1 )
!
!.....Print stop conditons
!
      if ( verbose ) THEN
         print "(3x,2A)", "  o Stop range type          : ", sc%stype
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
            & "  o Distance range      (km) : ", sc%d 
         print "(A)", ''
      end if
!
!
!

!
