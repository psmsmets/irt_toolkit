!*************************************************************************************************************************
!
!                                                 I R T _ C H E C K
!
!  INCLUDE:      irt_check.f90
!
!  Programmer:   Pieter S. M. Smets
!                Seismology Devision - Koninklijk Nederlands Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         October 30, 2014
!
!  Language:     Fortran-90
!
!  Description:  Include file with input checks for irt3d.f90, irtpl.f90 and ierf3d.f90. 
!
!*************************************************************************************************************************


! --------------------------------------------------------------------------------------------------
!
! CHECK INPUT
!
! --------------------------------------------------------------------------------------------------
!
!....check length input strings
!
      if ( len_trim( in_atmo ) .ge. 256 ) then
         print "(3A)", 'ERROR : atmo profile string contains too many ', &
           & 'characters (>= 256), PROGRAM terminated. Use symbolic ', &
           & 'links to avoid too long strings.'
         print "(A)", ''
         stop
      end if
      if ( len_trim( in_topo ) .ge. 256 ) then
         print "(3A)", 'ERROR : topography string contains too many ', &
           & 'characters (>= 256), PROGRAM terminated. Use symbolic ', &
           & 'links to avoid too long strings.'
         print "(A)", ''
         stop
      end if
      if ( len_trim( out_prefix ) + 10 .ge. 256 ) then
         print "(2A)", 'ERROR : output prefix contains too many ', &
            & 'characters (>= 256), PROGRAM terminated.'
         print "(A)", ''
         stop
      end if
!     
!....altitude
!
      if ( atmo%hmin .lt. 0.d0 ) then
         atmo%hmin = 0.d0
         print "(2A)", 'WARNING : minimal altitude cannot be smaller ', &
            & 'than zero! --> hMin = 0.0 km.'
         print "(A)", ''
      end if
      if ( atmo%hmax .lt. atmo%hmin ) then
         atmo%hmax = 60.d0
         print "(2A)", 'WARNING : maximal altitude cannot be smaller ', &
            & 'than minimal altitude! --> hMax = 60.0 km.'
         print "(A)", ''
      end if
      if ( atmo%hmax .gt. 80.d0 ) then
         atmo%hmax = 80.d0
         print "(2A)", 'WARNING : initial altitude cannot be larger ', &
            & 'than 80km! --> hmax=80km.'
         print "(A)", ''
      end if
      if ( h0 .lt. atmo%hmin ) then
         h0 = atmo%hmin
         print "(2A)", 'WARNING : initial altitude cannot be smaller ', &
            & 'than minimal altitude! --> h0 = hmin.'
         print "(A)", ''
      end if
      if ( h0 .gt. atmo%hmax ) then
         h0 = atmo%hmax
         print "(2A)", 'WARNING : initial altitude cannot be larger ', &
            & 'than maximal altitude! --> h0 = hmax.'
         print "(A)", ''
      end if
      if ( atmo%dh.lt.1.d-8 ) then
         atmo%dh = .5d0
         print "(2A)", 'WARNING : altitude stepsize cannot be zero! ', &
            & '--> dh = 0.5 km.'
         print "(A)", ''
      end if
      if ( atmo%ens .lt. 0 .or. atmo%step .lt. 0 ) then
         atmo%ens  = 0
         atmo%step = 0
         print "(2A)", 'WARNING : ensemble number or forecast step ' // &
            & ' cannot be smaller than zero! --> 0'
         print "(A)", ''
      end if
!
!.....elevation
!
      if ( elevMin .lt. -90.d0 .and. h0.gt.0.d0 ) then
         elevMin = -90.d0
         print "(2A)", 'WARNING : minumum elevation cannot be smaller ', &
            & 'than -90.0 deg when h0 > 0 ! --> elevMin = -90.0 deg.'
         print "(A)", ''
      end if
      if ( elevMin .lt. 0.d0 .and. h0.lt.1.d-8 ) then
         elevMin = 0.d0
         print "(2A)", 'WARNING : minumum elevation cannot be smaller ', &
            & 'than 0.0 deg when h0 = 0 ! --> elevMin = 0.0 deg.'
         print "(A)", ''
      end if
      if ( elevMax .gt. 90.d0 ) then
         elevMax = 90.d0
         print "(2A)", 'WARNING : maximum elevation cannot be larger ', &
            & 'than 90.0 deg! --> elevStep = 90.0 deg.'
         print "(A)", ''
      end if
      if ( abs(elevStep).lt.1.d-8 .and. abs(elevMin-elevMax).gt.1.d-8 ) then
         elevStep = 1.d0
         print "(2A)", 'WARNING : elevation stepsize cannot be zero! ', &
            & '--> elevStep = 1.0 deg.'
         print "(A)", ''
      end if
!
!.....azimuth
!
      if ( abs(aziStep).lt.1.d-8 .and. abs(aziMin-aziMax).gt.1.d-8 ) then
         aziStep = 1.d0
         print "(2A)", 'WARNING : azimuth stepsize cannot be zero! ', &
            & '--> aziStep = 1.0 deg.'
         print "(A)", ''
      end if
!
!.....rk stepsize
!
      if ( rk_stepsize.lt.1.d-8 ) then
         rk_stepsize = 0.1d0
         print "(2A)", 'WARNING : rk4 stepsize cannot be zero! ', &
            & '--> rk_stepsize = 0.1 km.'
         print "(A)", ''
      end if
!
!.....out_type
!
      if (out_type .ne. 1 .and. out_type .ne. 2 .and. out_type .ne. 3) then
         print "(2A)", 'ERROR : ray data output type should be ', &
            & 'equal to 1, 2, or 3 (1=reflections,2=rays,3=both).'
         print "(A)", ''
         stop
      end if
!
!.....Stop conditions
!
      if ( sc%d .le. 0.d0 ) then
         sc%d = 1000.d0
         print "(2A)", 'WARNING : stop distance should be > 0 km!', &
            & '--> sc%d = 1000 km.'
         print "(A)", ''
      end if
!
      if ( sc%d .gt. 40000.d0 ) then
         sc%d = 40000.d0
         print "(2A)", 'WARNING : max stop distance is 40,000 km!', &
            & '--> sc%d = 40,000 km.'
         print "(A)", ''
      end if
!     
