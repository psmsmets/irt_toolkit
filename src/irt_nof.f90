!*************************************************************************************************************************
!
!                                                         I R T _ N O F
!
!  INCLUDE:      irt_nof.f90
!
!  Programmer:   Pieter S. M. Smets
!                Seismology Devision - Koninklijk Nederlands Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         October 30, 2014
!
!  Language:     Fortran-90
!
!  Description:  Include file with procedures to set ranges and number of rays. 
!
!*************************************************************************************************************************


!
!.....check ranges and prepare atmo
!
      call atmo_set_ranges( trim( in_atmo ) )
!
!.....set number of rays
!
      if ( abs(elevStep).lt.1.d-8 ) then
         nofElev = 1
      else
         nofElev = int(  anint( ( elevMax - elevMin ) * 10**8 ) / &
            & anint( elevStep * 10**8 )  ) + 1
      end if
!
      if ( abs(aziStep).lt.1.d-8 ) then
         nofAzi = 1
      else
         nofAzi = int(  anint( ( aziMax - aziMin ) * 10**8 ) / &
            & anint( aziStep * 10**8 )  ) + 1
      end if
!
      nofRays = nofElev * nofAzi
!
