!*************************************************************************************************************************
!
!                                               I R T _ P A R A M S
!
!  INCLUDE:      irt_params.f90
!
!  Programmer:   Pieter S. M. Smets
!                Seismology Devision - Koninklijk Nederlands Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         October 30, 2014
!
!  Language:     Fortran-90
!
!  Description:  Include file with parameter print routine. 
!
!*************************************************************************************************************************

      if ( verbose ) then
!
         print "(A)",         'GENERAL PARAMETERS'
         print "(3x,A)",      "ECMWF atmospheric profile"
         print "(3x,2A)",     "  o File name                : ", &
            & trim( in_atmo )
         print "(3x,A,I4)",   "  o Levels                   : ", &
            & atmo%ranges%noflvl
         print "(3x,2(A,F6.2),A,F5.2,A,I4,A)", &
            & "  o Latitude           (deg) : ", &
            & atmo%latMin, " to ", atmo%latMax, ", ", &
            & atmo%dlat, " (#", atmo%noflat, ")"
         if ( atmo%l360 ) &
            & print "(3x,2(A,F6.2),A,F5.2,A,I4,A)", &
            & "  o Longitude          (deg) : ", 0.0d0, " to ", 360.0d0, &
            & ", ", atmo%ranges%dlon, " (#", atmo%ranges%noflon, ")"
         if ( .not. atmo%l360 ) &
            & print "(3x,2(A,F6.2),A,F5.2,A,I4,A)", &
            & "  o Longitude          (deg) : ", atmo%lonMin, &
            & " to ", atmo%lonMax, ", ", atmo%dlon, &
            & " (#", atmo%noflon, ")"
         print "(3x,2(A,F6.2),A,F5.2,A,I4,A)", &
            & "  o Altitude            (km) : ", atmo%hMin, &
            & " to ", atmo%hmax, ", ", atmo%dH, &
            & " (#", atmo%nofh, ")"
         print "(3x,2A)",     "  o Reverse horizontal wind  : ", &
            & trim( in_reverse )
         print "(3x,A)",      ""
         if ( ltopo ) then
            print "(3x,2A)", "Topography                   : ", trim( in_topo )
         else
            print "(3x,A)" , "Topography                   : no"
         end if
         print "(3x,A)",      ""
         print "(3x,A)",      "Source location"
         print "(3x,A,F8.4)", "  o Latitude           (deg) : ", lat0
         print "(3x,A,F8.4)", "  o Longitude          (deg) : ", lon0
         print "(3x,A,F8.4)", "  o Altitude            (km) : ", h0
         print "(3x,A)",      ""
         print "(3x,A)",      "Ray trace options"
         print "(3x,2(A,F7.3),A,F6.3,A,I4,A)", &
            & "  o Elevation          (deg) : ", elevMin, " to ", &
            & elevMax, ", ", elevStep, " (#", nofElev, ")"
         print "(3x,2(A,F7.3),A,F6.3,A,I4,A)", &
            & "  o Azimuth            (deg) : ", aziMin,   " to ", &
            & aziMax, ", ", aziStep, " (#", nofAzi, ")"
         print "(3x,A,F7.3)", "  o TLOSS frequency     (Hz) : ", frq
         select case ( re_unit )
            case ('km')
               print "(3x,A)","  o TLOSS unit               : dB Re 1km"
            case ('m')
               print "(3x,A)","  o TLOSS unit               : dB Re 1m" 
            case default
               print "(3x,A)","  o TLOSS unit               : normalize"
         end select
         print "(3x,A,F7.3)", "  o RK stepsize     (s) : ", &
            & rk_stepsize
         if ( jacobian ) then
            print "(3x,A)", "  o Jacobian                 : yes"
         else
            print "(3x,A)", "  o Jacobian                 : no"
         end if
         print "(3x,A)",      ""
         print "(3x,A)",      "Output"
         print "(3x,2A)",     "  o Output prefix            : ", &
            & trim( out_prefix )
         print "(3x,A,I4)", "  o Ray data output type     : ", out_type
         print "(3x,A)", ''
!
      end if
!
