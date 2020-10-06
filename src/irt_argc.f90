!*************************************************************************************************************************
!
!                                                 I R T _ A R G C
!
!  INCLUDE:      irt_argc.f90
!
!  Programmer:   Pieter S. M. Smets
!                Seismology Devision - Koninklijk Nederlands Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         August 29, 2014
!
!  Language:     Fortran-90
!
!  Description:  Include file with variable input procedures. 
!
!*************************************************************************************************************************

!
!.....get file_input from command-line argument 
!
      call getarg ( 1, in_inputfile )
!
!.....check length inputfile
!
      if ( len_trim( in_inputfile ) .ge. 256) then
         print "(2A)", 'Error : input file string contains too many', &
            & ' characters (>= 256), PROGRAM terminated.'
         stop
      end if
!
!.....open tmp input file
!
      inquire ( file = in_inputfile, exist = lexist ) 
!
      if ( lexist ) then
         open ( unit = 1, file = in_inputfile, &
            & status = 'old', action = 'read' )
      else
         print "(3A)", "Error : input file '", trim( in_inputfile ), &
            &  "' not found, PROGRAM terminated."
         stop ! stop PROGRAM due to lack of input profile
      end if
!
!.....read inputfile arguments 
!
!.....read in ECMWF profile filename
!
      read ( 1, "(a)", iostat = io ) in_atmo
      if ( io .ne. 0 ) stop 'Error: No ECMWF profile provided!'
!
!.....read in altitude parameters
!
      read ( 1, *, iostat = io ) atmo%hmin, atmo%hmax, atmo%dh, atmo%ens
      atmo%step = atmo%ens
      if ( io .ne. 0 ) stop 'Error: Not enough altitude ' // &
         & 'parameters listed! Provide hmin, hmax, dh (km), ' // &
         & 'and ensemble number (-).'
!
!.....read in reverse horizontal wind (yes/no)
!
      read ( 1, "(A3)", iostat = io ) in_reverse 
      if ( io .ne. 0 ) stop 'Error: No reverse parameter provided' // &
         & ' (yes/no)!'
      in_reverse = StrLowCase ( in_reverse )
!
!.....read in topography filename
!
      read ( 1, "(a)", iostat = io ) in_topo
      if ( io .ne. 0 ) stop 'Error: No topography option provided ' // &
         & '(yes/no/path)!'
!
!.....read in source parameters
!
      read ( 1, *, iostat = io ) lat0, lon0, h0
      if ( io .ne. 0 ) stop 'Error: Not correct source parameters' // &
         & ' listed! Provide lat0 and lon0 (deg), and h0 (km)'
!
      call LONRANGE ( lon0 )
!
!.....read in elevation parameters
!
      elevStep = 0.d0
!
      read ( 1, "(a)" ) dummy 
      i = scan( dummy, ' ', .FALSE. )
      read ( dummy( 1 : i ), *, iostat = io ) elevMin
      if ( io .ne. 0 ) stop 'Error: No elevation angle provided!'
      dummy = adjustl(  dummy( i + 1 : )  )
!
      i = scan( dummy, ' ', .FALSE. )
      read ( dummy( 1 : i ), *, iostat = io ) elevMax
!
      if ( io .eq. 0 ) then
         dummy = adjustl(  dummy( i + 1 : )  )
         i     = scan( dummy, ' ', .FALSE. )
         read ( dummy( 1 : i ), *, iostat = io ) elevStep
         if ( io.ne.0 .and. abs(elevMax).lt.1.d-8 ) stop &
            & 'Error: No elevation stepsize provided!'
      else
         elevMax = elevMin
      end if
!
!.....read in azimuth parameters
!
      aziStep = 0.d0
!
      read ( 1, "(a)" ) dummy 
      i = scan( dummy, ' ', .FALSE. )
      read ( dummy( 1 : i ), *, iostat = io ) aziMin
      if ( io .ne. 0 ) stop 'Error: No azimuth angle provided!'
      dummy = adjustl(  dummy( i + 1 : )  )
!
      i = scan( dummy, ' ', .FALSE. )
      read ( dummy( 1 : i ), *, iostat = io ) aziMax
!
      if ( io .eq. 0 ) then
         dummy = adjustl(  dummy( i + 1 : )  )
         i     = scan( dummy, ' ', .FALSE. )
         read ( dummy( 1 : i ), *, iostat = io ) aziStep
         if ( io .ne. 0 .and. abs(aziMax).lt.1.d-8 ) stop &
            & 'Error: No azimuth stepsize provided!'
      else
         aziMax = aziMin
      end if
!
!.....transmission loss frequency
!
      read ( 1, *, iostat = io ) frq
      if ( io .ne. 0 ) stop 'Error: No transmission loss frequency' // &
         & ' provided (Hz)!'
      if ( frq .LT. 0.002d0 ) frq = 0.002d0
!
!.....transmission loss unit
!
      read ( 1, "(a)", iostat = io ) dummy
      dummy = strlowcase ( dummy )
      if ( io .ne. 0 ) stop 'Error: No transmission loss re unit ' // &
         & 'provided (m or km)!'
      re_unit = dummy(1:2)
      re_unit_conversion = 1.d0
      if (re_unit.eq.'m') re_unit_conversion=1.d3
!
!.....read in integration stepsize (km^2/s)
!
      read ( 1, *, iostat = io ) rk_stepsize 
      if ( io .ne. 0 ) stop 'Error: No rk_stepsize provided!'
!
!.....read in output prefix
!
      read ( 1, *, iostat = io ) out_prefix
      if ( io .ne. 0 ) stop 'Error: No output prefix provided!'
!
!.....read in output type (1=refl,2=rays,3=both)
!
      read ( 1, "(i1)", iostat = io ) out_type
      if ( io .ne. 0 .OR. out_type .GT. 3 .OR. out_type .LT. 1 ) &
         & stop 'Error: No output type provided. Should be 1, 2, or 3!'
!
!.....read in output vertical cross-section
!
      read ( 1, "(a)", iostat = io ) cross_ver
      if ( io .ne. 0 ) stop 'Error: No output prefix provided!'
!
!.....read in verbose
!
      read ( 1, "(a)", iostat = io ) verb
      if ( io .ne. 0 ) stop 'Error: No verbose parameter provided!'
      verb = StrLowCase ( verb )
!
!.....stop condition
!
      read ( 1, "(a)", iostat = io ) dummy
      if ( io .ne. 0 ) stop 'Error: No stop condition provided!'
      sc%stype = strlowcase ( dummy( 1 : 2 ) )
      dummy    = adjustl(  dummy ( 3 : )  )
!
      sc%latmin =  -90.d0
      sc%latmax =   90.d0
      sc%lonmin =    0.d0
      sc%lonmax =  359.99d0
      sc%d      = 1000.d0
!
      select case ( sc%stype )
         case ( 'll' )
            read ( dummy, *, iostat = io ) sc%latmin, sc%latmax, sc%lonmin, sc%lonmax
         case ( 'ld' )
            read ( dummy, *, iostat = io ) sc%latmin, sc%latmax, sc%d
         case ( 'dl' )
            read ( dummy, *, iostat = io ) sc%d, sc%lonmin, sc%lonmax
         case ( 'd ' )
            read ( dummy, *, iostat = io ) sc%d
         case default
            sc%stype = 'd '
      end select
!
      call LONRANGE ( sc%lonmin )
      call LONRANGE ( sc%lonmax )

      if ( io .ne. 0 ) then
         print "(2a)", 'WARNING: illegal stop conditions provided, ', &
            & 'default settings will be used!'
         sc%stype = 'd '
         sc%latmin =  -90.d0
         sc%latmax =   90.d0
         sc%lonmin =    0.d0
         sc%lonmax =  360.d0
         sc%d      = 1000.d0
      end if
!
!.....read markers
!
      include 'irt_marker.f90'
!
!.....close input file
!
      close ( 1 )
!

