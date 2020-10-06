!*****************************************************************************80
!
!                   I E R F 3 D _ P A R S E _ C O M M A N D _ A R G S
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
!  Description:  Parse ierf3d command line arguments
!
!
!*****************************************************************************80
!
! init
!
  arg_cnt = 0_int32
  arg_parsed = .false.
!
! parse arguments
!
  do while ( arg_cnt.lt.command_argument_count() )
!
!   read argument type and change to lower case
!
    arg_cnt   = arg_cnt + 1_int32
    arg_value = ''
    call get_command_argument( arg_cnt, arg_type )
!
!   overwrite default values
!
    select case ( strlowcase(arg_type) )
!
!.....debug options
!
      case ('--debug') ! verbose all
        debug=.true.
!
      case ('--version')
        print "(a)", trim(irt_toolkit_version()) 
        stop
!
!.....help
!
      case ('-?','--help')
        call ierf3d_show_help( trim(irt_toolkit_version()) )
!
!.....atmospheric conditions
!
      case ('-a','--atmosphere')
        arg_type='-a'
        arg_cnt=arg_cnt+1; call get_command_argument( arg_cnt, arg_value )
        atmo_path=arg_value
!
!.....parameter arguments
!
      case ('-i','--input')
        call ierf3d_show_arguments( trim(irt_toolkit_version()) )
!
!.....overwrite netcdf output files
!
      case ('-o','--overwrite')
        arg_type='-o'
        overwrite = .true.
!
!.....output prefix
!
      case ('-p','--prefix')
        arg_type='-p'
        arg_cnt=arg_cnt+1; call get_command_argument( arg_cnt, arg_value )
        prefix=arg_value
!
!.....output suffix
!
      case ('-u','--suffix')
        arg_type='-u'
        arg_cnt=arg_cnt+1; call get_command_argument( arg_cnt, arg_value )
        suffix=arg_value
!
!.....output filename
!
      case ('-f','--filename')
        arg_type='-f'
        arg_cnt=arg_cnt+1; call get_command_argument( arg_cnt, arg_value )
        filename=arg_value
!
!.....verbose and level
!
      case ('-v','--verbose')
        arg_type='-v'
        verb = .true.
        call get_command_argument( arg_cnt + 1_int32, arg_value )
        if ( arg_value(1:1).eq.'-' .or. len_trim(arg_value).eq.0 ) then
          verblvl = 1_int8
          arg_value=''
        else
          arg_cnt=arg_cnt+1_int32
          read ( arg_value ,"(i2)", iostat=error ) verblvl
          if ( error.ne.0_int8 .or. verblvl.lt.0_int8 ) verblvl = 1_int8
          if ( verblvl.gt.2_int8 ) verblvl=2_int8
          if ( verblvl.eq.0_int8 ) verb=.false.
        end if
      case ('-n','--silent')
        arg_type='-n'
        verb = .false.
!
!.....openmp threads
!
      case ('-t','--threads')
        arg_type='-t'
        parallel = .true.
        call get_command_argument( arg_cnt + 1_int32, arg_value )
        if ( arg_value(1:1).eq.'-' .or. len_trim(arg_value).eq.0 ) then
          threads = max_threads
          arg_value=''
        else
          arg_cnt=arg_cnt+1_int32
          read ( arg_value ,"(i2)", iostat=error ) threads
          if ( error.ne.0 .or. threads.lt.-1 .or. threads.eq.0 .or. threads.gt.max_threads ) then
            threads = max_threads
            print "(3a,i2,a)", 'WARNING: defined threads for',&
              ' openmp (-t) is out of range.', &
              ' Number of threads is put to max_threads = ', &
              threads, '!'
          end if
        end if
!
!.....otherwise put in input file/folder list
!
      case default
!
        if (arg_cnt.eq.1) then
          include 'ierf3d_parse_file_args.f90'
        end if
!
    end select
!
  end do
!
! check
!
  if ( .not.arg_parsed ) then
    print "(a,1x,a,i3,a)", trim(errpref), 'no input argument file provided! See ierf3d --help for usage.'
    stop
  end if
