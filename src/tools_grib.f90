! =========================================================================================
!
! MODULE GRIBTOOLS
!
! --> usefull grib programs.
!
! =========================================================================================
!
Module tools_grib
!
   use grib_api
   implicit none
!
!.....types
!
      type type_gribtools_ranges
         logical           :: elda, ow, ltopo, isOctahedral, isECMWF
         integer(kind=4)   :: noflon, noflat, noflvl, nofstp, nofvar, nofpoints, nofadded, noffrq, nofdir, &
                              nofdims, nofens
         integer(kind=8)   :: epoch
         character(len=23) :: tstr
         character(len=20) :: centre, editionnumber, stream, gridType, typeoflevel, datatype, &
                              gridName, typeOfProcessedData
         double precision  :: latmin, latmax, lonmin, lonmax, dlon, dlat, dirscale, frqscale
         logical          , dimension(  3)  :: dims
         integer(kind=4)  , dimension(150)  :: lvllist, stplist, dirlist, frqlist, enslist
         character(len=20), dimension( 20)  :: varlist, varadded
      end type type_gribtools_ranges

Contains

! =========================================================================================
!
!.......function ISGRIB 
!
!       Check if file is of grib type
!
!****************************
   Logical function isgrib ( file )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character(len=*), intent(in) :: file
!
!.....Local variables
!
      integer(kind=4)  :: ifile, ierr
!
!  ---
!
!.....Open the grib file.
!
      CALL grib_open_file ( ifile, file, 'R', ierr )
!
!.....If no error, file is grib
!
      isgrib = ierr.eq.0
!
!.....If open succesfull, close file again
!
      if (isgrib) call grib_close_file ( ifile )
!
      return
!
   End function isgrib
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GRIBTOOLS_RANGES
!
!       Get ranges from grib file.
!
!****************************
   SUBROUTINE GRIBTOOLS_GET_TIME ( file1, epoch, tstr )
!****************************
!
      USE TIME
      implicit none
!
!.....Dummy variables
!
      character(len=*) , intent(in)             :: file1
      integer(kind=8)  , intent(out)            :: epoch
      character(len=23), intent(out), optional  :: tstr
!
!.....Local variables
!
      integer(kind=4)                               :: ifile, igrib, ierr, d, t
      character                                 :: d_str * 8, t_str * 4
!
!  ---
!
!.....Open file.
!
      CALL grib_open_file ( ifile, file1, 'R', ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_RANGES grib_open_file'
!
      CALL grib_new_from_file ( ifile, igrib, ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_RANGES grib_new_from_file'
!
!.....Get date and time from grib file.
!
      CALL grib_get( igrib, 'dataDate', d, ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_RANGES grib_get(dataDate)'
      CALL grib_get( igrib, 'dataTime', t, ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_RANGES grib_get(dataTime)'
!
      WRITE ( d_str, "(I8)" ) d
      WRITE ( t_str, "(I4.4)" ) t
!
      CALL STR2EPOCH ( d_str( 1 : 4 ) // '/' // d_str( 5 : 6 ) // '/' // d_str( 7 : ) // ' ' // &
         & t_str( 1 : 2 ) // ':' // t_str( 3 : 4 ), epoch )
      IF ( PRESENT ( tstr ) ) CALL EPOCH2STR ( epoch, tstr )
!
!.....Close file.
!
      CALL grib_release ( igrib )
!
      CALL grib_close_file ( ifile )
!
      RETURN
!
   END SUBROUTINE GRIBTOOLS_GET_TIME
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GRIBTOOLS_RANGES
!
!       Get ranges from grib file.
!
!****************************
   SUBROUTINE GRIBTOOLS_RANGES ( file1, ranges )
!****************************
!
      USE TIME
      implicit none
!
!.....Dummy variables
!
      character(len=*), intent(in)                 :: file1
      type(type_gribtools_ranges), intent(out)      :: ranges
!
!.....Local variables
!
      logical                     :: lexist
      integer(kind=4)             :: ifile, igrib, ierr, idx, d, t
      character                   :: dstr*8, tstr*4
!
!  ---
!
!.....Initialize
!
      ifile              = 50
      igrib              = 51
      ierr               = 0
!
      ranges%elda        = .FALSE.
      ranges%ow          = .FALSE.
!
      ranges%epoch       = 0
      ranges%tstr        = ''
      ranges%centre      = ''
      ranges%dataType    = ''
      ranges%typeoflevel = ''
      ranges%gridtype    = ''
      ranges%nofvar      = 0
      ranges%varlist     = ''
      ranges%nofadded    = 0
      ranges%varadded    = ''
      ranges%noflvl      = 0
      ranges%nofstp      = 0
      ranges%nofens      = 0
      ranges%latMin      = 0.d0
      ranges%latMax      = 0.d0
      ranges%dlat        = 0.d0
      ranges%noflat      = 0
      ranges%lonMin      = 0.d0
      ranges%lonMax      = 0.d0
      ranges%dlon        = 0.d0
      ranges%noflon      = 0
      ranges%nofdir      = 0
      ranges%dirscale    = 0.d0
      ranges%noffrq      = 0
      ranges%frqscale    = 0.d0
!
!.....Open profile, if exists, stop otherwise.
!
      inquire( file = file1, exist = lexist )
      if (.not.lexist) stop 'Error @ gribtools_ranges grib file does not exist'
!
!.....First get the type of level and data type
!
      call grib_open_file( ifile, file1, 'R', ierr )
      if (ierr.ne.0) stop 'Error @ gribtools_ranges grib_open_file'
!
      call grib_new_from_file( ifile, igrib, ierr )
      if (ierr.ne.0) stop 'Error @ gribtools_ranges grib_new_from_file'
!
!.....Get info
!
      call grib_get ( igrib, 'centre', ranges%centre, ierr )
      call grib_get ( igrib, 'typeOfProcessedData', ranges%typeOfProcessedData, ierr )
!
      call grib_get ( igrib, 'dataType', ranges%dataType, ierr )
!      if (ierr.ne.0) ranges%dataType='' !! fix for harmonie data
!
      call grib_get ( igrib, 'typeOfLevel', ranges%typeOfLevel, ierr )
      if (ierr.ne.0) stop 'Error @ gribtools_ranges grib_get(typeOfLevel)'
!
      call grib_get ( igrib, 'gridType', ranges%gridType, ierr )
      if (ierr.ne.0) stop 'Error @ gribtools_ranges grib_get(gridType)'
!
      call grib_get ( igrib, 'stream', ranges%stream, ierr )
      if (ierr.ne.0) ranges%stream='oper' !! fix for harmonie data
!! fix for harmonie
if (trim(ranges%centre).eq.'99') ranges%typeOfLevel='hybrid'
!
!.....Release and close
!
      call grib_release ( igrib )
      call grib_close_file ( ifile )
!
!.....create an index from a grib file using some keys
!
      if ( ranges%typeoflevel .ne. 'meanSea' ) then
         select case ( trim(ranges%stream) )
            case ('elda')
              call grib_index_create ( idx, file1, 'shortName,level,step,number', ierr )
              ranges%elda = .true.
            case default !('oper')
              call grib_index_create ( idx, file1, 'shortName,level,step', ierr )
         end select 
      else
         call grib_index_create ( idx, file1, 'shortName,level,step,direction,frequency', ierr )
      end if
      if ( ierr .ne. 0 ) then
         print *, file1
         stop 'ERROR @ GRIBTOOLS_RANGES grib_index_create'
      end if
!
!.....get the number of distinct values of shortName in the index
!
      CALL grib_index_get_size ( idx, 'shortName', ranges%nofvar, ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_RANGES grib_index_get_size(shortName)'
!
      CALL grib_index_get ( idx, 'shortName', ranges%varlist( 1 : ranges%nofvar ), ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_RANGES grib_index_get(shortName)'
!
      ranges%nofadded = 0
      ranges%ltopo    = .true.
!
!.....add parameters if not ocean wave model
!
      IF ( ranges%typeoflevel .NE. 'meanSea' ) THEN
!
!........add z?
!
         IF (  .NOT. GRIBTOOLS_PRESENT( ranges, 'z' )  ) THEN
            ranges%ltopo                       = .false.
            ranges%nofvar                      = ranges%nofvar + 1
            ranges%varlist( ranges%nofvar )    = 'z'
            ranges%nofadded                    = ranges%nofadded + 1
            ranges%varadded( ranges%nofadded ) = 'z'
         END IF
!
!........add pressure
!
         IF (  .NOT. GRIBTOOLS_PRESENT( ranges, 'p' )  ) THEN
            ranges%nofvar                      = ranges%nofvar + 1
            ranges%varlist(ranges%nofvar)      = 'p'
            ranges%nofadded                    = ranges%nofadded + 1
            ranges%varadded( ranges%nofadded ) = 'p'
         END IF
!
      ELSE
!
         ranges%ow = .TRUE.
!
      END IF
!   
!.....get the number of distinct values of level in the index
!
      CALL grib_index_get_size( idx, 'level', ranges%noflvl )
!
      CALL grib_index_get( idx, 'level', ranges%lvllist( 1 : ranges%noflvl ) )
!   
!.....get the number of distinct values of ensembles
!
      IF ( ranges%elda ) THEN
         CALL grib_index_get_size( idx, 'number', ranges%nofens )
         CALL grib_index_get( idx, 'number', ranges%enslist( 1 : ranges%nofens ) )
      END IF
!
!.....get the number of distinct values of step in the index
!
      CALL grib_index_get_size( idx, 'step', ranges%nofstp )
!
      CALL grib_index_get( idx, 'step', ranges%stplist( 1 : ranges%nofstp ) )
!   
!.....unique data for WAM model
!
      IF ( ranges%ow ) THEN
!   
!........get the number of distinct values of direction in the index
!
         CALL grib_index_get_size( idx, 'direction', ranges%nofdir )
         CALL grib_index_get( idx, 'direction', ranges%dirlist( 1 : ranges%nofdir ) )
!   
!........get the number of distinct values of frequency in the index
!
         CALL grib_index_get_size( idx, 'frequency', ranges%noffrq )
         CALL grib_index_get( idx, 'frequency', ranges%frqlist( 1 : ranges%noffrq ) )
!
         CALL grib_index_select( idx, 'direction', ranges%dirlist( 1 ) )
         CALL grib_index_select( idx, 'frequency', ranges%frqlist( 1 ) )
!
      END IF
!
!.....Get grid properties: select index! (first elements!)
!
      CALL grib_index_select( idx, 'shortName', ranges%varlist( 1 ) )
      CALL grib_index_select( idx, 'level', ranges%lvllist( 1 ) )
      CALL grib_index_select( idx, 'step', ranges%stplist( 1 ) )
!
      IF ( ranges%elda ) CALL grib_index_select( idx, 'number', ranges%enslist( 1 ) )
!
!.....get grib block from file opened by index
!
      CALL grib_new_from_index( idx, igrib)
!
!.....Get date and time from grib file.
!
      CALL grib_get( igrib, 'dataDate', d )
      CALL grib_get( igrib, 'dataTime', t )
!
      WRITE ( dstr, "(I8)" ) d
      WRITE ( tstr, "(I4.4)" ) t
!
      CALL STR2EPOCH ( dstr( 1 : 4 ) // '/' // dstr( 5 : 6 ) // '/' // dstr( 7 : ) // ' ' // &
         & tstr( 1 : 2 ) // ':' // tstr( 3 : 4 ),ranges%epoch )
      CALL EPOCH2STR ( ranges%epoch, ranges%tstr )
!
!.....get latitude first point
!
      call grib_get( igrib,'latitudeOfFirstGridPointInDegrees', ranges%latMax ) 
!
!.....get longitude first point
!
      call grib_get( igrib, 'longitudeOfFirstGridPointInDegrees',ranges%lonMin ) 
!
!.....get latitude last point
!
      call grib_get( igrib, 'latitudeOfLastGridPointInDegrees', ranges%latMin ) 
!
!.....get longitude last point
!
      call grib_get( igrib, 'longitudeOfLastGridPointInDegrees', ranges%lonMax )
!
!.....check longitude
!
      IF ( ranges%lonMin .GT. ranges%lonMax ) ranges%lonMin = ranges%lonMin - 360
      IF ( ranges%lonMax .LT. ranges%lonMin ) ranges%lonMax = ranges%lonMax + 360
!
!.....get longitude grid increment
!
      call grib_get( igrib, 'iDirectionIncrementInDegrees', ranges%dLon )
!
!.....get latitude grid increment
!
      call grib_get( igrib, 'jDirectionIncrementInDegrees', ranges%dLat )
!
!.....number of coordinates lat and lon
!
      CALL grib_get_int  ( igrib, 'Nx', ranges%noflon )
      CALL grib_get_int  ( igrib, 'Ny', ranges%noflat )
!
!.....get size of values of array
!
      call grib_get( igrib, 'numberOfPoints', ranges%nofpoints )
!
!.....ocean wave model direction and frequency scale
!
      IF ( ranges%ow ) THEN
!
!........direction
!
         call grib_get( igrib, 'directionScalingFactor', ranges%dirscale )
!
!........frequency
!
         call grib_get( igrib, 'frequencyScalingFactor', ranges%frqscale )
!
      END IF
!
!.....release igrib
!
      CALL grib_release ( igrib, ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_RANGES grib_release'
!
!.....Close grib file
!
      CALL grib_index_release ( idx, ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_RANGES grib_index_release'
!
!.....Dimension
!
      ranges%dims( 1 ) = ranges%noflon .GT. 1
      ranges%dims( 2 ) = ranges%noflat .GT. 1
      ranges%dims( 3 ) = ranges%noflvl .GT. 1
!
      ranges%nofdims = COUNT ( ranges%dims )
!
      RETURN
!
   END SUBROUTINE GRIBTOOLS_RANGES
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GRIBTOOLS_PRINT_RANGES
!
!       Print ranges from structure.
!
!****************************
   Function GRIBTOOLS_COMPARE_RANGES ( ranges1, ranges2 ) RESULT ( match )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      type(type_gribtools_ranges), intent(in)      :: ranges1, ranges2
      LOGICAL                                        :: match
!
!  ---
!
!.....Open profile, if exists, stop otherwise.
!
      match = ranges1%dataType .EQ. ranges2%dataType
      IF ( .NOT. match ) RETURN
!
      match = ranges1%typeoflevel .EQ. ranges2%typeoflevel
      IF ( .NOT. match ) RETURN
!
      match = ranges1%gridtype .EQ. ranges2%gridtype
      IF ( .NOT. match ) RETURN
!
      match = ranges1%nofvar .EQ. ranges2%nofvar
      IF ( .NOT. match ) RETURN
!
      match = ranges1%noflvl .EQ. ranges2%noflvl
      IF ( .NOT. match ) RETURN
!
      match = ranges1%nofstp .EQ. ranges2%nofstp
      IF ( .NOT. match ) RETURN
!
      match = abs(ranges1%latMin-ranges2%latMin).lt.1.d-8
      IF ( .NOT. match ) RETURN
!
      match = abs(ranges1%latMax-ranges2%latMax).lt.1.d-8
      IF ( .NOT. match ) RETURN
!
      match = abs(ranges1%dlat-ranges2%dlat).lt.1.d-8
      IF ( .NOT. match ) RETURN
!
      match = ranges1%noflat .EQ. ranges2%noflat
      IF ( .NOT. match ) RETURN
!
      match = abs(ranges1%lonMin-ranges2%lonMin).lt.1.d-8
      IF ( .NOT. match ) RETURN
!
      match = abs(ranges1%lonMax-ranges2%lonMax).lt.1.d-8
      IF ( .NOT. match ) RETURN
!
      match = abs(ranges1%dlon-ranges2%dlon).lt.1.d-8
      IF ( .NOT. match ) RETURN
!
      match = ranges1%noflon .EQ. ranges2%noflon
      IF ( .NOT. match ) RETURN
!
      match = ranges1%nofdir .EQ. ranges2%nofdir
      IF ( .NOT. match ) RETURN
!
      match = abs(ranges1%dirscale-ranges2%dirscale).lt.1.d-8
      IF ( .NOT. match ) RETURN
!
      match = ranges1%noffrq .EQ. ranges2%noffrq
      IF ( .NOT. match ) RETURN
!
      match = abs(ranges1%frqscale-ranges2%frqscale).lt.1.d-8
      IF ( .NOT. match ) RETURN
!
      RETURN
!
   END FUNCTION GRIBTOOLS_COMPARE_RANGES
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GRIBTOOLS_PRINT_RANGES
!
!       Print ranges from structure.
!
!****************************
   SUBROUTINE GRIBTOOLS_PRINT_RANGES ( ranges )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      type(type_gribtools_ranges), intent(in)       :: ranges
!
!.....Local variables
!
      integer(kind=4) :: i
!
!  ---
!
!.....Open profile, if exists, stop otherwise.
!
      PRINT "(2A)",     'date_time       = ', ranges%tstr
      PRINT "(2A)",     'data_type       = ', ranges%dataType
      PRINT "(2A)",     'data_stream     = ', ranges%stream
      PRINT "(2A)",     'type_of_level   = ', ranges%typeoflevel
      PRINT "(2A)",     'grid_type       = ', ranges%gridtype
      PRINT "(A,I4)",   'nofvar          = ', ranges%nofvar
      WRITE ( *, "(A)", ADVANCE = 'no' ) 'varlist         ='
      DO i = 1, ranges%nofvar
         WRITE ( *, "(1x,A)", ADVANCE = 'no' ) TRIM ( ranges%varlist( i ) )
      END DO
      WRITE ( *, "(A)" ) ''
      IF ( ranges%nofadded .GT. 0 ) THEN
        PRINT "(A,I4)",   'nofvar_added    = ', ranges%nofadded
        WRITE ( *, "(A)", ADVANCE = 'no' ) 'varadded        ='
        DO i = 1, ranges%nofadded
           WRITE ( *, "(1x,A)", ADVANCE = 'no' ) TRIM ( ranges%varadded( i ) )
        END DO
        WRITE ( *, "(A)" ) ''
      END IF
      PRINT "(A,I4)",   'nofstp          = ', ranges%nofstp
      WRITE ( *, "(A)", ADVANCE = 'no' ) 'stplist         ='
      DO i = 1, ranges%nofstp
         WRITE ( *, "(1x,I4)", ADVANCE = 'no' ) ranges%stplist( i )
      END DO
      WRITE ( *, "(A)" ) ''
      PRINT "(A,I4)",   'noflvl          = ', ranges%noflvl
      IF ( ranges%elda ) PRINT "(A,I4)", 'nofens          = ', ranges%nofens
      PRINT "(A,F7.2)", 'lat_min         = ', ranges%latMin
      PRINT "(A,F7.2)", 'lat_max         = ', ranges%latMax
      PRINT "(A,F7.2)", 'dlat            = ', ranges%dlat
      PRINT "(A,I4)",   'noflat          = ', ranges%noflat
      PRINT "(A,F7.2)", 'lon_min         = ', ranges%lonMin
      PRINT "(A,F7.2)", 'lon_max         = ', ranges%lonMax
      PRINT "(A,F7.2)", 'dlon            = ', ranges%dlon
      PRINT "(A,I4)",   'noflon          = ', ranges%noflon
      IF ( ranges%ow ) THEN
         PRINT "(A,I4)",   'nofdir          = ', ranges%nofdir
         PRINT "(A,F7.2)", 'dirscale        = ', ranges%dirscale
         PRINT "(A,I4)",   'noffrq          = ', ranges%noffrq
         PRINT "(A,F10.2)",'frqscale        = ', ranges%frqscale
      END IF
!
   END SUBROUTINE GRIBTOOLS_PRINT_RANGES
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GRIBTOOLS_READ_WAM_2DFD
!
!       read ECMWF wam 2dfd data from grib file
!
!****************************
   SUBROUTINE GRIBTOOLS_READ_WAM_2DFD ( file1, ranges, dat )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character(len=*), INTENT( IN )                                :: file1
      type(type_gribtools_ranges), INTENT( IN )                      :: ranges
      REAL, INTENT( OUT ), DIMENSION( ranges%noflon,  & 
         & ranges%noflat, ranges%noffrq, ranges%nofdir )                :: dat
!
!.....Local variables
!
      LOGICAL            :: lexist
!
      integer(kind=4)          :: ierr, idx, igrib, idir, ifrq
!
      DOUBLE PRECISION, PARAMETER  :: missingValue = 9999.d0
      REAL, DIMENSION ( ranges%nofpoints )              :: values
      REAL, DIMENSION ( ranges%noflon, ranges%noflat )  :: ll
!
!  ---
!
!.....initialize
!
      dat = missingValue
!
!.....Open profile, if exists, stop otherwise stop
!
      INQUIRE( FILE = file1, EXIST = lexist )
      IF ( .NOT. lexist ) STOP 'ERROR @ GRIBTOOLS_READ Grib file does not exist'
!
!.....create an index from a grib file using some keys
!
      CALL grib_index_create ( idx, file1, 'shortName,step,directionNumber,frequencyNumber', ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_READ Could not open file'
!
!.....set variable to 2dfd
!
      CALL grib_index_select (  idx, 'shortName', '2dfd'  )
!
!.....set step to first step
!
      CALL grib_index_select (  idx, 'step', ranges%stplist( 1 )  )
!
!$OMP DO
!
!.....loop on frequency
!
      DO ifrq = 1, ranges%noffrq
!
         CALL grib_index_select (  idx, 'frequencyNumber', ranges%frqlist( ifrq )  )
!
!........loop on direction
!
         DO idir = 1, ranges%nofdir
!
            CALL grib_index_select (  idx, 'directionNumber', ranges%dirlist( idir )  )
!
!...........get grib block for specified step, level, and variable
!
            CALL grib_new_from_index ( idx, igrib, ierr )
!
            IF ( ierr .EQ. 0 ) THEN
!
!..............get the data
!
               call grib_get( igrib, 'values', values )
!
!..............store to dat
!
               ll = RESHAPE ( values, (/ ranges%noflon, ranges%noflat /) )
               dat( :, :, ifrq, idir ) = ll( :, ranges%noflat : 1 : - 1 )
!
!..............release grib data
!
               CALL grib_release ( igrib, ierr )
!
            END IF
!
         END DO
!
      END DO
!
!$OMP END DO
!
!.....release grib file index
!
      CALL grib_index_release ( idx )
!
      RETURN
!
  END SUBROUTINE GRIBTOOLS_READ_WAM_2DFD
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GRIBTOOLS_READ_ALL_WAM
!
!       read ECMWF wam data from grib file
!
!****************************
   SUBROUTINE GRIBTOOLS_READ_ALL_WAM ( file1, ranges, dat )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character(len=*), INTENT( IN )                                :: file1
      type(type_gribtools_ranges), INTENT( IN )                      :: ranges
      REAL, INTENT( OUT ), DIMENSION( ranges%noflon,  & 
         & ranges%noflat, ranges%noffrq, ranges%nofdir, &
         & ranges%nofstp, ranges%nofvar )                               :: dat
!
!.....Local variables
!
      LOGICAL            :: lexist
!
      integer(kind=4)        :: ierr, idx, igrib
      integer(kind=4)        :: idir, ifrq, istp, ivar
!
      integer(kind=4), PARAMETER  :: missingValue = 9999
      REAL, DIMENSION ( ranges%nofpoints )              :: values
      REAL, DIMENSION ( ranges%noflon, ranges%noflat )  :: ll
!
!  ---
!
!.....initialize
!
      dat = REAL ( missingValue )
!
!.....Open profile, if exists, stop otherwise stop
!
      INQUIRE( FILE = file1, EXIST = lexist )
      IF ( .NOT. lexist ) STOP 'ERROR @ GRIBTOOLS_READ Grib file does not exist'
!
!.....create an index from a grib file using some keys
!
      CALL grib_index_create ( idx, file1, 'shortName,step,directionNumber,frequencyNumber', ierr )
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_READ Could not open file'
!
!.....loop on step
!
      DO istp = 1, ranges%nofstp
!
         CALL grib_index_select ( idx, 'step', ranges%stplist( istp ) )
!
!........loop on frequency
!
         DO ifrq = 1, ranges%noffrq
!
            CALL grib_index_select (  idx, 'frequencyNumber', ranges%frqlist( ifrq )  )
!
!...........loop on frequency
!
            DO idir = 1, ranges%nofdir
!
               CALL grib_index_select (  idx, 'directionNumber', ranges%dirlist( idir )  )
!
!..............loop on variables
!
               DO ivar = 1, ranges%nofvar
!
                  CALL grib_index_select (  idx, 'shortName', ranges%varlist( ivar )  )
!
!.................get grib block for specified step, level, and variable
!
                  CALL grib_new_from_index ( idx, igrib, ierr )
!
                  IF ( ierr .EQ. 0 ) THEN
!
!....................get the data
!
                     call grib_get( igrib, 'values', values )
!
!....................apply corrections / changes
!
                     IF ( ranges%varlist( ivar ) .EQ. '2dfd' ) THEN
                        WHERE ( int(values,kind=4).ne.missingValue ) 
                           values = 10 ** values
                        END WHERE
                     END IF
!
!....................store to dat
!
                     ll = RESHAPE ( values, (/ ranges%noflon, ranges%noflat /) )
                     dat( :, :, ifrq, idir, istp, ivar ) = ll( :, ranges%noflat : 1 : - 1 )
!
!....................release grib data
!
                     CALL grib_release ( igrib, ierr )
!
                  END IF
!
               END DO
!
            END DO
!
         END DO
!
      END DO
!
!.....release grib file index
!
      CALL grib_index_release ( idx )
!
      RETURN
!
  END SUBROUTINE GRIBTOOLS_READ_ALL_WAM
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GRIBTOOLS_READ_ALL_ISOBARIC
!
!       read ECMWF atmospheric profile from grib file
!
!****************************
   SUBROUTINE GRIBTOOLS_READ_ALL_ISOBARIC ( file1, ranges, dat, topo, ens )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character(len=*), INTENT( IN )                                :: file1
      type(type_gribtools_ranges), INTENT( IN )                      :: ranges
      DOUBLE PRECISION, INTENT( OUT ), DIMENSION( ranges%noflon,  & 
         & ranges%noflat, ranges%noflvl, ranges%nofstp, ranges%nofvar ) :: dat
      DOUBLE PRECISION, INTENT( OUT ), OPTIONAL, &
         & DIMENSION( ranges%noflon, ranges%noflat )                    :: topo
      integer(kind=4), INTENT( IN ), OPTIONAL                               :: ens
!
!.....Local variables
!
      LOGICAL            :: lexist, lmatch
!
      integer(kind=4)          :: ierr, idx, igrib
      integer(kind=4)          :: i, ilon, ilat, ilvl, istp, ivar
!
      DOUBLE PRECISION   :: kelvin, missingValue, grav
      DOUBLE PRECISION, DIMENSION( ranges%noflon * ranges%noflat ) :: values
!
!.....default values
!
      DATA    kelvin        /    273.15d0     /
      DATA    grav          /      9.80665d0  /
      DATA    missingValue  /   9999.d0       /
!
!  ---
!
!.....initialize
!
      lmatch=.true.
      dat = missingValue
!
!.....Open profile, if exists, stop otherwise stop
!
      INQUIRE( FILE = file1, EXIST = lexist )
      IF ( .NOT. lexist ) STOP 'ERROR @ GRIBTOOLS_READ Grib file does not exist'
!
!.....create an index from a grib file using some keys
!
      IF ( ranges%elda ) THEN
         CALL grib_index_create( idx, file1, 'shortName,level,step,perturbationNumber', ierr )
         IF ( PRESENT ( ens ) ) THEN
            DO i = 1, ranges%nofens
               lmatch = ens .EQ. ranges%enslist( i )
               IF ( lmatch ) EXIT
            END DO
            IF ( lmatch ) THEN
               CALL grib_index_select( idx, 'perturbationNumber', ens  )
            ELSE
               STOP 'ERROR @ GRIBTOOLS_READ Illegal ensemble number for elda provided'
            END IF
         ELSE
            CALL grib_index_select( idx, 'perturbationNumber', ranges%enslist( 1 )  )
         END IF
      ELSE
         CALL grib_index_create( idx, file1, 'shortName,level,step', ierr )
      END IF
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_READ Could not open file'
!
!.....loop on step
!
      DO istp = 1, ranges%nofstp
!
         CALL grib_index_select( idx, 'step', ranges%stplist( istp ) )
!
!........loop on level
!
         DO ilvl = ranges%noflvl, 1, - 1
!
            CALL grib_index_select(  idx, 'level', ranges%lvllist( ranges%noflvl - ilvl + 1 )  )
!
!...........loop on variables
!
            DO ivar = 1, ranges%nofvar - 1
!
               CALL grib_index_select(  idx, 'shortName', ranges%varlist( ivar )  )
!
!..............get grib block for specified step, level, and variable
!
               CALL grib_new_from_index( idx, igrib, ierr )
!
               IF ( ierr .EQ. 0 ) THEN
!
!.................get the data
!
                  CALL grib_get( igrib, 'values', values )
!
!.................apply corrections
!
                  IF ( ranges%varlist( ivar ) .EQ. 'z' )  values = values / grav
                  IF ( ranges%varlist( ivar ) .EQ. 'pv' ) values = values * 10.d0**6
!
!.................store to dat
!
                  i = 0
                  DO ilat = ranges%noflat , 1 , -1
                     DO ilon = 1 , ranges%noflon
                        i = i + 1 
                        dat( ilon, ilat, ilvl, istp, ivar ) = values( i )
                     END DO
                  END DO
!
!.................release grib data
!
                  CALL grib_release ( igrib, ierr )
!
               END IF
!
            END DO
!
!...........set pressure
!
            dat( :, :, ilvl, istp, ranges%nofvar ) = ranges%lvllist( ranges%noflvl - ilvl + 1 ) * 100
!
         END DO
!
      END DO
!
!.....release grib file index
!
      CALL grib_index_release ( idx )
!
      IF ( PRESENT ( topo ) ) topo = dat( :, :, 1, 1, GRIBTOOLS_INDEX ( ranges, 'z' ) )
!
      RETURN
!
  END SUBROUTINE GRIBTOOLS_READ_ALL_ISOBARIC
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GRIBTOOLS_READ_ALL_HYBRID
!
!       read ECMWF atmospheric profile from grib file
!
!****************************
   SUBROUTINE GRIBTOOLS_READ_ALL_HYBRID ( file1, ranges, dat, topo, ens )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character(len=*), INTENT( IN )                                :: file1
      TYPE( type_gribtools_ranges ), INTENT( IN )                       :: ranges
      DOUBLE PRECISION, INTENT( OUT ), DIMENSION( ranges%noflon,  & 
         & ranges%noflat, ranges%noflvl, ranges%nofstp, ranges%nofvar ) :: dat
      DOUBLE PRECISION, INTENT( OUT ), OPTIONAL, &
         & DIMENSION( ranges%noflon, ranges%noflat )                    :: topo
      integer(kind=4), INTENT( IN ), OPTIONAL                               :: ens
!
!.....Local variables
!
      LOGICAL            :: lexist, lset, lmatch, flip
!
      integer(kind=4)          :: i,idx,ilon,ilat,ilvl,istp,ivar
      integer(kind=4)          :: ixp, ixq, ixt, ixz
      integer(kind=4)          :: ierr, igrib 
!
      DOUBLE PRECISION   :: pnul, grav, kelvin, runiv
      DOUBLE PRECISION   :: pold, pstd, wsat, wprs, ppmv
      DOUBLE PRECISION   :: rn2, ro2, rar, rco2, rtot
      DOUBLE PRECISION   :: qh, qc, qn, qo, qar, qn2, qo2, qco2, qh2o, qdry
      DOUBLE PRECISION   :: qwet, shm, scale, sold, anew, pnew
      DOUBLE PRECISION   :: missingValue
!
      DOUBLE PRECISION, DIMENSION ( ranges%noflvl )                               :: a, b
      DOUBLE PRECISION, DIMENSION ( ranges%noflvl * 2 + 2 )                       :: pv
      DOUBLE PRECISION, DIMENSION ( ranges%nofLon * ranges%nofLat )               :: values
      DOUBLE PRECISION, DIMENSION ( ranges%nofLon, ranges%nofLat, ranges%nofstp ) :: psur, asur
!
!.....default values
!
      DATA   missingValue  /   9999.d0       /
      DATA      kelvin  /    273.15d0     /
      DATA      grav    /      9.80665d0  /  ! WMO Technical Regulations
      DATA      runiv   /      8.314472d0 /
      DATA      pnul    / 100000.0d0      /
      DATA      pstd    / 101325.0d0      / ! Pa (wiki)
!
!.....source:
!     http://en.wikipedia.org/wiki/Earth%27s_atmosphere
!     http://en.wikipedia.org/wiki/Periodic_table_%28detailed%29
!
      DATA      rn2     /      0.78084d0  /
      DATA      ro2     /      0.20946d0  /
      DATA      rar     /      0.00934d0  /
      DATA      rco2    /      0.00038d0  /
!
      DATA      qh      /      1.00794d0  /
      DATA      qc      /     12.0107d0   /
      DATA      qn      /     14.00674d0  /
      DATA      qo      /     15.9994d0   /
      DATA      qar     /     39.948d0    /
!
!  ---
!
!.....Check data properties
!
      IF ( ranges%typeOfLevel .NE. 'hybrid' ) STOP '@GRIBTOOLS_READ_HYBRID_ALL typeOfLevel is not hybrid'
      IF ( ranges%gridType .NE. 'regular_ll' ) STOP '@GRIBTOOLS_READ_HYBRID_ALL gridType is not regular_ll'

!!call gribtools_print_ranges( ranges )
!
!.....initialize
!
      lset = .true.
      lmatch=.true.
!
      asur = 0.d0
      psur = 1013.25d0
      dat  = missingValue
!
      rtot = rn2 + ro2 + rar + rco2
      qn2  = qn + qn
      qo2  = qo + qo
      qco2 = qc + qo + qo
      qh2o = qh + qh + qo
      qdry = qn2 * rn2 + qo2 * ro2 + qar * rar + qco2 * rco2
      qdry = qdry / rtot
!
!.....Open profile, if exists, stop otherwise stop
!
      INQUIRE ( FILE = file1, EXIST = lexist )
      IF ( .NOT. lexist ) STOP '@GRIBTOOLS_READ Grib file does not exist'
!
!.....create an index from a grib file using some keys
!
      IF ( ranges%elda ) THEN
         CALL grib_index_create( idx, file1, 'shortName,level,step,perturbationNumber', ierr )
         IF ( PRESENT ( ens ) ) THEN
            DO i = 1, ranges%nofens
               lmatch = ens .EQ. ranges%enslist( i )
               IF ( lmatch ) EXIT
            END DO
            IF ( lmatch ) THEN
               CALL grib_index_select( idx, 'perturbationNumber', ens  )
            ELSE
               STOP 'ERROR @ GRIBTOOLS_READ Illegal ensemble number for elda provided'
            END IF
         ELSE
            CALL grib_index_select( idx, 'perturbationNumber', ranges%enslist( 1 )  )
         END IF
      ELSE
         CALL grib_index_create( idx, file1, 'shortName,level,step', ierr )
      END IF
      IF ( ierr .NE. 0 ) STOP 'ERROR @ GRIBTOOLS_READ Could not open file'
!
!.....get indexes
!
      ixp    = GRIBTOOLS_INDEX ( ranges, 'p' )
      ixq    = GRIBTOOLS_INDEX ( ranges, 'q' )
      ixt    = GRIBTOOLS_INDEX ( ranges, 't' )
      ixz    = GRIBTOOLS_INDEX ( ranges, 'z' )
!!print *, ixp, ixq, ixt, ixz
!
!.....loop on step
!
      DO istp = 1, ranges%nofstp
!
         CALL grib_index_select(  idx, 'step', ranges%stplist( istp )  )
!!print *, ranges%stplist( istp )
!
!........loop on level
!
         DO ilvl = 1, ranges%noflvl
!
            CALL grib_index_select( idx, 'level', ranges%lvllist( ilvl ) )
!!print *, ranges%lvllist( ilvl )
!
!...........loop on variables
!
            DO ivar = 1, ranges%nofvar

               if (gribtools_added(ranges,ranges%varlist(ivar))) cycle
!
               CALL grib_index_select( idx, 'shortName', ranges%varlist( ivar ) )
!!print *, ranges%varlist( ivar )
!
!..............get grib block for specified step, level, and variable
!
               CALL grib_new_from_index( idx, igrib, ierr )
!
               IF ( ierr .EQ. 0 ) THEN
!
!.................get the data
!
                  CALL grib_get( igrib, 'values', values, ierr )
!
!.................apply corrections
!
                  flip=.true.
!
                  if (ranges%varlist(ivar).eq.'lnsp' ) then
                     values = DEXP( values )
                     flip=.false.
                  elseif (ranges%varlist(ivar).eq.'z' ) then
                     values = values / grav
                     flip=.false.
                  elseif (ranges%varlist(ivar).eq.'h' ) then
                     flip=.false.
                  elseif (ranges%varlist(ivar).eq.'pres' ) then
                     flip=.false.
                  elseif (ranges%varlist(ivar).eq.'pv' ) then
                     values = values * 10.d0**6
                  end if
!
!.................store to dat
!
                  i = 0
                  do ilat = ranges%noflat, 1, -1
                     do ilon = 1, ranges%noflon
                        i = i + 1
                        if ( flip ) then
                           dat( ilon, ilat, ranges%noflvl + 1 - ilvl, istp, ivar ) = values( i )
                        else
                           dat( ilon, ilat, ilvl, istp, ivar ) = values( i )
                        end if
                     end do
                  end do
!
!.................get pv values
!
                  IF ( lset ) THEN
                     call grib_get( igrib, 'pv', pv, ierr )
                     if (ierr.eq.0) then
                       do i=1, ranges%noflvl
                          a(ranges%noflvl+1-i)=( pv(i                ) + pv(i+1              ) )*0.5d0
                          b(ranges%noflvl+1-i)=( pv(i+ranges%noflvl+1) + pv(i+ranges%noflvl+2) )*0.5d0
                        end do
                     elseif (ierr.eq.-10.and.ranges%centre.eq.'99') then !harmonie
                          a=(/ 0.00000, 0.33262, 2.14204, 7.05859, 16.86313, 33.62249, &
                               59.68656, 97.68404, 150.51797, 221.36063, 313.64705, 431.06597, &
                               577.54740, 757.24521, 974.51338, 1233.87397, 1539.97471, 1897.53354, &
                               2307.24196, 2765.18335, 3270.49745, 3825.79219, 4433.29007, 5094.67211, &
                               5810.94701, 6582.34202, 7408.21146, 8286.95679, 9215.94957, 10191.44472, &
                               11208.46771, 12260.65515, 13340.02463, 14436.64712, 15538.19496, 16629.34171, &
                               17690.99893, 18699.39245, 19625.01060, 20431.50353, 21074.68088, 21501.84615, &
                               21654.72408, 21495.66550, 21021.57423, 20246.73475, 19194.25311, 17894.90172, &
                               16385.82649, 14709.27011, 12911.49495, 11042.22017, 9155.86698, 7324.44164, &
                               5614.33906, 4057.31113, 2674.80515, 1501.78627, 622.50823, 135.91409 /)
 
                          b=(/ 0.9988404, 0.9962437, 0.9928833, 0.9885790, 0.9833674, 0.9772583, &
                               0.9702498, 0.9623322, 0.9534903, 0.9437037, 0.9329483, 0.9211960, &
                               0.9084156, 0.8945727, 0.8796307, 0.8635503, 0.8462910, 0.8278109, &
                               0.8082468, 0.7878898, 0.7668276, 0.7449789, 0.7222823, 0.6986942, &
                               0.6741870, 0.6487477, 0.6223755, 0.5950806, 0.5668823, 0.5378080, &
                               0.5078923, 0.4771773, 0.4457134, 0.4135616, 0.3807978, 0.3475181, &
                               0.3138473, 0.2799486, 0.2460356, 0.2123856, 0.1793520, 0.1473753, &
                               0.1174521, 0.0909084, 0.0682819, 0.0494592, 0.0342425, 0.0223625, &
                               0.0134910, 0.0072549, 0.0032485, 0.0010406, 0.0001655, 0.0000000, &
                               0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000 /)
                     end if
                     lset = .false.
                  END IF
!
!.................release grib data
!
                  CALL grib_release ( igrib, ierr )
!
               END IF
!
            END DO
!
         END DO
!
      END DO
!
!.....release grib file index
!
      CALL grib_index_release ( idx )
!
!.....set ground pressure
!
      if ( gribtools_present(ranges,'lnsp') ) then
         psur = dat(:,:,1,:,gribtools_index(ranges,'lnsp') )
      elseif ( gribtools_present(ranges,'pres') ) then
         psur = dat(:,:,1,:,gribtools_index(ranges,'pres') )
      end if
!
!.....set ground level
!
      if ( gribtools_present(ranges,'z') .and. .not.gribtools_added(ranges,'z') ) then
         asur = dat( :, :, 1, :, ixz )
      elseif ( gribtools_present(ranges,'h') .and. .not.gribtools_added(ranges,'h') ) then
         asur = dat(:,:,1,:,gribtools_index(ranges,'h') )
      end if
      if ( PRESENT ( topo ) ) then
         if ( GRIBTOOLS_PRESENT ( ranges, 'z' ) ) then
            topo = asur( :, :, 1 )
         else
            topo = 0.d0
         end if
      end if
!
!.....set pressure
!
      do istp = 1, ranges%nofstp
         do ilvl = 1, ranges%noflvl
            dat( :, :, ilvl, istp, ixp ) = a( ilvl ) + b( ilvl ) * psur( :, :, istp )
         end do
      end do
!
!.....compute relative humidity and altitude
!     ( specific humidity = kg(water)/kg(wet air) )
!
      DO istp = 1, ranges%nofstp
         DO ilat = 1, ranges%noflat
            DO ilon = 1, ranges%noflon
!
               pold = psur( ilon, ilat, istp )
               anew = asur( ilon, ilat, istp )
!
               shm = MIN (  dat( ilon, ilat, 1, istp, ixq ), 0.9d0  )
               shm = shm / ( 1.0d0 - shm )
               qwet = ( 1.0d0 + shm ) / ( 1.0d0 + shm * qdry / qh2o ) * qdry
               sold = 1000.0d0 * runiv * dat( ilon, ilat, 1, istp, ixt ) / qwet / grav
!
               DO ilvl = 1 , ranges%noflvl
!
                  pnew = dat( ilon, ilat, ilvl, istp, ixp )
!
                  shm = MIN ( dat( ilon, ilat, ilvl, istp, ixq ), 0.9d0 )
                  shm = shm / ( 1.0d0 - shm )
                  qwet = ( 1.0d0 + shm ) / ( 1.0d0 + shm * qdry / qh2o ) * qdry
!
                  scale = 1000.0d0 * runiv * dat( ilon, ilat, ilvl, istp, ixt ) / qwet / grav
                  anew = anew + 0.5d0*( scale + sold ) * DLOG( pold / pnew )
                  dat( ilon, ilat, ilvl, istp, ixz ) = anew
                  pold = pnew
                  sold = scale
!
                  ppmv = dat( ilon, ilat, ilvl, istp, ixq) * qwet / qh2o
                  wprs = dat( ilon, ilat, ilvl, istp, ixp) * ppmv
                  wsat = WEXLER (  dat( ilon, ilat, ilvl, istp, ixt ), 0  )
                  dat( ilon, ilat, ilvl, istp, ixq ) = wprs / wsat 
!
               END DO
            END DO
         END DO
      END DO 
!
      RETURN
!
  END SUBROUTINE GRIBTOOLS_READ_ALL_HYBRID
! =========================================================================================


! =========================================================================================
!
!.......FUNCTION WEXLER
!
!       Hyland, R. W. and A. Wexler,
!       Formulations for the Thermodynamic Properties of the saturated
!       Phases of H2O from 173.15K to 473.15K,
!       ASHRAE Trans, 89(2A), 500-519, 1983. 
!       with T in [K] and pw in [Pa]
!
!       http://cires.colorado.edu/~voemel/vp.html
!
!****************************
      DOUBLE PRECISION FUNCTION wexler ( T, ice )
!****************************

      implicit none
!
!.....Dummy variables
!
      INTEGER         , intent(in)      :: ice
      DOUBLE PRECISION, intent(in)      :: T
!
!.....Local variables
!
      DOUBLE PRECISION                 :: wvpl
!
! ---
!
      IF ( ice .EQ. 0 ) THEN
         wvpl = -0.58002206E+04 / T         &
            &   + 0.13914993E+01            &
            &   - 0.48640239E-01 * T        &
            &   + 0.41764768E-04 * T ** 2   &
            &   - 0.14452093E-07 * T ** 3   &
            &   + 0.65459673E+01 * LOG ( T )
      ELSE
         wvpl = -0.56745359E+04 / T         &
            &   + 0.63925247E+01            &
            &   - 0.96778430E-02 * T        &
            &   + 0.62215701E-06 * T ** 2   &
            &   + 0.20747825E-08 * T ** 3   &
            &   - 0.94840240E-12 * T ** 4   &
            &   + 0.41635019E+01 * LOG ( T ) 
      END IF
!
      wexler = DEXP ( wvpl )
!
      RETURN
!
      END FUNCTION WEXLER
!
! =========================================================================================


! =========================================================================================
!
!.......INTEGER FUNCTION GRIBTOOLS_INDEX
!
!       get index of grib variable in list
!
!****************************
   INTEGER FUNCTION GRIBTOOLS_INDEX ( ranges, var ) RESULT ( ix )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      TYPE( type_gribtools_ranges ), INTENT( IN )              :: ranges
      character(len=*), INTENT( IN )                       :: var
!
!.....Local variables
!
      integer(kind=4)                                                :: ivar
!
!  ---
!
!.....initialize
!
      ix = - 1
!
!.....loop on variables
!
      DO ivar = 1, ranges%nofvar
         IF ( ranges%varlist( ivar ) .EQ. var ) THEN
            ix = ivar
            EXIT
         END IF
      END DO
!
!.....index
!
      IF ( ix .EQ. - 1 ) THEN
         PRINT "(A)", 'ERROR: GRIBTOOLS_INDEX variable "' // var // '" not found, no index returned.'
         STOP
      END IF
!
      RETURN
!
  END FUNCTION GRIBTOOLS_INDEX
!
! =========================================================================================


! =========================================================================================
!
!.......LOGICAL FUNCTION GRIBTOOLS_PRESENT
!
!       variable present in grib variable list
!
!****************************
   LOGICAL FUNCTION GRIBTOOLS_PRESENT ( ranges, var ) RESULT ( lvar )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      TYPE( type_gribtools_ranges ), INTENT( IN )             :: ranges
      character(len=*), INTENT( IN )                       :: var
!
!.....Local variables
!
      integer(kind=4)                                                :: ivar
!
!  ---
!
!.....initialize
!
      lvar = .TRUE.
!
!.....loop on variables
!
      DO ivar = 1, ranges%nofvar
         IF (  ranges%varlist( ivar ) .EQ. var  ) RETURN
      END DO
!
!.....no match
!
      lvar = .FALSE.
!
      RETURN
!
  END FUNCTION GRIBTOOLS_PRESENT
!
! =========================================================================================


! =========================================================================================
!
!.......LOGICAL FUNCTION GRIBTOOLS_PRESENT
!
!       variable present in grib variable list
!
!****************************
   LOGICAL FUNCTION GRIBTOOLS_ADDED ( ranges, var ) RESULT ( lvar )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      TYPE( type_gribtools_ranges ), intent(in)      :: ranges
      character(len=*)         , intent(in)      :: var
!
!.....Local variables
!
      integer(kind=4)                                   :: ivar
!
!  ---
!
!.....initialize
!
      lvar = .TRUE.
!
!.....loop on variables
!
      DO ivar = 1, ranges%nofadded
         IF (  ranges%varadded( ivar ) .EQ. var  ) RETURN
      END DO
!
!.....no match
!
      lvar = .FALSE.
!
      RETURN
!
  END FUNCTION GRIBTOOLS_ADDED
!
! =========================================================================================
!
END MODULE TOOLS_GRIB
