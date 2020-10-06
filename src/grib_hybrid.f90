!*************************************************************************************************************************
!
!                                                   G R I B _ H Y B R I D 
!
!  Module:       GRIB_HYBRID
!
!  Programmer:   Pieter S. M. Smets
!                R&D depart. of Seismology and Acoustics - Koninklijk Nederlands Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         September 9, 2016
!
!  Language:     Fortran-90
!
!  Description:  This module includes several subroutines to read 2dfd ocean wave spectra data and integrate over
!                direction (full direction range -> significant wave height, opposite direction -> Hasselmann )
!
!
!*************************************************************************************************************************

Module grib_hybrid 

   use io_types
   use syms
   use grib_api
   use iso_c_binding, only: C_LOC, C_F_POINTER
   use omp_lib

   implicit none

!
!..types
!
   type type_grib_grid_ll
     integer(int32)           :: count
     real(double)             :: first, last, increment, range
     real(double), dimension(:), allocatable  :: dat
   end type type_grib_grid_ll

   type type_grib_grid
     character(len=40)        :: gridName, gridType, packingType
     logical                  :: isOctahedral, iScansNegatively, jScansPositively, &
                                 jPointsAreConsecutive, distinctGridValues, earthIsOblate
     integer(int32)           :: nofpoints
     integer(int32), dimension(:), allocatable  :: pl !  number of points along a full parallel
     type(type_grib_grid_ll)  :: lon, lat
   end type type_grib_grid

   type type_grib_num
     logical                  :: defined 
     integer(int32)           :: count, sel
     integer(int32), dimension(:), allocatable :: list
   end type type_grib_num

   type type_grib_step
     character(len=40)        :: units
     integer(int32)           :: count, sel, startStep, endStep, stepRange
     integer(int32), dimension(:), allocatable :: list
   end type type_grib_step

   type type_grib_level
     character(len=40)        :: typeOfLevel
     integer(int32)           :: count
     real(double)  , dimension(:), allocatable :: pv, a, b
     integer(int32), dimension(:), allocatable :: list
     integer(int32), dimension(:), pointer     :: list_flip
   end type type_grib_level

   type type_grib_varSurf
     logical                                      :: defined, surface
     character(len=40)                            :: shortName, longName, units, typeOfLevel
     integer(int32)                               :: index, paramId, missingValue
     real(single), dimension(:)  , allocatable    :: dat
     real(single), dimension(:,:), pointer        :: regular_ll
     real(single), dimension(:,:), pointer        :: regular_ll_flip
   end type type_grib_varSurf

   type type_grib_varLevel
     logical                                      :: defined, surface 
     character(len=40)                            :: shortName, longName, units, typeOfLevel
     integer(int32)                               :: index, paramId, missingValue
     real(single), dimension(:,:)  , allocatable  :: dat
     real(single), dimension(:,:,:), pointer      :: regular_ll
     real(single), dimension(:,:,:), pointer      :: regular_ll_flip
   end type type_grib_varLevel

   type type_grib_hybrid
      type(type_grib_grid)     :: grid
      type(type_grib_num)      :: number
      type(type_grib_step)     :: step
      type(type_grib_level)    :: level
      logical                  :: isECMWF
      integer(int32)           :: idx, experimentVersionNumber, editionNumber, params
      integer(int64)           :: epoch
      character(len=26)        :: time
      character(len=40)        :: centre, dataClass, dataType, dataStream
      type(type_grib_varLevel) :: q,t,u,v,w,z,p
      type(type_grib_varSurf)  :: lnsp
      character(len=40),dimension(:),allocatable  :: paramList
   end type type_grib_hybrid

Contains

! =========================================================================================
!
!.......function isgrib 
!
!       Check if file is of grib type
!
!****************************
   Logical function isgrib ( gribFile )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character(len=*), intent(in) :: gribFile
!
!.....Local variables
!
      integer(int32)  :: ifile, iret
!
!  ---
!
!.....Open the grib file
!
      call grib_open_file ( ifile, gribFile, 'r', iret )
      isgrib = iret.eq.0_int32
      if (isgrib) call grib_close_file( ifile )
!
      return
!
   End function isgrib
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_referenceTime
!
!       Get ranges from grib file.
!
!****************************
   Subroutine grib_hybrid_referenceTime ( gribFile, epoch, timestring, mask )
!****************************
!
      use time
      implicit none
!
!.....Dummy variables
!
      character(len=*), intent(in)             :: gribFile
      integer(int64)  , intent(out)            :: epoch
      character(len=*), intent(out), optional  :: timestring 
      logical         , intent(out), optional  :: mask 
!
!.....Local variables
!
      logical         :: stopOnError, lexist
      integer(int32)  :: ifile, igrib, iret, d, t
      character       :: d_str*8, t_str*4, tstr*23
!
!  ---
!
!.....Init
!
      if (present(mask)) then 
         stopOnError=.false.
         mask=.false.
      else
         stopOnError=.true.
      end if
      epoch = 0_int64
      if (present(timestring)) timestring=''
!
!.....Open profile, if exists, stop otherwise.
!
      inquire( file = gribFile, exist = lexist )
      if (.not.lexist.and.stopOnError) then
         stop 'Error @ grib_hybrid_referenceTime : grib file does not exist'
      else
         return
      end if
!
!.....Open file and load first data
!
      call grib_open_file ( ifile, gribFile, 'r', iret )
      if (iret.ne.0_int32.and.stopOnError) then
         stop 'Error @ grib_hybrid_referenceTime : could not open file'
      else
         return
      end if
!
      call grib_new_from_file ( ifile, igrib, iret )
      if (iret.ne.0_int32.and.stopOnError) then
         stop 'Error @ grib_hybrid_referenceTime : could not get new from file'
      else
         return
      end if
!
!.....Get date and time from grib file.
!
      call grib_get( igrib, 'dataDate', d, iret )
      if (iret.ne.0_int32.and.stopOnError) then
         stop 'Error @ grib_hybrid_referenceTime : could not get dataDate'
      else
         return
      end if
      call grib_get( igrib, 'dataTime', t, iret )
      if (iret.ne.0_int32.and.stopOnError) then
         stop 'Error @ grib_hybrid_referenceTime : could not get dataTime'
      else
         return
      end if
!
      write (d_str,"(i8)"  ) d
      write (t_str,"(i4.4)") t
!
      tstr = d_str(1:4)//'/'//d_str(5:6)//'/'//d_str(7:)//' '//t_str(1:2)//':'//t_str(3:)//':00.000'
      call str2epoch ( tstr, epoch )
!
!.....Release grib and close file.
!
      call grib_release( igrib )
      call grib_close_file( ifile )
!
      if (present(timestring)) call epoch2str( epoch, timestring )
      if (present(mask)) mask=.true.
!
      return
!
   End subroutine grib_hybrid_referenceTime
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_clear_varLevel
!
!       Initialize varLevel structure.
!
!****************************
   Subroutine grib_hybrid_clear_varLevel ( varLevel )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      type(type_grib_varLevel), intent(inout), target  :: varLevel
!
!  ---
!
     varLevel%defined = .false.
     varLevel%shortName = ''
     varLevel%longName = ''
     varLevel%units = ''
     varLevel%typeOfLevel = ''
     varLevel%index = 0_int32
     varLevel%paramId = 0_int32
     varLevel%missingValue = 9999_int32
     if (allocated(varLevel%dat)) deallocate(varLevel%dat)
     if (associated(varLevel%regular_ll)) nullify(varLevel%regular_ll)
     if (associated(varLevel%regular_ll_flip)) nullify(varLevel%regular_ll_flip)
!
      return
!
   End subroutine grib_hybrid_clear_varLevel
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_compare_varLevel
!
!       Initialize varLevel structure.
!
!****************************
   function grib_hybrid_compare_varLevel ( varLevelA, varLevelB ) result ( match )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      type(type_grib_varLevel), intent(in)  :: varLevelA, varLevelB
      logical :: match
!
!  ---
!
     match = varLevelA%defined .eqv. varLevelB%defined
     if (.not.match) return
     match = varLevelA%shortName .eq. varLevelB%shortName
     if (.not.match) return
     match = varLevelA%longName .eq. varLevelB%longName
     if (.not.match) return
     match = varLevelA%units .eq. varLevelB%units
     if (.not.match) return
     match = varLevelA%typeOfLevel .eq. varLevelB%typeOfLevel
     if (.not.match) return
     match = varLevelA%paramId .eq. varLevelB%paramId
     if (.not.match) return
     match = varLevelA%missingValue .eq. varLevelB%missingValue
     if (.not.match) return
!
      return
!
   End function grib_hybrid_compare_varLevel
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_set_varIndex
!
!       
!
!****************************
   function grib_hybrid_set_varIndex ( varList, varName ) result ( varIndex )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character(len=*), intent(in), dimension(:)  :: varList
      character(len=*), intent(in)                :: varName
      integer(int32) :: varIndex
!
!.....Local variables
!
      logical, dimension(size(varList)) :: compare
!
!  ---
!
      compare=varList.eq.varName
      varIndex=-1_int32
      if (count(compare).ne.1_int32) return
      do varIndex=1_int32,size(varList)
         if (compare(varIndex)) exit
      end do
!
      return
!
   End function grib_hybrid_set_varIndex
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_list_varLevel
!
!       Initialize varLevel structure.
!
!****************************
   Subroutine grib_hybrid_list_varLevel ( idx, name, varLevel, iostat )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      integer(int32)          , intent(in)     :: idx
      character(len=*)        , intent(in)     :: name
      type(type_grib_varLevel), intent(inout)  :: varLevel
      integer(int32)          , intent(out), optional  :: iostat
!
!.....Local variables
!
      integer(int32)  :: iret, igrib
      logical         :: stopOnError
!
!  ---
!
      if (present(iostat)) then 
         stopOnError=.false.
         iostat=0_int32
      else
         stopOnError=.true.
      end if
!
      varLevel%shortName = name
      if (varLevel%index.eq.-1_int32) then
         varLevel%defined=.false.
         return
      else
         varLevel%defined=.true.
      end if
!
!.....Select variable index
!
      call grib_index_select( idx, 'shortName', name, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not select index for shortName'
         return
      end if
      call grib_new_from_index( idx, igrib, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not prepare data for selected index'
         return
      end if
!
!.....Examine
!
      call grib_get( igrib, 'paramId', varLevel%paramId, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varLevel : grib_get(paramId)'
         iostat=iret
         varLevel%defined=.false.
         return
      else
         varLevel%defined = .true.
      end if
      call grib_get( igrib, 'name', varLevel%longName, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varLevel : grib_get(name)'
         iostat=iret
         return
      end if
      call grib_get( igrib, 'units', varLevel%units, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varLevel : grib_get(units)'
         iostat=iret
         return
      end if
      call grib_get( igrib, 'missingValue', varLevel%missingValue, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varLevel : grib_get(missingValue)'
         iostat=iret
         return
      end if
      call grib_get( igrib, 'typeOfLevel', varLevel%typeOfLevel, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varLevel : grib_get(typeOfLevel)'
         iostat=iret
         return
      end if
!
!.....Release grib message
!
      call grib_release( igrib, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not release selected grib message'
         return
      end if
!
      return
!
   End subroutine grib_hybrid_list_varLevel
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_list_varLevel
!
!       Initialize varLevel structure.
!
!****************************
   Subroutine grib_hybrid_list_varSurf ( idx, name, varSurf, iostat )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      integer(int32)         , intent(in)     :: idx
      character(len=*)       , intent(in)     :: name
      type(type_grib_varSurf), intent(inout)  :: varSurf
      integer(int32)         , intent(out), optional  :: iostat
!
!.....Local variables
!
      integer(int32)  :: iret, igrib
      logical         :: stopOnError
!
!  ---
!
      if (present(iostat)) then 
         stopOnError=.false.
         iostat=0_int32
      else
         stopOnError=.true.
      end if
!
!  ---
!
      varSurf%shortName = name
      if (varSurf%index.eq.-1_int32) then
         varSurf%defined=.false.
         return
      else
         varSurf%defined=.true.
      end if
!
!.....Select variable index
!
      call grib_index_select( idx, 'shortName', name, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not select index for shortName'
         return
      end if
      call grib_index_select( idx, 'level', 1, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not select index for level'
         return
      end if
      call grib_new_from_index( idx, igrib, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not prepare data for selected index'
         return
      end if
!
!.....Examine variable
!
      call grib_get( igrib, 'paramId', varSurf%paramId, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varSurf : grib_get(paramId)'
         varSurf%defined = .false.
         iostat=iret
         return
      else
         varSurf%defined = .true.
      end if
      call grib_get( igrib, 'name', varSurf%longName, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varSurf : grib_get(name)'
         iostat=iret
         return
      end if
      call grib_get( igrib, 'units', varSurf%units, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varSurf : grib_get(units)'
         iostat=iret
         return
      end if
      call grib_get( igrib, 'missingValue', varSurf%missingValue, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varSurf : grib_get(missingValue)'
         iostat=iret
         return
      end if
      call grib_get( igrib, 'typeOfLevel', varSurf%typeOfLevel, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list_varSurf : grib_get(typeOfLevel)'
         iostat=iret
         return
      end if
!
!.....Release grib message
!
      call grib_release( igrib, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not release selected grib message'
         return
      end if
!
!.....Check if only 1 level
!
      call grib_index_select( idx, 'level', 2, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not select index for level'
         return
      end if
      call grib_new_from_index( idx, igrib, iret )
      varSurf%surface = iret.ne.0_int32
!
      call grib_release( igrib, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not release selected grib message'
         return
      end if
!
      return
!
   End subroutine grib_hybrid_list_varSurf
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_clear_varLevel
!
!       Initialize varLevel structure.
!
!****************************
   Subroutine grib_hybrid_clear_varSurf ( varSurf )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      type(type_grib_varSurf), intent(inout), target  :: varSurf
!
!  ---
!
     varSurf%defined = .false.
     varSurf%shortName = ''
     varSurf%longName = ''
     varSurf%units = ''
     varSurf%typeOfLevel = ''
     varSurf%index = 0_int32
     varSurf%paramId = 0_int32
     varSurf%missingValue = 9999_int32
     if (allocated(varSurf%dat)) deallocate(varSurf%dat)
     if (associated(varSurf%regular_ll)) nullify(varSurf%regular_ll)
     if (associated(varSurf%regular_ll_flip)) nullify(varSurf%regular_ll_flip)
!
      return
!
   End subroutine grib_hybrid_clear_varSurf
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_compare_varSurf
!
!       Initialize varLevel structure.
!
!****************************
   function grib_hybrid_compare_varSurf ( varSurfA, varSurfB ) result ( match )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      type(type_grib_varSurf), intent(in)  :: varSurfA, varSurfB
      logical :: match
!
!  ---
!
     match = varSurfA%defined .eqv. varSurfB%defined
     if (.not.match) return
     match = varSurfA%shortName .eq. varSurfB%shortName
     if (.not.match) return
     match = varSurfA%longName .eq. varSurfB%longName
     if (.not.match) return
     match = varSurfA%units .eq. varSurfB%units
     if (.not.match) return
     match = varSurfA%typeOfLevel .eq. varSurfB%typeOfLevel
     if (.not.match) return
     match = varSurfA%paramId .eq. varSurfB%paramId
     if (.not.match) return
     match = varSurfA%missingValue .eq. varSurfB%missingValue
     if (.not.match) return
!
      return
!
   End function grib_hybrid_compare_varSurf
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_clear
!
!       Initialize 2dfd structure.
!
!****************************
   Subroutine grib_hybrid_clear ( gribData )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      type(type_grib_hybrid), intent(inout), target  :: gribData
!
!  ---
!
!.....Grid
!
      gribData%grid%lon%count = 0_int32
      gribData%grid%lon%first = 0._double
      gribData%grid%lon%last = 0._double
      gribData%grid%lon%increment = 0._double
      gribData%grid%lon%range = 0._double
      if (allocated(gribData%grid%lon%dat)) deallocate(gribData%grid%lon%dat)
!
      gribData%grid%lat%count = 0_int32
      gribData%grid%lat%first = 0._double
      gribData%grid%lat%last = 0._double
      gribData%grid%lat%increment = 0._double
      gribData%grid%lat%range = 0._double
      if (allocated(gribData%grid%lat%dat)) deallocate(gribData%grid%lat%dat)
!
      gribData%grid%gridName = ''
      gribData%grid%gridType = ''
      gribData%grid%packingType = ''
      gribData%grid%isOctahedral = .false.
      gribData%grid%iScansNegatively = .false.
      gribData%grid%jScansPositively = .false.
      gribData%grid%jPointsAreConsecutive = .false.
      gribData%grid%distinctGridValues = .false.
      gribData%grid%nofpoints = 0_int32
      if (allocated(gribData%grid%pl)) deallocate(gribData%grid%pl)
!
!.....Ensemble numbers
!
      gribData%number%defined = .false.
      gribData%number%count = 0_int32
      gribData%number%sel = 1_int32
      if (allocated(gribData%number%list)) deallocate(gribData%number%list)
!
!.....Forecast steps
!
      gribData%step%count = 0_int32
      gribData%step%sel = 1_int32
      gribData%step%startStep = 0_int32
      gribData%step%endStep = 0_int32
      gribData%step%stepRange = 0_int32
      if (allocated(gribData%step%list)) deallocate(gribData%step%list)
!
!.....Level
!
      gribData%level%count = 0_int32
      if (allocated(gribData%level%list)) deallocate(gribData%level%list)
      if (associated(gribData%level%list_flip)) nullify(gribData%level%list_flip)
!
!.....General
!
      gribData%params = 0_int32
      if (allocated(gribData%paramList)) deallocate(gribData%paramList)
!
      gribData%isECMWF = .false.
      gribData%idx = 0_int32
      gribData%editionNumber = 0_int32
      gribData%experimentVersionNumber = 0_int32
      gribData%epoch = 0_int64
      gribData%time = ''
      gribData%centre = ''
      gribData%dataClass = ''
      gribData%dataType = ''
      gribData%dataStream = ''
!
!.....Surface variables
!
      call grib_hybrid_clear_varSurf( gribData%lnsp )
!
!.....Level variables
!
      call grib_hybrid_clear_varLevel( gribData%q )
      call grib_hybrid_clear_varLevel( gribData%t )
      call grib_hybrid_clear_varLevel( gribData%u )
      call grib_hybrid_clear_varLevel( gribData%v )
      call grib_hybrid_clear_varLevel( gribData%w )
      call grib_hybrid_clear_varLevel( gribData%z )
      call grib_hybrid_clear_varLevel( gribData%p )
!
      return
!
   End subroutine grib_hybrid_clear
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_list
!
!       List 2dfd variable information from grib file.
!
!****************************
   Subroutine grib_hybrid_list ( gribFile, gribData, release_index, mask )
!****************************
!
      use time
      implicit none
!
!.....Dummy variables
!
      character(len=*)      , intent(in)             :: gribFile
      type(type_grib_hybrid), intent(out), target    :: gribData
      logical               , intent(in) , optional  :: release_index
      logical               , intent(out), optional  :: mask 
!
!.....Local variables
!
      logical             :: lexist, lrelease, stopOnError
      integer(int32)      :: idx, ifile, igrib, iret, i, d, t, dummy
      character           :: dstr*8, tstr*4, dummy_str*40
!
!  ---
!
!.....Initialize
!
      if (present(mask)) then 
         stopOnError=.false.
         mask=.false.
      else
         stopOnError=.true.
      end if
!
      lrelease = .true.
      if (present(release_index)) lrelease=release_index
!
      if (.not.allocated(gribData%level%list).and..not.allocated(gribData%step%list).and. &
        & .not.allocated(gribData%number%list)) call grib_hybrid_clear(gribData)
!
!.....Open profile, if exists, stop otherwise.
!
      inquire( file=gribFile, exist=lexist )
      if (.not.lexist) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib file does not exist'
         call grib_hybrid_clear(gribData)
         return
      end if
!
!.....Open grib file and read first message
!
      call grib_open_file(ifile, gribFile,'r')
      call grib_new_from_file(ifile,igrib, iret)
!
!.....Get general info
!
      call grib_get( igrib, 'centre', gribData%centre, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(centre)'
         return
      end if
      if (iret.eq.0_int32) gribData%isECMWF=gribData%centre.eq.'ecmf'
      call grib_get( igrib, 'class', gribData%dataClass, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(class)'
         return
      end if
      call grib_get( igrib, 'type', gribData%dataType, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(type)'
         return
      end if
      call grib_get( igrib, 'stream', gribData%dataStream, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(stream)'
         return
      end if
      call grib_get( igrib, 'editionNumber', gribData%editionNumber, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(editionNumber)'
         return
      end if
      call grib_get( igrib, 'experimentVersionNumber', gribData%experimentVersionNumber, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(experimentVersionNumber)'
         return
      end if
!
      gribData%number%defined = gribData%dataType.eq.'ens'.or.gribData%dataStream.eq.'elda'
!
!.....Close grib file
!
      call grib_release(igrib)
      call grib_close_file(ifile)
!
!.....Open file and create index
!
      if (gribData%number%defined)then
         call grib_index_create ( idx, gribFile, 'shortName,number,level,step', iret )
      else
         call grib_index_create ( idx, gribFile, 'shortName,level,step', iret )
      end if
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not create index for grib file'
         return
      end if
!
!.....Params
!
      call grib_index_get_size( idx, 'shortName', gribData%params, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not get index size for shortName'
         return
      end if
      if (allocated(gribData%paramList).and.size(gribData%paramList,1).ne.gribData%params) &
         deallocate(gribData%paramList)
      if (.not.allocated(gribData%paramList)) allocate(gribData%paramList(gribData%params))
      call grib_index_get( idx, 'shortName', gribData%paramList, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not get list of shortNames'
         return
      end if
!
      gribData%lnsp%index=grib_hybrid_set_varIndex ( gribData%paramList, 'lnsp' )
      gribData%q%index=grib_hybrid_set_varIndex ( gribData%paramList, 'q' )
      gribData%t%index=grib_hybrid_set_varIndex ( gribData%paramList, 't' )
      gribData%u%index=grib_hybrid_set_varIndex ( gribData%paramList, 'u' )
      gribData%v%index=grib_hybrid_set_varIndex ( gribData%paramList, 'v' )
      gribData%w%index=grib_hybrid_set_varIndex ( gribData%paramList, 'w' )
      gribData%z%index=grib_hybrid_set_varIndex ( gribData%paramList, 'z' )
      gribData%p%index=grib_hybrid_set_varIndex ( gribData%paramList, 'p' )
!   
!.....Ensemble numbers
!
      if (gribData%number%defined) then
         call grib_index_get_size( idx, 'number', gribData%number%count, iret )
         if (iret.ne.0_int32) then
            if (stopOnError) stop 'Error @ grib_hybrid_list : could not get index size for ensemble number'
            return
         end if
         if (allocated(gribData%number%list).and.size(gribData%number%list,1).ne.gribData%number%count) &
            deallocate(gribData%number%list)
         if (.not.allocated(gribData%number%list)) allocate(gribData%number%list(gribData%number%count))
         call grib_index_get( idx, 'number', gribData%number%list, iret )
         if (iret.ne.0_int32) then
            if (stopOnError) stop 'Error @ grib_hybrid_list : could not get list of ensemble numbers'
            return
         end if
      end if
!
!.....Level
!
      call grib_index_get_size( idx, 'level', gribData%level%count, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not get index size for level'
         return
      end if
      if (allocated(gribData%level%list).and.size(gribData%level%list,1).ne.gribData%level%count) &
         deallocate(gribData%level%list)
      if (.not.allocated(gribData%level%list)) allocate(gribData%level%list(gribData%level%count))
      call grib_index_get( idx, 'level', gribData%level%list, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not get list of levels'
         return
      end if
      gribData%level%list_flip=>gribData%level%list(gribData%level%count:1_int32:-1_int32)
!
!.....Forecast steps
!
      call grib_index_get_size( idx, 'step', gribData%step%count, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not get index size for forecast step'
         return
      end if
      if (allocated(gribData%step%list).and.size(gribData%step%list,1).ne.gribData%step%count) &
         deallocate(gribData%step%list)
      if (.not.allocated(gribData%step%list)) allocate(gribData%step%list(gribData%step%count))
      call grib_index_get( idx, 'step', gribData%step%list, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not get list of forecast steps'
         return
      end if
!
!.....Select first step and ensemble
!
      if (gribData%number%defined) then
         call grib_index_select( idx, 'number', gribData%number%list(1), iret )
         if (iret.ne.0_int32) then
            if (stopOnError) stop 'Error @ grib_hybrid_list : could not select index 1 for ensemble number'
            return
         end if
      end if
      call grib_index_select( idx, 'step', gribData%step%list(1), iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not select index 1 for forecast step'
         return
      end if
      call grib_index_select( idx, 'level', gribData%level%list(1), iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not select index for level'
         return
      end if
      call grib_index_select( idx, 'shortName', gribData%paramList(gribData%t%index), iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not select index for shortName'
         return
      end if
!
!.....get grib block from file opened by index
!
      call grib_new_from_index( idx, igrib, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not prepare data for selected index'
         return
      end if
!
!.....Get date and time from grib file.
!
      call grib_get( igrib, 'dataDate', d, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not get dataDate'
         return
      end if
      call grib_get( igrib, 'dataTime', t, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not get dataTime'
         return
      end if
!
      write (dstr,"(i8)") d
      write (tstr,"(i4.4)") t
!
      gribData%time = dstr(1:4)//'/'//dstr(5:6)//'/'//dstr(7:)//' '//tstr(1:2)//':'//tstr(3:)//':00.000'
      call str2epoch( gribData%time, gribData%epoch )
!
!.....Get grid info
!
      call grib_get( igrib, 'gridType', gribData%grid%gridType, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(gridType)'
         return
      end if
      call grib_get( igrib, 'numberOfPoints', gribData%grid%nofpoints, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(numberOfPoints)'
         return
      end if
      call grib_get( igrib, 'gridName', gribData%grid%gridName, iret )
      call grib_get( igrib, 'gridPackingType', gribData%grid%gridName, iret )
      call grib_get( igrib, 'isOctahedral', dummy, iret )
      if (iret.eq.0_int32) gribData%grid%isOctahedral=dummy.eq.1_int32
      call grib_get( igrib, 'earthIsOblate', dummy, iret )
      if (iret.eq.0_int32) gribData%grid%earthIsOblate=dummy.eq.64_int32
!
!.....Step parameters
!
      call grib_get( igrib, 'stepUnits', gribData%step%units, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(stepUnits)'
         return
      end if
      call grib_get( igrib, 'stepRange', gribData%step%stepRange, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(stepRange)'
         return
      end if
      call grib_get( igrib, 'startStep', gribData%step%startStep, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(startStep)'
         return
      end if
      call grib_get( igrib, 'endStep', gribData%step%endStep, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(endStep)'
         return
      end if
!
!.....Get latitude info
!
      call grib_get( igrib,'latitudeOfFirstGridPointInDegrees', gribData%grid%lat%first, iret ) 
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(latitudeOfFirstGridPointInDegrees)'
         return
      end if
      call grib_get( igrib, 'latitudeOfLastGridPointInDegrees', gribData%grid%lat%last, iret ) 
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(latitudeOfLastGridPointInDegrees)'
         return
      end if
      call grib_get( igrib, 'jDirectionIncrementInDegrees', gribData%grid%lat%increment, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(jDirectionIncrementInDegrees)'
         return
      end if
      call grib_get_int ( igrib, 'Ny', gribData%grid%lat%count, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(Ny)'
         return
      end if
!
!.....Get longitude first point
!
      call grib_get( igrib, 'longitudeOfFirstGridPointInDegrees',gribData%grid%lon%first, iret ) 
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(longitudeOfFirstGridPointInDegrees)'
         return
      end if
      call grib_get( igrib, 'longitudeOfLastGridPointInDegrees', gribData%grid%lon%last, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(longitudeOfLastGridPointInDegrees)'
         return
      end if
      call grib_get( igrib, 'iDirectionIncrementInDegrees', gribData%grid%lon%increment, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(iDirectionIncrementInDegrees)'
         return
      end if
      call grib_get_int  ( igrib, 'Nx', gribData%grid%lon%count, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(Nx)'
         return
      end if
!
!.....Get latlon values
!
      gribData%grid%distinctGridValues = gribData%grid%gridType.eq.'regular_ll'
      if (gribData%grid%distinctGridValues) then
         if (allocated(gribData%grid%lat%dat).and.size(gribData%grid%lat%dat,1).ne.gribData%grid%lat%count) &
            deallocate(gribData%grid%lat%dat)
         if (.not.allocated(gribData%grid%lat%dat)) allocate(gribData%grid%lat%dat(gribData%grid%lat%count))
         call grib_get( igrib, 'distinctLatitudes' , gribData%grid%lat%dat, iret )
         if (iret.ne.0_int32) then
            if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(distinctLatitudes)'
            return
         end if
         if (allocated(gribData%grid%lon%dat).and.size(gribData%grid%lon%dat,1).ne.gribData%grid%lon%count) &
            deallocate(gribData%grid%lon%dat)
         if (.not.allocated(gribData%grid%lon%dat)) allocate(gribData%grid%lon%dat(gribData%grid%lon%count))
         call grib_get( igrib, 'distinctLongitudes', gribData%grid%lon%dat, iret ) 
         if (iret.ne.0_int32) then
            if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(distinctLongitudes)'
            return
         end if
      else
!         if (allocated(gribData%grid%lat%dat).and.size(gribData%grid%lat%dat,1).ne.gribData%grid%nofpoints) &
!            deallocate(gribData%grid%lat%dat)
!         if (.not.allocated(gribData%grid%lat%dat)) allocate(gribData%grid%lat%dat(gribData%grid%nofpoints))
!         call grib_get( igrib, 'latitudes' , gribData%grid%lat%dat, iret ) 
!         if (iret.ne.0_int32) then
!            if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(latitudes)'
!            return
!         end if
         if (allocated(gribData%grid%lat%dat).and.size(gribData%grid%lat%dat,1).ne.gribData%grid%lat%count) &
            deallocate(gribData%grid%lat%dat)
         if (.not.allocated(gribData%grid%lat%dat)) allocate(gribData%grid%lat%dat(gribData%grid%lat%count))
         call grib_get( igrib, 'distinctLatitudes' , gribData%grid%lat%dat, iret )
         if (iret.ne.0_int32) then
            if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(distinctLatitudes)'
            return
         end if
         if (allocated(gribData%grid%lon%dat).and.size(gribData%grid%lon%dat,1).ne.gribData%grid%nofpoints) &
            deallocate(gribData%grid%lon%dat)
         if (.not.allocated(gribData%grid%lon%dat)) allocate(gribData%grid%lon%dat(gribData%grid%nofpoints))
         call grib_get( igrib, 'longitudes', gribData%grid%lon%dat, iret ) 
         if (iret.ne.0_int32) then
            if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(longitudes)'
            return
         end if
         if (allocated(gribData%grid%pl).and.size(gribData%grid%pl,1).ne.gribData%grid%lat%count) &
            deallocate(gribData%grid%pl)
         if (.not.allocated(gribData%grid%pl)) allocate(gribData%grid%pl(gribData%grid%lat%count))
         call grib_get( igrib, 'pl', gribData%grid%pl, iret ) 
         if (iret.ne.0_int32) then
            if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(pl)'
            return
         end if
      end if
!
!.....Get pressure level params
!
      if (gribData%level%typeOfLevel.eq.'hybrid') then
         if (allocated(gribData%level%pv).and.size(gribData%level%pv,1).ne.gribData%level%count*2+2) &
            deallocate(gribData%level%pv)
         if (.not.allocated(gribData%level%pv)) allocate(gribData%level%pv(gribData%level%count*2+2))
         call grib_get( igrib, 'pv', gribData%level%pv, iret )
         if (iret.ne.0_int32) then
            if (stopOnError) stop 'Error @ grib_hybrid_list : grib_get(pv)'
            return
         end if
         if (allocated(gribData%level%a).and.size(gribData%level%a,1).ne.gribData%level%count) &
            deallocate(gribData%level%a)
         if (.not.allocated(gribData%level%a)) allocate(gribData%level%a(gribData%level%count))
         if (allocated(gribData%level%b).and.size(gribData%level%b,1).ne.gribData%level%count) &
            deallocate(gribData%level%b)
         if (.not.allocated(gribData%level%b)) allocate(gribData%level%b(gribData%level%count))
         do i=1_int32,gribData%level%count
            gribData%level%a(i)=( gribData%level%pv(i)+gribData%level%pv(i+1) )*0.5_double
            gribData%level%b(i)=( gribData%level%pv(gribData%level%count+i+1) + &
               & gribData%level%pv(gribData%level%count+i+2) )*0.5_double
         end do
      end if
!
!.....release select grib message
!
      call grib_release( igrib, iret )
      if (iret.ne.0_int32) then
         if (stopOnError) stop 'Error @ grib_hybrid_list : could not release selected grib message'
         return
      end if
!
!.....Get variable info
!
      call grib_hybrid_list_varSurf ( idx, gribData%paramList(gribData%lnsp%index), gribData%lnsp )
print *, gribData%lnsp%longName, gribData%lnsp%typeOfLevel
      call grib_hybrid_list_varLevel ( idx, gribData%paramList(gribData%q%index), gribData%q )
print *, gribData%q%longName, gribData%q%typeOfLevel
      call grib_hybrid_list_varLevel ( idx, gribData%paramList(gribData%t%index), gribData%t )
print *, gribData%t%longName, gribData%t%typeOfLevel
      call grib_hybrid_list_varLevel ( idx, gribData%paramList(gribData%u%index), gribData%u )
print *, gribData%u%longName, gribData%u%typeOfLevel
      call grib_hybrid_list_varLevel ( idx, gribData%paramList(gribData%v%index), gribData%v )
print *, gribData%v%longName, gribData%v%typeOfLevel
      call grib_hybrid_list_varLevel ( idx, gribData%paramList(gribData%w%index), gribData%w )
print *, gribData%w%longName, gribData%w%typeOfLevel
      call grib_hybrid_list_varLevel ( idx, gribData%paramList(gribData%z%index), gribData%z )
print *, gribData%z%longName, gribData%z%typeOfLevel
      call grib_hybrid_list_varLevel ( idx, gribData%paramList(gribData%p%index), gribData%p )
print *, gribData%p%longName, gribData%p%typeOfLevel
!
!.....Close grib file
!
      if (lrelease) then
         call grib_index_release( idx, iret )
         if (iret.ne.0_int32) then
            if (stopOnError) then
               stop 'Error @ grib_hybrid_list : could not release grib index'
            else
               gribData%idx = 0_int32
               return
            end if
         end if
      else
         gribData%idx = idx
      end if
!
      if (present(mask)) mask=.true.
!
      return
!
   End subroutine grib_hybrid_list
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine grib_hybrid_compare
!
!       Compare 2dfd structure.
!
!****************************
   Function grib_hybrid_compare ( gribDataA, gribDataB ) result ( match )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      type(type_grib_hybrid), intent(in)  :: gribDataA, gribDataB
      logical                           :: match
!
!  ---
!
!.....Grid
!
      match = gribDataA%grid%lon%count.eq.gribDataB%grid%lon%count
      if (.not.match) return
      match = abs(gribDataA%grid%lon%first-gribDataB%grid%lon%first).lt.1e-6
      if (.not.match) return
      match = abs(gribDataA%grid%lon%last-gribDataB%grid%lon%last).lt.1e-6
      if (.not.match) return
      match = abs(gribDataA%grid%lon%increment-gribDataB%grid%lon%increment).lt.1e-6
      if (.not.match) return
      match = abs(gribDataA%grid%lon%range-gribDataB%grid%lon%range).lt.1e-6
      if (.not.match) return
!
      match = gribDataA%grid%lat%count.eq.gribDataB%grid%lat%count
      if (.not.match) return
      match = abs(gribDataA%grid%lat%first-gribDataB%grid%lat%first).lt.1e-6
      if (.not.match) return
      match = abs(gribDataA%grid%lat%last-gribDataB%grid%lat%last).lt.1e-6
      if (.not.match) return
      match = abs(gribDataA%grid%lat%increment-gribDataB%grid%lat%increment).lt.1e-6
      if (.not.match) return
      match = abs(gribDataA%grid%lat%range-gribDataB%grid%lat%range).lt.1e-6
      if (.not.match) return
!
      match = gribDataA%grid%isOctahedral .eqv. gribDataB%grid%isOctahedral
      if (.not.match) return
      match = gribDataA%grid%iScansNegatively .eqv. gribDataB%grid%iScansNegatively
      if (.not.match) return
      match = gribDataA%grid%jScansPositively .eqv. gribDataB%grid%jScansPositively
      if (.not.match) return
      match = gribDataA%grid%jPointsAreConsecutive .eqv. gribDataB%grid%jPointsAreConsecutive
      if (.not.match) return
      match = gribDataA%grid%distinctGridValues .eqv. gribDataB%grid%distinctGridValues
      if (.not.match) return
      match = gribDataA%grid%gridName .eq. gribDataB%grid%gridName
      if (.not.match) return
      match = gribDataA%grid%gridType .eq. gribDataB%grid%gridType
      if (.not.match) return
      match = gribDataA%grid%packingType .eq. gribDataB%grid%packingType
      if (.not.match) return
      match = gribDataA%grid%nofpoints .eq. gribDataB%grid%nofpoints
      if (.not.match) return
!
!.....Ensemble numbers
!
      match = gribDataA%number%defined .eqv. gribDataB%number%defined
      if (.not.match) return
      match = gribDataA%number%count .eq. gribDataB%number%count
      if (.not.match) return
      match = gribDataA%number%sel .eq. gribDataB%number%sel
      if (.not.match) return
!
!.....Forecast steps
!
      match = gribDataA%step%count .eq. gribDataB%step%count
      if (.not.match) return
      match = gribDataA%step%startStep .eq. gribDataB%step%startStep
      if (.not.match) return
      match = gribDataA%step%endStep .eq. gribDataB%step%endStep
      if (.not.match) return
      match = gribDataA%step%stepRange .eq. gribDataB%step%stepRange
      if (.not.match) return
      match = gribDataA%step%sel .eq. gribDataB%step%sel
      if (.not.match) return
!      match = gribDataA%step%list .eq. gribDataB%step%list
!      if (.not.match) return
!
!.....General
!
      match = gribDataA%centre .eq. gribDataB%centre
      if (.not.match) return
      match = gribDataA%isECMWF .eqv. gribDataB%isECMWF
      if (.not.match) return
      match = gribDataA%editionNumber.eq.gribDataB%editionNumber
      if (.not.match) return
      match = gribDataA%experimentVersionNumber.eq.gribDataB%experimentVersionNumber
      if (.not.match) return
      match = gribDataA%dataClass .eq. gribDataB%dataClass
      if (.not.match) return
      match = gribDataA%dataType .eq. gribDataB%dataType
      if (.not.match) return
      match = gribDataA%dataStream .eq. gribDataB%dataStream
      if (.not.match) return
!
      return
!
   End function grib_hybrid_compare
!
! =========================================================================================


! =========================================================================================
!
!.......Function grib_hybrid_compare_gribFile_with_type
!
!       Compare 2dfd structure.
!
!****************************
   Function grib_hybrid_compare_gribFile_with_type ( gribFile, gribData ) result ( match )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character(len=*)      , intent(in)  :: gribFile
      type(type_grib_hybrid), intent(in)  :: gribData
      logical                             :: match
!
!.....Local variables
!
      type(type_grib_hybrid)  :: dummy
!
!  ---
!
      call grib_hybrid_list ( gribFile=trim(gribFile), gribData=dummy, mask=match )
      if (.not.match) return
      match = grib_hybrid_compare( gribData, dummy )
!
      return
!
   End function grib_hybrid_compare_gribFile_with_type
!
! =========================================================================================


! =========================================================================================
!
!.......subroutine grib_hybrid_print
!
!       print 2dfd from structure.
!
!****************************
   subroutine grib_hybrid_print ( gribData )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      type(type_grib_hybrid), intent(in)    :: gribData
!
!.....Local variables
!
      integer(int32) :: i
!
!  ---
!
      print "(3a)", '# // ', 'TEST', ' //'
!      print "(3a,i6,a)", '# shortname (id)    = ', trim(gribData%shortName), ' (', gribData%paramId, ')'
!      print "(2a)",      '# units             = ', trim(gribData%units)
      print "(2a)",      '# date_time         = ', gribData%time
      print "(2a)",      '# centre            = ', gribData%centre
      print "(2a)",      '# class             = ', gribData%dataClass
      print "(2a)",      '# type              = ', gribData%dataType
      print "(2a)",      '# stream            = ', gribData%dataStream
      print "(a,i10)",   '# editionNumber     = ', gribData%editionNumber
      print "(a,i10)",   '# ExpversionNumber  = ', gribData%experimentVersionNumber
!      print "(2a)",      '# type_of_level     = ', gribData%typeoflevel
!
      if (gribData%step%count.gt.1_int32) then
         print "(a,i4)", '# forecast steps    = ', gribData%step%count
         write ( *, "(a)", advance = 'no' ) '# step list         ='
         do i = 1, gribData%step%count
            write ( *, "(1x,i4)", advance = 'no' ) gribData%step%list( i )
         end do
         write ( *, "(a)" ) ''
      end if
      if (gribData%number%count.gt.1_int32) then
         print "(a,i4)", '# ensemble numbers  = ', gribData%number%count
         write ( *, "(a)", advance = 'no' ) '# ensemble list     ='
         do i = 1, gribData%number%count
            write ( *, "(1x,i4)", advance = 'no' ) gribData%number%list( i )
         end do
         write ( *, "(a)" ) ''
      end if
!
      if (len_trim(gribData%grid%gridName).gt.0_int32) &
         print "(2a)",  '# grid_gridName     = ', gribData%grid%gridName
      if (len_trim(gribData%grid%packingType).gt.0_int32) &
         print "(2a)",  '# grid_packingType  = ', gribData%grid%packingType
      print "(2a)",     '# grid_gridType     = ', gribData%grid%gridType
      print "(a,i8)",   '# grid_nofpoints    = ', gribData%grid%nofpoints
!
      print "(a,f9.4)", '# lat_first         = ', gribData%grid%lat%first
      print "(a,f9.4)", '# lat_last          = ', gribData%grid%lat%last
      print "(a,f9.4)", '# lat_increment     = ', gribData%grid%lat%increment
      print "(a,i5)"  , '# lat_N             = ', gribData%grid%lat%count
!
      print "(a,f9.4)", '# lon_first         = ', gribData%grid%lon%first
      print "(a,f9.4)", '# lon_last          = ', gribData%grid%lon%last
      print "(a,f9.4)", '# lon_increment     = ', gribData%grid%lon%increment
      print "(a,i5)"  , '# lon_N             = ', gribData%grid%lon%count
!
      return
!
   End subroutine grib_hybrid_print
!
! =========================================================================================


! =========================================================================================
!
!.......Reshape using pointers
!
!****************************
  Function grib_hybrid_reshape_ll(array, shape_) result(aptr)
!****************************
!
    use iso_c_binding, only: C_LOC, C_F_POINTER
    implicit none
!
    ! Pass in the array as an array of fixed size so that there
    ! is no array descriptor associated with it. This means we
    ! can get a pointer to the location of the data using C_LOC
    real(single), dimension(:), intent(in), target :: array
    integer(int32), intent(in), dimension(2) :: shape_
    real(single), dimension(:,:), pointer :: aptr

    ! Use C_LOC to get the start location of the array data, and
    ! use C_F_POINTER to turn this into a fortran pointer (aptr).
    ! Note that we need to specify the shape of the pointer using an
    ! integer array.
    call C_F_POINTER(C_LOC(array), aptr, shape_)
!
  End function grib_hybrid_reshape_ll
!
! =========================================================================================


! =========================================================================================
!
!.......Reshape using pointers
!
!****************************
  Function grib_hybrid_reshape_lll(array, shape_) result(aptr)
!****************************
!
    use iso_c_binding, only: C_LOC, C_F_POINTER
    implicit none
!
    ! Pass in the array as an array of fixed size so that there
    ! is no array descriptor associated with it. This means we
    ! can get a pointer to the location of the data using C_LOC
    real(single), dimension(:,:), intent(in), target :: array
    integer(int32), intent(in), dimension(3) :: shape_
    real(single), dimension(:,:,:), pointer :: aptr

    ! Use C_LOC to get the start location of the array data, and
    ! use C_F_POINTER to turn this into a fortran pointer (aptr).
    ! Note that we need to specify the shape of the pointer using an
    ! integer array.
    call C_F_POINTER(C_LOC(array), aptr, shape_)
!
  End function grib_hybrid_reshape_lll
!
! =========================================================================================
!
End module grib_hybrid
