!****************************
Program test 
!****************************
 
   use io_types
   use grib_hybrid
!
   implicit none
! 
!..declare varibles
!
   logical              :: verb
   type(type_grib_hybrid), target :: grib
   character(len=999)   :: file1
   integer(int32)  :: i, i0, j, stp

!   file1='/Users/psmets/Documents/KNMI/research/_studies/DPRK/raw_grib_data/fc_dprk_20160106_00.grib'
   file1='/Users/psmets/Documents/KNMI/research/_studies/DPRK/raw_grib_data/an_dprk_20160909_00.grib'
!   file1='/Users/psmets/Documents/KNMI/research/_studies/DPRK/raw_grib_data/eda_dprk_20160106_00.grib'
   call grib_hybrid_list ( gribFile=trim(file1), gribData=grib )
   call grib_hybrid_print( gribData=grib )
   call grib_hybrid_clear( gribData=grib )

!****************************
End program test
!****************************
