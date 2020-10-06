! =========================================================================================
!
! MODULE TOOLS_NETCDF 
!
! --> subroutines and functions for netcdf grd files.
!
! =========================================================================================
!
MODULE TOOLS_NETCDF
!
   USE NETCDF
   IMPLICIT NONE
!
   character (len = *), parameter, private :: UNITS      = "units"
   character (len = *), parameter, private :: TIME_UNITS = "-"
   character (len = *), parameter, private :: EPOC_UNITS = "milliseconds"
   character (len = *), parameter, private :: STEP_UNITS = "hours"
   character (len = *), parameter, private :: PRES_UNITS = "Pa"
   character (len = *), parameter, private :: TEMP_UNITS = "K"
   character (len = *), parameter, private :: DENS_UNITS = "kg/m^3"
   character (len = *), parameter, private :: WIND_UNITS = "m/s"
   character (len = *), parameter, private :: ALT__UNITS  = "m"
   character (len = *), parameter, private :: TOPO_UNITS = "m"
   character (len = *), parameter, private :: LAT__UNITS = "degrees north"
   character (len = *), parameter, private :: LON__UNITS = "degrees east"
   character (len = *), parameter, private :: LVL__UNITS = "-"
   character (len = *), parameter, private :: FREQ_UNITS = "Hz"
   character (len = *), parameter, private :: MBAP_UNITS = "Pa"
   character (len = *), parameter, private :: MBAD_UNITS = "dB"
   character (len = *), parameter, private :: MBSP_UNITS = "Pa^2/Hz"
   character (len = *), parameter, private :: MBSD_UNITS = "dB/Hz"
   character (len = *), parameter, private :: SWH__UNITS = "m"
   character (len = *), parameter, private :: TIME_SHORT_NAME  = "time"
   character (len = *), parameter, private :: EPOC_SHORT_NAME  = "epoch"
   character (len = *), parameter, private :: STEP_SHORT_NAME  = "step"
   character (len = *), parameter, private :: PRES_SHORT_NAME  = "p"
   character (len = *), parameter, private :: TEMP_SHORT_NAME  = "t"
   character (len = *), parameter, private :: DENS_SHORT_NAME  = "d"
   character (len = *), parameter, private :: LAT__SHORT_NAME  = "lat"
   character (len = *), parameter, private :: LON__SHORT_NAME  = "lon"
   character (len = *), parameter, private :: LVL__SHORT_NAME  = "lvl"
   character (len = *), parameter, private :: ALT__SHORT_NAME  = "z"
   character (len = *), parameter, private :: TOPO_SHORT_NAME  = "o"
   character (len = *), parameter, private :: W_X__SHORT_NAME  = "u"
   character (len = *), parameter, private :: W_Y__SHORT_NAME  = "v"
   character (len = *), parameter, private :: W_Z__SHORT_NAME  = "w"
   character (len = *), parameter, private :: FREQ_SHORT_NAME  = "f"
   character (len = *), parameter, private :: MB_A_SHORT_NAME  = "mb_ampl"
   character (len = *), parameter, private :: MB_S_SHORT_NAME  = "mb_spec"
   character (len = *), parameter, private :: SWH__SHORT_NAME  = "swh"
   character (len = *), parameter, private :: LONG_NAME       = "long_name"
   character (len = *), parameter, private :: TIME_LONG_NAME  = "Time"
   character (len = *), parameter, private :: EPOC_LONG_NAME  = "Epoch"
   character (len = *), parameter, private :: STEP_LONG_NAME  = "Forecast step hours"
   character (len = *), parameter, private :: PRES_LONG_NAME  = "Pressure"
   character (len = *), parameter, private :: TEMP_LONG_NAME  = "Temperature"
   character (len = *), parameter, private :: DENS_LONG_NAME  = "Density"
   character (len = *), parameter, private :: LAT__LONG_NAME  = "Latitude"
   character (len = *), parameter, private :: LON__LONG_NAME  = "Longitude"
   character (len = *), parameter, private :: LVL__LONG_NAME  = "Level"
   character (len = *), parameter, private :: ALT__LONG_NAME  = "Altitude"
   character (len = *), parameter, private :: TOPO_LONG_NAME  = "Orography"
   character (len = *), parameter, private :: W_X__LONG_NAME  = "Wind_zonal"
   character (len = *), parameter, private :: W_Y__LONG_NAME  = "Wind_meridional"
   character (len = *), parameter, private :: W_Z__LONG_NAME  = "Wind_vertical"
   character (len = *), parameter, private :: FREQ_LONG_NAME  = "Frequency"
   character (len = *), parameter, private :: MB_A_LONG_NAME  = "Microbarom source strength amplitude"
   character (len = *), parameter, private :: MB_S_LONG_NAME  = "Microbarom source strength spectrum"
   character (len = *), parameter, private :: SWH__LONG_NAME  = "Significant wave height"
!
   CONTAINS
!
! =========================================================================================
!
!.......FUNCTION ISNETCDF 
!
!       Check if file is of netcdf type
!
!****************************
   LOGICAL FUNCTION ISNETCDF ( file1  )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN )  :: file1
!
!.....Local variables
!
      integer(kind=4)  :: ncid, ierr
!
!  ---
!
!.....Open the netCDF file.
!
      ierr = nf90_open ( file1, 0, ncid )
!
!.....If no error, file is netcdf
!
      ISNETCDF = ierr .EQ. NF90_noerr
!
!.....If open succesfull, close file again
!
      IF ( ISNETCDF ) CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END FUNCTION ISNETCDF
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE NETCDF_WRITE_MB_AMPL
!
!       Write 2d grid to netcdf grd file
!
!****************************
   SUBROUTINE NETCDF_WRITE_MB_AMPL ( file1, xpos, ypos, dat, epoch, ldb  )
!****************************
!
      USE TIME
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * )            , INTENT ( IN )  :: file1
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN )  :: xpos
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN )  :: ypos
      real, DIMENSION ( SIZE ( xpos ), &
       & SIZE ( ypos ) )               , INTENT ( IN )  :: dat
      integer(kind=8)                      , INTENT ( IN )  :: epoch
      LOGICAL                          , INTENT ( IN )  :: ldb
!
!.....Local variables
!
      CHARACTER                    :: tstr*23
      integer(kind=4)                  :: ncid, x_dimid, y_dimid, timeDimID
      integer(kind=4)                  :: x_varid, y_varid, varid, timeVarID
      integer(kind=4), DIMENSION ( 2 ) :: dimids
!
!  ---
!
!.....Create the netCDF file.
!
      CALL CHECK ( nf90_create ( file1, NF90_CLOBBER, ncid ) )
!
!.....Define the dimensions.
!
      CALL CHECK ( nf90_def_dim ( ncid, LON__SHORT_NAME, SIZE ( xpos ), x_dimid ) )
      CALL CHECK ( nf90_def_dim ( ncid, LAT__SHORT_NAME, SIZE ( ypos ), y_dimid ) )
!
!.....Define coordinate variables
!
      CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, x_dimid, x_varid ) )
      CALL CHECK ( nf90_put_att ( ncid, x_varid, LONG_NAME, LON__LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, x_varid, UNITS    , LON__UNITS ) )
!
      CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, y_dimid, y_varid ) )
      CALL CHECK ( nf90_put_att ( ncid, y_varid, LONG_NAME, LAT__LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, y_varid, UNITS    , LAT__UNITS ) )
!
      dimids = (/ x_dimid, y_dimid /)
!
!.....Define time information
!
      CALL CHECK ( nf90_def_dim ( ncid, TIME_SHORT_NAME, len = 23, dimid = timeDimID ) )
      CALL CHECK ( nf90_def_var ( ncid, TIME_SHORT_NAME, NF90_CHAR, timeDimID, timeVarID ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, LONG_NAME, TIME_LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, UNITS, TIME_UNITS ) )
!
      CALL EPOCH2STR ( epoch, tstr )
!
!.....Define variable
!
      CALL CHECK ( nf90_def_var ( ncid, MB_A_SHORT_NAME, NF90_REAL, dimids, varid ) )
      CALL CHECK ( nf90_put_att ( ncid, varid, LONG_NAME, MB_A_LONG_NAME ) )
      IF ( .NOT. ldb ) CALL CHECK ( nf90_put_att ( ncid, varid, UNITS, MBAP_UNITS ) )
      IF ( ldb )       CALL CHECK ( nf90_put_att ( ncid, varid, UNITS, MBAD_UNITS ) )
!
!.....End definitions
!
      CALL CHECK ( nf90_enddef ( ncid ) )
!
!.....Write Data
!
      CALL CHECK ( nf90_put_var ( ncid, timeVarID, tstr ) )
!
      CALL CHECK ( nf90_put_var ( ncid, x_varid, REAL ( xpos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, y_varid, REAL ( ypos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid  , dat ) )
!
!.....Close netcdf file
!
      CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END SUBROUTINE NETCDF_WRITE_MB_AMPL
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE NETCDF_WRITE_SWH
!
!       Write 2d grid to netcdf grd file
!
!****************************
   SUBROUTINE NETCDF_WRITE_SWH ( file1, xpos, ypos, dat, epoch  )
!****************************
!
      USE TIME
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * )            , INTENT ( IN )  :: file1
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN )  :: xpos
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN )  :: ypos
      real, DIMENSION ( SIZE ( xpos ), &
       & SIZE ( ypos ) )               , INTENT ( IN )  :: dat
      integer(kind=8)                      , INTENT ( IN )  :: epoch
!
!.....Local variables
!
      CHARACTER                    :: tstr*23
      integer(kind=4)                  :: ncid, x_dimid, y_dimid, timeDimID
      integer(kind=4)                  :: x_varid, y_varid, varid, timeVarID
      integer(kind=4), DIMENSION ( 2 ) :: dimids
!
!  ---
!
!.....Create the netCDF file.
!
      CALL CHECK ( nf90_create ( file1, NF90_CLOBBER, ncid ) )
!
!.....Define the dimensions.
!
      CALL CHECK ( nf90_def_dim ( ncid, LON__SHORT_NAME, SIZE ( xpos ), x_dimid ) )
      CALL CHECK ( nf90_def_dim ( ncid, LAT__SHORT_NAME, SIZE ( ypos ), y_dimid ) )
!
!.....Define coordinate variables
!
      CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, x_dimid, x_varid ) )
      CALL CHECK ( nf90_put_att ( ncid, x_varid, LONG_NAME, LON__LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, x_varid, UNITS    , LON__UNITS ) )
!
      CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, y_dimid, y_varid ) )
      CALL CHECK ( nf90_put_att ( ncid, y_varid, LONG_NAME, LAT__LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, y_varid, UNITS    , LAT__UNITS ) )
!
      dimids = (/ x_dimid, y_dimid /)
!
!.....Define time information
!
      CALL CHECK ( nf90_def_dim ( ncid, TIME_SHORT_NAME, len = 23, dimid = timeDimID ) )
      CALL CHECK ( nf90_def_var ( ncid, TIME_SHORT_NAME, NF90_CHAR, timeDimID, timeVarID ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, LONG_NAME, TIME_LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, UNITS, TIME_UNITS ) )
!
      CALL EPOCH2STR ( epoch, tstr )
!
!.....Define variable
!
      CALL CHECK ( nf90_def_var ( ncid, SWH__SHORT_NAME, NF90_REAL, dimids, varid ) )
      CALL CHECK ( nf90_put_att ( ncid, varid, LONG_NAME, SWH__LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, varid, UNITS, SWH__UNITS ) )
!
!.....End definitions
!
      CALL CHECK ( nf90_enddef ( ncid ) )
!
!.....Write Data
!
      CALL CHECK ( nf90_put_var ( ncid, timeVarID, tstr ) )
!
      CALL CHECK ( nf90_put_var ( ncid, x_varid, REAL ( xpos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, y_varid, REAL ( ypos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid  , dat ) )
!
!.....Close netcdf file
!
      CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END SUBROUTINE NETCDF_WRITE_SWH
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE NETCDF_WRITE_MB_SPEC
!
!       Write 3d grid to netcdf grd file
!
!****************************
   SUBROUTINE NETCDF_WRITE_MB_SPEC ( file1, xpos, ypos, zpos, dat, epoch, ldb  )
!****************************
!
      USE TIME
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * )            , INTENT ( IN )  :: file1
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN )  :: xpos
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN )  :: ypos
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN )  :: zpos
      real, DIMENSION ( SIZE ( xpos ), &
       & SIZE ( ypos ), SIZE ( zpos ) ), INTENT ( IN )  :: dat
      integer(kind=8)                      , INTENT ( IN )  :: epoch
      LOGICAL                          , INTENT ( IN )  :: ldb
!
!.....Local variables
!
      CHARACTER                    :: tstr*23
      integer(kind=4)                  :: ncid, x_dimid, y_dimid, z_dimid, timeDimID
      integer(kind=4)                  :: x_varid, y_varid, z_varid, varid, timeVarID
      integer(kind=4), DIMENSION ( 3 ) :: dimids
!
!  ---
!
!.....Create the netCDF file.
!
      CALL CHECK ( nf90_create ( file1, NF90_CLOBBER, ncid ) )
!
!.....Define the dimensions.
!
      CALL CHECK ( nf90_def_dim ( ncid, LON__SHORT_NAME, SIZE ( xpos ), x_dimid ) )
      CALL CHECK ( nf90_def_dim ( ncid, LAT__SHORT_NAME, SIZE ( ypos ), y_dimid ) )
      CALL CHECK ( nf90_def_dim ( ncid, FREQ_SHORT_NAME, SIZE ( zpos ), z_dimid ) )
!
!.....Define coordinate variables
!
      CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, x_dimid, x_varid ) )
      CALL CHECK ( nf90_put_att ( ncid, x_varid, LONG_NAME, LON__LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, x_varid, UNITS    , LON__UNITS ) )
!
      CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, y_dimid, y_varid ) )
      CALL CHECK ( nf90_put_att ( ncid, y_varid, LONG_NAME, LAT__LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, y_varid, UNITS    , LAT__UNITS ) )
!
      CALL CHECK ( nf90_def_var ( ncid, FREQ_SHORT_NAME, NF90_REAL, z_dimid, z_varid ) )
      CALL CHECK ( nf90_put_att ( ncid, z_varid, LONG_NAME, FREQ_LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, z_varid, UNITS    , FREQ_UNITS ) )
!
      dimids = (/ x_dimid, y_dimid, z_dimid /)
!
!.....Define time information
!
      CALL CHECK ( nf90_def_dim ( ncid, TIME_SHORT_NAME, len = 23, dimid = timeDimID ) )
      CALL CHECK ( nf90_def_var ( ncid, TIME_SHORT_NAME, NF90_CHAR, timeDimID, timeVarID ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, LONG_NAME, TIME_LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, UNITS, TIME_UNITS ) )
!
      CALL EPOCH2STR ( epoch, tstr )
!
!.....Define variable
!
      CALL CHECK ( nf90_def_var ( ncid, MB_S_SHORT_NAME, NF90_REAL, dimids, varid ) )
      CALL CHECK ( nf90_put_att ( ncid, varid, LONG_NAME, MB_S_LONG_NAME ) )
      IF ( .NOT. ldb ) CALL CHECK ( nf90_put_att ( ncid, varid, UNITS, MBSP_UNITS ) )
      IF ( ldb )       CALL CHECK ( nf90_put_att ( ncid, varid, UNITS, MBSD_UNITS ) )
!
!.....End definitions
!
      CALL CHECK ( nf90_enddef ( ncid ) )
!
!.....Write Data
!
      CALL CHECK ( nf90_put_var ( ncid, timeVarID, tstr ) )
!
      CALL CHECK ( nf90_put_var ( ncid, x_varid, REAL ( xpos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, y_varid, REAL ( ypos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, z_varid, REAL ( zpos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid  , dat ) )
!
!.....Close netcdf file
!
      CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END SUBROUTINE NETCDF_WRITE_MB_SPEC
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE WRITE_NETCDF_ATMO_RAW
!
!       Write 3d grid to netcdf grd file
!
!****************************
   SUBROUTINE NETCDF_WRITE_ATMO_RAW_LIM ( file1, xv, yv, zv, h, t, u, v, o, epoch, fcstep, &
      & dataType, typeOfLevel, gridType, stream  )
!****************************
!
      USE TIME
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * )            , INTENT ( IN )  :: file1
      CHARACTER ( LEN = 20 )           , INTENT ( IN )  :: dataType, typeOfLevel, gridType, stream
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN )  :: xv, yv, zv
      DOUBLE PRECISION, DIMENSION ( SIZE ( xv ), SIZE ( yv ), SIZE ( zv ) ), &
         &                               INTENT ( IN )  :: t, u, v, h
      DOUBLE PRECISION, DIMENSION ( SIZE ( xv ), SIZE ( yv ) ), &
         &                               INTENT ( IN )  :: o
      integer(kind=8)                      , INTENT ( IN )  :: epoch
      integer(kind=4)                      , INTENT ( IN )  :: fcstep
!
!.....Local variables
!
      LOGICAL, DIMENSION ( 3 )     :: dims
!
      CHARACTER                    :: tstr*23
!
      integer(kind=4)                  :: nofdims, dimid1, ncid, nofx, nofy, nofz, &
                                   &  x_dimid, y_dimid, z_dimid, timeDimID, timeVarID, headerMaxLen, &
                                   & varid_x, varid_y, varid_z, varid_t, varid_u, varid_v, &
                                   & varid_h, varid_o, stepVarID, varid_dataType, varid_typeOfLevel, &
                                   & varid_gridType, varid_stream, varid_source
      integer(kind=4), DIMENSION ( 2 ) :: dimid2
      integer(kind=4), DIMENSION ( 3 ) :: dimid3
!
!  ---
!
!.....Create the netCDF file.
!
      CALL CHECK ( nf90_create ( file1, NF90_CLOBBER, ncid ) )
!
!.....Write time information
!
      CALL CHECK ( nf90_def_dim ( ncid, 'timeMaxLen', len = 23, dimid = timeDimID ) )
      CALL CHECK ( nf90_def_var ( ncid, TIME_SHORT_NAME, NF90_CHAR, timeDimID, timeVarID ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, LONG_NAME, TIME_LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, UNITS, TIME_UNITS ) )
!
      CALL EPOCH2STR ( epoch, tstr )
!
      CALL CHECK ( nf90_def_var ( ncid, STEP_SHORT_NAME, NF90_INT, stepVarID ) )
      CALL CHECK ( nf90_put_att ( ncid, stepVarID, LONG_NAME, STEP_LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, stepVarID, UNITS, STEP_UNITS ) )
!
!.....Write grib information
!
      CALL CHECK ( nf90_def_dim ( ncid, 'headerMaxLen', len = 20, dimid = headerMaxLen ) )
      CALL CHECK ( nf90_def_var ( ncid, 'dataType', NF90_CHAR, headerMaxLen, varid_dataType ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_dataType, LONG_NAME, 'Data type' ) )
!
      CALL CHECK ( nf90_def_var ( ncid, 'typeOfLevel', NF90_CHAR, headerMaxLen, varid_typeOfLevel ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_typeOfLevel, LONG_NAME, 'Type of level' ) )
!
      CALL CHECK ( nf90_def_var ( ncid, 'gridType', NF90_CHAR, headerMaxLen, varid_gridType ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_gridType, LONG_NAME, 'Grid type' ) )
!
      CALL CHECK ( nf90_def_var ( ncid, 'stream', NF90_CHAR, headerMaxLen, varid_stream ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_stream, LONG_NAME, 'Stream' ) )
!
      CALL CHECK ( nf90_def_var ( ncid, 'source', NF90_CHAR, headerMaxLen, varid_source ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_stream, LONG_NAME, 'source' ) )
!
!.....Get dimensions
!
      nofx = SIZE ( xv, 1 )
      nofy = SIZE ( yv, 1 )
      nofz = SIZE ( zv, 1 )
!
      dims( 1 ) = nofx .GT. 1
      dims( 2 ) = nofy .GT. 1
      dims( 3 ) = nofz .GT. 1
!
      nofdims = COUNT ( dims )
!
!.....3D
!
      IF ( nofdims .EQ. 3 ) THEN
!
!........Define the dimensions.
!
         CALL CHECK ( nf90_def_dim ( ncid, LON__SHORT_NAME, nofx, x_dimid ) )
         CALL CHECK ( nf90_def_dim ( ncid, LAT__SHORT_NAME, nofy, y_dimid ) )
         CALL CHECK ( nf90_def_dim ( ncid, LVL__SHORT_NAME, nofz, z_dimid ) )
!
!........Define coordinate variables
!
         CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, x_dimid, varid_x ) )
         CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, y_dimid, varid_y ) )
         CALL CHECK ( nf90_def_var ( ncid, LVL__SHORT_NAME, NF90_REAL, z_dimid, varid_z ) )
         dimid2 = (/ x_dimid, y_dimid /)
         dimid3 = (/ x_dimid, y_dimid, z_dimid /)
!
!........Define variable
!
         CALL CHECK ( nf90_def_var ( ncid, TEMP_SHORT_NAME, NF90_REAL, dimid3, varid_t ) )
         CALL CHECK ( nf90_def_var ( ncid, W_X__SHORT_NAME, NF90_REAL, dimid3, varid_u ) )
         CALL CHECK ( nf90_def_var ( ncid, W_Y__SHORT_NAME, NF90_REAL, dimid3, varid_v ) )
         CALL CHECK ( nf90_def_var ( ncid, ALT__SHORT_NAME, NF90_REAL, dimid3, varid_h ) )
         CALL CHECK ( nf90_def_var ( ncid, TOPO_SHORT_NAME, NF90_REAL, dimid2, varid_o ) )
!
!.....2D
!
      ELSE IF ( nofdims .EQ. 2 ) THEN
!
!........Define the dimensions.
!
         CALL CHECK ( nf90_def_dim ( ncid, LVL__SHORT_NAME, nofz, z_dimid ) )
         CALL CHECK ( nf90_def_var ( ncid, LVL__SHORT_NAME, NF90_REAL, z_dimid, varid_z ) )
         IF ( dims( 1 ) ) THEN
            CALL CHECK ( nf90_def_dim ( ncid, LON__SHORT_NAME, nofx, x_dimid ) )
            CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, x_dimid, varid_x ) )
            CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, varid_y ) )
            dimid1 = x_dimid
            dimid2 = (/ x_dimid, z_dimid /)
         ELSE
            CALL CHECK ( nf90_def_dim ( ncid, LAT__SHORT_NAME, nofy, y_dimid ) )
            CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, y_dimid, varid_y ) )
            CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, varid_x ) )
            dimid1 = y_dimid
            dimid2 = (/ y_dimid, z_dimid /)
         END IF
!
!........Define variable
!
         CALL CHECK ( nf90_def_var ( ncid, TEMP_SHORT_NAME, NF90_REAL, dimid2, varid_t ) )
         CALL CHECK ( nf90_def_var ( ncid, W_X__SHORT_NAME, NF90_REAL, dimid2, varid_u ) )
         CALL CHECK ( nf90_def_var ( ncid, W_Y__SHORT_NAME, NF90_REAL, dimid2, varid_v ) )
         CALL CHECK ( nf90_def_var ( ncid, ALT__SHORT_NAME, NF90_REAL, dimid2, varid_h ) )
         CALL CHECK ( nf90_def_var ( ncid, TOPO_SHORT_NAME, NF90_REAL, dimid1, varid_o ) )
!
!.....1D
!
      ELSE IF ( nofdims .EQ. 1 ) THEN
!
!........Define the dimensions.
!
         CALL CHECK ( nf90_def_dim ( ncid, LVL__SHORT_NAME, nofz, z_dimid ) )
!
!........Define coordinate variables
!
         CALL CHECK ( nf90_def_var ( ncid, LVL__SHORT_NAME, NF90_REAL, z_dimid, varid_z ) )
         CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, varid_x ) )
         CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, varid_y ) )
!
!........Define variable
!
         CALL CHECK ( nf90_def_var ( ncid, TEMP_SHORT_NAME, NF90_REAL, z_dimid, varid_t ) )
         CALL CHECK ( nf90_def_var ( ncid, W_X__SHORT_NAME, NF90_REAL, z_dimid, varid_u ) )
         CALL CHECK ( nf90_def_var ( ncid, W_Y__SHORT_NAME, NF90_REAL, z_dimid, varid_v ) )
         CALL CHECK ( nf90_def_var ( ncid, ALT__SHORT_NAME, NF90_REAL, z_dimid, varid_h ) )
         CALL CHECK ( nf90_def_var ( ncid, TOPO_SHORT_NAME, NF90_REAL, varid_o ) )
!
      END IF
!
!.....Assign units to variables
!
      CALL CHECK ( nf90_put_att(ncid, varid_x, UNITS, LON__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_y, UNITS, LAT__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_z, UNITS, LVL__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_t, UNITS, TEMP_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_u, UNITS, WIND_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_v, UNITS, WIND_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_h, UNITS, ALT__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_o, UNITS, TOPO_UNITS ) )
!
!.....Assign long_names to variables
!
      CALL CHECK ( nf90_put_att(ncid, varid_x, LONG_NAME, LON__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_y, LONG_NAME, LAT__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_z, LONG_NAME, LVL__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_t, LONG_NAME, TEMP_LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_u, LONG_NAME, W_X__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_v, LONG_NAME, W_Y__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_h, LONG_NAME, ALT__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_o, LONG_NAME, TOPO_LONG_NAME ) )
!
!.....end define mode
!
      CALL CHECK ( nf90_enddef ( ncid ) )
!
!.....Write Data
!
      CALL CHECK ( nf90_put_var ( ncid, varid_x, REAL ( xv ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_y, REAL ( yv ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_z, REAL ( zv ) ) )
!
      CALL CHECK ( nf90_put_var ( ncid, timeVarID, tstr ) )
      CALL CHECK ( nf90_put_var ( ncid, stepVarID, fcstep ) )
!
      CALL CHECK ( nf90_put_var ( ncid, varid_dataType, dataType ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_typeOfLevel, typeOfLevel ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_gridType, gridType ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_stream, stream ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_source, 'ECMWF' ) )
!
      IF ( nofdims .EQ. 3 ) THEN
         CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_h, REAL ( h ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_o, REAL ( o ) ) )
      ELSE IF ( nofdims .EQ. 2 ) THEN
         IF ( dims( 1 ) ) THEN
            CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_h, REAL ( h(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_o, REAL ( o(:,1) ) ) )
         ELSE
            CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_h, REAL ( h(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_o, REAL ( o(1,:) ) ) )
         END IF
      ELSE
         CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_h, REAL ( h(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_o, REAL ( o(1,1) ) ) )
      END IF
!
!.....Close netcdf file
!
      CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END SUBROUTINE NETCDF_WRITE_ATMO_RAW_LIM
!
! =========================================================================================



! =========================================================================================
!
!.......SUBROUTINE WRITE_NETCDF_ATMO_RAW
!
!       Write 3d grid to netcdf grd file
!
!****************************
   SUBROUTINE NETCDF_WRITE_ATMO_RAW ( file1, xv, yv, zv, h, t, u, v, w, p, d, o, epoch, fcstep, &
      & dataType, typeOfLevel, gridType, stream  )
!****************************
!
      USE TIME
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * )            , INTENT ( IN )  :: file1
      CHARACTER ( LEN = 20 )           , INTENT ( IN )  :: dataType, typeOfLevel, gridType, stream
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN )  :: xv, yv, zv
      DOUBLE PRECISION, DIMENSION ( SIZE ( xv ), SIZE ( yv ), SIZE ( zv ) ), &
         &                               INTENT ( IN )  :: t, u, v, w, p, d, h
      DOUBLE PRECISION, DIMENSION ( SIZE ( xv ), SIZE ( yv ) ), &
         &                               INTENT ( IN )  :: o
      integer(kind=8)                      , INTENT ( IN )  :: epoch
      integer(kind=4)                      , INTENT ( IN )  :: fcstep
!
!.....Local variables
!
      LOGICAL, DIMENSION ( 3 )     :: dims
!
      CHARACTER                    :: tstr*23
!
      integer(kind=4)                  :: nofdims, dimid1, ncid, nofx, nofy, nofz, headerMaxLen, &
                                   &  x_dimid, y_dimid, z_dimid, timeDimID, timeVarID, &
                                   & varid_x, varid_y, varid_z, varid_t, varid_u, varid_v, &
                                   & varid_w, varid_p, varid_d, varid_h, varid_o, stepVarID, varid_source, &
                                   & varid_dataType, varid_typeOfLevel, varid_gridType, varid_stream
      integer(kind=4), DIMENSION ( 2 ) :: dimid2
      integer(kind=4), DIMENSION ( 3 ) :: dimid3
!
!  ---
!
!.....Create the netCDF file.
!
      CALL CHECK ( nf90_create ( file1, NF90_CLOBBER, ncid ) )
!
!.....Write time information
!
      CALL CHECK ( nf90_def_dim ( ncid, 'timeMaxLen', len = 23, dimid = timeDimID ) )
      CALL CHECK ( nf90_def_var ( ncid, TIME_SHORT_NAME, NF90_CHAR, timeDimID, timeVarID ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, LONG_NAME, TIME_LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, timeVarID, UNITS, TIME_UNITS ) )
!
      CALL EPOCH2STR ( epoch, tstr )
!
      CALL CHECK ( nf90_def_var ( ncid, STEP_SHORT_NAME, NF90_INT, stepVarID ) )
      CALL CHECK ( nf90_put_att ( ncid, stepVarID, LONG_NAME, STEP_LONG_NAME ) )
      CALL CHECK ( nf90_put_att ( ncid, stepVarID, UNITS, STEP_UNITS ) )
!
!.....Write grib information
!
      CALL CHECK ( nf90_def_dim ( ncid, 'headerMaxLen', len = 20, dimid = headerMaxLen ) )
      CALL CHECK ( nf90_def_var ( ncid, 'dataType', NF90_CHAR, headerMaxLen, varid_dataType ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_dataType, LONG_NAME, 'Data type' ) )
!
      CALL CHECK ( nf90_def_var ( ncid, 'typeOfLevel', NF90_CHAR, headerMaxLen, varid_typeOfLevel ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_typeOfLevel, LONG_NAME, 'Type of level' ) )
!
      CALL CHECK ( nf90_def_var ( ncid, 'gridType', NF90_CHAR, headerMaxLen, varid_gridType ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_gridType, LONG_NAME, 'Grid type' ) )
!
      CALL CHECK ( nf90_def_var ( ncid, 'stream', NF90_CHAR, headerMaxLen, varid_stream ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_stream, LONG_NAME, 'Stream' ) )
!
      CALL CHECK ( nf90_def_var ( ncid, 'source', NF90_CHAR, headerMaxLen, varid_source ) )
      CALL CHECK ( nf90_put_att ( ncid, varid_stream, LONG_NAME, 'source' ) )
!
!.....Get dimensions
!
      nofx = SIZE ( xv, 1 )
      nofy = SIZE ( yv, 1 )
      nofz = SIZE ( zv, 1 )
!
      dims( 1 ) = nofx .GT. 1
      dims( 2 ) = nofy .GT. 1
      dims( 3 ) = nofz .GT. 1
!
      nofdims = COUNT ( dims )
!
!.....3D
!
      IF ( nofdims .EQ. 3 ) THEN
!
!........Define the dimensions.
!
         CALL CHECK ( nf90_def_dim ( ncid, LON__SHORT_NAME, nofx, x_dimid ) )
         CALL CHECK ( nf90_def_dim ( ncid, LAT__SHORT_NAME, nofy, y_dimid ) )
         CALL CHECK ( nf90_def_dim ( ncid, LVL__SHORT_NAME, nofz, z_dimid ) )
!
!........Define coordinate variables
!
         CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, x_dimid, varid_x ) )
         CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, y_dimid, varid_y ) )
         CALL CHECK ( nf90_def_var ( ncid, LVL__SHORT_NAME, NF90_REAL, z_dimid, varid_z ) )
         dimid2 = (/ x_dimid, y_dimid /)
         dimid3 = (/ x_dimid, y_dimid, z_dimid /)
!
!........Define variable
!
         CALL CHECK ( nf90_def_var ( ncid, TEMP_SHORT_NAME, NF90_REAL, dimid3, varid_t ) )
         CALL CHECK ( nf90_def_var ( ncid, W_X__SHORT_NAME, NF90_REAL, dimid3, varid_u ) )
         CALL CHECK ( nf90_def_var ( ncid, W_Y__SHORT_NAME, NF90_REAL, dimid3, varid_v ) )
         CALL CHECK ( nf90_def_var ( ncid, W_Z__SHORT_NAME, NF90_REAL, dimid3, varid_w ) )
         CALL CHECK ( nf90_def_var ( ncid, PRES_SHORT_NAME, NF90_REAL, dimid3, varid_p ) )
         CALL CHECK ( nf90_def_var ( ncid, DENS_SHORT_NAME, NF90_REAL, dimid3, varid_d ) )
         CALL CHECK ( nf90_def_var ( ncid, ALT__SHORT_NAME, NF90_REAL, dimid3, varid_h ) )
         CALL CHECK ( nf90_def_var ( ncid, TOPO_SHORT_NAME, NF90_REAL, dimid2, varid_o ) )
!
!.....2D
!
      ELSE IF ( nofdims .EQ. 2 ) THEN
!
!........Define the dimensions.
!
         CALL CHECK ( nf90_def_dim ( ncid, LVL__SHORT_NAME, nofz, z_dimid ) )
         CALL CHECK ( nf90_def_var ( ncid, LVL__SHORT_NAME, NF90_REAL, z_dimid, varid_z ) )
         IF ( dims( 1 ) ) THEN
            CALL CHECK ( nf90_def_dim ( ncid, LON__SHORT_NAME, nofx, x_dimid ) )
            CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, x_dimid, varid_x ) )
            CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, varid_y ) )
            dimid1 = x_dimid
            dimid2 = (/ x_dimid, z_dimid /)
         ELSE
            CALL CHECK ( nf90_def_dim ( ncid, LAT__SHORT_NAME, nofy, y_dimid ) )
            CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, y_dimid, varid_y ) )
            CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, varid_x ) )
            dimid1 = y_dimid
            dimid2 = (/ y_dimid, z_dimid /)
         END IF
!
!........Define variable
!
         CALL CHECK ( nf90_def_var ( ncid, TEMP_SHORT_NAME, NF90_REAL, dimid2, varid_t ) )
         CALL CHECK ( nf90_def_var ( ncid, W_X__SHORT_NAME, NF90_REAL, dimid2, varid_u ) )
         CALL CHECK ( nf90_def_var ( ncid, W_Y__SHORT_NAME, NF90_REAL, dimid2, varid_v ) )
         CALL CHECK ( nf90_def_var ( ncid, W_Z__SHORT_NAME, NF90_REAL, dimid2, varid_w ) )
         CALL CHECK ( nf90_def_var ( ncid, PRES_SHORT_NAME, NF90_REAL, dimid2, varid_p ) )
         CALL CHECK ( nf90_def_var ( ncid, DENS_SHORT_NAME, NF90_REAL, dimid2, varid_d ) )
         CALL CHECK ( nf90_def_var ( ncid, ALT__SHORT_NAME, NF90_REAL, dimid2, varid_h ) )
         CALL CHECK ( nf90_def_var ( ncid, TOPO_SHORT_NAME, NF90_REAL, dimid1, varid_o ) )
!
!.....1D
!
      ELSE IF ( nofdims .EQ. 1 ) THEN
!
!........Define the dimensions.
!
         CALL CHECK ( nf90_def_dim ( ncid, LVL__SHORT_NAME, nofz, z_dimid ) )
!
!........Define coordinate variables
!
         CALL CHECK ( nf90_def_var ( ncid, LVL__SHORT_NAME, NF90_REAL, z_dimid, varid_z ) )
         CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, varid_x ) )
         CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, varid_y ) )
!
!........Define variable
!
         CALL CHECK ( nf90_def_var ( ncid, TEMP_SHORT_NAME, NF90_REAL, z_dimid, varid_t ) )
         CALL CHECK ( nf90_def_var ( ncid, W_X__SHORT_NAME, NF90_REAL, z_dimid, varid_u ) )
         CALL CHECK ( nf90_def_var ( ncid, W_Y__SHORT_NAME, NF90_REAL, z_dimid, varid_v ) )
         CALL CHECK ( nf90_def_var ( ncid, W_Z__SHORT_NAME, NF90_REAL, z_dimid, varid_w ) )
         CALL CHECK ( nf90_def_var ( ncid, PRES_SHORT_NAME, NF90_REAL, z_dimid, varid_p ) )
         CALL CHECK ( nf90_def_var ( ncid, DENS_SHORT_NAME, NF90_REAL, z_dimid, varid_d ) )
         CALL CHECK ( nf90_def_var ( ncid, ALT__SHORT_NAME, NF90_REAL, z_dimid, varid_h ) )
         CALL CHECK ( nf90_def_var ( ncid, TOPO_SHORT_NAME, NF90_REAL, varid_o ) )
!
      END IF
!
!.....Assign units to variables
!
      CALL CHECK ( nf90_put_att(ncid, varid_x, UNITS, LON__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_y, UNITS, LAT__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_z, UNITS, LVL__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_t, UNITS, TEMP_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_u, UNITS, WIND_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_v, UNITS, WIND_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_w, UNITS, WIND_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_p, UNITS, PRES_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_d, UNITS, DENS_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_h, UNITS, ALT__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_o, UNITS, TOPO_UNITS ) )
!
!.....Assign long_names to variables
!
      CALL CHECK ( nf90_put_att(ncid, varid_x, LONG_NAME, LON__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_y, LONG_NAME, LAT__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_z, LONG_NAME, LVL__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_t, LONG_NAME, TEMP_LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_u, LONG_NAME, W_X__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_v, LONG_NAME, W_Y__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_w, LONG_NAME, W_Z__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_p, LONG_NAME, PRES_LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_d, LONG_NAME, DENS_LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_h, LONG_NAME, ALT__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_o, LONG_NAME, TOPO_LONG_NAME ) )
!
!.....end define mode
!
      CALL CHECK ( nf90_enddef ( ncid ) )
!
!.....Write Data
!
      CALL CHECK ( nf90_put_var ( ncid, varid_x, REAL ( xv ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_y, REAL ( yv ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_z, REAL ( zv ) ) )
!
      CALL CHECK ( nf90_put_var ( ncid, timeVarID, tstr ) )
      CALL CHECK ( nf90_put_var ( ncid, stepVarID, fcstep ) )
!
      CALL CHECK ( nf90_put_var ( ncid, varid_dataType, dataType ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_typeOfLevel, typeOfLevel ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_gridType, gridType ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_stream, stream ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_source, 'ECMWF' ) )
!
      IF ( nofdims .EQ. 3 ) THEN
         CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_w, REAL ( w ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_p, REAL ( p ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_d, REAL ( d ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_h, REAL ( h ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_o, REAL ( o ) ) )
      ELSE IF ( nofdims .EQ. 2 ) THEN
         IF ( dims( 1 ) ) THEN
            CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_w, REAL ( w(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_p, REAL ( p(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_d, REAL ( d(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_h, REAL ( h(:,1,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_o, REAL ( o(:,1) ) ) )
         ELSE
            CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_w, REAL ( w(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_p, REAL ( p(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_d, REAL ( d(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_h, REAL ( h(1,:,:) ) ) )
            CALL CHECK ( nf90_put_var ( ncid, varid_o, REAL ( o(1,:) ) ) )
         END IF
      ELSE
         CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_w, REAL ( w(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_p, REAL ( p(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_d, REAL ( d(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_h, REAL ( h(1,1,:) ) ) )
         CALL CHECK ( nf90_put_var ( ncid, varid_o, REAL ( o(1,1) ) ) )
      END IF
!
!.....Close netcdf file
!
      CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END SUBROUTINE NETCDF_WRITE_ATMO_RAW
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE WRITEGRID3_ATMO
!
!       Write 3d grid to netcdf grd file
!
!****************************
   SUBROUTINE WRITE_NETCDF_ATMO_GRD ( file1, nofx, nofy, nofz, t, u, v, w, p, d, o, xmin, dx, ymin, dy, zmin, dz  )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN )                       :: file1
      integer(kind=4), INTENT ( IN )                                 :: nofx, nofy, nofz
      DOUBLE PRECISION                    , INTENT ( IN )        :: xmin, ymin, zmin, dx, dy, dz
      DOUBLE PRECISION, DIMENSION ( nofx, nofy, nofz ), INTENT ( IN )  :: t, u, v, w, p, d
      DOUBLE PRECISION, DIMENSION ( nofx, nofy ), INTENT ( IN )  :: o
!
!.....Local variables
!
      integer(kind=4)                  :: ncid, x_dimid, y_dimid, z_dimid, i
      integer(kind=4)                  :: varid_x, varid_y, varid_z, varid_t, varid_u, varid_v, &
                                   & varid_w, varid_p, varid_d, varid_o
      integer(kind=4), DIMENSION ( 2 ) :: dimid2
      integer(kind=4), DIMENSION ( 3 ) :: dimid3
!
      DOUBLE PRECISION, DIMENSION ( nofx ) :: xpos
      DOUBLE PRECISION, DIMENSION ( nofy ) :: ypos
      DOUBLE PRECISION, DIMENSION ( nofz ) :: zpos
!
!  ---
!
!.....Make position vectors
!
      xpos = (/ ( xmin + i * dx, i = 0, nofx - 1 ) /)
      ypos = (/ ( ymin + i * dy, i = 0, nofy - 1 ) /)
      zpos = (/ ( zmin + i * dz, i = 0, nofz - 1 ) /)
!
!.....Create the netCDF file.
!
      CALL CHECK ( nf90_create ( file1, NF90_CLOBBER, ncid ) )
!
!.....Define the dimensions.
!
      CALL CHECK ( nf90_def_dim ( ncid, LON__SHORT_NAME, nofx, x_dimid ) )
      CALL CHECK ( nf90_def_dim ( ncid, LAT__SHORT_NAME, nofy, y_dimid ) )
      CALL CHECK ( nf90_def_dim ( ncid, ALT__SHORT_NAME, nofz, z_dimid ) )
!
!.....Define coordinate variables
!
      CALL CHECK ( nf90_def_var ( ncid, LON__SHORT_NAME, NF90_REAL, x_dimid, varid_x ) )
      CALL CHECK ( nf90_def_var ( ncid, LAT__SHORT_NAME, NF90_REAL, y_dimid, varid_y ) )
      CALL CHECK ( nf90_def_var ( ncid, ALT__SHORT_NAME, NF90_REAL, z_dimid, varid_z ) )
      dimid2 = (/ x_dimid, y_dimid /)
      dimid3 = (/ x_dimid, y_dimid, z_dimid /)
!
!.....Define variable
!
      CALL CHECK ( nf90_def_var ( ncid, TEMP_SHORT_NAME, NF90_REAL, dimid3, varid_t ) )
      CALL CHECK ( nf90_def_var ( ncid, W_X__SHORT_NAME, NF90_REAL, dimid3, varid_u ) )
      CALL CHECK ( nf90_def_var ( ncid, W_Y__SHORT_NAME, NF90_REAL, dimid3, varid_v ) )
      CALL CHECK ( nf90_def_var ( ncid, W_Z__SHORT_NAME, NF90_REAL, dimid3, varid_w ) )
      CALL CHECK ( nf90_def_var ( ncid, PRES_SHORT_NAME, NF90_REAL, dimid3, varid_p ) )
      CALL CHECK ( nf90_def_var ( ncid, DENS_SHORT_NAME, NF90_REAL, dimid3, varid_d ) )
      CALL CHECK ( nf90_def_var ( ncid, TOPO_SHORT_NAME, NF90_REAL, dimid2, varid_o ) )
!
!.....Assign units to variables
!
      CALL CHECK ( nf90_put_att(ncid, varid_x, UNITS, LON__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_y, UNITS, LAT__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_z, UNITS, ALT__UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_t, UNITS, TEMP_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_u, UNITS, WIND_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_v, UNITS, WIND_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_w, UNITS, WIND_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_p, UNITS, PRES_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_d, UNITS, DENS_UNITS ) )
      CALL CHECK ( nf90_put_att(ncid, varid_o, UNITS, TOPO_UNITS ) )
!
!.....Assign long_names to variables
!
      CALL CHECK ( nf90_put_att(ncid, varid_x, LONG_NAME, LON__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_y, LONG_NAME, LAT__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_z, LONG_NAME, ALT__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_t, LONG_NAME, TEMP_LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_u, LONG_NAME, W_X__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_v, LONG_NAME, W_Y__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_w, LONG_NAME, W_Z__LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_p, LONG_NAME, PRES_LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_d, LONG_NAME, DENS_LONG_NAME ) )
      CALL CHECK ( nf90_put_att(ncid, varid_o, LONG_NAME, TOPO_LONG_NAME ) )
!
!.....end define mode
!
      CALL CHECK ( nf90_enddef ( ncid ) )
!
!.....Write Data
!
      CALL CHECK ( nf90_put_var ( ncid, varid_x, REAL ( xpos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_y, REAL ( ypos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_z, REAL ( zpos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_w, REAL ( w ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_p, REAL ( p ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_d, REAL ( d ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_o, REAL ( o ) ) )
!
!.....Close netcdf file
!
      CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END SUBROUTINE WRITE_NETCDF_ATMO_GRD
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE WRITEGRID3_HYDRO
!
!       Write 3d grid to netcdf grd file
!
!****************************
   SUBROUTINE WRITE_NETCDF_HYDRO ( file1, nofx, nofy, nofz, t, u, v, w, p, d, s, xmin, dx, ymin, dy, zmin, dz  )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN )                       :: file1
      integer(kind=4), INTENT ( IN )                                 :: nofx, nofy, nofz
      DOUBLE PRECISION                    , INTENT ( IN )        :: xmin, ymin, zmin, dx, dy, dz
      DOUBLE PRECISION, DIMENSION ( nofx, nofy, nofz ), INTENT ( IN )  :: t, u, v, w, p, d
      DOUBLE PRECISION, DIMENSION ( nofx, nofy )      , INTENT ( IN )  :: s
!
!.....Local variables
!
      integer(kind=4)                  :: ncid, x_dimid, y_dimid, z_dimid, i
      integer(kind=4)                  :: varid_x, varid_y, varid_z, varid_t, varid_u, varid_v, &
                                   & varid_w, varid_p, varid_d, varid_s
      integer(kind=4), DIMENSION ( 2 ) :: dimid2
      integer(kind=4), DIMENSION ( 3 ) :: dimid3
!
      DOUBLE PRECISION, DIMENSION ( nofx ) :: xpos
      DOUBLE PRECISION, DIMENSION ( nofy ) :: ypos
      DOUBLE PRECISION, DIMENSION ( nofz ) :: zpos
!
!  ---
!
!.....Make position vectors
!
      xpos = (/ ( xmin + i * dx, i = 0, nofx - 1 ) /)
      ypos = (/ ( ymin + i * dy, i = 0, nofy - 1 ) /)
      zpos = (/ ( zmin + i * dz, i = 0, nofz - 1 ) /)
!
!.....Create the netCDF file.
!
      CALL CHECK ( nf90_create ( file1, NF90_CLOBBER, ncid ) )
!
!.....Define the dimensions.
!
      CALL CHECK ( nf90_def_dim ( ncid, "lon", nofx, x_dimid ) )
      CALL CHECK ( nf90_def_dim ( ncid, "lat", nofy, y_dimid ) )
      CALL CHECK ( nf90_def_dim ( ncid, "depth", nofz, z_dimid ) )
!
!.....Define coordinate variables
!
      CALL CHECK ( nf90_def_var ( ncid, "lon", NF90_REAL, x_dimid, varid_x ) )
      CALL CHECK ( nf90_def_var ( ncid, "lat", NF90_REAL, y_dimid, varid_y ) )
      CALL CHECK ( nf90_def_var ( ncid, "depth", NF90_REAL, z_dimid, varid_z ) )
      dimid2 = (/ x_dimid, y_dimid /)
      dimid3 = (/ x_dimid, y_dimid, z_dimid /)
!
!.....Define variable
!
      CALL CHECK ( nf90_def_var ( ncid, "water_temp", NF90_REAL, dimid3, varid_t ) )
      CALL CHECK ( nf90_def_var ( ncid, "water_u"   , NF90_REAL, dimid3, varid_u ) )
      CALL CHECK ( nf90_def_var ( ncid, "water_v"   , NF90_REAL, dimid3, varid_v ) )
      CALL CHECK ( nf90_def_var ( ncid, "water_w"   , NF90_REAL, dimid3, varid_w ) )
      CALL CHECK ( nf90_def_var ( ncid, "pressure"  , NF90_REAL, dimid3, varid_p ) )
      CALL CHECK ( nf90_def_var ( ncid, "density"   , NF90_REAL, dimid3, varid_d ) )
      CALL CHECK ( nf90_def_var ( ncid, "topography", NF90_REAL, dimid2, varid_s ) )
      CALL CHECK ( nf90_enddef ( ncid ) ) !End Definitions
!
!.....Write Data
!
      CALL CHECK ( nf90_put_var ( ncid, varid_x, REAL ( xpos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_y, REAL ( ypos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_z, REAL ( zpos ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_t, REAL ( t ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_u, REAL ( u ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_v, REAL ( v ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_w, REAL ( w ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_p, REAL ( p ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_d, REAL ( d ) ) )
      CALL CHECK ( nf90_put_var ( ncid, varid_s, REAL ( s ) ) )
!
!.....Close netcdf file
!
      CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END SUBROUTINE WRITE_NETCDF_HYDRO
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GRDDIMS 
!
!       Get dimensions of a netCDF 2D gridfile
!
!****************************
   SUBROUTINE GRDDIMS ( file1, nofx, nofy )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN )                       :: file1
      integer(kind=4), INTENT ( OUT )                                :: nofx, nofy
!
!.....Local variables
!
      integer(kind=4)                  :: ncid
      CHARACTER ( LEN = 50 )       :: xname, yname
!
!  ---
!
!.....Open netCDF file
!
      CALL CHECK ( nf90_open ( file1, nf90_nowrite, ncid ) )
!
!.....Inquire about the dimensions
!
      CALL CHECK ( nf90_inquire_dimension ( ncid, 1, xname, nofx ) )
      CALL CHECK ( nf90_inquire_dimension ( ncid, 2, yname, nofy ) )
!
!.....Close netCDF file
!
      CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END SUBROUTINE GRDDIMS 
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE READGRD 
!
!       Read 2d grid to netcdf grd file
!
!****************************
   SUBROUTINE READGRD ( file1, nofx, nofy, dat, xpos, ypos )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN )                           :: file1
      integer(kind=4), INTENT ( IN )                                     :: nofx, nofy
      DOUBLE PRECISION, DIMENSION ( nofx, nofy ), INTENT ( OUT )     :: dat
      DOUBLE PRECISION, DIMENSION ( nofx ), INTENT ( OUT ), OPTIONAL :: xpos
      DOUBLE PRECISION, DIMENSION ( nofy ), INTENT ( OUT ), OPTIONAL :: ypos
!
!.....Local variables
!
      integer(kind=4)                  :: ncid, xtype, ndims, varid
      integer(kind=4), DIMENSION ( 2 ) :: dimids
      CHARACTER ( LEN = 50 )       :: xname, yname, vname
!
!  ---
!
!.....Open netCDF file
!
      CALL CHECK ( nf90_open( file1, nf90_nowrite, ncid ) )
!
!.....Get the values of the coordinates and put them in xpos & ypos
!
      IF ( PRESENT ( xpos ) ) THEN
         CALL CHECK ( nf90_inquire_variable ( ncid, 1, xname, xtype, ndims, dimids ) )
         CALL CHECK ( nf90_inq_varid ( ncid, xname, varid ) )
         CALL CHECK ( nf90_get_var ( ncid, varid, xpos ) )
      END IF
!
      IF ( PRESENT ( ypos ) ) THEN
         CALL CHECK ( nf90_inquire_variable ( ncid, 2, yname, xtype, ndims, dimids ) )
         CALL CHECK ( nf90_inq_varid ( ncid, yname, varid ) )
         CALL CHECK ( nf90_get_var ( ncid, varid, ypos ) )
      END IF
!
!.....Get the values in z direction and put them in dat
!
      CALL CHECK ( nf90_inquire_variable ( ncid, 3, vname, xtype, ndims, dimids ) )
      CALL CHECK ( nf90_inq_varid ( ncid, vname, varid ) )
      CALL CHECK ( nf90_get_var ( ncid, varid, dat ) )
!
!.....Close netcdf file
!
      CALL CHECK ( nf90_close ( ncid ) )
!
      RETURN
!
   END SUBROUTINE READGRD
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE CHECK 
!
!       Check ierr of NETCDF routines and print error message and terminate program.
!
!****************************
   SUBROUTINE CHECK ( istatus )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      INTEGER, INTENT (IN) :: istatus
!
!  ---
!
      IF ( istatus .NE. nf90_noerr ) THEN
         PRINT *, TRIM (  ADJUSTL ( nf90_strerror ( istatus ) )  )
         STOP
      END IF
!
      RETURN
!
   END SUBROUTINE CHECK
!
! =========================================================================================
!
END MODULE TOOLS_NETCDF
