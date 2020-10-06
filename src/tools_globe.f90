! =========================================================================================
!
! MODULE GLOBE 
!
! --> subroutines and functions for globe bathymetry and topography.
!
! =========================================================================================
!
MODULE TOOLS_GLOBE
!
   IMPLICIT NONE
!
!.....types
!
      TYPE type_globetools_ranges
         integer(kind=4)                        :: noflon, noflat, nofpoints
         DOUBLE PRECISION                 :: latMin, latMax, lonMin, lonMax, dLon, dLat
      END TYPE type_globetools_ranges
!
   CONTAINS
!
! =========================================================================================
!
!.......SUBROUTINE SET_GLOBERANGES 
!
!       Calculate the frequency based on the start frequency and the 2DFD index
!
!****************************
   SUBROUTINE SET_GLOBERANGES ( latMin, latMax, lonMin, lonMax, dLon, dLat, ranges )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      TYPE ( type_globetools_ranges ), INTENT ( INOUT ) :: ranges
      DOUBLE PRECISION, INTENT ( IN )                   :: latMin, latMax, lonMin, lonMax, dLon, dLat
!
!  ---
!
!.....copy
!
      ranges%latmin = latmin
      ranges%lonmin = lonmin
      ranges%latmax = latmax
      ranges%lonmax = lonmax
      ranges%dlat   = dlat
      ranges%dlon   = dlon
!
!.....calculate nof
!
      ranges%noflat    = INT (  ( latmax - latmin ) / dlat + 1  )
      ranges%noflon    = INT (  ( lonmax - lonmin ) / dlon + 1  )
      ranges%nofpoints = ranges%noflat * ranges%noflon
!
      RETURN
!
   END SUBROUTINE SET_GLOBERANGES
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GET_GLOBE 
!
!       Get globe data from gribfile
!
!****************************
   SUBROUTINE GET_GLOBE ( file1, ranges, dat, xpos, ypos )
!****************************
!
      USE TOOLS_NETCDF
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN )                                         :: file1
      TYPE ( type_globetools_ranges ), INTENT ( IN )                               :: ranges
      DOUBLE PRECISION, DIMENSION ( ranges%noflon, ranges%noflat ), INTENT ( OUT ) :: dat
      DOUBLE PRECISION, DIMENSION ( ranges%noflon ), INTENT ( OUT ), OPTIONAL      :: xpos
      DOUBLE PRECISION, DIMENSION ( ranges%noflat ), INTENT ( OUT ), OPTIONAL      :: ypos
!
!.....Local variables
!
      CHARACTER   :: str_noflon*9, str_noflat*9, str_lon1*10, str_lon2*10, str_lat1*10, str_lat2*10, &
                   & str_i*30, str_r*50
      LOGICAL     :: lexist
      integer(kind=4) ::  nofx, nofy
!
!  ---
!
!.....Write parameters to string
!
      WRITE ( str_noflon, "(I9)" ) ranges%noflon
      WRITE ( str_noflat, "(I9)" ) ranges%noflat
      WRITE ( str_lon1  , "(F10.5)" ) ranges%lonmin
      WRITE ( str_lon2  , "(F10.5)" ) ranges%lonmax
      WRITE ( str_lat1  , "(F10.5)" ) ranges%latmin
      WRITE ( str_lat2  , "(F10.5)" ) ranges%latmax
!
!.....Adjust and combine strings
!
      str_i = TRIM ( ADJUSTL ( str_noflon ) ) // '+/' // TRIM ( ADJUSTL ( str_noflat ) ) // '+'
      str_r = TRIM ( ADJUSTL ( str_lon1 ) ) // '/' // TRIM ( ADJUSTL ( str_lon2 ) ) // '/' &
         & // TRIM ( ADJUSTL ( str_lat1 ) ) // '/' // TRIM ( ADJUSTL ( str_lat2 ) )
!
!.....Resample global gebco bathymetry gridone.nc
!
      CALL SYSTEM ( 'grdsample ' // TRIM ( file1 ) // ' -G.globe.grd -I' &
         & // TRIM ( str_i ) // ' -R' // TRIM ( str_r ) // ' -fg' )
!
!.....Check if tmp.grd is created
!
      INQUIRE ( file = '.globe.grd', EXIST = lexist )
!
      IF ( lexist ) THEN
!
!........Get grid dimensions
!
         CALL GRDDIMS ( '.globe.grd', nofx, nofy )
!
!........Compare with ranges
!
         IF ( ranges%noflon .NE. nofx ) THEN
            PRINT "(A)", 'ERROR: grd bathymetry lon dimensions do not match!'
            STOP
         END IF
!
         IF ( ranges%noflat .NE. nofy ) THEN
            PRINT "(A)", 'ERROR: grd bathymetry lat dimensions do not match!'
            STOP
         END IF
!
!........Get grd values
!
         IF ( PRESENT ( xpos ) .AND. PRESENT ( ypos ) ) THEN
            CALL READGRD ( '.globe.grd', nofx, nofy, dat, xpos, ypos )
         ELSE IF ( PRESENT ( xpos ) ) THEN
            CALL READGRD ( '.globe.grd', nofx, nofy, dat, xpos )
         ELSE IF ( PRESENT ( ypos ) ) THEN
            CALL READGRD ( '.globe.grd', nofx, nofy, dat, ypos )
         ELSE
            CALL READGRD ( '.globe.grd', nofx, nofy, dat )
         END IF
!
!........Remove temporary grd
!
!!         CALL SYSTEM ( 'rm -rf .globe.grd' )
!
      ELSE
!
         WRITE ( *, "(2A)" )'ERROR: .globe.grd could not be created from "' // TRIM ( file1 ) // '".', &
         & 'Please check the bash variable, file location or grdsample (GMT).'
         STOP
!
      END IF
!
      RETURN
!
   END SUBROUTINE GET_GLOBE
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GET_BATHYMETRY  
!
!       Get bathymetry from $TOOLS_GLOBE_BATH
!
!****************************
   SUBROUTINE GET_BATHYMETRY ( ranges, dat, xpos, ypos )
!****************************
!
      USE TOOLS_NETCDF
      IMPLICIT NONE
!
!.....Dummy variables
!
      TYPE ( type_globetools_ranges ), INTENT ( IN )                               :: ranges
      DOUBLE PRECISION, DIMENSION ( ranges%noflon, ranges%noflat ), INTENT ( OUT ) :: dat
      DOUBLE PRECISION, DIMENSION ( ranges%noflon ), INTENT ( OUT ), OPTIONAL      :: xpos
      DOUBLE PRECISION, DIMENSION ( ranges%noflat ), INTENT ( OUT ), OPTIONAL      :: ypos
!
!  ---
!
      IF ( PRESENT ( xpos ) .AND. PRESENT ( ypos ) ) THEN
         CALL GET_GLOBE ( '$TOOLS_GLOBE_BATH', ranges, dat, xpos, ypos )
      ELSE IF ( PRESENT ( xpos ) .AND. .NOT. PRESENT ( ypos ) ) THEN
         CALL GET_GLOBE ( '$TOOLS_GLOBE_BATH', ranges, dat, xpos )
      ELSE IF ( .NOT. PRESENT ( xpos ) .AND. PRESENT ( ypos ) ) THEN
         CALL GET_GLOBE ( '$TOOLS_GLOBE_BATH', ranges, dat, ypos )
      ELSE
         CALL GET_GLOBE ( '$TOOLS_GLOBE_BATH', ranges, dat )
      END IF
!
      RETURN
!
   END SUBROUTINE GET_BATHYMETRY
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GET_TOPOGRAPHY  
!
!       Get topography from $TOOLS_GLOBE_TOPO
!
!****************************
   SUBROUTINE GET_TOPOGRAPHY ( ranges, dat, xpos, ypos )
!****************************
!
      USE TOOLS_NETCDF
      IMPLICIT NONE
!
!.....Dummy variables
!
      TYPE ( type_globetools_ranges ), INTENT ( IN )                               :: ranges
      DOUBLE PRECISION, DIMENSION ( ranges%noflon, ranges%noflat ), INTENT ( OUT ) :: dat
      DOUBLE PRECISION, DIMENSION ( ranges%noflon ), INTENT ( OUT ), OPTIONAL      :: xpos
      DOUBLE PRECISION, DIMENSION ( ranges%noflat ), INTENT ( OUT ), OPTIONAL      :: ypos
!
!  ---
!
      IF ( PRESENT ( xpos ) .AND. PRESENT ( ypos ) ) THEN
         CALL GET_GLOBE ( '$TOOLS_GLOBE_TOPO', ranges, dat, xpos, ypos )
      ELSE IF ( PRESENT ( xpos ) .AND. .NOT. PRESENT ( ypos ) ) THEN
         CALL GET_GLOBE ( '$TOOLS_GLOBE_TOPO', ranges, dat, xpos )
      ELSE IF ( .NOT. PRESENT ( xpos ) .AND. PRESENT ( ypos ) ) THEN
         CALL GET_GLOBE ( '$TOOLS_GLOBE_TOPO', ranges, dat, ypos )
      ELSE
         CALL GET_GLOBE ( '$TOOLS_GLOBE_TOPO', ranges, dat )
      END IF
!
      RETURN
!
   END SUBROUTINE GET_TOPOGRAPHY
!
! =========================================================================================
!
END MODULE TOOLS_GLOBE
