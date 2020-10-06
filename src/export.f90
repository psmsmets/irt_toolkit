! =========================================================================================
!
! MODULE EXPORT
!
! --> subroutines to export cross section through atmospheric profile.
!
! =========================================================================================
!
MODULE EXPORT
!
   USE SYMS
   USE TOOLS_ATMO
!
CONTAINS
! =========================================================================================
!
!.......SUBROUTINE EXPORT_VERTICAL_CROSSECTION 
!
!       export vertical cross-section of the atmospheric atmo% for a specific bearing
!
!****************************
  SUBROUTINE EXPORT_VERTICAL_CROSSECTION ( prefix, lat, lon, bearing, dx, dh, dmax )
!****************************
!  
      USE TOOLS_TOPO
      USE MATH
      IMPLICIT NONE
!          
!.....dummy variables
!
      DOUBLE PRECISION     , INTENT ( IN  )  :: lon, lat, bearing, dx, dh
      CHARACTER ( LEN = * ), INTENT ( IN  )  :: prefix
      DOUBLE PRECISION     , INTENT ( OUT )  :: dmax
!
!.....local variables
!
      integer(kind=4)                    :: i, nofh, hskip, ierr
      DOUBLE PRECISION               :: lon1, lat1, bearing1, h, x, r, dr, Ceff, htopo, w_at, w_ct
      DOUBLE PRECISION, PARAMETER    :: ex = -999.9d0
      DOUBLE PRECISION, DIMENSION ( 0:2, 0:2, 0:2 ) :: C, Wlon, Wlat, Wv
!
! ----
!          
!.....initialize
!
      lon1        = lon
      lat1        = lat
      bearing1    = bearing
      hskip       = -1
      x           = 0.d0
      r           = 0.d0
      dr          = dx / RE * R2D
      dmax        = floor( sc%d / dx ) * dx
!
!.....set atmo options
!
      nofh  = INT (  ( sc%hmax - sc%hmin ) / dh, KIND = 4  )
!
!.....open output file
!
      OPEN ( UNIT = 13, STATUS = 'replace', FILE = prefix // '_ceff.xyz', FORM = 'formatted' )
!
!.....loop from source to end of profile
!
      DO WHILE ( lonrangef(lon1+atmo%lonshift).LT.sc%lonrange .AND. lonrangef(lon1+atmo%lonshift).GT.0.d0 &
         & .AND. lat1 .LT. sc%latmax .AND. lat1 .GT. sc%latmin .AND. x .LE. dmax )
!
!........WITH TOPOGRAPHY
!              
         IF ( ltopo ) THEN
!
!...........get topography for this location
!
            CALL TOPO_GET ( lon1, lat1, htopo )
!
!...........values below surface of the earth
!
            hskip = FLOOR (  ( htopo - sc%hmin ) / dh  )
!
!...........skip points below the surface, WRITE 1001 as value
!
            DO i = 0, hskip
               h = sc%hmin + DBLE ( i ) * dh
               WRITE (13,"(F9.2,1x,F9.5,1x,F10.5)") x, h, ex
            END DO
!
         END IF
!
!........loop for altitude
!
         DO i = hskip + 1, nofh
!
            h = sc%hmin + DBLE ( i ) * dh
!
            CALL ATMO_GET ( lon1 + atmo%lonshift, lat1, h, .FALSE., ierr, C, Wlon, Wlat, Wv )
!
            IF ( ierr .EQ. 0 ) THEN
!
!..............Wind components
!
               w_at = DSIN ( bearing1 * D2R ) * Wlon( 0, 0, 0 ) &
                  & + DCOS ( bearing1 * D2R ) * ( - Wlat( 0, 0, 0 ) )
               w_ct = DCOS ( bearing1 * D2R ) * Wlon( 0, 0, 0 ) &
                  & + DSIN ( bearing1 * D2R ) * ( - Wlat( 0, 0, 0 ) )
!
!..............Calculate effective speed of sound
!              = sum of speed of sound and horizontal wind
!
               Ceff = C( 0, 0, 0 ) + w_at
!
!..............Write distance, altitude and ceff to file
!
               WRITE ( 13, "(F9.2,1x,F9.5,1x,F10.5)" ) x, h, Ceff 
!
            END IF
!
         END DO
!
!........next point
!
         x = x + dx
         r = r + dr
         CALL RADIAL_DISTANCE ( lat, lon, r, bearing, lat1, lon1, bearing1 )
!
!........correct longitude range
!
         CALL LONRANGE ( lon1 )
!
!.....end of loop along profile
!
      END DO
!
!.....close output file
!
      CLOSE ( 13 )
!
!.....maximum distance (required for plot ranges)
!
      dmax = x - dx
!
      RETURN
!
   END SUBROUTINE EXPORT_VERTICAL_CROSSECTION 
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE exportHorizontalProfile
!
!       export horizontal cross-section of the atmospheric profile for a specific altitude
!
!****************************
   SUBROUTINE exportHorizontalProfile ( altitude, res, out_cross_hor )
!****************************
!
      USE MATH, ONLY: LONRANGEF
      IMPLICIT NONE
!          
! ....dummy variables
!
      DOUBLE PRECISION      , INTENT ( IN )  :: altitude,res
      CHARACTER ( LEN = 32 ), INTENT ( IN )  :: out_cross_hor
!
!.....local variables
!
      integer(kind=4)                                   :: i, j, nofLon_cross, nofLat_cross, ierr
      DOUBLE PRECISION                              :: lon, lat, lon0, lat0
      DOUBLE PRECISION, DIMENSION ( 0:2, 0:2, 0:2 ) :: C, Wlon, Wlat, Wv
!
! ----
!
!.....profile options
!
      nofLon_cross = INT (  ( atmo%ranges%dlon * ( atmo%ranges%noflon - 1 ) ) / res + 1, KIND = 4  )
      nofLat_cross = INT (  ( atmo%ranges%dlat * ( atmo%ranges%noflat - 1 ) ) / res + 1, KIND = 4  )
!
!.....initialize
!
      lon0 = atmo%ranges%lonmin
      lat0 = atmo%ranges%latmin
!
!.....output file
!
      OPEN ( UNIT = 14, STATUS = 'replace', FILE = out_cross_hor, FORM = 'formatted' )
!
!.....Write header with altitude, and number of points in lon and lat direction
!
      WRITE ( 14, "(F5.2,2I8)" ) altitude, nofLat_cross, nofLon_cross
!
!.....Loop to export horizontal plane
!
      DO i = 0, nofLon_cross - 1
!
         lon = LONRANGEF( lon0 + res * i )
!
         DO j = 0, nofLat_cross - 1
!
            lat = lat0 + res * j
            CALL ATMO_GET (  lon + atmo%lonshift, lat, altitude, .FALSE., ierr, C, Wlon, Wlat, Wv  )
            IF ( ierr .EQ. 0 ) WRITE ( 14, "(2F8.2,4F9.5)" ) lon, lat, C( 0, 0, 0 ), &
               & Wlon( 0, 0, 0 ), Wlat( 0, 0, 0 ), Wv( 0, 0, 0 )
!
         END DO
!
      END DO
!
      CLOSE( 14 )
!
      RETURN
!
   END SUBROUTINE exportHorizontalProfile
!
! =========================================================================================



! =========================================================================================
!
!.......SUBROUTINE exportOffsetTime
!
!       export offset and time with azimuth plane
!
!****************************
   SUBROUTINE EXPORT_OFFSET_TIME ( prefix, rayMin, rayMax, azi, lat0, lon0, tmax, omin, omax )
!****************************
!
      USE IRTPL_READ
      USE MATH, ONLY: COURSE
      IMPLICIT NONE
!          
! ....dummy variables
!
      DOUBLE PRECISION     , INTENT ( IN  )  :: azi, lat0, lon0
      CHARACTER ( LEN = * ), INTENT ( IN  )  :: prefix
      integer(kind=4)          , INTENT ( IN  )  :: rayMin, rayMax
      DOUBLE PRECISION     , INTENT ( OUT )  :: tmax, omin, omax
!
!.....local variables
!
      DOUBLE PRECISION                :: d, a, offset, dangle
      integer(kind=4)                     :: i, nofrec 
!
! ----
!
!.....initialize
!
      nofrec = SIZE ( store_refl, 1 )
      tmax   = 0
      omax   = 0
      omin   = 99999.d0
!
!.....open output file
!
      OPEN ( UNIT = 14, STATUS = 'replace', FILE = prefix // '_offset.xy', FORM = 'formatted' )
      OPEN ( UNIT = 15, STATUS = 'replace', FILE = prefix // '_time.xy'  , FORM = 'formatted' )
!
!.....reflections?
!
      IF (  ALLOCATED ( store_refl )  ) THEN
!
!........Loop to export reflections offset and time
!
         DO i = 1, nofrec
!
            IF (  store_refl( i, 1 ) .GE. rayMin .AND. store_refl( i, 1 ) .LE. rayMax  ) THEN
!
               CALL COURSE ( lat0, lon0, store_refl( i, 4 ), store_refl( i, 3 ), d, a )
!
               dangle = azi - a
               IF (  DABS ( dangle ) .GT. 180.d0  ) dangle = SIGN ( 1.d0, dangle ) * ( dangle - 360.d0 )
               offset = dangle * D2R * d
!
               WRITE( 14, "(F9.2,1x,F9.2)" ) d, offset
               WRITE( 15, "(F9.2,1x,F9.2)" ) d, store_refl( i, 2 )
!
               IF ( offset .GT. omax ) omax = offset
               IF ( offset .LT. omin ) omin = offset
               IF ( store_refl(i,2) .GT. tmax ) tmax = store_refl(i,2)
!
            END IF
!
         END DO
!
      ELSE
!
         omin =  -10.d0
         omax =   10.d0
         tmax = 1000.d0
!
      END IF
!
      CLOSE ( 14 )
      CLOSE ( 15 )
!
      RETURN
!
   END SUBROUTINE EXPORT_OFFSET_TIME
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE EXPORT_PROFILE
!
!       export profile for a specified coordinate
!
!****************************
  SUBROUTINE EXPORT_PROFILE ( prefix, lat, lon, bearing, dh )
!****************************
!  
      USE TOOLS_TOPO
      USE MATH, ONLY: LONRANGEF
      IMPLICIT NONE
!          
!.....dummy variables
!
      DOUBLE PRECISION     , INTENT ( IN )   :: lon, lat, bearing, dh
      CHARACTER ( LEN = * ), INTENT ( IN )   :: prefix
!
!.....local variables
!
      integer(kind=4)                                   :: i, nofh, hskip, ierr
      DOUBLE PRECISION                              :: h1, htopo, Ceff, w_at, w_ct
      DOUBLE PRECISION, DIMENSION ( 0:2, 0:2, 0:2 ) :: C, Wlon, Wlat, Wv
!
! ----
!
!.....initialize 
!
      nofh  = INT (  ( sc%hmax - sc%hmin ) / dh, KIND = 4  )
      hskip = 0
!
!.....open output file
!
      OPEN ( UNIT = 16, STATUS = 'replace', FILE = prefix // '_c.xy', FORM = 'formatted' )
      OPEN ( UNIT = 17, STATUS = 'replace', FILE = prefix // '_ceff.xy', FORM = 'formatted' )
!
!.....Topography?
!              
      IF ( ltopo ) THEN
!
!........get topography for this location
!
         CALL TOPO_GET ( lon, lat, htopo )
!
!........values below surface of the earth
!
         hskip = FLOOR (  ( htopo - sc%hmin ) / dh  ) + 1
!
      END IF
!
!.....Export vertical profile
!
      DO i = hskip, nofh
!
         h1 = sc%hmin + DBLE ( i ) * dh
!
         CALL ATMO_GET (  lon + atmo%lonshift, lat, h1, .FALSE., ierr, C, Wlon, Wlat, Wv  )
!
         IF ( ierr .EQ. 0 ) THEN
!
!...........Wind components
!
            w_at = DSIN ( bearing * D2R ) * Wlon( 0, 0, 0 ) &
               & + DCOS ( bearing * D2R ) * ( - Wlat( 0, 0, 0 ) )
            w_ct = DCOS ( bearing * D2R ) * Wlon( 0, 0, 0 ) &
               & + DSIN ( bearing * D2R ) * ( - Wlat( 0, 0, 0 ) )
!
!...........Calculate effective speed of sound
!             = sum of speed of sound and horizontal wind
!
            Ceff = C( 0, 0, 0 ) + w_at
!
            WRITE ( 16, "(F10.5,1x,F9.5)" ) C( 0, 0, 0 ), h1
            WRITE ( 17, "(F10.5,1x,F9.5)" ) Ceff, h1
!
         END IF
!
      END DO
!
      CLOSE ( 16 )
      CLOSE ( 17 )
!
      RETURN
!
   END SUBROUTINE EXPORT_PROFILE
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine EXPORT_MARKERS 
!
!****************************
   SUBROUTINE EXPORT_MARKERS ( prefix, marker, marker_sym )
!****************************
!
      IMPLICIT NONE
!          
!.....dummy variables
!
      CHARACTER ( LEN = * )                     , INTENT ( IN  )  :: prefix
      CHARACTER ( LEN = * ), DIMENSION ( : )    , INTENT ( IN  )  :: marker_sym
      DOUBLE PRECISION     ,  DIMENSION ( :, : ), INTENT ( IN  )  :: marker
!          
!.....local variables
!
      integer(kind=4) :: marker_cnt, i, munit = 21
!
! ----
!
!.....get markers
!
      marker_cnt = SIZE ( marker, 1 )
!
      OPEN ( UNIT = munit, STATUS = 'replace', FILE = prefix // '_markers.xy', FORM = 'formatted' )
!
      DO i = 1, marker_cnt
!
         WRITE ( munit, "(F10.5,1x,F9.5,1x,A)" ) marker( i, 2 ), marker( i, 1 ), TRIM (  marker_sym( i )  )
!
      END DO
!
      CLOSE ( munit )
!
      RETURN
!
   END SUBROUTINE EXPORT_MARKERS
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine EXPORT_REFL 
!
!****************************
   SUBROUTINE EXPORT_REFL ( prefix )
!****************************
!
      USE IRTPL_READ
      USE SORT, ONLY : cocktailSortIx_DBLE
!
      IMPLICIT NONE
!          
!.....dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN  ) :: prefix
!          
!.....local variables
!
      INTEGER :: i, nofrefl, munit1 = 31, munit2 = 32, munit3 = 33
      INTEGER, DIMENSION ( SIZE ( store_refl, 1 ) )  :: ix1, ix2, ix3 
!
! ----
!
!.....open output files
!
      OPEN ( UNIT = munit1, STATUS = 'replace', FILE = prefix // '_refl_h.xy', FORM = 'formatted' )
      OPEN ( UNIT = munit2, STATUS = 'replace', FILE = prefix // '_refl_t.xy', FORM = 'formatted' )
      OPEN ( UNIT = munit3, STATUS = 'replace', FILE = prefix // '_refl_tloss.xy', FORM = 'formatted' )
!
!.....reflections?
!
      IF (  ALLOCATED ( store_refl )  ) THEN
!
!........initialize
!
         nofrefl = SIZE ( store_refl, 1 )
!
!........Sort output
!
         CALL cocktailSortIx_DBLE (  store_refl( :, 7 ), 0, ix1 )
         CALL cocktailSortIx_DBLE (  store_refl( :, 2 ), 0, ix2 )
         CALL cocktailSortIx_DBLE ( -store_refl( :, 6 ), 0, ix3 )
!
!........write to file
!
         DO i = 1, nofrefl
!
            WRITE ( munit1, "(F10.5,1x,F9.5,1x,F9.5)" ) store_refl( ix1( i ), 3 : 4 ), store_refl( ix1( i ), 7 )
            WRITE ( munit2, "(F10.5,1x,F9.5,1x,F9.2)" ) store_refl( ix2( i ), 3 : 4 ), store_refl( ix2( i ), 2 )
            WRITE ( munit3, "(F10.5,1x,F9.5,1x,F12.5)") store_refl( ix3( i ), 3 : 4 ), store_refl( ix3( i ), 6 )
!
         END DO
!
      END IF
!
      CLOSE ( munit1 )
      CLOSE ( munit2 )
      CLOSE ( munit3 )
!
      RETURN
!
   END SUBROUTINE EXPORT_REFL
!
! =========================================================================================
!
END MODULE EXPORT
