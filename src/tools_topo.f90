! =========================================================================================
!
! MODULE TOPO
!
! --> topo subroutines.
!
! =========================================================================================
!
Module tools_topo
!
   use tools_netcdf
   use tools_globe
!
   IMPLICIT NONE
!
!.....topography type
!
      TYPE type_topo
         LOGICAL         , DIMENSION ( 2 )                 :: circ
         DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE    :: h
         TYPE ( type_globetools_ranges )                   :: ranges
      END TYPE type_topo
!
!.....global variables
!
      TYPE( type_topo ), PUBLIC  :: topo
      LOGICAL, PUBLIC            :: ltopo

Contains

! =========================================================================================
!
!.......SUBROUTINE TOPOGET
!
!       get topography by use of 2D cubic Hermite Splines
!
!****************************
   SUBROUTINE TOPO_GET ( lon, lat, h, hdLon, hdLat )
!****************************
!
      USE CUBICWEIGHT
      use math, only: lonrangef
      use syms
      implicit none
!  
!.....dummy variables
!
      double precision, intent ( in  )           :: lon, lat
      double precision, intent ( out )           :: h
      double precision, intent ( out ), optional :: hdLon, hdLat
!
!.....local variables
!
      integer                                       :: ierr
      integer         , dimension ( 16, 2 )         :: ix
      double precision, dimension ( 16, 0:1, 0:1 )  :: cw
      double precision, dimension ( 0:1, 0:1 )      :: val
      double precision                              :: z
      double precision, dimension ( 2 )             :: scalar
!
! ----
!
!.....correct for horizontal stepsize in km!!
!
      z  = RE + h
      scalar( 1 ) = D2R * z * DCOS( lat * D2R )
      scalar( 2 ) = D2R * z
!
!.....Calculate cubic weights for 2-D
!
      CALL CALC_CW2 ( lon - topo%ranges%lonmin, lat - topo%ranges%latmin, cw, ix, &
         & topo%ranges%noflon, topo%ranges%dlon, topo%ranges%noflat, topo%ranges%dlat, 1, ierr, scalar  )
!
!.....Optain value and first order derivatives
!
      IF ( ierr .EQ. 0 ) THEN
!
         CALL EVAL_CW2 ( val, topo%h, cw, ix, 1, topo%circ )
!
         h = val( 0, 0 )
         IF (  PRESENT ( hdLat )  ) hdLat = val( 0, 1 )
         IF (  PRESENT ( hdLon )  ) hdLon = val( 1, 0 )
!
      ELSE
!
         h = 0.d0
         IF (  PRESENT ( hdLat )  ) hdLat = 0.d0
         IF (  PRESENT ( hdLon )  ) hdLon = 0.d0
!
      END IF
!
      RETURN
!
   END SUBROUTINE TOPO_GET 
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE TOPO_DEALLOCATE
!
!       Deallocate topo arrays
!
!****************************
   SUBROUTINE TOPO_DEALLOCATE ()
!****************************
!
      IMPLICIT NONE
!
      IF ( ALLOCATED ( topo%h ) ) DEALLOCATE ( topo%h   )
!
   END SUBROUTINE TOPO_DEALLOCATE
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE TOPO_SET_NEW
!
!       Initialize topo data.
!
!****************************
   SUBROUTINE TOPO_SET_NEW ( latMin, latMax, lonMin, lonMax, frq, l360 )
!****************************
!
      USE SYMS
      USE CUBICWEIGHT, ONLY : CW_SET_BOUND2
      IMPLICIT NONE
!
!.....dummy arguments
!
      DOUBLE PRECISION, INTENT ( IN ) :: latMin, latMax, lonMin, lonMax, frq
      LOGICAL         , INTENT ( IN ) :: l360
!
! ----
!
!.....copy parameters
!
      topo%ranges%latmin = latmin
      topo%ranges%lonmin = lonmin
      topo%ranges%latmax = latmax
      topo%ranges%lonmax = lonmax
!
      topo%ranges%dlat   = 2.d0 / frq / RE * R2D
      topo%ranges%dlon   = topo%ranges%dlat
!
      topo%circ = (/ l360, .FALSE. /)
!
!.....calculate nof
!
      topo%ranges%noflat    = INT (  ( latmax - latmin ) / topo%ranges%dlat + 1  )
      topo%ranges%noflon    = INT (  ( lonmax - lonmin ) / topo%ranges%dlon + 1  )
      topo%ranges%nofpoints = topo%ranges%noflat * topo%ranges%noflon
!
      topo%ranges%dlat = ( latmax - latmin ) / ( topo%ranges%noflat - 1 ) 
      topo%ranges%dlon = ( lonmax - lonmin ) / ( topo%ranges%noflon - 1 ) 
!
!.....Allocate
!
      ALLOCATE ( topo%h( 0:topo%ranges%nofLon+1, 0:topo%ranges%nofLat+1 ) )
!
!.....Get data from grd file
!
      CALL GET_TOPOGRAPHY ( topo%ranges, topo%h( 1:topo%ranges%nofLon, 1:topo%ranges%nofLat ) )
!
!.....Fix edges
!
      CALL CW_SET_BOUND2 ( topo%ranges%nofLon, topo%ranges%nofLat, topo%h, topo%circ )
!
!.....Remove ocean
!
      WHERE ( topo%h .LT. 0.d0 )
         topo%h = 0.d0
      END WHERE
!
      topo%h = topo%h / 1000.d0
!
      RETURN
!
   END SUBROUTINE TOPO_SET_NEW
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE TOPO_READ
!
!       Initialize topo data.
!
!****************************
   SUBROUTINE TOPO_READ ( file1, latMin, latMax, lonMin, lonMax, l360 )
!****************************
!
      USE SYMS
      USE CUBICWEIGHT, ONLY : CW_SET_BOUND2
      USE MATH, ONLY : LONRANGEF
      IMPLICIT NONE
!
!.....dummy arguments
!
      CHARACTER ( LEN = * ), INTENT ( IN )  :: file1
      DOUBLE PRECISION     , INTENT ( IN )  :: latMin, latMax, lonMin, lonMax
      LOGICAL              , INTENT ( IN )  :: l360
!
!.....local arguments
!
      DOUBLE PRECISION, DIMENSION ( : ), ALLOCATABLE :: lon, lat 
!
! ----
!
!.....get topo netcdf dimensions
!
!!print *, 'input: ', latmin, latmax, lonmin, lonmax
      CALL GRDDIMS ( file1, topo%ranges%noflon, topo%ranges%noflat )
!
!.....allocate arrays
!
      ALLOCATE ( lon( topo%ranges%noflon ), lat( topo%ranges%noflat ), &
         & topo%h( 0 : topo%ranges%nofLon + 1, 0 : topo%ranges%nofLat + 1 )  )
!
!.....Get data from grd file
!
      CALL READGRD ( file1, topo%ranges%nofLon, topo%ranges%nofLat, &
         & topo%h( 1:topo%ranges%nofLon, 1:topo%ranges%nofLat ), lon, lat )
!
!.....get parameters
!
!!print *, 'file: ', lat(1), lat(ubound(lat)), lon(1), lon(ubound(lon))
      topo%ranges%dlat   = lat( 2 ) - lat( 1 )
      topo%ranges%dlon   = lon( 2 ) - lon( 1 )
!
      topo%ranges%latmin = lat( 1 ) - topo%ranges%dlat
      topo%ranges%lonmin = lon( 1 ) - topo%ranges%dlon
      topo%ranges%latmax = lat( topo%ranges%noflat ) + topo%ranges%dlat
      topo%ranges%lonmax = lon( topo%ranges%noflon ) + topo%ranges%dlon
!
!!print *, 'topo: ', topo%ranges%latmin, topo%ranges%latmax, topo%ranges%lonmin, topo%ranges%lonmax
      IF (  topo%ranges%latmin .GT. latmin .OR. topo%ranges%lonmin .GT. lonmin &
         & .OR. topo%ranges%latmax .LT. latmax .OR. ( .NOT. l360 .AND. topo%ranges%lonmax .LT. lonmax )  ) &
         & STOP 'ERROR : ranges of topography does not match atmospheric ranges!'
!
      topo%circ = (/ l360, .FALSE. /)
!
!.....Fix edges
!
      CALL CW_SET_BOUND2 ( topo%ranges%nofLon, topo%ranges%nofLat, topo%h, topo%circ )
!
!.....deallocate
!
      DEALLOCATE ( lon, lat )
!
!.....Remove ocean
!
      WHERE ( topo%h .LT. 0.d0 )
         topo%h = 0.d0
      END WHERE
!
      topo%h = topo%h / 1000
!
      RETURN
!
   END SUBROUTINE TOPO_READ
!
! =========================================================================================
!
END MODULE TOOLS_TOPO
