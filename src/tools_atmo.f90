! =========================================================================================
!
! MODULE ATMO
!
!        Get ECMWF atmospheric profiles from grib files, regularize profiles, 3D grid
!        cubic interpolation routine
!
! =========================================================================================
!
MODULE TOOLS_ATMO
!
   use tools_grib
   implicit none
!
!.....custom data types
!
      TYPE type_atmo_raw
         DOUBLE PRECISION, DIMENSION ( :, :, : ), ALLOCATABLE  :: C, W1, W2, W3, P, Rh, H
         DOUBLE PRECISION, DIMENSION ( :, : )   , ALLOCATABLE  :: topo
         DOUBLE PRECISION, DIMENSION ( : )      , ALLOCATABLE  :: lon, lat, lvl
      END TYPE type_atmo_raw
!
      TYPE type_atmo_grd
         DOUBLE PRECISION, DIMENSION ( :, :, : ), ALLOCATABLE  :: C, W1, W2, W3, P, Rh
         DOUBLE PRECISION, DIMENSION ( :, : )   , ALLOCATABLE  :: H
         DOUBLE PRECISION, DIMENSION ( : )      , ALLOCATABLE  :: lon, lat, alt
      END TYPE type_atmo_grd
!
      TYPE type_atmo
         LOGICAL                            :: l360, lreverse
         integer(kind=4)                    :: nofh, noflon, noflat, step, ens
         DOUBLE PRECISION                   :: lonshift, hmin, hmax, dh, lonmin, lonmax, dlon, lonrange, &
                                             & latmin, latmax, dlat, latrange
         TYPE ( type_gribtools_ranges )                       :: ranges
         TYPE ( type_atmo_grd )                               :: grd
         TYPE ( type_atmo_raw )                               :: raw
      END TYPE type_atmo
!
!.....stop condition type
!
      TYPE type_sc
         CHARACTER                                            :: stype * 2
         DOUBLE PRECISION                                     :: latmin, latmax, lonmin, lonmax, lonrange, &
                                                               & latrange, d, hmin, hmax 
      END TYPE type_sc
!
!.....global variables
!
      TYPE( type_atmo ), PUBLIC, target  :: atmo
      TYPE( type_sc )  , PUBLIC, target  :: sc
!
!.....molecular fraction contants
!
      double precision, dimension(6), private, parameter  :: O_a1   = &
         (/ -1.1195d1, 1.5408d-1, -1.4348d-3, 1.0166d-5, 0.d0, 0.d0 /)
      double precision, dimension(6), private, parameter  :: O_a2   = &
         (/ -3.2456d0, 4.6642d-2, -2.6894d-2, 5.264d-7, 0.d0, 0.d0 /)
      double precision, dimension(6), private, parameter  :: O2_a1  = &
         (/ -0.67887d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
      double precision, dimension(6), private, parameter  :: O2_a2  = &
         (/ 49.296d0, -1.5524d0, 1.8714d-2, -1.1069d-4, 3.1990d-7, -3.6211d-10 /)
      double precision, dimension(6), private, parameter  :: O3_a1  = &
         (/ -1.9027d1, 1.3093d0, -4.6496d-2, 7.8543d-4, -6.5169d-6, 2.1343d-8 /)
      double precision, dimension(6), private, parameter  :: O3_a2  = &
         (/ -4.234d0, -3.0975d-2, 0.d0, 0.d0, 0.d0, 0.d0 /)
      double precision, dimension(6), private, parameter  :: N_a1   = &
         (/ -5.3746d1, 1.5439d0, -1.8824d-2, 1.1587d-4, -3.5399d-7, 4.2609d-10 /)
      double precision, dimension(6), private, parameter  :: N_a2   = &
         (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
      double precision, dimension(6), private, parameter  :: N2_a1  = &
         (/ -0.10744d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
      double precision, dimension(6), private, parameter  :: N2_a2  = &
         (/ 1.3972d-1, -5.6269d-3, 3.9407d-5, -1.0737d-7, 0.d0, 0.d0 /)
      double precision, dimension(6), private, parameter  :: CO2_a1 = &
         (/ -3.3979d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
      double precision, dimension(6), private, parameter  :: CO2_a2 = &
         (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
!
Contains

! =========================================================================================
!
!.......SUBROUTINE ATMO_GET 
!
!       get atmospheric specifications of wind and temperature and certain position,
!       provided in longitude, latitude, and altitude by use of atmo splines.
!
!       Note that longitude should be provided as dlon, the angular distance with
!       respect to lon_min ( dlon = lon - lon_min ).
!
!****************************
   SUBROUTINE ATMO_GET ( dlon, lat, h, unit_km, ierr, C, W1, W2, W3, P, Rh )
!****************************
!
      USE CUBICWEIGHT
      USE MATH, ONLY: LONRANGEF
      USE SYMS
      IMPLICIT NONE
!  
!.....dummy variables
!
      LOGICAL                                      , INTENT ( IN  )           :: unit_km
      DOUBLE PRECISION                             , INTENT ( IN  )           :: dlon, lat, h
      integer(kind=4)                              , INTENT ( OUT )           :: ierr
      DOUBLE PRECISION, DIMENSION ( 0:2, 0:2, 0:2 ), INTENT ( OUT )           :: C, W1, W2, W3
      DOUBLE PRECISION, DIMENSION ( 0:2, 0:2, 0:2 ), INTENT ( OUT ), OPTIONAL :: P, Rh
!
!.....local variables
!
      INTEGER         , DIMENSION ( 64, 3 )              :: ix
      DOUBLE PRECISION, DIMENSION ( 3 )                  :: scalar
      DOUBLE PRECISION, DIMENSION ( 64, 0:2, 0:2, 0:2 )  :: cw
      DOUBLE PRECISION, DIMENSION ( 0:2, 0:2, 0:2 )      :: P_, Rh_
!
! --------------------
!
!.....correct for <longitude,latitude,h> to <phi,theta,r>
!
      if ( atmo%lreverse ) then
         scalar = (/ -D2R*cos(lat*D2R), D2R, 1.d0 /)
      else
         scalar = (/ D2R*cos(lat*D2R), -D2R, 1.d0 /)
      end if
!
!.....Calculate cubic weights for 3-D
!
      CALL CALC_CW3 (   LONRANGEF ( dlon ), lat - atmo%latmin, h - atmo%hmin, cw, ix, &
         & atmo%noflon, atmo%dlon, atmo%noflat, atmo%dlat, atmo%nofh, atmo%dh, 2, ierr, scalar  )
!
      IF ( ierr .EQ. 0 ) THEN
!
!........Optain value and first order derivatives
!
         CALL ATMO_EVAL_CW3 ( C, W1, W2, W3, P_, Rh_, cw, ix )
!
         IF ( PRESENT ( P  ) ) P  = P_
         IF ( PRESENT ( Rh ) ) Rh = Rh_
!
!........m/s to km/s
!
         IF ( unit_km ) THEN
            C  = C  / 1000.d0
            W1 = W1 / 1000.d0
            W2 = W2 / 1000.d0
            W3 = W3 / 1000.d0
         END IF
!
      END IF
!
      RETURN
!
   END SUBROUTINE ATMO_GET
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CALC_CW3
!
!       Get cubic weight up and data indices for 3-D data.
!
!****************************
   SUBROUTINE ATMO_EVAL_CW3 ( C, W1, W2, W3, P, Rh, cw, ix )
!****************************
!
      USE CUBICWEIGHT
      IMPLICIT NONE
!
!.....Dummy variables
!
      double precision, dimension( 0:2, 0:2, 0:2 )    , intent( OUT )  :: C, W1, W2, W3, P, Rh
      double precision, dimension( 64, 0:2, 0:2, 0:2 ), intent( IN  ), target  :: cw
      integer         , dimension( 64, 3 )            , intent( IN  ), target  :: ix
!
!.....Local variables
!
      integer           :: i
      integer, pointer :: ix1, ix2, ix3
      double precision, dimension(:,:,:), pointer :: cwi
!
!  ---
!    
      C  = 0.d0
      W1 = 0.d0
      W2 = 0.d0
      W3 = 0.d0
      P  = 0.d0
      Rh = 0.d0
!
      DO i = 1, 64
         ix1=>ix(i,1)
         ix2=>ix(i,2)
         ix3=>ix(i,3)
         cwi=>cw( i, 0:2, 0:2, 0:2 )
         C  = C  + atmo%grd%C ( ix1, ix2, ix3 ) * cwi
         W1 = W1 + atmo%grd%W1( ix1, ix2, ix3 ) * cwi
         W2 = W2 + atmo%grd%W2( ix1, ix2, ix3 ) * cwi
         W3 = W3 + atmo%grd%W3( ix1, ix2, ix3 ) * cwi
         P  = P  + atmo%grd%P ( ix1, ix2, ix3 ) * cwi
         Rh = Rh + atmo%grd%Rh( ix1, ix2, ix3 ) * cwi
      END DO
!
      RETURN
!
   END SUBROUTINE ATMO_EVAL_CW3
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE ATMO_DEALLOCATE
!
!       deallocate atmospheric arrays.
!
!****************************
   SUBROUTINE ATMO_DEALLOCATE ( what )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = 3 ), INTENT ( IN ) :: what
!
! ---
!
      IF ( what .EQ. 'raw' .OR. what .EQ. 'all' ) THEN
         IF (  ALLOCATED ( atmo%raw%H    )  ) DEALLOCATE ( atmo%raw%H    )
         IF (  ALLOCATED ( atmo%raw%C    )  ) DEALLOCATE ( atmo%raw%C    )
         IF (  ALLOCATED ( atmo%raw%W1   )  ) DEALLOCATE ( atmo%raw%W1   )
         IF (  ALLOCATED ( atmo%raw%W2   )  ) DEALLOCATE ( atmo%raw%W2   )
         IF (  ALLOCATED ( atmo%raw%W3   )  ) DEALLOCATE ( atmo%raw%W3   )
         IF (  ALLOCATED ( atmo%raw%P    )  ) DEALLOCATE ( atmo%raw%P    )
         IF (  ALLOCATED ( atmo%raw%Rh   )  ) DEALLOCATE ( atmo%raw%Rh   )
         IF (  ALLOCATED ( atmo%raw%H    )  ) DEALLOCATE ( atmo%raw%H    )
         IF (  ALLOCATED ( atmo%raw%TOPO )  ) DEALLOCATE ( atmo%raw%TOPO )
         IF (  ALLOCATED ( atmo%raw%lon  )  ) DEALLOCATE ( atmo%raw%lon  )
         IF (  ALLOCATED ( atmo%raw%lat  )  ) DEALLOCATE ( atmo%raw%lat  )
         IF (  ALLOCATED ( atmo%raw%lvl  )  ) DEALLOCATE ( atmo%raw%lvl  )
      END IF
!
      IF ( what .EQ. 'grd' .OR. what .EQ. 'all' ) THEN
         IF (  ALLOCATED ( atmo%grd%C   )  ) DEALLOCATE ( atmo%grd%C   )
         IF (  ALLOCATED ( atmo%grd%W1  )  ) DEALLOCATE ( atmo%grd%W1  )
         IF (  ALLOCATED ( atmo%grd%W2  )  ) DEALLOCATE ( atmo%grd%W2  )
         IF (  ALLOCATED ( atmo%grd%W3  )  ) DEALLOCATE ( atmo%grd%W3  )
         IF (  ALLOCATED ( atmo%grd%P   )  ) DEALLOCATE ( atmo%grd%P   )
         IF (  ALLOCATED ( atmo%grd%Rh  )  ) DEALLOCATE ( atmo%grd%Rh  )
         IF (  ALLOCATED ( atmo%grd%H   )  ) DEALLOCATE ( atmo%grd%H   )
         IF (  ALLOCATED ( atmo%grd%lon )  ) DEALLOCATE ( atmo%grd%lon )
         IF (  ALLOCATED ( atmo%grd%lat )  ) DEALLOCATE ( atmo%grd%lat )
         IF (  ALLOCATED ( atmo%grd%alt )  ) DEALLOCATE ( atmo%grd%alt )
      END IF
!
      RETURN
!
   END SUBROUTINE ATMO_DEALLOCATE
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE ATMO_ALLOCATE
!
!       allocate atmospheric arrays.
!
!****************************
   SUBROUTINE ATMO_ALLOCATE ( what )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = 3 ), INTENT ( IN ) :: what
!
! ---
!
      CALL ATMO_DEALLOCATE ( what )
!
      SELECT CASE ( what )
         CASE ( 'grd' )
            ALLOCATE ( &
               & atmo%grd%C   ( 0:atmo%noflon+1, 0:atmo%noflat+1, 0:atmo%nofh+1 ), &
               & atmo%grd%W1  ( 0:atmo%noflon+1, 0:atmo%noflat+1, 0:atmo%nofh+1 ), &
               & atmo%grd%W2  ( 0:atmo%noflon+1, 0:atmo%noflat+1, 0:atmo%nofh+1 ), &
               & atmo%grd%W3  ( 0:atmo%noflon+1, 0:atmo%noflat+1, 0:atmo%nofh+1 ), &
               & atmo%grd%P   ( 0:atmo%noflon+1, 0:atmo%noflat+1, 0:atmo%nofh+1 ), &
               & atmo%grd%Rh  ( 0:atmo%noflon+1, 0:atmo%noflat+1, 0:atmo%nofh+1 ), &
               & atmo%grd%H   ( 0:atmo%noflon+1, 0:atmo%noflat+1 )               , &
               & atmo%grd%lon ( 0:atmo%noflon+1 )                                , &
               & atmo%grd%lat ( 0:atmo%noflat+1 )                                , &
               & atmo%grd%alt ( 0:atmo%nofh  +1 )                                  &
               & )
               atmo%grd%C=0.d0
               atmo%grd%W1=0.d0
               atmo%grd%W2=0.d0
               atmo%grd%W3=0.d0
               atmo%grd%P=0.d0
               atmo%grd%Rh=0.d0
               atmo%grd%H=0.d0
               atmo%grd%lon=0.d0
               atmo%grd%lat=0.d0
               atmo%grd%alt=0.d0
         CASE ( 'raw' )
            ALLOCATE ( &
               & atmo%raw%H ( atmo%ranges%noflon, atmo%ranges%noflat, atmo%ranges%noflvl ), &
               & atmo%raw%C ( atmo%ranges%noflon, atmo%ranges%noflat, atmo%ranges%noflvl ), &
               & atmo%raw%W1( atmo%ranges%noflon, atmo%ranges%noflat, atmo%ranges%noflvl ), &
               & atmo%raw%W2( atmo%ranges%noflon, atmo%ranges%noflat, atmo%ranges%noflvl ), &
               & atmo%raw%W3( atmo%ranges%noflon, atmo%ranges%noflat, atmo%ranges%noflvl ), &
               & atmo%raw%P ( atmo%ranges%noflon, atmo%ranges%noflat, atmo%ranges%noflvl ), &
               & atmo%raw%Rh( atmo%ranges%noflon, atmo%ranges%noflat, atmo%ranges%noflvl ), &
               & atmo%raw%TOPO( atmo%noflon, atmo%noflat )                                , &
               & atmo%raw%lon ( atmo%ranges%noflon )                                      , &
               & atmo%raw%lat ( atmo%ranges%noflat )                                      , &
               & atmo%raw%lvl ( atmo%ranges%noflvl )                                        &
               & )
               atmo%raw%H=0.d0
               atmo%raw%C=0.d0
               atmo%raw%W1=0.d0
               atmo%raw%W2=0.d0
               atmo%raw%W3=0.d0
               atmo%raw%P=0.d0
               atmo%raw%Rh=0.d0
               atmo%raw%topo=0.d0
               atmo%raw%lat=0.d0
               atmo%raw%lon=0.d0
               atmo%raw%lvl=0.d0
      END SELECT
!
   END SUBROUTINE ATMO_ALLOCATE
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE ATMO_SET_RANGES
!
!       scan ECMWF atmospheric profile from grib file and set dimensions
!
!****************************
   Subroutine atmo_set_ranges ( file1 )
!****************************
!
      use math, only: lonrangef 
      IMPLICIT NONE
!
!.....Dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN )  :: file1
!
!  ---
!
!.....read ranges from grib file
!
      call gribtools_ranges ( file1, atmo%ranges )
!
!.....full zonal atmospheric profile?
!
      atmo%l360 =  ( atmo%ranges%noflon * atmo%ranges%dlon .ge. 360.d0 - atmo%ranges%dlon )
!
!.....copy user input ranges to atmospheric profile ranges
!
      atmo%nofh = INT (  ( atmo%hmax - atmo%hmin ) / atmo%dh  ) + 1
!
!.....longitude ranges
!
      if ( atmo%l360 ) then
         atmo%lonmin =   0.d0
         atmo%lonmax = 360.d0
         atmo%noflon = int ( 360.d0 / atmo%ranges%dlon, kind = 4 ) + 1
      else
         atmo%lonmin = LONRANGEF( atmo%ranges%lonmin )
         atmo%lonmax = atmo%ranges%lonmax
         if ( atmo%lonmax .lt. atmo%lonmin ) atmo%lonmax = LONRANGEF( atmo%lonmax )
         atmo%noflon = atmo%ranges%noflon
      end if
!
! to do: change to lon%(first,last,range,step,shift,count,size)
!
!
      atmo%lonshift  = - atmo%lonmin
      atmo%dlon      = atmo%ranges%dlon
      atmo%lonrange  = (atmo%noflon-1)*atmo%ranges%dlon
!
!.....latitude ranges
!
      atmo%latmin    = atmo%ranges%latmin
      atmo%latmax    = atmo%ranges%latmax
      atmo%dlat      = atmo%ranges%dlat
      atmo%noflat    = atmo%ranges%noflat
      atmo%latrange  = (atmo%noflat-1)*atmo%ranges%dlat
!
      RETURN
!
  End subroutine atmo_set_ranges
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE ATMO_READ_RAW
!
!       read ECMWF atmospheric profile from grib file, grid and regularize.
!
!****************************
   Subroutine atmo_read_raw ( file1, ltopo )
!****************************
!
      implicit none
!
!.....Dummy variables
!
      character ( len = * ), intent ( in )  :: file1
      logical              , intent ( in )  :: ltopo
!
!.....Local variables
!
      integer :: i, j, k, step
      double precision, dimension ( atmo%ranges%noflon, atmo%ranges%noflat, &
         & atmo%ranges%noflvl, atmo%ranges%nofstp, atmo%ranges%nofvar ) :: atmo_tmp_raw
!
!  ---
!
!.....Read raw grib data
!
      if ( atmo%ranges%typeOfLevel .eq. 'hybrid' ) then
         call GRIBTOOLS_READ_ALL_HYBRID ( file1, atmo%ranges, atmo_tmp_raw, atmo%raw%topo, atmo%ens )
      else
         call GRIBTOOLS_READ_ALL_ISOBARIC ( file1, atmo%ranges, atmo_tmp_raw, atmo%raw%topo, atmo%ens )
      endif
!
!.....get step index
!
      step = 1
!
      do i = 1, atmo%ranges%nofstp
         if ( atmo%step .eq. atmo%ranges%stplist( i ) ) then
            step = i
            exit
         end if
      end do
!
      atmo%step = atmo%ranges%stplist ( step )
!
!.....Units in m/s and reverse if required
!
      atmo%raw%h    = atmo_tmp_raw( :, :, :, step, GRIBTOOLS_INDEX( atmo%ranges, 'z' ) ) / 1000
      atmo%raw%topo = atmo%raw%topo / 1000
      where ( atmo%raw%topo .lt. 0.d0 )
         atmo%raw%topo = 0.d0
      end where
!
!.....orography based instead of sea-level based (if topography provided)
!
      if ( .not. ltopo ) then
         do i = 1, atmo%ranges%noflvl
            atmo%raw%h( :, :, i ) = atmo%raw%h( :, :, i ) - atmo%raw%topo
         end do
      end if
!
      atmo%raw%c  = atmo_tmp_raw( :, :, :, step, GRIBTOOLS_INDEX( atmo%ranges, 't' ) )
!
      atmo%raw%W1 = atmo_tmp_raw( :, :, :, step, GRIBTOOLS_INDEX( atmo%ranges, 'u' ) )
      atmo%raw%W2 = atmo_tmp_raw( :, :, :, step, GRIBTOOLS_INDEX( atmo%ranges, 'v' ) )
      if (gribtools_present(atmo%ranges,'w') ) &
         atmo%raw%W3 = atmo_tmp_raw( :, :, :, step, GRIBTOOLS_INDEX( atmo%ranges, 'w' ) )
!
      atmo%raw%p  = atmo_tmp_raw( :, :, :, step, GRIBTOOLS_INDEX( atmo%ranges, 'p' ) )
      atmo%raw%rh = atmo_tmp_raw( :, :, :, step, GRIBTOOLS_INDEX( atmo%ranges, 'q' ) )
!
!.....calculate density from temperature, pressure and humidity
!
      do i = 1, atmo%ranges%noflon
         do j = 1, atmo%ranges%noflat
            do k = 1, atmo%ranges%noflvl
               atmo%raw%rh( i, j, k ) = ATMO_DENSITY( atmo%raw%c( i, j, k ), &
                  & atmo%raw%p( i, j, k ), atmo%raw%rh( i, j, k ) )
            end do
         end do
      end do
!
!.....temperature to speed of sound (dry air, sqrt(1.4*287.05) )
!
      atmo%raw%c = 20.0466954883d0 * dsqrt( atmo%raw%c )
!
!.....reverse wind?
!
      IF ( atmo%lreverse ) THEN
         atmo%raw%W1 = - atmo%raw%W1
         atmo%raw%W2 =   atmo%raw%W2
      ELSE
         atmo%raw%W1 =   atmo%raw%W1
         atmo%raw%W2 = - atmo%raw%W2
      ENDIF
!
!.....Construct vectors
!
      if ( atmo%ranges%dims(1) ) then 
         atmo%raw%lon = (/ ( atmo%ranges%lonmin + i * atmo%ranges%dlon, i = 0, atmo%ranges%noflon - 1 ) /)
      else
         atmo%raw%lon( 1 ) = atmo%ranges%lonmin
      end if
      if ( atmo%ranges%dims(2) ) then 
         atmo%raw%lat = (/ ( atmo%ranges%latmin + i * atmo%ranges%dlat, i = 0, atmo%ranges%noflat - 1 ) /)
      else
         atmo%raw%lat( 1 ) = atmo%ranges%latmin
      end if
      atmo%raw%lvl = (/ ( i, i = 1, atmo%ranges%noflvl ) /)
!
      return
!
  End subroutine atmo_read_raw
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE atmo_raw2grd
!
!       Make regular spaced vectors from irregular input
!
!****************************
   SUBROUTINE atmo_raw2grd ()
!****************************
!
      use cubicweight, only : cw_set_bound2
      use cubicweight, only : cw_set_bound3
      IMPLICIT NONE
!
!.....local variables
!
      integer(kind=4) :: i, ix, iy, ihMin
      integer(kind=4) :: ni, nalt 
!
      logical, dimension(3)  :: circ
!
      integer(kind=4), pointer  :: nd
      double precision, pointer               :: h0
      double precision, pointer, dimension(:) :: hg, cg, w1g, w2g, w3g, pg, rhg ! grd
      double precision, pointer, dimension(:) :: hr, cr, w1r, w2r, w3r, pr, rhr ! raw
      double precision, pointer, dimension(:) :: hi, ci, w1i, w2i, w3i, pi, rhi ! raw
!
! ---
!
!.....take log of pressure to prevent negative values
!
      atmo%raw%p  = dlog( atmo%raw%p  )
      atmo%raw%rh = dlog( atmo%raw%rh )
!
!.....set altitude
!
      nd=>atmo%ranges%noflvl
      nalt=atmo%nofh+1
      atmo%grd%alt=(/ (atmo%hmin-atmo%dh+i*atmo%dh,i=0,nalt)  /)
      hg=>atmo%grd%alt
!
!.....loop over grid
!
      do ix=1, atmo%ranges%noflon
         do iy=1, atmo%ranges%noflat
!
!...........point fun: raw
!
            Hr=>atmo%raw%h(ix,iy,:)
            Cr=>atmo%raw%c(ix,iy,:)
            W1r=>atmo%raw%w1(ix,iy,:)
            W2r=>atmo%raw%w2(ix,iy,:)
            W3r=>atmo%raw%w3(ix,iy,:)
            Pr=>atmo%raw%p(ix,iy,:)
            Rhr=>atmo%raw%rh(ix,iy,:)
!
!...........point fun: grd
!
            Cg=>atmo%grd%c(ix,iy,1:nalt)
            W1g=>atmo%grd%w1(ix,iy,1:nalt)
            W2g=>atmo%grd%w2(ix,iy,1:nalt)
            W3g=>atmo%grd%w3(ix,iy,1:nalt)
            Pg=>atmo%grd%p(ix,iy,1:nalt)
            Rhg=>atmo%grd%rh(ix,iy,1:nalt)
!
!...........find lower index to prevent extrapolation
!
            h0=>hr(1)
            do ihmin=1,nalt
              if (atmo%grd%alt(ihmin).gt.h0) exit
            end do
!
!...........pchip interpolate
!
            ni=nalt-ihmin+1
            hi=>hg(ihmin:nalt)
            ci=>cg(ihmin:nalt)
            w1i=>w1g(ihmin:nalt)
            w2i=>w2g(ihmin:nalt)
            w3i=>w3g(ihmin:nalt)
            pi=>pg(ihmin:nalt)
            rhi=>rhg(ihmin:nalt)
!
!...........Compute the interpolant and evaluate
!
            call atmo_regularize(nd,hr,cr ,ni,hi,ci )
            call atmo_regularize(nd,hr,w1r,ni,hi,w1i)
            call atmo_regularize(nd,hr,w2r,ni,hi,w2i)
            call atmo_regularize(nd,hr,w3r,ni,hi,w3i)
            call atmo_regularize(nd,hr,pr ,ni,hi,pi )
            call atmo_regularize(nd,hr,rhr,ni,hi,rhi)
!
!...........set boundery layer (and copy below)
!
            Cg (1:ihmin-1) = Cg (ihmin+2) - 3*Cg (ihmin+1) + 3*Cg (ihmin)
!            W1g(1:ihmin-1) = W1g(ihmin+2) - 3*W1g(ihmin+1) + 3*W1g(ihmin)
!            W2g(1:ihmin-1) = W2g(ihmin+2) - 3*W2g(ihmin+1) + 3*W2g(ihmin)
!            W3g(1:ihmin-1) = W3g(ihmin+2) - 3*W3g(ihmin+1) + 3*W3g(ihmin)
            Pg (1:ihmin-1) = Pg (ihmin+2) - 3*Pg (ihmin+1) + 3*Pg (ihmin)
            Rhg(1:ihmin-1) = Rhg(ihmin+2) - 3*Rhg(ihmin+1) + 3*Rhg(ihmin)
W1g(1:ihmin-1)=0.d0
W2g(1:ihmin-1)=0.d0
W3g(1:ihmin-1)=0.d0
!
         enddo
!
!........Exit loop early for continuous grid (prevent repeating data)
!
         if ( atmo%l360 .and. abs((ix-1)*atmo%dlon-360.d0+atmo%dlon).lt.1.d-8 ) exit
!
      enddo
!
!........remove log of pressure
!
         atmo%grd%P  = dexp( atmo%grd%P  )
         atmo%grd%Rh = dexp( atmo%grd%Rh )
!
!.....correct for continuous longitude circles
!
      IF ( atmo%l360 ) THEN
!
!........lon = 360 degrees
!
         atmo%grd%C ( atmo%noflon, :, : ) = atmo%grd%C ( 1, :, : )
         atmo%grd%W1( atmo%noflon, :, : ) = atmo%grd%W1( 1, :, : )
         atmo%grd%W2( atmo%noflon, :, : ) = atmo%grd%W2( 1, :, : )
         atmo%grd%W3( atmo%noflon, :, : ) = atmo%grd%W3( 1, :, : )
         atmo%grd%P ( atmo%noflon, :, : ) = atmo%grd%P ( 1, :, : )
         atmo%grd%Rh( atmo%noflon, :, : ) = atmo%grd%Rh( 1, :, : )
!
      ENDIF
!
!.....Fill gaps at side
!
      circ = (/atmo%l360,.false.,.false./)
!
      CALL CW_SET_BOUND3 ( atmo%noflon, atmo%noflat, atmo%nofh, atmo%grd%C , circ )
      CALL CW_SET_BOUND3 ( atmo%noflon, atmo%noflat, atmo%nofh, atmo%grd%w1, circ )
      CALL CW_SET_BOUND3 ( atmo%noflon, atmo%noflat, atmo%nofh, atmo%grd%w2, circ )
      CALL CW_SET_BOUND3 ( atmo%noflon, atmo%noflat, atmo%nofh, atmo%grd%w3, circ )
      CALL CW_SET_BOUND3 ( atmo%noflon, atmo%noflat, atmo%nofh, atmo%grd%P , circ )
      CALL CW_SET_BOUND3 ( atmo%noflon, atmo%noflat, atmo%nofh, atmo%grd%Rh, circ )
!
      atmo%grd%H( 1:atmo%noflon, 1:atmo%noflat ) = atmo%raw%topo
      CALL CW_SET_BOUND2 ( atmo%noflon, atmo%noflat, atmo%grd%H, (/atmo%l360,.FALSE./) )
!
      RETURN
! 
   END SUBROUTINE ATMO_RAW2GRD 
!
! =========================================================================================


subroutine atmo_regularize(nx,x,y,nxi,xi,yi)
!
  integer(kind=4), intent(in)  :: nx, nxi
  real(kind=8), intent(in) , dimension(nx)  :: x, y
  real(kind=8), intent(in) , dimension(nxi) :: xi
  real(kind=8), intent(out), dimension(nxi) :: yi
!
  logical, parameter  :: spline=.false.
  integer(kind=4)     :: ierr, nwk, i, nxi2, lb, ub
  real(kind=8)        :: dxi
  real(kind=8), dimension(nx*2)    :: wk
  real(kind=8), dimension(nx)      :: dx
  real(kind=8), dimension(nxi*2-1), target :: xi2, yi2, dxi2
  real(kind=8), dimension(:), pointer :: yul

  nwk=nx*2
  call dpchez ( nx, x, y, dx, spline, wk, nwk, ierr )
  if ( ierr .lt. 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DPCHEZ - Fatal error!'
    write ( *, '(a,i6)' ) '  DPCHEZ returned error flag IERR= ', ierr
    stop 1
  end if
  nxi2=nxi*2-1
  dxi=(xi(2)-xi(1))/2
  xi2=(/ (xi(1)+i*dxi,i=0,nxi2-1) /)
  call dpchev ( nx, x, y, dx, nxi2, xi2, yi2, dxi2, ierr )
  if ( ierr .lt. 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DPCHEV - Fatal error!'
    write ( *, '(a,i6)' ) '  DPCHEV returned error flag IERR= ', ierr
    stop 1
  end if
! smooth
  yi(1)=yi2(1)
  lb=1; ub=5;
  do i=2,nxi-1
    yul=>yi2(lb:ub)
    yi(i)=sum(yul)/5
    lb=lb+2; ub=ub+2;
  end do
  yi(nxi)=yi2(nxi2)

  return

end subroutine atmo_regularize


! =========================================================================================
!
!.......FUNCTION T2C
!
!       temperature to speed of sound
!
!****************************
   DOUBLE PRECISION FUNCTION T2C ( t )
!****************************
!
      IMPLICIT NONE
!          
!.....dummy variables
!
      DOUBLE PRECISION, INTENT ( IN )         :: t
!
! ----
!          
      t2c = DSQRT ( 1.4d0 * 288.15d0 * ( t + 273.15d0 ) )
!
      RETURN
!
   END FUNCTION T2C
!
! =========================================================================================


! =========================================================================================
!
!.......FUNCTION C2T
!
!       speed of sound to temperature (degrees Kelvin)
!
!****************************
   DOUBLE PRECISION FUNCTION C2K ( c )
!****************************
!
      IMPLICIT NONE
!          
!.....dummy variables
!
      DOUBLE PRECISION, INTENT ( IN )         :: c
!
! ----
!
      c2k = ( c ** 2 ) / 1.4d0 / 288.15d0
!
      RETURN
!
   END FUNCTION C2K
!
! =========================================================================================


! =========================================================================================
!
!.......FUNCTION C2T
!
!       speed of sound to temperature (degrees Celcius)
!
!****************************
   DOUBLE PRECISION FUNCTION C2T ( c )
!****************************
!
      IMPLICIT NONE
!          
!.....dummy variables
!
      DOUBLE PRECISION, INTENT ( IN )         :: c
!
! ----
!
      c2t = ( c ** 2 ) / 1.4d0 / 288.15d0 - 273.15d0
!
      RETURN
!
   END FUNCTION C2T
!
! =========================================================================================


! =========================================================================================
!
!.......FUNCTION ATMO_DENSITY
!
!       Calculate atmosphere density based on temperature, pressure and humidity
!
!****************************
   DOUBLE PRECISION FUNCTION ATMO_DENSITY ( T, P, RH ) RESULT ( Rho )
!****************************
!
      IMPLICIT NONE
!          
!.....dummy variables
!
      DOUBLE PRECISION, INTENT ( IN )     :: T, P, RH
!
!.....local variables
!
      DOUBLE PRECISION                    :: Es, Pv
      DOUBLE PRECISION, DIMENSION ( 10 )  :: Tv
      DOUBLE PRECISION, DIMENSION ( 10 ), PARAMETER :: polyn = &
         & (/ 0.99999683E0, -0.90826951E-2, 0.78736169E-4, -0.61117958E-6, 0.43884187E-8, &
         & -0.29883885E-10, 0.21874425E-12, -0.17892321E-14, 0.11112018E-16, -0.30994571E-19 /)
!
! ----
!
!.....Saturation pressure of water vapour
!
      Tv = (/ 1.d0, T, T ** 2, T ** 3, T ** 4, T ** 5, T ** 6, T ** 7, T ** 8, T ** 9 /)
      Es = 6.1078d0 / (  SUM ( Tv * polyn )  ) ** 8
!
!.....Water vapour pressure
!
      Pv = RH * Es * 100
!
!.....Density of a mixture of dry air molecules and water vapor moleculs
!
      Rho = ABS (  P / 287.05d0 / T * ( 1.d0 - 0.378d0 * Pv / P )  )
!
      RETURN
!
   END FUNCTION ATMO_DENSITY
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE ATMO_ABSORPTION
!
!       Calculate atmospheric absorption based on altitude
!
!****************************
   SUBROUTINE ATMO_ABSORPTION ( h, c, p, rho, frq, unit_km, alpha )
!****************************
!
      USE syms
      IMPLICIT NONE
!          
!.....dummy variables
!
      DOUBLE PRECISION, INTENT ( IN  )  :: h, c, p, rho, frq
      LOGICAL         , INTENT ( IN  )  :: unit_km
      DOUBLE PRECISION, INTENT ( OUT )  :: alpha
!
!.....local parameters
!
      DOUBLE PRECISION, PARAMETER  :: mu0 = 18.192d-6
      DOUBLE PRECISION, PARAMETER  :: p0  = 1.d5
      DOUBLE PRECISION, PARAMETER  :: t0  = 293.15d0
      DOUBLE PRECISION, PARAMETER  :: gam = 1.4d0
      DOUBLE PRECISION, PARAMETER  :: R   = 8.31432d0
!
!.....local variables
!
      DOUBLE PRECISION  :: t, xo, xo2, xo3, xn, xn2, xco2, xon, w, mu, v, zrot, s, n, x
      DOUBLE PRECISION  :: alpha_cl, alpha_diff, alpha_rot, tr, tr0, rpm
      DOUBLE PRECISION, DIMENSION ( 4 )  :: theta, Xvec, Cpinf, Cvinf, fvib, Cr, alpha_vib, amax
!
! ----
!
!.....get atmospheric concentrations
!
      CALL FRAC_MOL ( h, xo, xo2, xo3, xn, xn2, xco2 )
      xon  = ( xo2 + xn2 ) / 0.9903d0
!
      t = c2k ( c )
!
!.....prepare
!
      w    = TWOPI * frq / c
      mu   = 1.4866d-6 * t ** 1.5d0 / ( t + 117.d0 )
      v    = 8.d0 * PI * frq * mu / 3.d0 / p
      zrot = 1.d0 / (   xn2 / (  63.3d0 * DEXP (  - 16.7d0 * t ** ( - 1.d0 / 3.d0 )  )  ) + &
         &   xo2 / 54.1d0 * DEXP ( - 17.3d0 * t ** ( - 1.d0 / 3.d0 )  )   )
      s    = 5.d0 / DSQRT ( 21.d0 )
      n    = 4.d0 / 5.d0 * DSQRT ( 3.d0 / 7.d0 ) * zrot
      x    = .75d0 * n * v
!
!.....viscosity
!
      alpha_cl = w * DSQRT (   &
         & .5d0 * ( DSQRT ( 1 + v ** 2 ) - 1 ) * ( 1 + ( 2.36d0 * x ) ** 2  ) / &
         & (  ( 1 + v ** 2 ) * ( 1 + ( s * 2.36 * x ) ** 2 )  )   )
!
!.....thermal diffusion
!
      alpha_diff = 0.003d0 * alpha_cl
!
!.....rotation 
!
      alpha_rot = w * (  ( s ** 2 - 1 ) * x / 2 / s * DSQRT (  .5d0 * ( DSQRT ( 1 + v ** 2 ) + 1 ) / &
         & (  ( 1 + v ** 2 ) * ( 1 + ( 2.36d0 * x ) ** 2 ) * ( 1 + ( s * 2.36d0 * x ) ** 2 )  )  )   ) * xon
!
!.....relaxation loss ( vibration )
!
      theta = (/ 2239.1d0, 3352.d0, 915.d0, 1037.d0 /)
      xvec  = (/ 0.20948d0, 0.78084d0, 0.000628d0, 7.d-8 /)
      Cpinf = (/ 3.5d0, 3.5d0, 4.d0, 4.d0 /)
      Cvinf = (/ 2.5d0, 2.5d0, 3.d0, 3.d0 /)
!
      rpm = p / p0 * mu0 / mu
      tr0 = ( t0 / t ) ** (   1 / 3.d0 ) - 1
      tr  = ( t / t0 ) ** ( - 1 / 3.d0 ) - 1
!
!........O2
!
         fvib( 1 ) = rpm * (   ( xo2 + xn2 ) * 24 * DEXP ( - 9.16d0 * tr0 ) + &
            & 40400.d0 * DEXP ( 10.d0 * tr ) * ( h + xo3 * 100 ) * &
            & (  0.02d0 * DEXP ( - 11.2d0 * tr ) + ( h + xo3 * 100 )  ) / &
            & (  0.391d0 * DEXP ( 8.41d0 * tr ) + ( h + xo3 * 100 )  ) + &
            & ( xo + xn ) * 2400   )
!
!........N2
!
         fvib( 2 ) = rpm * (  9 * DEXP ( -19.9d0 * tr ) + 28000 * DEXP ( -4.17d0 * tr ) * h / 100 &
            & + 60000 * xo3  )
!
!........CO2
!
         fvib( 3 ) = rpm * (   22000 * DEXP ( -7.68d0 * tr ) * xco2 + &
            & 15100 * DEXP ( -10.4d0 * tr ) * ( xo2 + 0.5d0 * xo ) + &
            & 11500 * DEXP ( -9.17d0 * tr ) * ( xn2 + 0.5d0 * xn ) + &
            & 8.48d8 * DEXP ( 9.17d0 * tr ) * ( h / 100 + xo3 )   )
!
!........O3
!
         fvib( 4 ) = rpm * 1.2d5 * DEXP ( -7.72d0 * tr )
!
!.....Calculate relaxtion loss attenuation coefficients
!
      cr        = ( theta / t ) ** 2 * DEXP ( - theta / t ) / ( 1 - DEXP ( - theta / t )  ) ** 2
      amax      = xvec * PIO2 * cr / ( Cpinf * ( Cvinf + cr ) )
      alpha_vib = amax * 2 * frq ** 2 / (  c * fvib * ( 1 + ( frq / fvib ) ** 2 )  )
!
!.....total absorption coefficient
!
      alpha = alpha_cl + alpha_diff + alpha_rot !+ SUM ( alpha_vib )
!
!.....N/m to N/km
!
      IF ( unit_km ) alpha = alpha * 1000
!
      RETURN
!
   END SUBROUTINE ATMO_ABSORPTION
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE FRAC_MOL
!
!       Get molecular fractions of O, O2, O3, N, N2, CO2
!
!****************************
   SUBROUTINE FRAC_MOL ( h, o, o2, o3, n, n2, co2 )
!****************************
!
      IMPLICIT NONE
!
!.....dummy variables
!
      DOUBLE PRECISION, INTENT ( IN  )  :: h
      DOUBLE PRECISION, INTENT ( OUT )  :: O, O2, O3, N, N2, CO2
!
! ----
!
!.....O
!
      CALL CALC_FRAC (  h, O, 95.d0, O_a1, O_a2  )
!
!.....O2
!
      CALL CALC_FRAC (  h, O2, 90.d0, O2_a1, O2_a2  )
!
!.....O3
!
      CALL CALC_FRAC (  h, O3, 80.d0, O3_a1, O3_a2  )
!
!.....N
!
      CALL CALC_FRAC (  h, N, 180.d0, N_a1, N_a2  )
!
!.....N2
!
      CALL CALC_FRAC (  h, n2, 76.d0, N2_a1, N2_a2  )
!
!.....CO2
!
      CALL CALC_FRAC (  h, co2, 180.d0, CO2_a1, CO2_a2  )
!
      RETURN
!
   END SUBROUTINE FRAC_MOL
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE CALC_FRAC
!
!       Calculate molecular fractions from coefficients
!
!****************************
   SUBROUTINE CALC_FRAC ( h, x, href, a1, a2 )
!****************************
!
      IMPLICIT NONE
!
!.....dummy variables
!
      DOUBLE PRECISION                 , INTENT ( IN  )  :: h, href
      DOUBLE PRECISION, DIMENSION ( 6 ), INTENT ( IN  )  :: a1, a2
      DOUBLE PRECISION                 , INTENT ( OUT )  :: x
!          
!.....local variables
!
      DOUBLE PRECISION, DIMENSION ( 6 )  :: hs    
!
! ----
!
!.....construct vector
!
      hs = (/ 1.d0, h, h ** 2, h ** 3, h ** 4, h ** 5 /)
!
!.....calculate molecular fractions from polynomial coefficients
!
      IF ( h .LE. href ) THEN
         x = 10 ** (  SUM ( a1 * hs )  )
      ELSE
         x = 10 ** (  SUM ( a2 * hs )  )
      END IF
!
      RETURN
!
   END SUBROUTINE CALC_FRAC
!
! =========================================================================================
!
END MODULE TOOLS_ATMO
