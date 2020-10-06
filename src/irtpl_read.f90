! =========================================================================================
!
! MODULE IRTPL_READ
!
! --> subroutines to read IRT3D binary data.
!
! =========================================================================================
!
MODULE IRTPL_READ
!
   IMPLICIT NONE
!
!.....public variable raytrace reflections
!
      DOUBLE PRECISION, DIMENSION ( :, : ), ALLOCATABLE, PUBLIC :: store_refl
!
CONTAINS
!
! =========================================================================================
!
!.......Subroutine read reflections
!
!
!****************************
   SUBROUTINE READ_REFL ( prefix, ierr )
!****************************
!
      IMPLICIT NONE
!          
!.....dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN  )  :: prefix
      integer(kind=4)          , INTENT ( OUT )  :: ierr
!          
!.....local variables
!
      LOGICAL           :: lexist
      integer(kind=2)       :: munit = 10
      integer(kind=4)       :: i, j, io, nrec
      DOUBLE PRECISION  :: dummy
!
! ----
!
!.....does the data file exist?
!
      INQUIRE ( FILE = prefix // '_refl.out', EXIST = lexist )
      IF ( .NOT. lexist ) THEN
         PRINT("(A)"), 'WARNING: File ' // prefix // '_refl.out not found.'
         ierr=-1
         RETURN
      END IF
!
!.....open file
!
      OPEN ( UNIT = munit, FILE = prefix // '_refl.out', ACTION = 'READ', STATUS = 'OLD', &
         & FORM = 'UNFORMATTED', ACCESS = 'STREAM' )
!
!.....loop over records
!
      nrec = 0
!
      DO
!
!........read record
!
         READ ( UNIT = munit, IOSTAT = io ) dummy
!
!........error reading file
!
         IF ( io .GT. 0 )  THEN
            PRINT ("(2A)"), 'ERROR: Something wrong with the file ', prefix // '_refl.out'
            ierr = 1
            RETURN
!
!........end of file reached
!
         ELSE IF ( io .LT. 0 ) THEN
!
            EXIT
!
!........record ok
!
         ELSE
!
            nrec = nrec + 1
!
         END IF
!
      END DO
!
!.....check nofrec
!
      IF ( nrec / 7 .GT. 0 ) THEN
!
!........allocate store_refl
!
         ALLOCATE (  store_refl( nrec / 7, 7 )  )
!
!........rewind file
!
         REWIND (munit)
!
!........read data to store_refl
!
         DO i = 1, nrec / 7
            READ ( UNIT = munit ) ( store_refl(i,j), j = 1, 7 )
         END DO
!
      END IF
!
!.....close file
!
      CLOSE ( munit )
!
!.....no error
!
      ierr = 0
!
      RETURN
!
   END SUBROUTINE READ_REFL
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine read rays, select and export 
!
!
!****************************
   SUBROUTINE READ_RAYS ( prefix, rayMin, rayMax, skipRays, lat0, lon0, rk_stepsize, hmax1, angles, verbose1, ierr )
!****************************
!
      USE MATH, ONLY: DISTANCE, SMOOTH
!
      IMPLICIT NONE
!          
!.....dummy variables
!
      CHARACTER ( LEN = * )                    , INTENT ( IN  )  :: prefix
      integer(kind=4)                          , INTENT ( IN  )  :: rayMin, rayMax, skipRays
      DOUBLE PRECISION                         , INTENT ( IN  )  :: lat0, lon0, rk_stepsize, hmax1
      DOUBLE PRECISION     , DIMENSION ( :, : ), INTENT ( IN  )  :: angles
      LOGICAL                                  , INTENT ( IN  )  :: verbose1
      integer(kind=4)                          , INTENT ( OUT )  :: ierr
!          
!.....local variables
!
      LOGICAL                                  :: lexist, lskip, lprint, verbose
      integer(kind=2)                          :: munit, munit2, munit3, munit4
      integer(kind=4)                          :: i, io, rayPrev, cnt, skip, s, nofs, nofrefl
      DOUBLE PRECISION                         :: d, altRefl, hmax, meanRefl
      DOUBLE PRECISION, DIMENSION( 7 )         :: dummy
      DOUBLE PRECISION, DIMENSION( 6, 9999 )   :: store_ray
!
      data munit / 10 /, munit2 / 11 /, munit3 / 12 /, munit4 / 13 /
!
! ----
!
!.....open ray output file
!
      OPEN ( UNIT = munit2, STATUS = 'replace', FILE = prefix // '_rays_v.xy', FORM = 'formatted' )
      OPEN ( UNIT = munit3, STATUS = 'replace', FILE = prefix // '_rays_v_tloss.xy', FORM = 'formatted' )
      OPEN ( UNIT = munit4, STATUS = 'replace', FILE = prefix // '_rays_h.xy', FORM = 'formatted' )
!
!.....does the data file exist?
!
      INQUIRE ( FILE = prefix // '_rays.out', EXIST = lexist)
      IF ( .NOT.lexist ) THEN
         PRINT("(A)"), 'WARNING: File ' // prefix // '_rays.out not found.'
         ierr = - 1
         CLOSE( munit2 )
         RETURN
      ENDIF
!
!.....open file
!
      OPEN ( UNIT = munit, FILE = prefix // '_rays.out', ACTION = 'READ', STATUS = 'OLD', &
         & FORM = 'UNFORMATTED', ACCESS = 'STREAM' )
!
!.....initialize
!
      i       = 0
      rayPrev = rayMin
      cnt     = - 1                ! start with printing first ray (mod(0,skip)=0)
      skip    = skipRays + 1
      lskip   = ( skipRays .GT. 0 ) 
      lprint  = .TRUE.
      nofs    = 0
      store_ray = 0.d0
      nofrefl   = 0
      altRefl   = 0.d0
      meanRefl  = 0.d0

!
!.....Clip rays?
!
      IF ( hmax1 .LE. 0.d0 ) THEN
         verbose = .FALSE.
         hmax    = 9999.d0
      ELSE
         verbose = verbose1
         hmax    = hmax1
      END IF
!
!.....loop over records
!
      DO WHILE ( rayPrev .LE. rayMax )
!
!........read block of 7 records
!
         READ ( UNIT = munit, IOSTAT = io ) ( dummy( i ), i = 1, 7 )
!
!........error reading file
!
         IF ( io .GT. 0 )  THEN
            PRINT ( "(2A)" ), 'ERROR: Something wrong with the file ', TRIM ( prefix ) // '_rays.out'
            ierr = 1
            RETURN
!
!........end of file reached
!
         ELSE IF ( io .LT. 0 ) THEN
!
            EXIT
!
!........record ok
!
         ELSE
!
!...........Check if ray number is ok
!
            IF ( dummy( 1 ) .GE. rayMin .AND. dummy( 1 ) .LE. rayMax ) THEN
!
!..............new ray?
!
               IF ( dummy( 1 ) .GT. rayPrev ) THEN
!
                  IF ( nofs .GT. 0 ) THEN
                     IF ( store_ray( 3, nofs ) .LT. hmax ) THEN
!
                        CALL smooth ( store_ray( 4, 1:nofs ), int ( 25 / rk_stepsize ) )
!
                        if (nofrefl.gt.0) meanRefl=altRefl/nofrefl
                        IF ( verbose ) PRINT "(5x,A,F8.3,A,F7.3,A,F7.3,A,I3,A)", &
                           & 'o plot ray : azi = ', angles( INT ( dummy( 1 ) ), 2 ), &
                           & ', elev = ', angles( INT ( dummy( 1 ) ), 1 ), &
                           & ', avg refl alt = ', meanRefl ,' km (', nofrefl,' bounces)'
!
                        do s = 1, nofs
                           write (munit2,"(f10.3,1x,f9.4,1x,f10.4)") store_ray(2,s), store_ray(3,s), store_ray(4,s) 
                           write (munit3,"(f10.3,1x,f10.4)")         store_ray(2,s), store_ray(4,s)
                           write (munit4,"(f9.4,1x,f8.4,1x,f10.4)")  store_ray(5,s), store_ray(6,s), store_ray(4,s) 
                        end do
!
                       write (munit2,"(a)") '>'
                       write (munit3,"(a)") '>'
                       write (munit4,"(a)") '>'
!
                     END IF
                  END IF
!
                  rayPrev = rayPrev + 1
!
                  IF ( lskip ) THEN
                     cnt    = cnt + 1
                     lprint = ( MOD( cnt, skip ) .EQ. 0 )
                  END IF
!
!.................reset
!
                  nofs      = 0
                  store_ray = 0.d0
                  nofrefl   = 0
                  altRefl   = 0.d0
                  meanRefl  = 0.d0
!
               END IF
!
!..............distance to source and print
!
               IF ( lprint ) THEN
!
                  CALL DISTANCE( lat0, lon0, dummy( 4 ), dummy( 3 ), d )
!
                  nofs = nofs + 1
                  store_ray( :, nofs ) = (/ dummy( 2 ), d, dummy( 5 ), dummy( 6 ), dummy( 3 ), dummy( 4 ) /)
!
                  IF ( abs(dummy(7)).gt.1.d-8 ) THEN
                     nofrefl = nofrefl + 1
                     altRefl = altRefl + dummy( 7 )
                  END IF
!
               END IF
!
            END IF
!
         END IF
!
      END DO
!
!.....close file
!
      CLOSE ( munit  )
      CLOSE ( munit2 )
      CLOSE ( munit3 )
      CLOSE ( munit4 )
!
!.....no error
!
      ierr = 0
!
      RETURN
!
   END SUBROUTINE READ_RAYS
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine read reflections
!
!
!****************************
   SUBROUTINE GET_REFL ( prefix, ierr )
!****************************
!
      IMPLICIT NONE
!          
!.....dummy variables
!
      CHARACTER ( LEN = * ), INTENT ( IN  )  :: prefix
      integer(kind=4)          , INTENT ( OUT )  :: ierr
!          
!.....local variables
!
      LOGICAL                            :: lexist
      integer(kind=2)                    :: munit = 10
      integer(kind=4)                    :: i, io, nofrefl
      DOUBLE PRECISION, DIMENSION ( 7 )  :: dummy
!
! ----
!
!.....does the data file exist?
!
      INQUIRE ( FILE = prefix // '_rays.out', EXIST = lexist )
!
      IF ( lexist ) THEN
         PRINT("(A)"), 'WARNING: File ' // prefix // '_rays.out not found.'
         ierr = 1
         RETURN
      END IF
!
!.....open file
!
      OPEN ( UNIT = munit, FILE = prefix // '_rays.out', ACTION = 'READ', STATUS = 'OLD', &
         & FORM = 'UNFORMATTED', ACCESS = 'STREAM' )
!
!.....initialize
!
      nofrefl = 0
      i       = 0
!
!.....loop over records
!
      DO 
!
!........read block of 7 records
!
         READ ( UNIT = munit, IOSTAT = io) ( dummy( i ), i = 1, 7 )
!
!........error reading file
!
         IF ( io .GT. 0 )  THEN
            PRINT ("(2A)"), 'ERROR: Something wrong with the file ', prefix // '_rays.out'
            ierr = 1
            RETURN
!
!........end of file reached
!
         ELSE IF ( io .LT. 0 ) THEN
!
            EXIT
!
!........record ok
!
         ELSE
!
!...........Check if reflection
!
            IF ( abs(dummy(7)).gt.1.d-8 ) nofrefl = nofrefl + 1
!
         END IF
!
      END DO
!
!.....check nofrefl
!
      IF ( nofrefl .GT. 0 ) THEN
!
!........rewind file
!
         REWIND ( munit )
!
!........allocate reflections
!
         ALLOCATE ( store_refl( nofrefl, 7 ) )
!
!........loop over records
!
         DO WHILE ( i .LT. nofrefl )
!
!...........read block of 9 records
!
            READ ( UNIT = munit, IOSTAT = io ) ( dummy( i ), i = 1, 7 )
!
!...........error reading file
!
            IF ( io .GT. 0 )  THEN
               PRINT ("(2A)"), 'ERROR: Something wrong with the file ', prefix // '_rays.out'
               ierr = 1
               RETURN
!
!...........end of file reached
!
            ELSE IF ( io .LT. 0 ) THEN
!
               EXIT
!
!...........record ok
!
            ELSE
!
!..............Check if reflection
!
               IF ( abs(dummy(7)).gt.1.d-8 ) THEN
                  i = i + 1
                  store_refl( i, : ) = dummy
               END IF
!
            END IF
!
         END DO
!
       ENDIF
!
!.....close file
!
      CLOSE ( munit )
!
!.....no error
!
      ierr = 0
!
      RETURN
!
   END SUBROUTINE GET_REFL
!
! =========================================================================================
!
END MODULE IRTPL_READ
