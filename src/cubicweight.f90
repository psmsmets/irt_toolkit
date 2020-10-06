! =========================================================================================
!
! MODULE CUBICSPLINES 
!
! --> subroutines and functions for globe bathymetry and topography.
!
! =========================================================================================
!
MODULE CUBICWEIGHT 
!
   IMPLICIT NONE
!
!.....Define irrational numbers
!
      DOUBLE PRECISION, PRIVATE, PARAMETER :: ONE_THIRDS  = 1.d0 / 3.d0
      DOUBLE PRECISION, PRIVATE, PARAMETER :: TWO_THIRDS  = 2.d0 / 3.d0
      DOUBLE PRECISION, PRIVATE, PARAMETER :: FOUR_THIRDS = 4.d0 / 3.d0
      DOUBLE PRECISION, PRIVATE, PARAMETER :: ONE_SIXTH   = 1.d0 / 6.d0
!
!.....Set spline coefficients
!
!.....Keys (1981), little sharpening
!
!!      DOUBLE PRECISION, DIMENSION ( 0:3 ), PRIVATE, PARAMETER :: CW_A = (/ 2.d0,  4.d0,  2.5d0,  0.5d0 /)
!!      DOUBLE PRECISION, DIMENSION ( 0:3 ), PRIVATE, PARAMETER :: CW_B = (/ 1.d0,  0.d0, -2.5d0, -1.5d0 /)
!!      DOUBLE PRECISION, DIMENSION ( 0:3 ), PRIVATE, PARAMETER :: CW_C = (/ 1.d0,  0.d0, -2.5d0,  1.5d0 /)
!!      DOUBLE PRECISION, DIMENSION ( 0:3 ), PRIVATE, PARAMETER :: CW_D = (/ 2.d0, -4.d0,  2.5d0, -0.5d0 /)
!
!.....Dirk, smoothening
!
      DOUBLE PRECISION, DIMENSION ( 0:3 ), PRIVATE, PARAMETER :: CW_A = (/ FOUR_THIRDS,  2.d0,  1.d0,  ONE_SIXTH /)
      DOUBLE PRECISION, DIMENSION ( 0:3 ), PRIVATE, PARAMETER :: CW_B = (/ TWO_THIRDS ,  0.d0, -1.d0, -.5d0 /)
      DOUBLE PRECISION, DIMENSION ( 0:3 ), PRIVATE, PARAMETER :: CW_C = (/ TWO_THIRDS ,  0.d0, -1.d0,  .5d0 /)
      DOUBLE PRECISION, DIMENSION ( 0:3 ), PRIVATE, PARAMETER :: CW_D = (/ FOUR_THIRDS, -2.d0,  1.d0, -ONE_SIXTH /)
!
!.....Shift vector
!
      INTEGER, DIMENSION (  4 ), PRIVATE, PARAMETER :: lshift = (/ -1, 0, 1, 2 /)
      INTEGER, DIMENSION (  4 ), PRIVATE, PARAMETER :: rshift = - lshift
      INTEGER, DIMENSION (  4 ), PRIVATE, PARAMETER :: ix1    = (/ 1,2,3,4 /)
      INTEGER, DIMENSION ( 16 ), PRIVATE, PARAMETER :: ix2i   = (/ ix1,ix1,ix1,ix1 /)
      INTEGER, DIMENSION ( 16 ), PRIVATE, PARAMETER :: ix2j   = (/ 1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4 /) 
      INTEGER, DIMENSION ( 64 ), PRIVATE, PARAMETER :: ix3i   = (/ ix2i,ix2i,ix2i,ix2i /)
      INTEGER, DIMENSION ( 64 ), PRIVATE, PARAMETER :: ix3j   = (/ ix2j,ix2j,ix2j,ix2j /)
      INTEGER, DIMENSION ( 64 ), PRIVATE, PARAMETER :: ix3k   = (/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
         & 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4 /)
!
   CONTAINS
!
! =========================================================================================
!
!.......Subroutine EVAL_CW1
!
!       Get cubic weight up and data indices for 1-D data.
!
!****************************
   SUBROUTINE EVAL_CW1 ( val, dat, cw, ix, order, circ )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      INTEGER                                   , INTENT ( IN  )  :: order
      DOUBLE PRECISION, DIMENSION ( 0:order )   , INTENT ( OUT )  :: val
      DOUBLE PRECISION, DIMENSION ( : )         , INTENT ( IN  )  :: dat
      DOUBLE PRECISION, DIMENSION ( 4, 0:order ), INTENT ( IN  )  :: cw
      INTEGER         , DIMENSION ( 4 )         , INTENT ( IN  )  :: ix
      LOGICAL                                   , INTENT ( IN  )  :: circ
!
!.....Local variables
!
      INTEGER           :: i
      DOUBLE PRECISION  :: val_
!
!  ---
!     
      val = 0.d0
!
      DO i = 1, 4
         CALL GETVAL_CW1 ( val_, dat, ix( i ), circ )
         val = val + val_ * cw( i, : )
      END DO
!
      RETURN
!
   END SUBROUTINE EVAL_CW1 
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine EVAL_CW2
!
!       Get cubic weight up and data indices for 2-D data.
!
!****************************
   SUBROUTINE EVAL_CW2 ( val, dat, cw, ix, order, circ )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      INTEGER                                             , INTENT ( IN  )  :: order
      DOUBLE PRECISION, DIMENSION ( 0:order, 0:order )    , INTENT ( OUT )  :: val
      DOUBLE PRECISION, DIMENSION ( :, : )                , INTENT ( IN  )  :: dat
      DOUBLE PRECISION, DIMENSION ( 16, 0:order, 0:order ), INTENT ( IN  )  :: cw
      INTEGER         , DIMENSION ( 16, 2 )               , INTENT ( IN  )  :: ix
      LOGICAL         , DIMENSION ( 2 )                   , INTENT ( IN  )  :: circ
!
!.....Local variables
!
      INTEGER                   :: i
      INTEGER, DIMENSION ( 2 )  :: ix_
      DOUBLE PRECISION          :: val_
!
!  ---
!     
      val = 0.d0
!
      DO i = 1, 16
         ix_ = ix( i, : )
         CALL GETVAL_CW2 ( val_, dat, ix_, circ )
         val = val + val_ * cw( i, :, : )
      END DO
!
      RETURN
!
   END SUBROUTINE EVAL_CW2
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine EVAL_CW3
!
!       Get cubic weight up and data indices for 3-D data.
!
!****************************
   SUBROUTINE EVAL_CW3 ( val, dat, cw, ix, order, circ )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      INTEGER                                                      , INTENT ( IN  )  :: order
      DOUBLE PRECISION, DIMENSION ( 0:order, 0:order, 0:order )    , INTENT ( OUT )  :: val
      DOUBLE PRECISION, DIMENSION ( :, :, : )                      , INTENT ( IN  )  :: dat
      DOUBLE PRECISION, DIMENSION ( 64, 0:order, 0:order, 0:order ), INTENT ( IN  )  :: cw
      INTEGER         , DIMENSION ( 64, 3 )                        , INTENT ( IN  )  :: ix
      LOGICAL         , DIMENSION ( 3 )                            , INTENT ( IN  )  :: circ
!
!.....Local variables
!
      INTEGER                   :: i
      INTEGER, DIMENSION ( 3 )  :: ix_
      DOUBLE PRECISION          :: val_
!
!  ---
!     
      val = 0.d0
!
      DO i = 1, 64
         ix_ = ix( i, : )
         CALL GETVAL_CW3 ( val_, dat, ix_, circ )
         val = val + val_ * cw( i, :, :, : )
      END DO
!
      RETURN
!
   END SUBROUTINE EVAL_CW3
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine GETVAL_CW1
!
!       Get cubic weight up and data indices for 1-D data.
!
!****************************
   SUBROUTINE GETVAL_CW1 ( val, dat, ix, circ )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                 , INTENT ( OUT )  :: val
      DOUBLE PRECISION, DIMENSION ( : ), INTENT ( IN  )  :: dat
      INTEGER                          , INTENT ( IN  )  :: ix
      LOGICAL                          , INTENT ( IN  )  :: circ
!
!.....Local variables
!
      INTEGER  :: i, l
!
!  ---
!
!.....Get dimensions
!
      l = SIZE ( dat, 1 )
!
!.....Get index
!
      i = ix
!
!.....Circular grid?
!
      IF ( circ ) THEN
         IF ( i .EQ. 0     ) i = l
         IF ( i .EQ. l + 1 ) i = 1
      END IF
!
!.....get value from dat based on index
!
      IF ( i .GE. 1 .AND. i .LE. l ) THEN
         val = dat( i )
!
!.....Left-Hand Boundary
!
      ELSE IF ( i .EQ. 0 ) THEN
         val = dat( 3 ) - 3 * dat( 2 ) + 3 * dat( 1 )
!
!.....Right-Hand Boundary
!
      ELSE IF ( i .EQ. l + 1 ) THEN
         val = 3 * dat( l ) - 3 * dat( l - 1 ) + dat( l - 2 )
!
      END IF
!
      RETURN
!
   END SUBROUTINE GETVAL_CW1 
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine GETVAL_CW2
!
!       Get cubic weight up and data indices for 2-D data.
!
!****************************
   SUBROUTINE GETVAL_CW2 ( val, dat, ix, circ )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                    , INTENT ( OUT )  :: val
      DOUBLE PRECISION, DIMENSION ( :, : ), INTENT ( IN  )  :: dat
      INTEGER         , DIMENSION ( 2 )   , INTENT ( IN  )  :: ix
      LOGICAL         , DIMENSION ( 2 )   , INTENT ( IN  )  :: circ
!
!.....Local variables
!
      INTEGER  :: i, j, l, m
!
!  ---
!
!.....Get dimensions
!
      l = SIZE ( dat, 1 )
      m = SIZE ( dat, 2 )
!
!.....Get index
!
      i = ix( 1 )
      j = ix( 2 )
!
!.....Circular grid?
!
      IF ( circ( 1 ) ) THEN
         IF ( i .EQ. 0     ) i = l
         IF ( i .EQ. l + 1 ) i = 1
      END IF
!
      IF ( circ( 2 ) ) THEN
         IF ( j .EQ. 0     ) j = m
         IF ( j .EQ. m + 1 ) j = 1
      END IF
!
!.....get value from dat based on index
!
      IF ( i .GE. 1 .AND. i .LE. l .AND. j .GE. 1 .AND. j .LE. m ) THEN
         val = dat( i, j )
!
!.....Fix corners!
!
      ELSE IF ( i .EQ. 0 .AND. j .EQ. 0 ) THEN
         val = (  dat( 1, 3 )  + dat( 3, 1 ) - 3 * dat( 1, 2 )  - 3 * dat( 2, 1 ) + 6 * dat( 1, 1 )  ) / 2
      ELSE IF ( i .EQ. l + 1 .AND. j .EQ. 0 ) THEN
         val = (  6 * dat( l, 1 ) - 3 * dat( l - 1, 1 ) + dat( l - 2, 1 ) + dat( l, 3 ) - 3 * dat( l, 2 )  ) / 2
      ELSE IF ( i .EQ. 0 .AND. j .EQ. m + 1 ) THEN
         val = (  6 * dat( 1, m ) - 3 * dat( 1, m - 1 ) + dat( 1, m - 2 ) + dat( 3, m ) - 3 * dat( 2, m )  ) / 2
      ELSE IF ( i .EQ. l + 1 .AND. j .EQ. m + 1 ) THEN
         val = (  6 * dat( l, m ) - 3 * dat( l - 1, m ) + dat( l - 2, m ) - 3 * dat( l, m - 1 ) + dat( l, m - 2 )  ) / 2
!
!.....Left-Hand Boundary
!
      ELSE IF ( i .EQ. 0 ) THEN
         val = dat( 3, j ) - 3 * dat( 2, j ) + 3 * dat( 1, j )
      ELSE IF ( j .EQ. 0 ) THEN
         val = dat( i, 3 ) - 3 * dat( i, 2 ) + 3 * dat( i, 1 )
!
!.....Right-Hand Boundary
!
      ELSE IF ( i .EQ. l + 1 ) THEN
         val = 3 * dat( l, j ) - 3 * dat( l - 1, j ) + dat( l - 2, j )
      ELSE IF ( j .EQ. m + 1 ) THEN
         val = 3 * dat( i, m ) - 3 * dat( i, m - 1 ) + dat( i, m - 2 )
!
      END IF
!
      RETURN
!
   END SUBROUTINE GETVAL_CW2
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine GETVAL_CW3
!
!       Get cubic weight up and data indices for 3-D data.
!
!****************************
   SUBROUTINE GETVAL_CW3 ( val, dat, ix, circ )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                       , INTENT ( OUT )  :: val
      DOUBLE PRECISION, DIMENSION ( :, :, : ), INTENT ( IN  )  :: dat
      INTEGER         , DIMENSION ( 3 )      , INTENT ( IN  )  :: ix
      LOGICAL         , DIMENSION ( 3 )      , INTENT ( IN  )  :: circ
!
!.....Local variables
!
      INTEGER  :: i, j, k, l, m, n
!
!  ---
!
!.....Get index
!
      i = ix( 1 )
      j = ix( 2 )
      k = ix( 3 )
!
!.....Get dimensions
!
      l = SIZE ( dat, 1 )
      m = SIZE ( dat, 2 )
      n = SIZE ( dat, 3 )
!
!.....Circular grid?
!
      IF ( circ( 1 ) ) THEN
         IF ( i .EQ. 0     ) i = l
         IF ( i .EQ. l + 1 ) i = 1
      END IF
!
      IF ( circ( 2 ) ) THEN
         IF ( j .EQ. 0     ) j = m
         IF ( j .EQ. m + 1 ) j = 1
      END IF
!
      IF ( circ( 3 ) ) THEN
         IF ( k .EQ. 0     ) k = n
         IF ( k .EQ. n + 1 ) k = 1
      END IF
!
!.....get value from dat based on index
!
      IF ( i .GE. 1 .AND. i .LE. l .AND. j .GE. 1 .AND. j .LE. m .AND. k .GE. 1 .AND. k .LE. n  ) THEN
         val = dat( i, j, k )
!
!.....Fix corners!
!
      ELSE IF ( i .EQ. 0 .AND. j .EQ. 0 .AND. k .EQ. 0 ) THEN
         val = (  dat( 3, 1, 1 ) + dat( 1, 3, 1 ) + dat( 1, 1, 3 ) - 3 * dat( 2, 1, 1 ) &
            & - 3 * dat( 1, 2, 1 )  - 3 * dat( 1, 1, 2 ) + 9 * dat( 1, 1, 1 )  ) / 3
      ELSE IF ( i .EQ. l + 1 .AND. j .EQ. 0 .AND. k .EQ. 0 ) THEN
         val = (  9 * dat( l, 1, 1 ) - 3 * dat( l - 1, 1, 1 ) + dat( l - 2, 1, 1 ) &
            & + dat( l, 3, 1 ) - 3 * dat( l, 2, 1 ) + dat( l, 1, 3 ) - 3 * dat( l, 1, 2 )  ) / 3
      ELSE IF ( i .EQ. 0 .AND. j .EQ. m + 1 .AND. k .EQ. 0 ) THEN
         val = (  dat( 3, m, 1 ) - 3 * dat( 2, m, 1 ) + 3 * dat( 1, m, 1 ) &
            & + 3 * dat( 1, m, 1 ) - 3 * dat( 1, m - 1, 1 ) + dat( 1, m - 2, 1 ) &
            & + dat( 1, m, 3 ) - 3 * dat( 1, m, 2 ) + 3 * dat( 1, m, 1 )  ) / 3
      ELSE IF ( i .EQ. l + 1 .AND. j .EQ. m + 1 .AND. k .EQ. 0 ) THEN
         val = (  3 * dat( l, m, 1 ) - 3 * dat( l - 1, m, 1 ) + dat( l - 2, m, 1 ) &
            & + 3 * dat( l, m, 1 ) - 3 * dat( l, m - 1, 1 ) + dat( l, m - 2, 1 ) &
            & + dat( l, m, 3 ) - 3 * dat( l, m, 2 ) + 3 * dat( l, m, 1 )  ) / 3
      ELSE IF ( i .EQ. 0 .AND. j .EQ. 0 .AND. k .EQ. n + 1 ) THEN
         val = (  dat( 3, 1, n ) - 3 * dat( 2, 1, n ) + 3 * dat( 1, 1, n ) &
            & + dat( 1, 3, n ) - 3 * dat( 1, 2, n ) + 3 * dat( 1, 1, n ) &
            & + 3 * dat( 1, 1, n ) - 3 * dat( 1, 1, n - 1 ) + dat( 1, 1, n - 2 )  ) / 3
      ELSE IF ( i .EQ. l + 1 .AND. j .EQ. 0 .AND. k .EQ. n + 1 ) THEN
         val = (  3 * dat( l, 1, n ) - 3 * dat( l - 1, 1, n ) + dat( l - 2, 1, n ) &
            & + dat( l, 3, n ) - 3 * dat( l, 2, n ) + 3 * dat( l, 1, n ) &
            & + 3 * dat( l, 1, n ) - 3 * dat( l, 1, n - 1 ) + dat( l, 1, n - 2 )  ) / 3
      ELSE IF ( i .EQ. 0 .AND. j .EQ. m + 1 .AND. k .EQ. n + 1 ) THEN
         val = (  dat( 3, m, n ) - 3 * dat( 2, m, n ) + 3 * dat( 1, m, n ) &
            & + 3 * dat( 1, m, n ) - 3 * dat( 1, m - 1, n ) + dat( 1, m - 2, n ) &
            & + 3 * dat( 1, m, n ) - 3 * dat( 1, m, n - 1 ) + dat( 1, m, n - 2 )  ) / 3
      ELSE IF ( i .EQ. l + 1 .AND. j .EQ. m + 1 .AND. k .EQ. n + 1 ) THEN
         val = (  9 * dat( l, m, n ) - 3 * dat( l - 1, m, n ) + dat( l - 2, m, n ) &
            & - 3 * dat( l, m - 1, n ) + dat( l, m - 2, n ) &
            & - 3 * dat( l, m, n - 1 ) + dat( l, m, n - 2 )  ) / 3
!
!.....Offside x
!
      ELSE IF ( j .EQ. 0 .AND. k .EQ. 0 ) THEN
         val = (  dat( i, 3, 1 ) - 3 * dat( i, 2, 1 ) + 3 * dat( i, 1, 1 ) &
            & + dat( i, 1, 3 ) - 3 * dat( i, 1, 2 ) + 3 * dat( i, 1, 1 )  ) / 2
      ELSE IF ( j .EQ. m + 1 .AND. k .EQ. 0 ) THEN
         val = (  3 * dat( i, m, 1 ) - 3 * dat( i, m - 1, 1 ) + dat( i, m - 2, 1 ) &
            & + dat( i, m, 3 ) - 3 * dat( i, m, 2 ) + 3 * dat( i, m, 1 )  ) / 2
      ELSE IF ( j .EQ. 0 .AND. k .EQ. n + 1 ) THEN
         val = (  dat( i, 3, n ) - 3 * dat( i, 2, n ) + 3 * dat( i, 1, n ) &
            & + 3 * dat( i, 1, n ) - 3 * dat( i, 1, n - 1 ) + dat( i, 1, n - 2 )  ) / 2
      ELSE IF ( j .EQ. m + 1 .AND. k .EQ. n + 1 ) THEN
         val = (  3 * dat( i, m, n ) - 3 * dat( i, m - 1, n ) + dat( i, m - 2, n ) &
            & + 3 * dat( i, m, n ) - 3 * dat( i, m, n - 1 ) + dat( i, m, n - 2 )  ) / 2
!
!.....Offside y
!
      ELSE IF ( i .EQ. 0 .AND. k .EQ. 0 ) THEN
         val = (  dat( 3, j, 1 ) - 3 * dat( 2, j, 1 ) + 3 * dat( 1, j, 1 ) &
            & + dat( 1, j, 3 ) - 3 * dat( 1, j, 2 ) + 3 * dat( 1, j, 1 )  ) / 2
      ELSE IF ( i .EQ. l + 1 .AND. k .EQ. 0 ) THEN
         val = (  3 * dat( l, j, 1 ) - 3 * dat( l - 1, j, 1 ) + dat( l - 2, j, 1 ) &
            & + dat( l, j, 3 ) - 3 * dat( l, j, 2 ) + 3 * dat( l, j, 1 )  ) / 2
      ELSE IF ( i .EQ. 0 .AND. k .EQ. n + 1 ) THEN
         val = (  dat( 3, j, n ) - 3 * dat( 2, j, n ) + 3 * dat( 1, j, n ) &
            & + 3 * dat( 1, j, n ) - 3 * dat( 1, j, n - 1 ) + dat( 1, j, n - 2 )  ) / 2
      ELSE IF ( i .EQ. l + 1 .AND. k .EQ. n + 1 ) THEN
         val = (  3 * dat( l, j, n ) - 3 * dat( l - 1, j, n ) + dat( l - 2, j, n ) &
            & + 3 * dat( l, j, n ) - 3 * dat( l, j, n - 1 ) + dat( l, j, n - 2 )  ) / 2
!
!.....Offside z
!
      ELSE IF ( i .EQ. 0 .AND. j .EQ. 0 ) THEN
         val = (  dat( 3, 1, k ) - 3 * dat( 2, 1, k ) + 3 * dat( 1, 1, k ) &
            & + dat( 1, 3, k ) - 3 * dat( 1, 2, k ) + 3 * dat( 1, 1, k )  ) / 2
      ELSE IF ( i .EQ. l + 1 .AND. j .EQ. 0 ) THEN
         val = (  3 * dat( l, 1, k ) - 3 * dat( l - 1, 1, k ) + dat( l - 2, 1, k ) &
            & + dat( l, 3, k ) - 3 * dat( l, 2, k ) + 3 * dat( l, 1, k )  ) / 2
      ELSE IF ( i .EQ. 0 .AND. j .EQ. m + 1 ) THEN
         val = (  dat( 3, 1, k ) - 3 * dat( 2, 1, k ) + 3 * dat( 1, 1, k ) &
            & + 3 * dat( 1, m, k ) - 3 * dat( 1, m - 1, k ) + dat( 1, m - 2, k )  ) / 2
      ELSE IF ( i .EQ. l + 1 .AND. j .EQ. m + 1 ) THEN
         val = (  3 * dat( l, m, k ) - 3 * dat( l - 1, m, k ) + dat( l - 2, m, k ) &
            & + 3 * dat( l, m, k ) - 3 * dat( l, m - 1, k ) + dat( l, m - 2, k )  ) / 2
!
!.....Left-Hand Boundary
!
      ELSE IF ( i .EQ. 0 ) THEN
         val = dat( 3, j, k ) - 3 * dat( 2, j, k ) + 3 * dat( 1, j, k )
      ELSE IF ( j .EQ. 0 ) THEN
         val = dat( i, 3, k ) - 3 * dat( i, 2, k ) + 3 * dat( i, 1, k )
      ELSE IF ( k .EQ. 0 ) THEN
         val = dat( i, j, 3 ) - 3 * dat( i, j, 2 ) + 3 * dat( i, j, 1 )
!
!.....Right-Hand Boundary
!
      ELSE IF ( i .EQ. l + 1 ) THEN
         val = 3 * dat( l, j, k ) - 3 * dat( l - 1, j, k ) + dat( l - 2, j, k )
      ELSE IF ( j .EQ. m + 1 ) THEN
         val = 3 * dat( i, m, k ) - 3 * dat( i, m - 1, k ) + dat( i, m - 2, k )
      ELSE IF ( k .EQ. n + 1 ) THEN
         val = 3 * dat( i, j, n ) - 3 * dat( i, j, n - 1 ) + dat( i, j, n - 2 )
!
      END IF
!
      RETURN
!
   END SUBROUTINE GETVAL_CW3
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CALC_CW1
!
!       Get cubic weight up and data indices for 1-D data.
!
!****************************
   SUBROUTINE CALC_CW1 ( x, cw, ix, nofx, dx, order, ierr, scalar )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                          , INTENT ( IN  )  :: x, dx
      INTEGER                                   , INTENT ( IN  )  :: order, nofx
      DOUBLE PRECISION, DIMENSION ( 4, 0:order ), INTENT ( OUT )  :: cw
      INTEGER         , DIMENSION ( 4 )         , INTENT ( OUT )  :: ix
      INTEGER                                   , INTENT ( OUT )  :: ierr
      DOUBLE PRECISION                          , INTENT ( IN  ), OPTIONAL  :: scalar

!
!.....Local variables
!
      INTEGER                                :: i, u_
      INTEGER         , DIMENSION ( 4 )      :: l
      DOUBLE PRECISION                       :: u
      DOUBLE PRECISION, DIMENSION ( 4, 0:2 ) :: cw_i
      DOUBLE PRECISION, DIMENSION ( 0:2 )    :: dx_
!
!  ---
!     
      u  = x / dx + 1
      u_ = FLOOR ( u )
!
!.....check boundaries
!
      IF ( u_ .LT. 1 .OR. u_ .GT. nofx ) THEN
         ierr = 1
         RETURN
      END IF
!
      IF ( u_ .GE. nofx - 1 ) THEN
         u_   = CEILING ( u )
         l    = u_ + rshift
      ELSE
         l    = u_ + lshift
      END IF
!
!.....copy index
!
      ix = l
!
!.....calculate stepsize for derivative
!
      dx_ = (/ 1.d0, dx, dx ** 2 /)
!
      IF (  PRESENT ( scalar )  ) THEN
         dx_ = dx_ * (/ 1.d0, scalar, scalar ** 2 /)
      END IF
!
!.....get coefficients
!
      DO i = 1, 4
         CALL CW_COEFF2 ( u - DBLE( l( i ) ), cw_i( i , : ) )
         cw_i( i, : ) = cw_i( i, : ) / dx_ 
      END DO
!
      cw = cw_i( :, 0:order )
!
      ierr = 0
!
      RETURN
!
   END SUBROUTINE CALC_CW1 
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CALC_CW2
!
!       Get cubic weight and data indices for 2-D data.
!
!****************************
   SUBROUTINE CALC_CW2 ( x, y, cw, ix, nofx, dx, nofy, dy, order, ierr, scalar )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                                    , INTENT ( IN  )  :: x, dx, y, dy
      INTEGER                                             , INTENT ( IN  )  :: order, nofx, nofy
      DOUBLE PRECISION, DIMENSION ( 16, 0:order, 0:order ), INTENT ( OUT )  :: cw
      INTEGER         , DIMENSION ( 16, 2 )               , INTENT ( OUT )  :: ix
      INTEGER                                             , INTENT ( OUT )  :: ierr
      DOUBLE PRECISION, DIMENSION ( 2 )                   , INTENT ( IN  ), OPTIONAL  :: scalar

!
!.....Local variables
!
      INTEGER                                 :: i, o1, o2, u_, v_
      INTEGER         , DIMENSION ( 4 )       :: l, m
      DOUBLE PRECISION                        :: u, v
      DOUBLE PRECISION, DIMENSION ( 4, 0:2 )  :: cw_i, cw_j
      DOUBLE PRECISION, DIMENSION ( 0:2 )     :: dx_, dy_, cw_tmp
      DOUBLE PRECISION, DIMENSION ( 16, 0:2, 0:2 ) :: cw_
!
!  ---
!     
      u  = x / dx + 1
      u_ = FLOOR ( u )
!     
      v  = y / dy + 1
      v_ = FLOOR ( v )
!
!.....check boundaries
!
      IF ( u_ .LT. 1 .OR. u_ .GT. nofx .OR. &
         & v_ .LT. 1 .OR. v_ .GT. nofy ) THEN
         ierr = 1
         RETURN
      END IF
!
      IF ( u_ .GE. nofx - 1 ) THEN
         u_   = CEILING ( u )
         l    = u_ + rshift
      ELSE
         l    = u_ + lshift
      END IF
!
      IF ( v_ .GE. nofy - 1 ) THEN
         v_   = CEILING ( v )
         m    = v_ + rshift
      ELSE
         m    = v_ + lshift
      END IF
!
!.....calculate stepsize for derivative
!
      dx_ = (/ 1.d0, dx, dx ** 2 /)
      dy_ = (/ 1.d0, dy, dy ** 2 /)
!
      IF (  PRESENT ( scalar )  ) THEN
         dx_ = dx_ * (/ 1.d0, scalar( 1 ), scalar( 1 ) ** 2 /)
         dy_ = dy_ * (/ 1.d0, scalar( 2 ), scalar( 2 ) ** 2 /)
      END IF
!
!.....get coefficients
!
      DO i = 1, 4
         CALL CW_COEFF2 ( u - DBLE( l( i ) ), cw_tmp ) 
         cw_i( i, : ) = cw_tmp / dx_ 
         CALL CW_COEFF2 ( v - DBLE( m( i ) ), cw_tmp ) 
         cw_j( i, : ) = cw_tmp / dy_ 
      END DO
!
!.....apply convolution
!
      ix( :, 1 ) = l( ix2i )
      ix( :, 2 ) = m( ix2j )
!
      cw_ = 0.d0
!
      DO o1 = 0, 1
         DO o2 = 0, 1
            cw_( :, o1, o2 ) =  cw_i( ix2i, o1 ) * cw_j( ix2j, o2 )
         END DO
      END DO
!
      cw_( :, 2, 0 ) = cw_i( ix2i, 2 ) * cw_j( ix2j, 0 )
      cw_( :, 0, 2 ) = cw_i( ix2i, 0 ) * cw_j( ix2j, 2 )
!
      cw = cw_( :, 0:order, 0:order )
!
      ierr = 0
!
      RETURN
!
   END SUBROUTINE CALC_CW2
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CALC_CW3
!
!       Get cubic weight up and data indices for 3-D data.
!
!****************************
   SUBROUTINE CALC_CW3 ( x, y, z, cw, ix, nofx, dx, nofy, dy, nofz, dz, order, ierr, scalar )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                                             , INTENT ( IN  )  :: x, dx, y, dy, z, dz
      INTEGER                                                      , INTENT ( IN  )  :: order, nofx, nofy, nofz
      DOUBLE PRECISION, DIMENSION ( 64, 0:order, 0:order, 0:order ), INTENT ( OUT )  :: cw
      INTEGER         , DIMENSION ( 64, 3 )                        , INTENT ( OUT )  :: ix
      INTEGER                                                      , INTENT ( OUT )  :: ierr
      DOUBLE PRECISION, DIMENSION ( 3 )                            , INTENT ( IN  ), OPTIONAL  :: scalar
!
!.....Local variables
!
      INTEGER                                 :: i, o1, o2, o3, u_, v_, w_
      INTEGER         , DIMENSION ( 4 )       :: l, m, n
      DOUBLE PRECISION                        :: u, v, w
      DOUBLE PRECISION, DIMENSION ( 4, 0:2 )  :: cw_i, cw_j, cw_k
      DOUBLE PRECISION, DIMENSION ( 0:2 )     :: dx_, dy_, dz_, cw_tmp
      DOUBLE PRECISION, DIMENSION ( 64, 0:2, 0:2, 0:2 ) :: cw_
!
!  ---
!
      u  = x / dx + 1
      u_ = FLOOR ( u )
!     
      v  = y / dy + 1
      v_ = FLOOR ( v )
!     
      w  = z / dz + 1
      w_ = FLOOR ( w )
!
!.....check boundaries
!
      IF ( u .LT. 1 .OR. u .GE. nofx .OR. &
         & v .LT. 1 .OR. v .GE. nofy .OR. &
         & w .GE. nofz ) THEN
         ierr = 1
         RETURN
      ELSE IF ( w_ .LT. 1 ) THEN
         ierr = -1
         RETURN
      END IF
!
      IF ( u_ .GE. nofx - 1 ) THEN
         u_   = CEILING ( u )
         l    = u_ + rshift
      ELSE
         l    = u_ + lshift
      END IF
!
      IF ( v_ .GE. nofy - 1 ) THEN
         v_   = CEILING ( v )
         m    = v_ + rshift
      ELSE
         m    = v_ + lshift
      END IF
!
      IF ( w_ .GE. nofz - 1 ) THEN
         w_   = CEILING ( w )
         n    = w_ + rshift
      ELSE
         n    = w_ + lshift
      END IF
!
!.....calculate stepsize for derivative
!
      dx_ = (/ 1.d0, dx, dx ** 2 /)
      dy_ = (/ 1.d0, dy, dy ** 2 /)
      dz_ = (/ 1.d0, dz, dz ** 2 /)
!
      IF (  PRESENT ( scalar )  ) THEN
         dx_ = dx_ * (/ 1.d0, scalar( 1 ), scalar( 1 ) ** 2 /)
         dy_ = dy_ * (/ 1.d0, scalar( 2 ), scalar( 2 ) ** 2 /)
         dz_ = dz_ * (/ 1.d0, scalar( 3 ), scalar( 3 ) ** 2 /)
      END IF
!
!.....get coefficients
!
      DO i = 1, 4
         CALL CW_COEFF2 ( u - DBLE( l( i ) ), cw_tmp ) 
         cw_i( i, : ) = cw_tmp / dx_ 
         CALL CW_COEFF2 ( v - DBLE( m( i ) ), cw_tmp ) 
         cw_j( i, : ) = cw_tmp / dy_ 
         CALL CW_COEFF2 ( w - DBLE( n( i ) ), cw_tmp ) 
         cw_k( i, : ) = cw_tmp / dz_ 
      END DO
!
!.....apply convolution
!
      ix( :, 1 ) = l( ix3i )
      ix( :, 2 ) = m( ix3j )
      ix( :, 3 ) = n( ix3k )
!
      cw_ = 0.d0
!
      DO o1 = 0, 1
         DO o2 = 0, 1
            DO o3 = 0, 1
               cw_( :, o1, o2, o3 ) =  cw_i( ix3i, o1 ) * cw_j( ix3j, o2 ) * cw_k( ix3k, o3 )
            END DO
         END DO
      END DO
!
      cw_( :, 2, 0, 0 ) = cw_i( ix3i, 2 ) * cw_j( ix3j, 0 ) * cw_k( ix3k, 0 )
      cw_( :, 0, 2, 0 ) = cw_i( ix3i, 0 ) * cw_j( ix3j, 2 ) * cw_k( ix3k, 0 )
      cw_( :, 0, 0, 2 ) = cw_i( ix3i, 0 ) * cw_j( ix3j, 0 ) * cw_k( ix3k, 2 )
!
      cw = cw_( :, 0:order, 0:order, 0:order )
!
      ierr = 0
!
      RETURN
!
   END SUBROUTINE CALC_CW3
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CW_COEFF0
!
!       Get cubic weight up to derivative order 0.
!
!****************************
   SUBROUTINE CW_COEFF0 ( t, cw )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION, INTENT ( IN  )  :: t
      DOUBLE PRECISION, INTENT ( OUT )  :: cw
!
!  ---
!
      IF ( t .GE. -2.d0 .AND. t .LT. -1.d0 ) THEN
         CALL cw_0 ( t, cw, CW_A )
      ELSE IF ( t .GE. -1.d0 .AND. t .LT. 0.d0 ) THEN
         CALL cw_0 ( t, cw, CW_B )
      ELSE IF ( t .GE. 0.d0 .AND. t .LT. 1.d0 ) THEN
         CALL cw_0 ( t, cw, CW_C )
      ELSE IF ( t .GE. 1.d0 .AND. t .LE. 2.d0 ) THEN
         CALL cw_0 ( t, cw, CW_D )
      ELSE
         cw = 0.d0
      ENDIF
!
      RETURN
!
   END SUBROUTINE CW_COEFF0
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CW_COEFF1
!
!       Get cubic weight up to derivative order 1.
!
!****************************
   SUBROUTINE CW_COEFF1 ( t, cw )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                   , INTENT ( IN  )  :: t
      DOUBLE PRECISION, DIMENSION ( 0:1 ), INTENT ( OUT )  :: cw
!
!  ---
!
      IF ( t .GE. -2.d0 .AND. t .LT. -1.d0 ) THEN
         CALL cw_1 ( t, cw, CW_A )
      ELSE IF ( t .GE. -1.d0 .AND. t .LT. 0.d0 ) THEN
         CALL cw_1 ( t, cw, CW_B )
      ELSE IF ( t .GE. 0.d0 .AND. t .LT. 1.d0 ) THEN
         CALL cw_1 ( t, cw, CW_C )
      ELSE IF ( t .GE. 1.d0 .AND. t .LE. 2.d0 ) THEN
         CALL cw_1 ( t, cw, CW_D )
      ELSE
         cw = 0.d0
      ENDIF

!
      RETURN
!
   END SUBROUTINE CW_COEFF1
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CW_COEFF2
!
!       Get cubic weight up to derivative order 2.
!
!****************************
   SUBROUTINE CW_COEFF2 ( t, cw )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                   , INTENT ( IN  )  :: t
      DOUBLE PRECISION, DIMENSION ( 0:2 ), INTENT ( OUT )  :: cw
!
!  ---
!
      IF ( t .GE. -2.d0 .AND. t .LT. -1.d0 ) THEN
         CALL cw_2 ( t, cw, CW_A )
      ELSE IF ( t .GE. -1.d0 .AND. t .LT. 0.d0 ) THEN
         CALL cw_2 ( t, cw, CW_B )
      ELSE IF ( t .GE. 0.d0 .AND. t .LT. 1.d0 ) THEN
         CALL cw_2 ( t, cw, CW_C )
      ELSE IF ( t .GE. 1.d0 .AND. t .LE. 2.d0 ) THEN
         CALL cw_2 ( t, cw, CW_D )
      ELSE
         cw = 0.d0
      ENDIF
!
      RETURN
!
   END SUBROUTINE CW_COEFF2
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CW_0
!
!       Get cubic spline weight from coefficients up to derivative order 0.
!
!****************************
   SUBROUTINE CW_0 ( t, cw, coeff )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                   , INTENT ( IN  )  :: t
      DOUBLE PRECISION, DIMENSION ( 0:3 ), INTENT ( IN  )  :: coeff
      DOUBLE PRECISION                   , INTENT ( OUT )  :: cw
!
!  ---
!
      cw = coeff( 0 ) + t * ( coeff( 1 ) + t * ( coeff( 2 ) + t * coeff( 3 ) ) )
!
      RETURN
!
   END SUBROUTINE CW_0
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CW_1
!
!       Get cubic spline weight from coefficients up to derivative order 1.
!
!****************************
   SUBROUTINE CW_1 ( t, cw, coeff )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                   , INTENT ( IN  )  :: t
      DOUBLE PRECISION, DIMENSION ( 0:3 ), INTENT ( IN  )  :: coeff
      DOUBLE PRECISION, DIMENSION ( 0:1 ), INTENT ( OUT )  :: cw
!
!  ---
!
      cw( 1 ) = coeff( 3 ); cw( 0 ) = coeff( 3 )
!
      cw( 0 ) = cw( 0 ) * t + coeff( 2 )
      cw( 1 ) = cw( 1 ) * t + cw( 0 )
      cw( 0 ) = cw( 0 ) * t + coeff( 1 )
      cw( 1 ) = cw( 1 ) * t + cw( 0 )
      cw( 0 ) = cw( 0 ) * t + coeff( 0 )
!
      RETURN
!
   END SUBROUTINE CW_1
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CW_2
!
!       Get cubic spline weight from coefficients up to derivative order 2.
!
!****************************
   SUBROUTINE CW_2 ( t, cw, coeff )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      DOUBLE PRECISION                   , INTENT ( IN  )  :: t
      DOUBLE PRECISION, DIMENSION ( 0:3 ), INTENT ( IN  )  :: coeff
      DOUBLE PRECISION, DIMENSION ( 0:2 ), INTENT ( OUT )  :: cw
!
!  ---
!
      cw( 2 ) = coeff( 3 ); cw( 1 ) = coeff( 3 ); cw( 0 ) = coeff( 3 )
!
      cw( 0 ) = cw( 0 ) * t + coeff( 2 )
      cw( 1 ) = cw( 1 ) * t + cw( 0 )
      cw( 0 ) = cw( 0 ) * t + coeff( 1 )
      cw( 2 ) = cw( 2 ) * t + cw( 1 )
      cw( 1 ) = cw( 1 ) * t + cw( 0 )
      cw( 0 ) = cw( 0 ) * t + coeff( 0 )
      cw( 2 ) = cw( 2 ) * 2.d0
!
      RETURN
!
   END SUBROUTINE CW_2
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CW_SET_BOUND1
!
!       Set boundary values for a 1-D data set.
!
!****************************
   SUBROUTINE CW_SET_BOUND1 ( l, dat, circ )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      INTEGER                              , INTENT ( IN    )  :: l
      DOUBLE PRECISION, DIMENSION ( 0:l+1 ), INTENT ( INOUT )  :: dat
      LOGICAL                              , INTENT ( IN    )  :: circ
!
!  ---
!
!.....Circular repetition? Copy! Otherwise fix left and right hand side boundary values.
!
      IF ( circ ) THEN
         dat( 1   ) = dat( l )
         dat( l+1 ) = dat( 1 )
      ELSE
         dat(   0 ) = dat( 3 ) - 3 * dat( 2 ) + 3 * dat( 1 )
         dat( l+1 ) = 3 * dat( l ) - 3 * dat( l - 1 ) + dat( l - 2 )
      END IF
!
      RETURN
!
   END SUBROUTINE CW_SET_BOUND1
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CW_SET_BOUND2
!
!       Set boundary values for a 2-D data set.
!
!****************************
   SUBROUTINE CW_SET_BOUND2 ( l, m, dat, circ )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      INTEGER                                     , INTENT ( IN    )  :: l, m
      DOUBLE PRECISION, DIMENSION ( 0:l+1, 0:m+1 ), INTENT ( INOUT )  :: dat
      LOGICAL         , DIMENSION ( 2 )           , INTENT ( IN    )  :: circ
!
!  ---
!
!.....Circular repetition? Copy! Otherwise fix left and right hand side boundary values.
!
      IF ( circ( 1 ) ) THEN
         dat( 1  , : ) = dat( l, : )
         dat( l+1, : ) = dat( 1, : )
      ELSE
         dat(   0, : ) = dat( 3, : ) - 3 * dat( 2, : ) + 3 * dat( 1, : )
         dat( l+1, : ) = 3 * dat( l, : ) - 3 * dat( l - 1, : ) + dat( l - 2, : )
      END IF
!
      IF ( circ( 2 ) ) THEN
         dat( :, 1   ) = dat( :, m )
         dat( :, m+1 ) = dat( :, 1 )
      ELSE
         dat( :,   0 ) = dat( :, 3 ) - 3 * dat( :, 2 ) + 3 * dat( :, 1 )
         dat( :, m+1 ) = 3 * dat( :, m ) - 3 * dat( :, m - 1 ) + dat( :, m - 2 )
      END IF
!
!.....Fix corners!
!
      dat( 0, 0 ) = (  dat( 1    , 0 ) + dat( 0, 1     )  ) / 2
      dat( l, 0 ) = (  dat( l - 1, 0 ) + dat( l, 1     )  ) / 2
      dat( 0, m ) = (  dat( 1    , m ) + dat( 0, m - 1 )  ) / 2
      dat( l, m ) = (  dat( l - 1, m ) + dat( l, m - 1 )  ) / 2
!
      RETURN
!
   END SUBROUTINE CW_SET_BOUND2
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine CW_SET_BOUND3
!
!       Set boundary values for a 3-D data set.
!
!****************************
   SUBROUTINE CW_SET_BOUND3 ( l, m, n, dat, circ )
!****************************
!
      IMPLICIT NONE
!
!.....Dummy variables
!
      INTEGER                                            , INTENT ( IN    )  :: l, m, n
      DOUBLE PRECISION, DIMENSION ( 0:l+1, 0:m+1, 0:n+1 ), INTENT ( INOUT )  :: dat
      LOGICAL         , DIMENSION ( 3 )                  , INTENT ( IN    )  :: circ
!
!  ---
!
!.....Circular repetition? Copy! Otherwise fix left and right hand side boundary values.
!
      IF ( circ( 1 ) ) THEN
         dat(   0, :, : ) = dat( l, :, : )
         dat( l+1, :, : ) = dat( 1, :, : )
      ELSE
         dat(   0, :, : ) = dat( 3, :, : ) - 3 * dat( 2, :, : ) + 3 * dat( 1, :, : )
         dat( l+1, :, : ) = 3 * dat( l, :, : ) - 3 * dat( l - 1, :, : ) + dat( l - 2, :, : )
      END IF
!
      IF ( circ( 2 ) ) THEN
         dat( :,   0, : ) = dat( :, m, : )
         dat( :, m+1, : ) = dat( :, 1, : )
      ELSE
         dat( :,   0, : ) = dat( :, 3, : ) - 3 * dat( :, 2, : ) + 3 * dat( :, 1, : )
         dat( :, m+1, : ) = 3 * dat( :, m, : ) - 3 * dat( :, m - 1, : ) + dat( :, m - 2, : )
      END IF
!
      IF ( circ( 3 ) ) THEN
         dat( :, :,   0 ) = dat( :, :, n )
         dat( :, :, n+1 ) = dat( :, :, 1 )
      ELSE
         dat( :, :,   0 ) = dat( :, :, 3 ) - 3 * dat( :, :, 2 ) + 3 * dat( :, :, 1 )
         dat( :, :, n+1 ) = 3 * dat( :, :, n ) - 3 * dat( :, :, n - 1 ) + dat( :, :, n - 2 )
      END IF
!
!.....Fix corners!
!
      dat( 0  , 0  , 0   ) = (  dat( 1    , 0, 0 ) + dat( 0, 1    , 0 ) + dat( 0, 0, 1     )  ) / 3
      dat( l+1, 0  , 0   ) = (  dat( l - 1, 0, 0 ) + dat( l, 1    , 0 ) + dat( l, 0, 1     )  ) / 3
      dat( 0  , m+1, 0   ) = (  dat( 1    , m, 0 ) + dat( 0, m - 1, 0 ) + dat( 0, m, 1     )  ) / 3
      dat( 0  , 0  , n+1 ) = (  dat( 1    , 0, n ) + dat( 0, 1    , n ) + dat( 0, 0, n - 1 )  ) / 3
      dat( l+1, m+1, 0   ) = (  dat( l - 1, m, 0 ) + dat( l, m - 1, 0 ) + dat( l, m, 1     )  ) / 3
      dat( l+1, 0  , n+1 ) = (  dat( l - 1, 0, n ) + dat( l, 1    , n ) + dat( l, 0, n - 1 )  ) / 3
      dat( 0  , m+1, n+1 ) = (  dat( 1    , m, n ) + dat( 0, m - 1, n ) + dat( 0, m, n - 1 )  ) / 3
      dat( l+1, m+1, n+1 ) = (  dat( l - 1, m, n ) + dat( l, m - 1, n ) + dat( l, m, n - 1 )  ) / 3
!
      RETURN
!
   END SUBROUTINE CW_SET_BOUND3
!
! =========================================================================================
!
!
END MODULE CUBICWEIGHT
