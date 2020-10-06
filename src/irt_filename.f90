!*****************************************************************************80
!
!                            I R T _ F I L E N A M E
!
!  Module:       / (include)
!
!  Programmer:   Pieter S. M. Smets
!                R&D depart. of Seismology and Acoustics - Koninklijk Nederlands
!                Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         May 20, 2016
!
!  Language:     Fortran-90
!
!  Description:  Substitute filename variables
!
!
!*****************************************************************************80


!*****************************************************************************80
Function irt_substitute_filename ( string, e, c, t, s, n ) result( filename )
!*****************************************************************************80
!
! Subroutine irt_substitute_filename
!
! Description
!
!   Substitute filename variables, e.g., filename='eigenrays_%T_%Y%M%D_%H'
!
!   Filename variables:
!     %C = centre
!     %T = datatype
!     %Y = year
!     %M = month
!     %D = day of month
!     %O = day of year
!     %H = hour
!     %S = forecast step
!     %N = ensemble number
!
!*****************************************************************************80
!
! Dummy variables
!
  character(len=*), intent(in)  :: string, c, t
  integer(int32)  , intent(in)  :: s, n 
  integer(int64)  , intent(in)  :: e
  character(len=max_len)        :: filename
!
! Local variables
!
  character(len=len_trim(string))  :: pattern
  character(len=2)   :: var 
  integer(int32)     :: i, y, m, d, o, h
!
! Convert epoch
!
  call epoch2time ( epoch=e, year=y, month=m, day=d, h=h )
  call monthday2doy ( month=m, day=d, doy=o, year=y )
!
! Evaluate filename pattern
!
  i=0_int32
  pattern=string
  filename=''
!
  do while (.true.)
    i=scan(pattern,'%',.false.)
    if (i.eq.0_int32) then
      filename=trim(filename)//pattern(1_int32:len_trim(pattern))
      exit
    end if
    filename=trim(filename)//pattern(1_int32:i-1_int32)
    if (len_trim(pattern).lt.i+1_int32) exit
    var=pattern(i:i+1_int32)
    select case (strlowcase(var))
      case ('%c') ! centre
        write(filename,'(2a)') trim(filename), trim(c)
      case ('%t') ! datatype
        write(filename,'(2a)') trim(filename), trim(t)
      case ('%y') ! year
        write(filename,'(a,i4.4)') trim(filename), y
      case ('%m') ! month
        write(filename,'(a,i2.2)') trim(filename), m
      case ('%d') ! day
        write(filename,'(a,i2.2)') trim(filename), d
      case ('%o') ! doy
        write(filename,'(a,i3.3)') trim(filename), o
      case ('%h') ! hour
        write(filename,'(a,i2.2)') trim(filename), h
      case ('%s') ! forecast step
        write(filename,'(a,i3.3)') trim(filename), s
      case ('%n') ! ensemble number
        write(filename,'(a,i2.2)') trim(filename), n
      case default
    end select
    if (i+2_int32.eq.len_trim(pattern)) exit
    pattern=pattern(i+2_int32:len_trim(pattern))
  end do
!
  return
!
!*****************************************************************************80
End function irt_substitute_filename
!*****************************************************************************80
