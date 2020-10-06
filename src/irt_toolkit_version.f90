Subroutine irt_toolkit_package_info ( package_version, package_fullname )
!*****************************************************************************80
!
!! IRT_TOOLKIT_PACKAGE_INFO
!
!  Description:
!
!    Get irt_toolkit package info  
!
!  Modified:
!
!    16 November 2016
!
!  Author:
!
!    Pieter Smets
!
!  Parameters:
!
!
!****************************
!
  implicit none
!
! Dummy variable
!
  character(len=5) , intent(out), optional  :: package_version
  character(len=25), intent(out), optional  :: package_fullname
!
! Local variables
!
  character(len=20)  :: package_name
  character(len=5)   :: irt_version
  character(len=2)   :: irt_major_version_number
  character(len=1)   :: irt_minor_version_number, irt_revision_version_number
  integer(int16)     :: irt_major_version, irt_minor_version, irt_revision_version

  include '../version.sh'

  write(irt_major_version_number,"(i2)") irt_major_version
  write(irt_minor_version_number,"(i1)") irt_minor_version
  write(irt_revision_version_number,"(i1)") irt_revision_version

  write(irt_version,"(a,'.',2a)") trim(adjustl(irt_major_version_number)), &
    irt_minor_version_number, irt_revision_version_number

  if (present(package_version)) package_version=irt_version
  if (present(package_fullname)) write(package_fullname,"(3a)") &
    trim(package_name), ' v', trim(irt_version)
!
  return
!
!*****************************************************************************80
End subroutine irt_toolkit_package_info


Function irt_toolkit_version() result(version)
!*****************************************************************************80
!
!! IRT_TOOLKIT_VERSION
!
!  Description:
!
!    Get irt_toolkit version number  
!
!  Modified:
!
!    20 March 2017
!
!  Author:
!
!    Pieter Smets
!
!  Parameters:
!
!
!****************************
!
  implicit none
  character(len=5) :: version
  call irt_toolkit_package_info(package_version=version)
!
!*****************************************************************************80
End function irt_toolkit_version


Function irt_toolkit_fullname() result(fullname)
!*****************************************************************************80
!
!! IRT_TOOLKIT_FULLNAME
!
!  Description:
!
!    Get irt_toolkit full version name  
!
!  Modified:
!
!    20 March 2017
!
!  Author:
!
!    Pieter Smets
!
!  Parameters:
!
!
!****************************
!
  implicit none
  character(len=25) :: fullname
  call irt_toolkit_package_info(package_fullname=fullname)
!
!*****************************************************************************80
End function irt_toolkit_fullname
