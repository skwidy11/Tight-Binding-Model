! Subroutine to swap two vectors without using auxiliary storage
!
! Input:  |b>, |a>
!
! Output: |a>, |b>
!
! ND    = actual number of sites
! veca  = array with first vector
! vecb  = array with second vector
!
! E. Mucciolo, UCF (Apr 01, 2014)
!
  subroutine swap ( ND, norb, veca, vecb )
!
  implicit none
!
  integer, parameter :: double = KIND(0.d0)
  integer, intent(in) :: ND, norb
  complex (kind=double), intent(out), dimension(ND, norb) :: veca, vecb
!
  veca = veca + vecb
  vecb = veca - vecb
  veca = veca - vecb
!
  end subroutine swap
