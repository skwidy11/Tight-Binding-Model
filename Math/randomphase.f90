! Subroutine to generate an array of random phases.
!
! E. Mucciolo, UCF (Apr 01, 2014)
! modified P. Schelling, UCF (June 8, 2018)
!
  subroutine randomphase ( ND, norb, theta )
!
  implicit none
!
  integer, parameter :: double = KIND(0.d0)
  real (kind=double), parameter :: pi = 3.14159265d0
!
  integer, intent(in) :: ND
  real (kind=double), intent(out), dimension(ND*norb) :: theta
!
  real (kind=double) :: phase
  integer  k, norb
!
! generate the random phases (uniformly distributed)
!
  do k = 1, ND*norb
     call random_number(phase)
     theta(k) = 2.d0*pi*(phase-.5d0)
  end do
!
  end subroutine randomphase
