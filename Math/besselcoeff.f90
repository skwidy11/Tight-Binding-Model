! Surbroutine to compute Bessel functions of order n of a fixed
! argument. It uses functions from Numerical Recipes.
!
! E. Mucciolo, UCF (Apr 01, 2014)
!
  subroutine besselcoeff ( NT, time, Jbessel )
!
  implicit none
!
  integer, parameter :: double = KIND(0.d0)
!
  integer, intent(in) :: NT
  real (kind=double), intent(out), dimension(0:NT-1) :: Jbessel
  real (kind=double), intent(in) :: time
!
  integer :: n
!
  real (kind=double) bessj0, bessj1, bessjn
  external bessj0, bessj1, bessjn
!
!!  Jbessel(0) = BESSEL_J0(time) 
!!  Jbessel(1) = BESSEL_J1(time)
  Jbessel(0) = bessj0(time) 
  Jbessel(1) = bessj1(time) 
 
!
  do n = 2, NT-1
!!     Jbessel(n) = BESSEL_JN(n,time) 
     Jbessel(n) = bessjn(n,time)
  end do
!
  end subroutine besselcoeff
