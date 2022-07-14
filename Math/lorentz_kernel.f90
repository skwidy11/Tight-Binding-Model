! Subroutine to generate the Lorentz kernel coefficients.
!
! Input: NT, lambda
!
! Output: g(n)
!
! NT      = actual number of kernel polynomials
! g       = array containing the kernel coefficients
! lambda  = parameter
!
! E. Mucciolo, UCF (Apr 01, 2014)
!
  subroutine lorentz_kernel ( NT, lambda, g )
!
  implicit none
!
  integer, parameter :: double = KIND(0.d0)
!
  integer, intent(in) :: NT
  real (kind=double), intent(out), dimension(0:NT-1) :: g
  real (kind=double), intent(in) :: lambda
!
  integer :: n
!
! Lorentz kernel coefficients
!
  do n = 0, NT-1
     g(n) = dsinh(lambda*(1-dfloat(n)/NT))/dsinh(lambda)
  end do
!
  end subroutine lorentz_kernel
