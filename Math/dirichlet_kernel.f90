! Subroutine to generate the Dirichlet kernel coefficients.
!
! Input: NT
!
! Output: g(n)
!
! NT      = actual number of kernel polynomials
! g       = array containing the kernel coefficients
!
! E. Mucciolo, UCF (Apr 01, 2014)
!
  subroutine dirichlet_kernel ( NT, g )
!
  implicit none
!
  integer, parameter :: double = KIND(0.d0)
  integer, intent(in) :: NT
  real (kind=double), intent(out), dimension(0:NT-1) :: g
!
! Dirichlet kernel coefficients
!
     g = 1.d0
!
  end subroutine dirichlet_kernel
