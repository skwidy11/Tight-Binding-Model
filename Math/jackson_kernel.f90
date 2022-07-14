! Subroutine to generate the Jackson kernel coefficients.
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
  subroutine jackson_kernel ( NT, g )
!
  implicit none
!
  integer, parameter :: double = KIND(0.d0)
  real (kind=double), parameter :: pi = 3.14159265358979323846264338327950288d0
!
  integer, intent(in) :: NT
  real (kind=double), intent(out), dimension(0:NT-1) :: g
!
  integer :: n
!
! Jackson kernel coefficients
!
  do n = 0, NT-1
     g(n) = ((NT-n+1)*dcos(pi*n/(NT+1)) + dsin(pi*n/(NT+1))/dtan(pi &
            /(NT+1)))/(NT+1)
  end do
!
  end subroutine jackson_kernel
