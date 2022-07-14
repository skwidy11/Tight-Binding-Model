! Subroutine to multiply a vector by the velocity operator, and output the new
! resulting vector
!
! Input:  |psi>
!
! Output: |psiv> = v |psi>
!
! ND     = actual number of sites
! NZ     = maximum coordination number per site
! neighb = array containing the number of nearest neighbors of a 
!          given site and their numbers
! hamilt = Hamiltonian matrix element between a given site and its
!          nearest neighbors
! rvec0  = array with the zeroth order vector components
! rvec1  = array with the first order vector components
!
! E. Mucciolo, UCF (Apr 01, 2014)
!
  subroutine vmult( ND, norb, neighb, vel, psi, psiv )
!
!
  include 'params.inc'
  integer, parameter :: double = KIND(0.d0)
  integer, intent(in) :: ND, norb
  real (kind=double), intent(in), dimension(ND,NZ+1,norb,norb) :: vel
  complex (kind=double), intent(in), dimension(ND,norb) :: psi
  complex (kind=double), intent(out), dimension(ND,norb) :: psiv
  integer, intent(in), dimension(ND,0:NZ+1) :: neighb
!
  integer :: k, nn, i, kk, mu, nu
!

  psiv=0.0d0
  do k = 1, ND
! determine how many neighbors the k-site has (including itself)
  nn = neighb(k,0)
  do mu=1,norb
!
!
     do i = 1, nn
     kk = neighb(k,i)
     do nu=1, norb
        psiv(k,mu) = psiv(k,mu) + vel(k,i,mu,nu)*psi(kk,nu)
     enddo
     enddo
  enddo
  enddo
!
  end subroutine vmult
