! Subroutine to generate a 1st order Chebyshev recurrent vector
!
! Input:  |r>_0
!
! Output: |r>_1 = H |r>_0
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
  subroutine chebvec1 ( ND, norb, neighb, hamilt, rvec0, rvec1 )
!
!
  include 'params.inc'
  integer, parameter :: double = KIND(0.d0)
  integer, intent(in) :: ND, norb
  real (kind=double), intent(in), dimension(ND,NZ+1,norb,norb) :: hamilt
  complex (kind=double), intent(in), dimension(ND,norb) :: rvec0
  complex (kind=double), intent(out), dimension(ND,norb) :: rvec1
  integer, intent(in), dimension(ND,0:NZ+1) :: neighb
!
  complex (kind=double) :: rvec1aux
  integer :: k, nn, i, kk, mu, nu
!
  do k = 1, ND
! determine how many neighbors the k-site has (including itself)
  nn = neighb(k,0)
  do mu=1,norb
     rvec1aux = (0.d0,0.d0)
!
! take the product between the k component and the components
! corresponding to its nearest-neighbors
!
     do i = 1, nn
     kk = neighb(k,i)
     do nu=1, norb
        rvec1aux = rvec1aux + hamilt(k,i,mu,nu)*rvec0(kk,nu)
     enddo
     enddo
     rvec1(k,mu) = rvec1aux
  enddo
  enddo
!
  end subroutine chebvec1
