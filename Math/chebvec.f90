! Subroutine to generate a Chebyshev recurrent vector (order n > 1).
!
! Input:  |r>_0, |r>_1 (for |r>_{n-2} and |r>_{n-1})
!
! Output: |r>_0 = 2H |r>_1 - |r>_0 (for |r>_n)
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
! E. Mucciolo, UCF (April 01, 2014)
!
  subroutine chebvec ( ND, norb, neighb, hamilt, rvec0, rvec1 )
!
!
  include 'params.inc'
  integer, parameter :: double = KIND(0.d0)
  integer, intent(in) :: ND, norb
  real (kind=double), intent(in), dimension(ND,NZ+1,norb,norb) :: hamilt
  complex (kind=double), intent(out), dimension(ND,norb) :: rvec0
  complex (kind=double), intent(in), dimension(ND,norb) :: rvec1
  integer, intent(in) :: neighb(ND,0:NZ+1)
  complex (kind=double) :: rvec0aux
  integer  k, i, nn, kk
  integer  mu, nu
!
  do k = 1, ND
! determine how many neighbors the k-site has (including itself)
  nn = neighb(k,0)
  do mu= 1, norb
     rvec0aux = (0.d0,0.d0)
!
!
! take the product between the k component and the components
! corresponding to its nearest-neighbors
!
     do i = 1, nn
        kk = neighb(k,i)
        do nu= 1, norb
           rvec0aux = rvec0aux + 2.d0*hamilt(k,i,mu,nu)*rvec1(kk,nu)
        enddo
     enddo
     rvec0(k,mu) = rvec0aux - rvec0(k,mu)
  enddo
  enddo
!
  end subroutine chebvec
