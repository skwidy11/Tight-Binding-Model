! Surbroutine to compute the updated state vector in an iterative time
! evolution using the Chebyshev expansion of the evolution operator
! with a time-dependent Hamiltonian.
!
! E. Mucciolo, UCF (Apr 01, 2014)
!
  subroutine updatevec ( ND, NT, Ntau, norb, gamma, g, Jbessel, psi, norm, anint)
!
  include 'params.inc'
!
  integer, parameter :: double = KIND(0.d0)
  integer, intent(in) :: ND, NT, norb, Ntau, anint
  complex (kind=double), intent(in), dimension(ND,norb,0:NT-1) :: gamma
  complex (kind=double), intent(out), dimension(ND,norb) :: psi
  real (kind=double), intent(in), dimension(0:NT-1) :: g, Jbessel
  real (kind=double), intent(out) :: norm
!
  complex (kind=double) :: factor
  integer :: n, k, mu
!
  psi = (0.d0,0.d0)
!
! evaluating the state vector
!
  norm = 0.d0

  do k = 1, ND
  do mu = 1,norb
     factor = (1.d0,0.d0)
     do n = 1, NT-1
        factor = factor*(0.d0,-1.d0)
        psi(k,mu) = psi(k,mu) + g(n)*Jbessel(n)*gamma(k,mu,n)*factor !
     end do
     psi(k,mu) = 2.d0*psi(k,mu) + g(0)*Jbessel(0)*gamma(k,mu,0)
     norm = norm + abs(dconjg(psi(k,mu))*psi(k,mu))
  end do
  end do

!  if(norm - ND > 0.01 .or. norm - ND < -0.01) then
    write(6,*) 'update, norm=',norm,"Delta = ", (1728d0 + 1152d0-norm)/dble(anint), anint
!  end if
!
  end subroutine updatevec
