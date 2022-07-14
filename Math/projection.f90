! Subroutine to evaluate the projection of the current state vector
! onto the initial state vector.
!
! E. Mucciolo (Apr 01, 2014)
!
  subroutine projection ( ND, norb, psi0, psi, ampl )
!
  implicit none
!
  integer, parameter :: double = KIND(0.d0)
  integer, intent(in) :: ND, norb
  complex (kind=double), intent(in), dimension(ND,norb) :: psi0, psi
  complex (kind=double), intent(out) :: ampl
!
  integer :: k, mu
!
  ampl = 0.d0
  do k = 1, ND
  do mu=1,norb
     ampl = ampl + dconjg(psi0(k,mu))*psi(k,mu)
  enddo
  enddo
!
 end subroutine projection
