! Subroutine to evaluate the coefficients of the Chebyshev expansion
! of time-dependent vectors that go into the calculation of a
! single-particle localization problem. It implements the calculation
! recursively.
!
! E. Mucciolo, UCF (Apr 01, 2014)
!
  subroutine gammacoeff_st ( ND, NT, norb, neighb, hamilt, psi0, gamma)
!
  include 'params.inc'
!
  integer, parameter :: double = KIND(0.d0)
!
  integer, intent(in) :: ND, NT, norb
  real (kind=double), intent(in), dimension(ND,NZ+1,norb,norb) :: hamilt
  complex (kind=double), intent(out), dimension(ND,norb,0:NT-1) :: gamma
  complex (kind=double), intent(in), dimension(ND,norb) :: psi0
  integer, intent(in), dimension(ND,0:NZ+1) :: neighb
!
  complex (kind=double), allocatable, dimension(:,:) :: psi1, psin
!
  integer :: n
!
! set array dimensions
!
  allocate (psi1(ND,norb),psin(ND,norb))
!
! recursion for |psi>
!
  gamma(:,:,0) = psi0(1:ND,1:norb)
!

  call chebvec1 ( ND, norb, neighb, hamilt, psi0, psi1 )
     gamma(:,:,1) = psi1(1:ND,1:norb)
!
  psin = psi0
!
  do n = 2, NT-1
     call chebvec ( ND, norb, neighb, hamilt, psin, psi1 )
        gamma(:,:,n) = psin(1:ND,1:norb)
     call swap ( ND, norb, psi1, psin )
  end do
!
! unset array dimensions
!
  deallocate (psi1,psin)
!
  end subroutine gammacoeff_st
