! Subroutine to set the initial state of a single-particle system
!
! E. Mucciolo, UCF (Apr 01, 2014)
! Modified P. Schelling, UCF (June 8, 2018)
!
! imode=0 initial state on a single site
! imode=1 random initial vector, using theta(k) vector
!
  subroutine setinitsta ( ND, norb, imode, theta, psi0, ntype)
!
  implicit none
!
  integer, parameter :: double = KIND(0.d0)
!
  integer, intent(in) :: ND
  integer :: imode,k,n,l,norb
  integer, dimension(ND), intent(in) :: ntype
  complex (kind=double), intent(out), dimension(ND,norb) :: psi0
  double precision, dimension(ND*norb) :: theta
  double precision norm
!

   psi0 = (0.d0,0.d0)
   if(imode.eq.1) then
   do k=1,ND
   do n=1,norb-1
    l=n+(k-1)*norb
    psi0(k,n) = (1.0d0,0.0d0)*dcos(theta(l))+(0.0d0,1.0d0)*dsin(theta(l))
    if(n.eq.2 .and. ntype(k).eq.1) psi0(k,n) = 0
    
!    if(n.eq.1 .and. ntype(k).eq.2) psi0(k,n) = 0 for potent4
!    if(n.eq.3 .and. ntype(k).eq.1) psi0(k,n) = 0
   enddo
   enddo
   endif

   if(imode.eq.2) then
        k=829
        write(6,*) 'Starting off k=',k
        psi0(k,1)=(0.0d0,2.0d0)
   endif

! Compute and output norm of wave function

  norm=0.0d0
  do k=1,ND
  do n=1,norb
     norm=norm+abs(psi0(k,n))**2
  enddo
  enddo
  write(6,*) 'norm=',norm

  return
  end
