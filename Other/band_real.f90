program bands
  implicit none
 double precision, dimension(5,5) :: matrix
 double precision, dimension(1000) :: doso
 double precision :: kx,ky,pi,eps,eps_p, tpd, tpp, tpp_p
 double precision :: emax,emin, delta, energy,U,nocc
 integer :: ien,ikx,iky
 double precision,allocatable :: eig(:)
 double precision,allocatable :: work(:)
 integer :: cells,n,i,j,l,info

 pi=4.0d0*datan(1.0d0)
 cells = 200

 eps   = -1.5d0
 eps_p = -1.8d0
 tpd   =  1.0d0
 tpp   =  0.3d0
 tpp_p =  0.0d0
 nocc = 0.25d0
 U = 3d0

 open(unit=11,file='bands.dat')
 l = 170
 allocate(eig(5),work(l))
 n=0

! Gamma to X
  ky = 0d0
  do j = 0, cells/2
    kx =pi*dble(j)/dble(cells)
    call update(kx,ky,matrix)
    n=n+1
    call dsyev('N','U',5,matrix,5,eig,work,l,info)
    write(11,*) n,kx,ky,eig(1),eig(2),eig(3),eig(4),eig(5)
  end do
! X to M
  kx = pi/2d0
  do j = 0, cells/2
    ky = -pi*dble(j)/dble(cells)
    n=n+1
    call update(kx,ky,matrix)
    call dsyev('N','U',5,matrix,5,eig,work,l,info)
    write(11,*) n,kx,ky,eig(1),eig(2),eig(3),eig(4),eig(5)
  end do

! M to Gamma
  do j = 0, cells/2
    kx = pi*(dble(cells)/2d0-dble(j))/dble(cells)
    ky = -pi*(dble(cells)/2d0-dble(j))/dble(cells)
    n=n+1
    call update(kx,ky,matrix)
    call dsyev('N','U',5,matrix,5,eig,work,l,info)
    write(11,*) n,kx,ky,eig(1),eig(2),eig(3),eig(4),eig(5)
  end do

  delta=12.d0/dfloat(cells)
  emin=-6.0d0-delta
  emax=-6.0d0
  doso=0.0d0

  do ikx=-cells/2,cells/2-1
    do iky=-cells/2,cells/2-1
      emin=-6.0d0-delta
      emax=-6.0d0

      kx = pi*dble(ikx)/dble(cells)
      ky = pi*dble(iky)/dble(cells)
      call update(kx,ky,matrix)
      call dsyev('N','U',5,matrix,5,eig,work,l,info)

      do ien=1,cells
        emin=emin+delta
        emax=emax+delta
        if(eig(1).gt.emin.and.eig(1).le.emax)then
          doso(ien)=doso(ien)+1.
        endif
        if(eig(2).gt.emin.and.eig(2).le.emax)then
          doso(ien)=doso(ien)+1.
        endif
        if(eig(3).gt.emin.and.eig(3).le.emax)then
          doso(ien)=doso(ien)+1.
        endif
        if(eig(4).gt.emin.and.eig(4).le.emax)then
          doso(ien)=doso(ien)+1.
        endif
        if(eig(5).gt.emin.and.eig(5).le.emax)then
          doso(ien)=doso(ien)+1.
        endif
    enddo ! ien
  enddo !iky
  enddo !ikx

  emin=-6.0d0-delta
  open(unit=3,file='banddos.dat')
  do ien=2,cells
    energy=emin+(dfloat(ien)+0.5)*delta
    write(3,*) energy,doso(ien)
  enddo

contains

subroutine update(kx,ky,matrix)
  double precision, dimension(5,5), intent(out) :: matrix
  double precision,intent(in) :: kx,ky
  matrix = 0d0
  matrix(1,1) = U*(0.5d0-nocc)
  matrix(2,2) = eps
  matrix(3,3) = eps_p
  matrix(4,4) = eps_p
  matrix(5,5) = eps
  matrix(1,2) = 2d0*tpd*dsin(kx)
  matrix(2,1) = matrix(1,2)
  matrix(1,5) = -2d0*tpd*dsin(ky)
  matrix(5,1) = matrix(1,5)
  matrix(2,4) = -4d0*tpp*dcos(kx)*dcos(ky)
  matrix(4,2) = matrix(2,4)
  matrix(3,5) = -4d0*tpp*dcos(kx)*dcos(ky)
  matrix(5,3) = matrix(3,5)
  matrix(2,5) = -4d0*tpp_p*dsin(kx)*dsin(ky)
  matrix(5,2) = matrix(2,5)
  matrix(3,4) = -4d0*tpp_p*dsin(kx)*dsin(ky)
  matrix(4,3) = matrix(3,4)
end subroutine



end program bands
