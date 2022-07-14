program bands
implicit none
 complex, dimension(5,5) :: matrix
 double precision :: kx,ky,pi,eps,eps_p, tpd, tpp, tpp_p
 real,allocatable :: eig(:),rwork(:)
 complex,allocatable :: work(:)
 integer :: cells,n,i,j,lwork,info

 pi=4.0d0*datan(1.0d0)
 cells = 200

 eps   = -1.5d0
 eps_p = -1.8d0
 tpd   =  1.0d0
 tpp   =  0.3d0
 tpp_p =  0.0d0

 open(unit=11,file='bands.dat')
 lwork = 165
 allocate(eig(5),work(lwork),rwork(12))
 n=0
! Cheev get's eigenvalues of complex hermitian matrix.
! Gamma to X
  ky = 0d0
  do j = 0, cells/2
    kx = pi*dble(j)/dble(cells)
    call update(kx,ky,matrix)
    n=n+1
    call cheev('N','U',5,matrix,5,eig,work,lwork,rwork,info)
    write(11,*) n,kx,ky,eig(1),eig(2),eig(3),eig(4),eig(5)
  end do

! X to M
  kx = pi/2d0
  do j = 0, cells/2
    ky = pi*dble(j)/dble(cells)
    n=n+1
    call update(kx,ky,matrix)
    call cheev('N','U',5,matrix,5,eig,work,lwork,rwork,info)
    write(11,*) n,kx,ky,eig(1),eig(2),eig(3),eig(4),eig(5)
  end do

! M to Gamma
  do j = 0, cells/2
    kx = pi*(dble(cells)/2d0-dble(j))/dble(cells)
    ky = pi*(dble(cells)/2d0-dble(j))/dble(cells)
    n=n+1
    call update(kx,ky,matrix)
    call cheev('N','U',5,matrix,5,eig,work,lwork,rwork,info)
    write(11,*) n,kx,ky,eig(1),eig(2),eig(3),eig(4),eig(5)
  end do

contains

subroutine update(kx,ky,matrix)
  complex, dimension(5,5), intent(out) :: matrix
  double precision,intent(in) :: kx,ky
  matrix = 0d0
  matrix(2,2) = eps
  matrix(3,3) = eps_p
  matrix(4,4) = eps_p
  matrix(5,5) = eps
  matrix(1,2) = (0,-1d0)*2d0*tpd*dsin(kx)*exp((0,1d0)*kx)
  matrix(2,1) = conjg(matrix(1,2))
  matrix(1,5) = (0,1d0)*2d0*tpd*dsin(ky)*exp((0,1d0)*ky)
  matrix(5,1) = conjg(matrix(1,5))
  matrix(2,4) = -4d0*tpp*dcos(kx)*dcos(ky)*exp((0,1d0)*(ky-kx))
  matrix(4,2) = conjg(matrix(2,4))
  matrix(3,5) = -4d0*tpp*dcos(kx)*dcos(ky)*exp((0,1d0)*(ky-kx))
  matrix(5,3) = conjg(matrix(3,5))
  matrix(2,5) = -4d0*tpp_p*dsin(kx)*dsin(ky)*exp((0,1d0)*(ky-kx))
  matrix(5,2) = conjg(matrix(2,5))
  matrix(3,4) = -4d0*tpp_p*dsin(kx)*dsin(ky)*exp((0,1d0)*(ky-kx))
  matrix(4,3) = conjg(matrix(3,4))
end subroutine


end program bands
