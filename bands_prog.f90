! The paramagnetic unit cell (not AF)
program bands
  implicit none
 complex, dimension(10,10) :: matrix
 double precision, dimension(1000) :: doso
 double precision :: kx,ky,pi,eps,eps_p, tpd, tpp, tpp_p, U,n_up, n_down
 double precision :: costrans,sintrans,wave,xpos,ypos,rx,ry,rz
 double precision :: ef,dum1,dum2,dum3, phasearg, overlap,rc,a
 double precision :: emax,emin, delta, energy, delta_sp,fkx,fky
 double precision ::sx(10),sy(10), vol, h(3,3)
 integer :: ntype(2880),nt1
 integer :: ien,hx,hy,iky
 integer :: ikx,iy,ie, il,icount, pos, m
 real,allocatable :: eig(:),rwork(:),eig_f(:)
 complex,allocatable :: work(:)
 double precision, dimension(:,:), allocatable :: weight
 integer :: cells,n,i,j,lwork,info,ND

 pi=4.0d0*datan(1.0d0)
 hx = 24
 hy = 24
 cells = 24
 allocate(weight(-2*hx:2*hx,-2*hy:2*hy))

 eps   = -1.5d0
 eps_p = -1.8d0
 tpd   =  1.0d0
 tpp   =  0.3d0
 tpp_p =  0.0d0
 U     =  3d0
 ef = 0.634694869312
! For fermi = 0.7 with above values.
 n_up  =  0.960d0
 n_down=  0.349d0
 !n_up = 0.25d0
 !n_down = 0.25d0
 open(unit=11,file='bands.dat')
 lwork = 330
 allocate(eig(10),work(lwork),rwork(28))
 allocate(eig_f(10))
 n=0

 sx(1)=0.0
 sy(1)=0.0
 sx(2)=0.5
 sy(2)=0.0
 sx(3)=sx(2)
 sy(3)=sy(2)
 sx(4)=0.0
 sy(4)=0.5
 sx(5)=sx(4)
 sy(5)=sy(4)
 do i=1,5
   sx(i+5)=sx(i)+1
   sy(i+5)=sy(i)
 end do


! Cheev get's eigenvalues of complex hermitian matrix.
! Gamma to X
  do j = 0, cells/2-1
    ky = 0d0
    kx = pi*dble(j)/dble(cells)
    call update(kx,ky,matrix)
    n=n+1
    call cheev('N','U',10,matrix,10,eig,work,lwork,rwork,info)
!   Folded part
    ky = pi*(cells/2d0 - dble(j))/dble(cells)
    kx = pi/2d0
    call update(kx,ky,matrix)
    call cheev('N','U',10,matrix,10,eig_f,work,lwork,rwork,info)
    write(11,*) n,kx,ky,eig(1),eig(2),eig(3),eig(4),eig(5),eig(6),eig(7),eig(8),eig(9),eig(10)!,&
      !eig_f(1),eig_f(2),eig_f(3),eig_f(4),eig_f(5),eig_f(6),eig_f(7),eig_f(8),eig_f(9),eig_f(10)
  end do

! X to M

  do j = 0, cells/2-1
    n=n+1
!    kx = pi/2d0
!    ky = pi*dble(j)/dble(cells)
!    call update(kx,ky,matrix)
!    call cheev('N','U',10,matrix,10,eig,work,lwork,rwork,info)
!   For folding
    ky = pi*dble(j)/(dble(cells)*2d0)
    kx = pi*(dble(cells) - dble(j))/(dble(cells)*2d0)
    call update(kx,ky,matrix)
    call cheev('N','U',10,matrix,10,eig_f,work,lwork,rwork,info)
    call update(kx,ky,matrix)
    call cheev('N','U',10,matrix,10,eig,work,lwork,rwork,info)
    write(11,*) n,kx,ky,eig(1),eig(2),eig(3),eig(4),eig(5),eig(6),eig(7),eig(8),eig(9),eig(10)!,&
      !eig_f(1),eig_f(2),eig_f(3),eig_f(4),eig_f(5),eig_f(6),eig_f(7),eig_f(8),eig_f(9),eig_f(10)
  end do

! M to Gamma
  do j = 0, cells/2
    kx = pi*(dble(cells)/2d0-dble(j))/(dble(cells)*2d0)
    ky = pi*(dble(cells)/2d0-dble(j))/(dble(cells)*2d0)
    n=n+1
    call update(kx,ky,matrix)
    call cheev('N','U',10,matrix,10,eig,work,lwork,rwork,info)
!   Folded Part
    kx = pi*(dble(cells)/2d0+dble(j))/(dble(cells)*2d0)
    ky = pi*(dble(cells)/2d0+dble(j))/(dble(cells)*2d0)
    call update(kx,ky,matrix)
    call cheev('N','U',10,matrix,10,eig_f,work,lwork,rwork,info)
    write(11,*) n,kx,ky,eig(1),eig(2),eig(3),eig(4),eig(5),eig(6),eig(7),eig(8),eig(9),eig(10)!,&
    !  eig_f(1),eig_f(2),eig_f(3),eig_f(4),eig_f(5),eig_f(6),eig_f(7),eig_f(8),eig_f(9),eig_f(10)
  end do

  delta=12.d0/dfloat(cells)
  delta_sp=4.0d0/48d0
  emin=-6.0d0-delta
  emax=-6.0d0
  doso=0.0d0
  open(unit=7,file='weightatfs_b')
  do ikx=-cells/2,cells/2-1
    do iky=-cells/2,cells/2-1
      emin=-6.0d0-delta
      emax=-6.0d0
      kx = pi*dble(ikx)/dble(cells)
      ky = pi*dble(iky)/dble(cells)
      call update(kx,ky,matrix)
      call cheev('V','U',10,matrix,10,eig,work,lwork,rwork,info)

      do ien=1,cells
        emin=emin+delta
        emax=emax+delta
        do i=1,10
          if(eig(i).gt.emin.and.eig(i).le.emax) doso(ien)=doso(ien)+1.
        enddo
      enddo ! ien

      emax=ef+delta_sp
      emin=ef-delta_sp
      do i=1,10
        if (eig(i).gt.emin.and.eig(i).le.emax) then
        !  print *, i,ikx,iky
          costrans=0.0d0
          sintrans=0.0d0
          do il=1,10
            phasearg=(sx(il)*kx+sy(il)*ky)
            costrans=costrans+real(matrix(il,i))*dcos(phasearg)
            sintrans=sintrans+real(matrix(il,i))*dsin(phasearg)
          end do
          overlap=(costrans**2+sintrans**2)
          weight(ikx,iky)=weight(ikx,iky)+overlap
        end if
      end do
      write(7,*)  kx,ky,weight(ikx,iky)

  enddo !iky
  enddo !ikx


  !emin=-6.0d0-delta
  !open(unit=3,file='banddos.dat')
  !do ien=2,cells
  !  energy=emin+(dfloat(ien)+0.5)*delta
  !  write(3,*) energy,doso(ien)
  !enddo

contains

subroutine update(kx,ky,matrix)
  complex, dimension(10,10), intent(out) :: matrix
  double precision,intent(in) :: kx,ky
  matrix = 0d0
  do i=0,1
    if(i.eq.0) then
      matrix(1,1) = U*(0.5d0-n_up)
    else
      matrix(6,6) = U*(0.5d0-n_down)
    end if
    matrix(2+5*i,2+5*i) = eps
    matrix(3+5*i,3+5*i) = eps_p
    matrix(4+5*i,4+5*i) = eps_p
    matrix(5+5*i,5+5*i) = eps
    matrix(1+5*i,2+5*i) = (0,-1d0)*2d0*tpd*dsin(kx)
    matrix(2+5*i,1+5*i) = conjg(matrix(1+5*i,2+5*i))
    matrix(1+5*i,5+5*i) = (0,1d0)*2d0*tpd*dsin(ky)
    matrix(5+5*i,1+5*i) = conjg(matrix(1+5*i,5+5*i))
    matrix(2+5*i,4+5*i) = -4d0*tpp*dcos(kx)*dcos(ky)
    matrix(4+5*i,2+5*i) = matrix(2+5*i,4+5*i)
    matrix(3+5*i,5+5*i) = -4d0*tpp*dcos(kx)*dcos(ky)
    matrix(5+5*i,3+5*i) = matrix(3+5*i,5+5*i)
    matrix(2+5*i,5+5*i) = -4d0*tpp_p*dsin(kx)*dsin(ky)
    matrix(5+5*i,2+5*i) = matrix(2+5*i,5+5*i)
    matrix(3+5*i,4+5*i) = -4d0*tpp_p*dsin(kx)*dsin(ky)
    matrix(4+5*i,3+5*i) = matrix(3+5*i,4+5*i)
  end do
end subroutine

include "update_band.f90"


end program bands
