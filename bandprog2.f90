program bands_prog
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
   integer :: ien,iky,num_elec
   integer :: ikx,iy,ie, il,icount, pos, m
   real,allocatable :: eig(:),rwork(:),eig_f(:)
   complex,allocatable :: work(:)
   double precision, dimension(:,:), allocatable :: weight
   integer :: cells,n,i,j,lwork,info,ND

   pi=4.0d0*datan(1.0d0)
   hx = 24
   hy = 24
   cells = 30
   num_elec = 18

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
   allocate(eig(20),work(lwork),rwork(28))
   allocate(eig_f(10))
   n=0


  ! Cheev get's eigenvalues of complex hermitian matrix.
  ! Gamma to X
    do j = 0, cells/2-1
      kx = pi*dble(j)/dble(cells)
      ky = 0d0
      call update_band_AF(kx,ky,band_mat,t,eps,nocc,site,U)
      n=n+1
      call cheev('N','U',20,band_mat,20,eig,work,lwork,rwork,info)
      write(11,*) n,kx,ky,eig(1)-fermi,eig(2)-fermi,eig(3)-fermi,eig(4)-fermi,eig(5)-fermi,&
            eig(6)-fermi,eig(7)-fermi,eig(8)-fermi,eig(9)-fermi,eig(10)-fermi,eig(11)-fermi,&
            eig(12)-fermi,eig(13)-fermi,eig(14)-fermi,eig(15)-fermi,eig(16)-fermi,eig(17)-fermi,&
            eig(18)-fermi,eig(19)-fermi,eig(20)-fermi
    end do

    ! X to M
      do j = 0, cells/2-1
        ! AF Path
        !kx = -pi*(dble(cells)-dble(j))/(dble(cells)*2d0)
        !ky = pi*dble(j)/(dble(cells)*2d0)
        ! Non-AF Path
        kx = pi/2
        ky = pi*dble(j)/dble(cells)
        n=n+1
        call update_band_AF(kx,ky,band_mat,t,eps,nocc,site,U)
        call cheev('N','U',20,band_mat,20,eig,work,lwork,rwork,info)
        write(11,*) n,kx,ky,eig(1)-fermi,eig(2)-fermi,eig(3)-fermi,eig(4)-fermi,eig(5)-fermi,&
              eig(6)-fermi,eig(7)-fermi,eig(8)-fermi,eig(9)-fermi,eig(10)-fermi,eig(11)-fermi,&
              eig(12)-fermi,eig(13)-fermi,eig(14)-fermi,eig(15)-fermi,eig(16)-fermi,eig(17)-fermi,&
              eig(18)-fermi,eig(19)-fermi,eig(20)-fermi
      end do

  ! M to Gamma
  do j = 0, cells/2
    ! AF Path
    !kx = pi*(dble(cells)/2-dble(j))/(dble(cells)*2d0)
    !ky = pi*(dble(cells)/2-dble(j))/(dble(cells)*2d0)
    ! Non AF Path
    kx = pi/2 - pi*dble(j)/dble(cells)
    ky = pi/2 - pi*dble(j)/dble(cells)
    n=n+1
    call update_band_AF(kx,ky,band_mat,t,eps,nocc,site,U)
    call cheev('N','U',20,band_mat,20,eig,work,lwork,rwork,info)
    write(11,*) n,kx,ky,eig(1)-fermi,eig(2)-fermi,eig(3)-fermi,eig(4)-fermi,eig(5)-fermi,&
          eig(6)-fermi,eig(7)-fermi,eig(8)-fermi,eig(9)-fermi,eig(10)-fermi,eig(11)-fermi,&
          eig(12)-fermi,eig(13)-fermi,eig(14)-fermi,eig(15)-fermi,eig(16)-fermi,eig(17)-fermi,&
          eig(18)-fermi,eig(19)-fermi,eig(20)-fermi
  end do

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

  include "update_band_af.f90"

end program
