! DOS calculation using KPM
! E. Mucciolo, UCF (August 05, 2011)
!
! modified from above version by P. Schelling, June 8, 2018
  program chbsh_dos
!
! ND      = Hilbert space size
! NT      = order of expansion for a given time increment
! NZ      = maximum coordination number in the Hilbert space
! Emax  = upper bound for the spectral rescaling
! Emin  = upper bound for the spectral rescaling
! t       = hopping matrix element amplitude
! iseed   = seed for random potential sequence
! neighb  = array containing the number of nearest neighbors of a
!           given site and their numbers
! hamilt  = array containing Hamiltonian matrix element between a given
!           site and its nearest neighbors
! g       = array containing the kernel coefficients
! lambda  = parameter for the Lorentz kernel
! kernel  = character defining the kernel choice
! boundx,y,z = choice of boundary conditions along different directions
! convfac = energy conversion factor
! psi0    = initial state
! norm    = norm of the current state vector
! theta   = array with random phases
! Rmax    = integer number of random vectors used for truncated basis approxima
! Anchovies
  include 'params.inc'
  integer :: imode, Rmax, ir, m, natoms, jj, kk, matz, nn, mu, nu, j
  integer :: ntypes,nt1
!
  double complex, allocatable :: vel(:,:,:,:), ampl(:), l_ampl(:,:)
!  double precision, allocatable :: mu_n(:,:)
  double complex, allocatable :: psi0(:,:), psi(:,:), psi_new(:,:), psi_old(:,:)
!  double complex, allocatable :: psi0_x(:,:), psi_x(:,:), psi_x_new(:,:), psi_x_old(:,:)
  double complex, allocatable :: gamma(:,:,:)
!  double complex, allocatable :: gamma_x(:,:,:)
!
  double precision, allocatable :: g(:),tn(:),Jbessel(:),hamilt(:,:,:,:,:)
  double precision, allocatable :: hmatup(:,:), hmatdown(:,:)
  double precision, allocatable :: theta(:), eigup(:), eigdown(:)
  double precision, allocatable :: sx(:), sy(:), sz(:)
  double precision, allocatable :: nocc(:,:), nocc_old(:,:)
  double complex, allocatable :: dos(:)
  double precision, dimension(3,3) :: h, hinv
   double precision, dimension(1000) :: doso
!
  integer, allocatable :: neighb(:,:)
  integer, allocatable :: ntype(:)
!
  double precision :: Emax, Emin, E, lambda, convfac, epsmax, epsmin, rx, ry, rz, vol, pi, efac, diff, sum, time, time_interval
  double precision :: norm, norm_x, fac, tact, time_max, w, p, Etot
  double precision :: delta_sc, energy
  double precision :: rxij,ryij,rzij,sxij,syij,szij, nocc_tot_up, nocc_tot_down
  double precision :: nocc_mean_up,nocc_mean_down,nocc_dev_up,nocc_dev_down
  double precision, dimension(2) :: n_cu, n_ox
  double precision :: n_tot, spin
  double complex :: mag, amplt
  double precision :: magr,ampr,ampi, U, fermi, dc_cor, hub_cor, threshold
  double precision, dimension(4,9) :: eps
  double precision, dimension(4,4,9,9) :: t, delta

  complex, dimension(20,20) :: band_mat
  real, allocatable :: eig(:), rwork(:)
  complex, allocatable :: work(:)
  double precision :: kx,ky
  integer :: cells, lwork, site, runs
  integer :: ien,ikx,iky, CupVac,OxyVac

  double precision :: dos_me1,dos_me2,dos_their1,dos_their2,rsquare, dummy

  integer ND, NT, L, KT, Ntau, nesteps, numstates
!
  integer iseed, ne, oxy, cu, info,a_cu_site
  integer n, k, r, i, it, norb, dum1, dum2, num_elec
!
  character(len=1) :: kernel, boundx, boundy, boundz
  character(len=4) :: ext
!
  pi=4.0d0*datan(1.0d0)

! Read structure file
  open(unit=12,file='structure')
  read(12,*) natoms,norb, p ! Number of atoms and number of orbitals per site
  allocate (sx(natoms),sy(natoms),sz(natoms))
  allocate (ntype(natoms))
  read(12,*) h(1,1),h(1,2),h(1,3)
  read(12,*) h(2,1),h(2,2),h(2,3)
  read(12,*) h(3,1),h(3,2),h(3,3)
  read(12,*) a,rc ! Lattice parameter in Angstroms
  call matinv(h,hinv,vol) ! Invert the h-matrix
  vol=vol*a**3 ! System volume, A**3
  print *, "p =", p
  print *, "vol", vol
  numstates = natoms + natoms*2/3





  do n=1,natoms ! Cycle over all atoms in structure file
   read(12,*) m,rx,ry,rz,nt1 ! Number, rx,ry,rz coordinates
   sx(n)=hinv(1,1)*rx+hinv(1,2)*ry+hinv(1,3)*rz
   sy(n)=hinv(2,1)*rx+hinv(2,2)*ry+hinv(2,3)*rz
   sz(n)=hinv(3,1)*rx+hinv(3,2)*ry+hinv(3,3)*rz
   ntype(n)=nt1
!   write(6,*) sx(n),sy(n),sz(n),ntype(n)
  enddo
  close(12)
!
! Total number of sites

  ND=natoms
! Set oxygen and copper atom to look at\


 call input_super(ND, norb, NT, Ntau, KT, time_max, iseed, &
                  ntypes, eps, t, delta,      &
                  kernel, lambda, Rmax, boundx, boundy, boundz,  &
                  Emax, Emin, epsmax, epsmin, nesteps,U)
  Ntau=NT ! Keep the same for now

  allocate (dos(nesteps+1))
  dos=0.0d0
  call initrand(iseed)

!
! * - maximum order of expansion
!
  print *,'NT = ',NT
!
! * - setting array dimensions
!
  allocate (hamilt(ND,NZ+1,norb,norb,2))
  allocate (vel(ND,NZ+1,norb,norb))
  allocate (psi0(ND,norb),psi(ND,norb),psi_new(ND,norb),psi_old(ND,norb))
!  allocate (psi0_x(ND,norb),psi_x(ND,norb),psi_x_new(ND,norb),psi_x_old(ND,norb))
  allocate (gamma(ND,norb,0:NT-1))
!  allocate (gamma_x(ND,norb,0:NT-1))
  allocate (neighb(ND,0:NZ+1))
  allocate (g(0:NT-1),Jbessel(0:NT-1))
  allocate (theta(ND*norb))
  allocate (tn(0:NT-1))
  allocate(ampl(0:KT-1))
  allocate(l_ampl(0:KT-1,2)) ! n = 1 Cu, n = 2 Ox

  allocate(hmatup(numstates,numstates),hmatdown(numstates,numstates))
  allocate(eigup(numstates),eigdown(numstates))
  allocate(nocc(numstates,2))

!  allocate  (mu_n(0:KT-1,0:NT-1))
  nocc = 0.0d0
  open(unit=43,file='Occup.dat')
  open(unit=47,file='AntiF.dat')
  a_cu_site = ND*2/3 + 50

! Find a Copper Vac site, oxygen vac site
  !do j=1,ND
  !  if (ntype(j)==3) then
  !    CupVac = j
  !    exit
  !  endif
  !enddo
  !do j=1,ND
  !  if (ntype(j)==4) then
  !    OxyVac = j
  !    exit
  !  endif
  !enddo

! ****************
! Diagonalization
! ****************
  !U = 3.0d0
  print *, "U = ", U
  threshold = 5*10d0**(-2)
  num_elec = 9*natoms/3
!  do while (fermi.le.3.26d0)
!  fermi = fermi + 0.25d0
  nocc_dev_up = 1
  do j = 1, numstates ! Initialize anti-ferro state
    !if (ntype(mod(j-1,numstates) + 1).eq.1.or.ntype(mod(j-1,numstates) + 1).eq.3) then
    if (ntype(j).eq.1.or.ntype(j).eq.3) then
      n = j-natoms*2/3
      k = (n-1) / h(1,1)
       if (mod(n+k,2).eq.0) then! May need to update if lattice size changes.
!        Antiferromagnetic
        !nocc(j,1) = 1.0d0
        !nocc(j,2) = 0.0d0
        nocc(j,1) = 0.95d0
        nocc(j,2) = 0.40d0
!       Paramagnetic
        !nocc(j,1) = 0.25d0
        !nocc(j,2) = 0.25d0
      else
!       Antiferromagnetic
        !nocc(j,1) = 0.0d0
        !nocc(j,2) = 1.0d0
        nocc(j,1) = 0.40d0
        nocc(j,2) = 0.95d0
!       Paramagnetic
        !nocc(j,1) = 0.25d0
        !nocc(j,2) = 0.25d0
      end if
    else ! oxygen
      nocc(j,1) = 0.91d0
      nocc(j,2) = 0.91d0
    end if
  end do
  nocc_old = nocc
  n=1

  do while (nocc_dev_up.ge.threshold) ! Go to self-consistency
    call ham_super ( ND, norb, hamilt, eps, t, delta, neighb, Emax, Emin, &
                    boundx, boundy, boundz, h, sx, sy, sz, ntype, nocc_old, U,numstates)

    call diagonalize(hamilt, neighb, ND, norb, ntype, fermi, eigup, eigdown, &
                      nocc,hmatup,hmatdown,num_elec,numstates)

    print *, "Cu Old: ", nocc_old(a_cu_site,1),nocc_old(a_cu_site,2), "New: ",nocc(a_cu_site,1),nocc(a_cu_site,2)

    print *, "Cu Old: ", nocc_old(289,1),nocc_old(289,2), "New: ",nocc(289,1),nocc(289,2)
    print *, "Cu Old: ", nocc_old(290,1),nocc_old(301,2), "New: ",nocc(290,1),nocc(301,2)
    !print *, "CuVac Old: ", nocc_old(CupVac,1),nocc_old(CupVac,2), "New: ",nocc(CupVac,1),nocc(CupVac,2)
    print *, "Oxy Old: ", nocc_old(112,1),nocc_old(112,2), "New: ",nocc(112,1),nocc(112,2)
    !print *, "OxyVac Old: ", nocc_old(OxyVac,1),nocc_old(OxyVac,2), "New: ",nocc(OxyVac,1),nocc(OxyVac,2)

!   Straight term by term deviation
    nocc_dev_up = 0.0d0
    nocc_dev_down = 0.0d0
    do j=1,ND
      if(ntype(j).eq.1.or.ntype(j).eq.3) nocc_dev_up = nocc_dev_up + (nocc(j,1)-nocc_old(j,1))**2
      if(ntype(j).eq.1.or.ntype(j).eq.3) nocc_dev_down = nocc_dev_down + (nocc(j,2)-nocc_old(j,2))**2
    end do
    print *, "Up Deviation: ", nocc_dev_up, "Down Deviation: ", nocc_dev_down,"Run:", n

    n_cu = 0.0d0
    n_tot = 0.0d0
    n_ox = 0.0d0
    do j=0,numstates-1
      if(ntype(mod(j,ND)+1).eq.1.or.ntype(mod(j,ND)+1).eq.3) then
         n_cu(1) = n_cu(1) + nocc(j+1,1)
         n_cu(2) = n_cu(2) + nocc(j+1,2)
      else
         n_ox(1) = n_ox(1) + nocc(j+1,1)
         n_ox(2) = n_ox(2) + nocc(j+1,2)
      end if
      n_tot = n_tot + nocc(j+1,1) + nocc(j+1,2)
    end do

    nocc_old = nocc
    n = n + 1
  end do ! End of Self-consistency

  open(unit=51,file='eigvec')
  open(unit=52,file='eigval')
  do i=1,numstates
    do j=1,numstates
      write(51,*) hmatup(i,j)
    end do
    write(52,*) eigup(i)
  end do
  fermi = eigup(num_elec/2)

  call dos2(eigup,eigdown,hmatup,ntype,num_elec,numstates,ND)
  ! Calculate square deviation
  open(unit=81,file='dos.dat')
  open(unit=82,file='dossy.dat')
  rsquare=0.0d0
  do i=1,160
    read(81,*) energy,dummy,dos_me1,dos_me2
    read(82,*) energy,dos_their1,dos_their2
    rsquare = rsquare + (dos_me1-2*dos_their1)**2 + (dos_me2-2*dos_their2)**2
  end do
  print *, "Rsquare = ", rsquare

  do i=1 , ND
    spin = 0.0d0
    spin = nocc(i,1) - nocc(i,2)
    if(ntype(i).eq.2) spin = spin + nocc(i+ND,1)-nocc(i+ND,2)
    write(47,*) sx(i)*h(1,1), sy(i)*h(2,2),ntype(i),spin
    sx(i) = h(1,1)*sx(i) ! Hopefully doesn't break anything lol.
    sy(i) = h(2,2)*sy(i)
  end do
  close(43)
  close(47)
!  end do

! *****************************************************************
! Band Structure
! *****************************************************************

  open(unit=11,file='bands.dat')
  lwork = 660
  cells = 30
  site = a_cu_site ! Chooses a copper site
  allocate(eig(20),work(lwork),rwork(3*20-2))
  n=0
! Cheev get's eigenvalues of complex hermitian matrix
! Update_band_AF updates the band matrix (using antiferromagnetic zone)
! Gamma to X
  do j = 0, cells/2-1
    kx = pi*dble(j)/dble(cells)
    ky = 0d0
    call update_band_AF(kx,ky,band_mat,t,eps,nocc,site,U,numstates)
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
    call update_band_AF(kx,ky,band_mat,t,eps,nocc,site,U,numstates)
    call cheev('N','U',20,band_mat,20,eig,work,lwork,rwork,info)
    write(11,*) n,kx,ky,eig(1)-fermi,eig(2)-fermi,eig(3)-fermi,eig(4)-fermi,eig(5)-fermi,&
          eig(6)-fermi,eig(7)-fermi,eig(8)-fermi,eig(9)-fermi,eig(10)-fermi,eig(11)-fermi,&
          eig(12)-fermi,eig(13)-fermi,eig(14)-fermi,eig(15)-fermi,eig(16)-fermi,eig(17)-fermi,&
          eig(18)-fermi,eig(19)-fermi,eig(20)-fermi
  end do

! M to Gamma
  do j = 0, cells
    ! AF Path
    !kx = pi*(dble(cells)/2-dble(j))/(dble(cells)*2d0)
    !ky = pi*(dble(cells)/2-dble(j))/(dble(cells)*2d0)
    ! Non AF Path
    kx = pi/2 - pi*dble(j)/dble(2*cells)
    ky = pi/2 - pi*dble(j)/dble(2*cells)
    n=n+1
    call update_band_AF(kx,ky,band_mat,t,eps,nocc,site,U,numstates)
    call cheev('N','U',20,band_mat,20,eig,work,lwork,rwork,info)
    write(11,*) n,kx,ky,eig(1)-fermi,eig(2)-fermi,eig(3)-fermi,eig(4)-fermi,eig(5)-fermi,&
          eig(6)-fermi,eig(7)-fermi,eig(8)-fermi,eig(9)-fermi,eig(10)-fermi,eig(11)-fermi,&
          eig(12)-fermi,eig(13)-fermi,eig(14)-fermi,eig(15)-fermi,eig(16)-fermi,eig(17)-fermi,&
          eig(18)-fermi,eig(19)-fermi,eig(20)-fermi
  end do
  ! *************************************************************************

  call spectral_weight(fermi,eigup,eigdown,hmatup,hmatdown,sx,sy,numstates,ND,ntype,num_elec,NINT(h(1,1)))


  ! %%%%%%%%%%%%% Band Structure's DOS %%%%%%%%%%%%%%%
  runs = 200
  delta_sc = 10.d0/runs
  emin=-5.0d0-delta_sc
  emax=-5.0d0
  doso=0.0d0
  n=0
  open(unit=34,file='energykvec.dat')
  do ikx=-cells/2,cells/2-1
    do iky=-cells/2,cells/2-1
      emin=-5.0d0-delta_sc
      emax=-5.0d0

      kx = pi*dble(ikx)/dble(cells) !this is k / 2 actually.
      ky = pi*dble(iky)/dble(cells)
      call update_band_AF(kx,ky,band_mat,t,eps,nocc,site,U,numstates)
      call cheev('N','U',20,band_mat,20,eig,work,lwork,rwork,info)
      write(34,*) n,kx,ky,eig(1)-fermi,eig(2)-fermi,eig(3)-fermi,eig(4)-fermi,eig(5)-fermi,&
            eig(6)-fermi,eig(7)-fermi,eig(8)-fermi,eig(9)-fermi,eig(10)-fermi,eig(11)-fermi,&
            eig(12)-fermi,eig(13)-fermi,eig(14)-fermi,eig(15)-fermi,eig(16)-fermi,eig(17)-fermi,&
            eig(18)-fermi,eig(19)-fermi,eig(20)-fermi
      n=n+1

      do ien=1,runs
        emin=emin+delta_sc
        emax=emax+delta_sc
        do i=1,20
          if(eig(i).gt.emin.and.eig(i).le.emax) doso(ien)=doso(ien)+1.
        enddo
      enddo ! ien
  enddo !iky
  enddo !ikx

  ! Dos from band struc
  emin=-5.0d0-delta_sc
  open(unit=3,file='banddos.dat')
  do ien=2,runs
    energy=emin+(dfloat(ien)+0.5)*delta_sc
    write(3,*) energy-fermi,doso(ien)
  enddo
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






! *******************************************************
! Chebyshev From on out
! *******************************************************



goto 45

!
! * - generate the kernel coefficients
!
ampl = 0
l_ampl = 0
if (kernel.eq.'j') then
   call jackson_kernel ( NT, g )
else
   if (kernel.eq.'l') then
      call lorentz_kernel ( NT, lambda, g )
   else
      call dirichlet_kernel ( NT, g )
   end if
end if
open(unit=80, file='evolvenew.dat')
open(unit=70, file='local_evolvenew.dat')

! Compute time
  time_interval = time_max/dfloat(KT)
  write(6,*) Emax,Emin
  write(6,*) 'Time step =',time_interval

  call besselcoeff(NT, time_interval, Jbessel)
call scale_ham ( ND, norb, hamilt, eps, t, delta, neighb, Emax, Emin, &
                boundx, boundy, boundz, h, sx, sy, sz, ntype, nocc_old, U)
!  mu_n=0.0d0
  do ir=0,Rmax-1   ! Loop on random state vectors for truncated basis approximation
     write(6,*) 'New random state, ir, Rmax-1=',ir,Rmax-1

     call randomphase ( ND, norb, theta ) ! Random phases

     imode=1 ! Random initial vector

     call setinitsta ( ND, norb, imode, theta, psi0, ntype ) ! Initial state

     psi_old=psi0

! Loop over time

     do it = 0, KT-1
       time = dfloat(it)*time_interval

! Evolve vectors forward in time by one time interval
! Expansion coefficients
       if(it.gt.0) then
        call gammacoeff_st(ND, NT, norb, neighb, hamilt, psi_old, gamma)
        ! Gamma factor can hit infinity, careful.
! Update state vectors

        call updatevec(ND, NT, NT, norb, gamma, g, Jbessel, psi_new, norm,it)

       else
        psi_new=psi_old
       endif

       call projection(ND, norb, psi0, psi_new, amplt)
       call local_projection(ND, KT, it, norb, psi0, psi_new, ntype, 1, l_ampl)
       call local_projection(ND, KT, it, norb, psi0, psi_new, ntype, 2, l_ampl)

       psi_old=psi_new
       ampl(it) = ampl(it) + amplt
     enddo ! End loop on time steps
  enddo ! End loop on random vectors

write(80,*) KT, time_max, Emax, Emin, vol, p
do it =0,KT-1
  time = dble(it) * time_interval
  ampl(it) = ampl(it)/dble(Rmax)
  l_ampl(it,1) = l_ampl(it,1)/dble(Rmax)
  l_ampl(it,2) = l_ampl(it,2)/dble(Rmax)
  ampr=real(ampl(it)); ampi=imag(ampl(it))
  magr=dsqrt(ampr**2+ampi**2)
  write(80,122) time, real(ampl(it)), imag(ampl(it)), magr
  write(70,*) time, real(l_ampl(it,1)), imag(l_ampl(it,1)), real(l_ampl(it,2)), imag(l_ampl(it,2))
enddo

45  print *, "The End."
100 format(i9,2x,3f20.7)
122     format(4f20.8)

!
  end program
!
!***********************************************************************
! Include subroutines
!
  include "Math/matinv.f90"
  include "dos2.f90"
  include "Math/cheb.f90"
  include "update_band.f90"
  include "update_band_af.f90"
  include "diagonalize.f90"
  include "Math/chebvec.f90"
  include "Math/chebvec1.f90"
  include "Math/dirichlet_kernel.f90"
  include "Input/initrand.f90"
  include "Input/input_super.f90"
  include "Math/jackson_kernel.f90"
  include "Math/lorentz_kernel.f90"
!  include "numchar4.f90"
  include "Math/projection.f90"
  include "Math/local_projection.f90"
  include "Math/randomphase.f90"
  include "Input/setinitsta.f90"
  include "Math/swap.f90"
  include "Input/ham_super.f90"
  include "Math/besselcoeff.f90"
  include "spectralweight5b.f90"
  include "Math/bessj0.f90"
  include "Math/bessj1.f90"
  include "Math/bessjn.f90"
  include "Math/gammacoeff_st.f90"
  include "updatevec.f90"
  include "scale_ham.f90"
