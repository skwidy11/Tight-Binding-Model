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
! Rmax    = integer number of random vectors used for truncated basis approximation
!

  include 'params.inc'


  integer :: imode, Rmax, ir, m, natoms, jj, kk, matz, nn, mu, nu, j
  integer :: ntypes, nt1
!
  double precision, allocatable :: hamilt(:,:,:,:)
  double precision, allocatable :: vel(:,:,:,:)
  double complex, allocatable :: psi0(:,:), psi(:,:), psi_new(:,:),psi_old(:,:)
!  double complex, allocatable :: psiv(:,:), psiv_old(:,:), psiv0(:,:)
!  double complex, allocatable :: psil(:,:,:),psir(:,:,:)
   double complex, allocatable :: mu_n(:)
!
  double precision, allocatable :: g(:),tn(:,:)
  double precision, allocatable :: theta(:)
  double precision, allocatable :: sx(:), sy(:), sz(:)
  integer, allocatable :: ntype(:)
  double complex, allocatable :: dos(:),doss(:)
  double precision, dimension(3,3) :: h, hinv
  double precision, dimension(3,9) :: eps
  double precision, dimension(3,3,9,9) :: t, delta
!
  integer, allocatable :: neighb(:,:)
!
  double precision :: Emax, Emin, E, lambda, convfac
  double precision :: eps0, epsmax, epsmin, rx, ry, rz, vol
  double precision :: pi, efac, diff, sum, time, time_interval, time_max
  double precision :: norm, norm_x, fac, tact, dE
  double precision rxij,ryij,rzij,sxij,syij,szij
  double complex ampl
                 

  integer ND, NT, L, KT, Ntau, nesteps
!
  integer iseed, ne
  integer n, k, r, i, it, norb, idr
!
  character(len=1) :: kernel, boundx, boundy, boundz
  character(len=4) :: ext
!

! Initialize MPI
!  call MPI_INIT(ierr)
!  call MPI_COMM_RANK(MPI_COMM_WORLD,mynod,ierr)
!  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)


  pi=4.0d0*datan(1.0d0)
! Read structure file, for mynod=0

        open(unit=12,file='structure')
        read(12,*) natoms,norb ! Number of atoms and number of orbitals per site

!  call MPI_BCAST(natoms,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(norb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate (sx(natoms),sy(natoms),sz(natoms))
  allocate (ntype(natoms))

        read(12,*) h(1,1),h(1,2),h(1,3)
        read(12,*) h(2,1),h(2,2),h(2,3)
        read(12,*) h(3,1),h(3,2),h(3,3)
        read(12,*) a,rc ! Lattice parameter in Angstroms
!  call MPI_BCAST(h,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(a,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(rc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call matinv(h,hinv,vol) ! Invert the h-matrix
  vol=vol*a**3 ! System volume, A**3

  do n=1,natoms ! Cycle over all atoms in structure file
    read(12,*) m,rx,ry,rz,nt1 ! Number, rx,ry,rz coordinates
!  call MPI_BCAST(rx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!   call MPI_BCAST(ry,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!   call MPI_BCAST(rz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!   call MPI_BCAST(nt1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   sx(n)=hinv(1,1)*rx+hinv(1,2)*ry+hinv(1,3)*rz
   sy(n)=hinv(2,1)*rx+hinv(2,2)*ry+hinv(2,3)*rz
   sz(n)=hinv(3,1)*rx+hinv(3,2)*ry+hinv(3,3)*rz
   ntype(n)=nt1
  enddo
!
! Total number of sites

  ND=natoms


!
! * - read input data
!
        
  call input_super (ND, norb, NT, Ntau, KT, time_max, iseed, &
                     ntypes, eps, t, delta, kernel,      &
                     lambda, Rmax, boundx, boundy, boundz, &
                     Emax, Emin, epsmax, epsmin, nesteps)
  write(6,*) 'done input_super NT=',NT
  write(6,*) 'epsmax,epsmin,nesteps=',epsmax,epsmin,nesteps

  allocate (dos(nesteps+1),doss(nesteps+1))
  allocate(tn(1:nesteps+1,0:NT-1))

  dos=0.0d0

  call initrand(iseed)

!
! * - maximum order of expansion
!
  print *,'NT = ',NT
!
! * - setting array dimensions
!
  allocate (hamilt(ND,NZ+1,norb,norb))
!  allocate (vel(ND,NZ+1,norb,norb))
  allocate (psi0(ND,norb),psi(ND,norb),psi_new(ND,norb),psi_old(ND,norb))
!  allocate (psiv(ND,norb),psiv_old(ND,norb),psiv0(ND,norb))
!  allocate (psil(nesteps+1,ND,norb),psir(nesteps+1,ND,norb))
  allocate (neighb(ND,0:NZ+1))
  allocate (g(0:NT-1))
  allocate (theta(ND*norb))
  allocate (mu_n(0:NT-1))
!
! * - generate the kernel coefficients
!

  if (kernel.eq.'j') then
     call jackson_kernel ( NT, g )
  else
     if (kernel.eq.'l') then
        call lorentz_kernel ( NT, lambda, g )
     else
        call dirichlet_kernel ( NT, g )
     end if
  end if
!
! * - save the Kernel coefficients
!
  open(unit=17,file='g.d',status='unknown')
  do n = 0, NT-1
     write(17,"(i4,1x,e10.3)")n,g(n)
  end do
  close(unit=17)


  call cheb(NT, nesteps, epsmin, epsmax, g, tn )


  dE=((Emax-Emin)/2.0d0)*(epsmax-epsmin)/dfloat(nesteps)

  write(6,*) 'calling hamiltonian'
  call ham_super ( ND, norb, hamilt, eps, t, delta, neighb, Emax, Emin, &
                    boundx, boundy, boundz, h, sx, sy, sz, ntype)
  write(6,*) 'done hamiltonian'

! Zero the mu array

  mu_n=0.0d0
  do ir=0,Rmax-1   ! Loop on random state vectors for truncated basis approximation
        write(6,*) 'ir=',ir

!     if(mod(ir,nproc).eq.mynod) then

     call randomphase ( ND, norb, theta ) ! Random phases

      imode=1 ! Random initial vector
!      imode=2 ! Start on a Cu next to a vacancy

        k=829
        write(6,*) 'Cu starting at site k, ntype=',k,ntype(k)
        k=900
        write(6,*) 'Other neighboring Cu site at k,ntype=',k,ntype(k)
     call setinitsta ( ND, norb, imode, theta, psi0 ) ! Initial state

! * - open output files

     open(unit=19,file='norm.d',status='unknown')

     if(imode.eq.1) call projection(ND, norb, psi0, psi0, ampl)
     if(imode.eq.2) ampl=psi0(k,2)
     mu_n(0)=mu_n(0)+ampl
!     write(6,*) 'starting amplitude=',ampl
     call chebvec1 (ND, norb, neighb, hamilt, psi0, psi)
     if(imode.eq.1) call projection(ND, norb, psi0, psi, ampl)
     if(imode.eq.2) ampl=psi(k,2)
!     write(6,*) 'ampl=',ampl
     mu_n(1)=mu_n(1)+ampl
     psi_old=psi0 

     do n=2,NT-1 
        call chebvec ( ND, norb, neighb, hamilt, psi_old, psi )
        psi_new=psi_old
        psi_old=psi
        psi=psi_new
        if(imode.eq.1) call projection(ND, norb, psi0, psi_new, ampl)
        if(imode.eq.2) ampl=psi_new(k,2)
!        write(6,*) 'n,ampl=',n,ampl
        mu_n(n)=mu_n(n)+ampl
     enddo

     do ne= 1, nesteps+1 ! Loop on energies for calcualtion of n(E)

        dos(ne)=dos(ne)+g(0)*mu_n(0)  ! First mu term is series is equal to 1
        if(ne.eq.1) write(6,*) mu_n(0)
        do n=1,NT-1 ! Summation on series for terms from n=1 to n=NT-1
          dos(ne)=dos(ne)+2.0d0*g(n)*mu_n(n)*tn(ne,n) 
        if(ne.eq.1) write(6,*) mu_n(n)
        enddo
     enddo   ! End of loop on energies for n(E) final calculation

  enddo ! End loop on random vectors

! MPI sum the dos array to get the density of states

!  call MPI_REDUCE(dos,doss,nesteps+1,MPI_DOUBLE_COMPLEX,   &
!              MPI_SUM,0,MPI_COMM_WORLD,ierr)

! Output averaged DOS
  fac=2.0d0/dfloat(nesteps) ! normalized energy interval for integration
! MPI sum the dos array to get the energy-dependent conductivity  
    open(unit=44,file='dos.out')
    do ne=1,nesteps+1
        eps0=epsmin+dfloat(ne-1)*(epsmax-epsmin)/dfloat(nesteps)
        efac=fac/(pi*dsqrt(1.0d0-eps0**2)) 
        E=eps0*(Emax-Emin)/2.0d0+(Emax+Emin)/2.0d0
        dos(ne)=dos(ne)*efac/dfloat(Rmax)
        write(44,100) ne,eps0,E,dreal(dos(ne)),aimag(dos(ne))
        if(E.le.0.0d0) sum=sum+aimag(dos(ne))*dE
     enddo
   close(44)
   write(6,*) 'normal termination'
   write(6,*) 'order parameter (summed to E=0)=',sum
!   call MPI_FINALIZE(ierr)
100 format(i9,2x,4f20.7)

!
  stop
  end program
!
!***********************************************************************
! Include subroutines
!
  include "matinv.f90"
  include "cheb.f90"
  include "chebvec.f90"
  include "chebvec1.f90"
  include "dirichlet_kernel.f90"
  include "initrand.f90"
  include "input_super.f90"
!  include "input_ch_diff.f90"
  include "jackson_kernel.f90"
  include "lorentz_kernel.f90"
  include "projection.f90"
  include "randomphase.f90"
  include "setinitsta.f90"
  include "ham_super.f90"
  include "updatevec.f90"
!  include "vmult.f90"
