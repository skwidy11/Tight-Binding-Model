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
!
  double complex, allocatable :: hamilt(:,:,:,:),mu_n(:,:),vel(:,:,:,:)
  double complex, allocatable :: psi0(:,:), psi(:,:), psi_new(:,:), psi_old(:,:)
  double complex, allocatable :: psi0_x(:,:), psi_x(:,:), psi_x_new(:,:), psi_x_old(:,:)
  double complex, allocatable :: gamma(:,:,:),gamma_x(:,:,:)
!
  double precision, allocatable :: g(:),tn(:),Jbessel(:)
  double precision, allocatable :: theta(:)
  double precision, allocatable :: sx(:), sy(:), sz(:)
  double complex, allocatable :: dos(:)
  double precision, dimension(3,3) :: h, hinv
!
  integer, allocatable :: neighb(:,:)
!
  double precision t, Emax, Emin, E, lambda, convfac, eps, epsmax, epsmin, rx, ry, rz, vol, pi, efac, diff, sum, time, time_interval, time_max, norm, norm_x, fac, tact
  double precision rxij,ryij,rzij,sxij,syij,szij
  double complex ampl
                 

  integer ND, NT, L, KT, Ntau, nesteps, ierr
!
  integer iseed, ne
  integer n, k, r, i, it, norb
!
  character(len=1) :: kernel, boundx, boundy, boundz
  character(len=4) :: ext
!

  pi=4.0d0*datan(1.0d0)
! Read structure file
  open(unit=12,file='structure')
  read(12,*) natoms,norb ! Number of atoms and number of orbitals per site
  allocate (sx(natoms),sy(natoms),sz(natoms))
  read(12,*) h(1,1),h(1,2),h(1,3)
  read(12,*) h(2,1),h(2,2),h(2,3)
  read(12,*) h(3,1),h(3,2),h(3,3)
  read(12,*) a,rc ! Lattice parameter in Angstroms
  call matinv(h,hinv,vol) ! Invert the h-matrix
  vol=vol*a**3 ! System volume, A**3


  do n=1,natoms ! Cycle over all atoms in structure file
   read(12,*) m,rx,ry,rz ! Number, rx,ry,rz coordinates
   sx(n)=hinv(1,1)*rx+hinv(1,2)*ry+hinv(1,3)*rz
   sy(n)=hinv(2,1)*rx+hinv(2,2)*ry+hinv(2,3)*rz
   sz(n)=hinv(3,1)*rx+hinv(3,2)*ry+hinv(3,3)*rz
  enddo
  close(12)
!
! Total number of sites

  ND=natoms


!
! * - read input data
!
  call input_ch_diff (ND, norb, NT, Ntau, KT, time_max, iseed, &
                     kernel, lambda, Rmax, boundx, boundy, boundz, &
                     Emax, Emin, epsmax, epsmin, nesteps)

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
  allocate (hamilt(ND,NZ+1,norb,norb),vel(ND,NZ+1,norb,norb))
  allocate (psi0(ND,norb),psi(ND,norb),psi_new(ND,norb),psi_old(ND,norb))
  allocate (psi0_x(ND,norb),psi_x(ND,norb),psi_x_new(ND,norb),psi_x_old(ND,norb))
  allocate (gamma(ND,norb,0:NT-1),gamma_x(ND,norb,0:NT-1))
  allocate (neighb(ND,0:NZ+1))
  allocate (g(0:NT-1),Jbessel(0:NT-1))
  allocate (theta(ND*norb))
  allocate (tn(0:NT-1))
  allocate  (mu_n(0:KT-1,0:NT-1))
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

  time_interval = time_max/dfloat(KT)

! Compute actual time
  
  tact=time_interval*4.135667d0/((Emax-Emin)*pi) ! Time in fs

  write(6,*) 'Time step (fs)=',tact

! Factor to obtain results in Ohm**-1 cm**1

  fac=(2.0d0*pi*1.602d0)**2/(vol*4.135667**2)

! Take into account time in fs for integration, factor of 2 for
! positive/negative time integration

  fac=2.0d0*fac*tact

! Account for normalization of Hamiltonian for velocity operators
! and the fact that the delta function defined in terms of the normalized
! energy and Hamiltonian
! Final result to be integrated units 10**-15 Ohm**-1 cm**-1

  fac=fac*((Emax-Emin)/2.0d0)

  call besselcoeff(NT, time_interval, Jbessel)

  call ham_3D ( ND, norb, hamilt, vel, neighb, Emax, Emin, &
                    boundx, boundy, boundz, h, sx, sy, sz)

! Zero the mu_n array

  mu_n=0.0d0

  do ir=0,Rmax-1   ! Loop on random state vectors for truncated basis approximation
     write(6,*) 'New random state, ir, Rmax-1=',ir,Rmax-1

     call randomphase ( ND, norb, theta ) ! Random phases

     imode=1 ! Random initial vector

     call setinitsta ( ND, norb, imode, theta, psi0 ) ! Initial state

! Determine psi0_x vector for initial time
     do k=1,ND
       nn=neighb(k,0)
       do n=2,nn
         j=neighb(k,n)
         sxij=sx(i)-sx(j)
         syij=sy(i)-sy(j)
         szij=sz(i)-sz(j)

         if(boundx.eq.'p') sxij=sxij-dble(dnint(sxij)) ! Periodic b.c. along x direction
         if(boundy.eq.'p') syij=syij-dble(dnint(syij)) ! Periodic b.c. along y direction
         if(boundz.eq.'p') szij=szij-dble(dnint(szij)) ! Periodic b.c. along z direction

         rxij=h(1,1)*sxij+h(1,2)*syij+h(1,3)*szij
         ryij=h(2,1)*sxij+h(2,2)*syij+h(2,3)*szij
         rzij=h(3,1)*sxij+h(3,2)*syij+h(3,3)*szij
         rxij=rxij*a; ryij=ryij*a; rzij=rzij*a

         do mu=1,norb
         do nu=1,norb
           psi0_x(k,mu)=psi0_x(k,mu)+vel(k,n,mu,nu)*psi0(j,nu)
         enddo
         enddo
       enddo
     enddo


     psi_old=psi0
     psi_x_old=psi0_x

! Loop over time
     
     do it = 0, KT-1
       time = dfloat(it)*time_max/dfloat(KT)
       write(6,*) 'time=',time

! Evolve vectors forward in time by one time interval
! Expansion coefficients
       if(it.gt.0) then
        write(6,*) 'calling gamma'
        call gammacoeff_st(ND, NT, norb, neighb, hamilt, psi_old, gamma)
        write(6,*) 'calling gamma'
        call gammacoeff_st(ND, NT, norb, neighb, hamilt, psi_x_old, gamma_x)

! Update state vectors
        write(6,*) 'calling update'
        call updatevec(ND, NT, NT, norb, gamma, g, Jbessel, psi_new, norm)
        write(6,*) 'calling update'
        call updatevec(ND, NT, NT, norb, gamma_x, g, Jbessel, psi_x_new, norm_x)
        write(6,*) 'norm,norm_x',norm,norm_x
        write(22,*) norm,norm_x
       else
        psi_new=psi_old
        psi_x_new=psi_x_old
       endif

       psi0=0.0d0

       do k=1,ND
         nn=neighb(k,0)
         do n=2,nn
           j=neighb(k,n)
           sxij=sx(i)-sx(j)
           syij=sy(i)-sy(j)
           szij=sz(i)-sz(j)

           if(boundx.eq.'p') sxij=sxij-dble(dnint(sxij)) ! Periodic b.c. along x direction
           if(boundy.eq.'p') syij=syij-dble(dnint(syij)) ! Periodic b.c. along y direction
           if(boundz.eq.'p') szij=szij-dble(dnint(szij)) ! Periodic b.c. along z direction

           rxij=h(1,1)*sxij+h(1,2)*syij+h(1,3)*szij
           ryij=h(2,1)*sxij+h(2,2)*syij+h(2,3)*szij
           rzij=h(3,1)*sxij+h(3,2)*syij+h(3,3)*szij
           rxij=rxij*a; ryij=ryij*a; rzij=rzij*a

           do mu=1,norb
           do nu=1,norb
             psi0(k,mu)=psi0(k,mu)+vel(k,n,mu,nu)*psi_new(j,nu)
           enddo
           enddo
         enddo
       enddo


        call projection(ND, norb, psi_x_new, psi0, ampl)
        mu_n(it,0)=mu_n(it,0)+ampl
	call chebvec1 (ND, norb, neighb, hamilt, psi0, psi)
        call projection(ND, norb, psi_x_new, psi, ampl)
	mu_n(it,1)=mu_n(it,1)+ampl
        psi_old=psi0 

       do n=2,NT-1 
          call chebvec ( ND, norb, neighb, hamilt, psi_old, psi )
          call projection(ND, norb, psi_x_new, psi_old, ampl)
          mu_n(it,n)=mu_n(it,n)+ampl
          call swap(ND,norb,psi_old,psi)
       enddo

! Update vectors for time evolution
       psi_old=psi_new
       psi_x_old=psi_x_new

     enddo ! End loop on time steps
  enddo ! End loop on random vectors


  do it=0,KT-1 ! Loop on time steps
     time = dfloat(it)*time_max/dfloat(KT)
     dos=0.0d0
     do ne= 1, nesteps+1 ! Loop on energies for calcualtion of n(E)
        eps=epsmin+dfloat(ne-1)*(epsmax-epsmin)/dfloat(nesteps)

        call cheb(NT, eps, tn )

        dos(ne)=dos(ne)+g(0)*mu_n(it,0)  ! First mu term is series is equal to 1
        do n=1,NT-1 ! Summation on series for terms from n=1 to n=NT-1
          dos(ne)=dos(ne)+2.0d0*g(n)*mu_n(it,n)*tn(n) ! Need the term related to Cheb expansion in terms of H
        enddo
     enddo   ! End of loop on energies for n(E) final calculation
  
! Output averaged DOS
    do ne=1,nesteps+1
        eps=epsmin+dfloat(ne-1)*(epsmax-epsmin)/dfloat(nesteps)
        efac=fac/(pi*dsqrt(1.0d0-eps**2)) 
        E=eps*(Emax-Emin)/2.0d0+(Emax+Emin)/2.0d0
        dos(ne)=dos(ne)*efac/dfloat(Rmax)
        write(6,100) ne,time,E,dreal(dos(ne))
   enddo
 enddo

100 format(i9,2x,3f20.7)

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
  include "input_ch_diff.f90"
  include "jackson_kernel.f90"
  include "lorentz_kernel.f90"
!  include "numchar4.f90"
  include "projection.f90"
  include "randomphase.f90"
  include "setinitsta.f90"
  include "swap.f90"
  include "ham_3D.f90"
  include "besselcoeff.f90"
  include "bessj0.f90"
  include "bessj1.f90"
  include "bessjn.f90"
  include "gammacoeff_st.f90"
  include "updatevec.f90"
