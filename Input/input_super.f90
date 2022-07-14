! Subroutine to read the input data for the evaluation of time
! evolution ing a tight-binding model using the (Chebyshev) kernel polynomial
! approximation.
!
! E. Mucciolo, UCF (August 05, 2011)
!
! modified by P. Schelling, UCF (June 8, 2018)

!
  subroutine input_super  (ND, norb, NT, Ntau, KT, time_max, iseed, &
                            ntypes, eps, t, delta,      &
                            kernel, lambda, Rmax, boundx, boundy, boundz,  &
                            Emax, Emin, epsmax, epsmin, nesteps,U)
  include 'params.inc'

!
  integer ND, NT, NR, L, KT, Ntau
  integer iseed
  double precision lambda, W, time_max,U
  double precision rx,ry,rz,Emax,Emin,epsmax,epsmin
  double precision, dimension(4,9) :: eps
  double precision, dimension(4,4,9,9) :: t, delta
  character(len=1) :: kernel, boundx, boundy, boundz
  character(len=15) :: label1,label2,label3
  integer :: m,n,natoms,i,j,norb,Rmax,nesteps, ntypes
  integer :: mu, nu

  t = 0.0
  open(unit=11,file='Input/input_super.dat') ! control parameters
  open(unit=14,file='Input/potent.txt') ! tight-binding parameters

  write(6,*) 'reading tight-binding potential parameters'
! Read tight-binding parameters
!
  read(14,*)
  write(6,*) 'through first line'
  read(14,*) ntypes,norb
  write(6,*) 'ntypes,norb=',ntypes,norb
  read(14,*)
  do n=1,ntypes
    read(14,*) (eps(n,mu),mu=1,norb)
  enddo
  read(14,*) U
  write(6,*) 'diagonal energies from input file'
  write(6,*) eps(1,1),eps(1,2)!,eps(1,3)
  write(6,*) eps(2,1),eps(2,2)!,eps(2,3)
  write(6,*) eps(3,1),eps(3,2)!,eps(3,3)
  write(6,*) eps(4,1),eps(4,2)!,eps(4,3)
  do n=1,ntypes
  do m=n,ntypes
    read(14,*)
    read(14,*)
    do mu=1,norb
      read(14,*) (t(n,m,mu,nu),nu=1,norb)
    enddo
    write(6,*) t(n,m,1,1),t(n,m,1,2)!,t(n,m,1,3)
    write(6,*) t(n,m,2,1),t(n,m,2,2)!,t(n,m,2,3)
    !write(6,*) t(n,m,3,1),t(n,m,3,2)!,t(n,m,3,3)
    read(14,*)
    do mu=1,norb
      read(14,*) (delta(n,m,mu,nu),nu=1,norb)
    enddo
    read(14,*)
  enddo
  enddo


  do n=1,ntypes
  do m=1,n-1
     do mu=1,norb
     do nu=1,norb
          t(n,m,mu,nu)=t(m,n,nu,mu)
          delta(n,m,mu,nu)=delta(m,n,nu,mu)
     enddo
     enddo
  enddo
  enddo

  !print*, "Oxy Oxy Hop: ", eps(1,1)

!  call MPI_BCAST(eps,ntypes*norb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(t,ntypes*ntypes*norb*norb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(delta,ntypes*ntypes*norb*norb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


        read(11,*) KT ! Number of time steps
        read(11,*) time_max ! max time

        read(11,*)NT ! order in the expansion

        read(11,*) iseed ! seed for the random potential sequence

        read(11,*)kernel ! choice of kernel (j=Jackson, l=Lorentz, d=Dirichlet)
        write(6,*) kernel

        read(11,*)lambda ! Lorentz kernel parameter

        if (lambda.eq.0.d0) return

        read(11,*)boundx, boundy, boundz ! choice of boundary conditions (p=periodic, h=hard wall)

        read(11,*) Rmax
        write(6,*) "Rmax=", Rmax

        read(11,*) Emax   ! Highest energy for normalization of Hamiltonian
        read(11,*) Emin   ! Highest energy for normalization of Hamiltonian

        read(11,*) nesteps ! Number of energy steps

        read(11,*) epsmin  ! Lowest normalized energy for DOS calculation

        read(11,*) epsmax  ! Highest energy for DOS calculation
        close(unit=11)
! call MPI_BCAST(NT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(Rmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(kernel,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(lambda,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  write(6,*) 'someways through bcast'
!  call MPI_BCAST(boundx,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(boundy,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(boundz,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(Emax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(Emin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  write(6,*) 'mostly through bcast'
!  call MPI_BCAST(nesteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(epsmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(epsmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!   write(6,*) 'mynod,boundx,boundy,boundz=',mynod,boundx,boundy,boundz
!  write(6,*) 'mynod,kernel=',mynod,kernel

  return
  end subroutine
