! Subroutine to construct the elements of a nearest-neighbor
! tight-binding Hamiltonian. To be used in
! conjunction with the Chebyshev recursive vector multiplications.
!
! ND      = number of sites
! norb    = 9 (number of orbitals/site) PKS
! NZ      =  (maximum number of nearest neighbors) PKS
! neighb  = array containing the number of nearest neighbors of a
!           given site and their numbers
! hamilt  = array containing Hamiltonian matrix element between a
!           given site and its nearest neighbors
! t       = hopping amplitude
! convfac = energy conversion factor
!
! E. Mucciolo (July 26, 2011)
!
!
! neighb(k,0)  = NZ+1 = number of neighboring sites of site k, including
!                      itself
! neighb(k,1)  = k (the k site itself)
! neighb(k,i)  = site number of the (i-1)-th neighbor, i=2,....,NZ+1
!
! hamilt(k,i,mu,nu)  = Hamiltonian matrix element between site k and its
!                (i-1)-th neighbor (i=1 means the k site itself)
!
  subroutine scale_ham ( ND, norb, hamilt, eps, t, delta, neighb, Emax, Emin, &
                    boundx, boundy, boundz, h, sx, sy, sz, ntype, nocc, U)
!
  include 'params.inc'
!
  integer ND, norb, mu, nu
  double precision, dimension(ND,NZ+1,norb,norb,2) :: hamilt
! double precision, dimension(ND,NZ+1,norb,norb) ::  vel
  double precision, dimension(3,3) :: h, hinv
!  double precision es,ep,ed,e,f,g,hop,xl,yl,zl,fc
!  double precision :: tm, dm
  double precision :: rxij,ryij,rzij,rij,sxij,syij,szij
  double precision, dimension(4,9) :: eps
  double precision, dimension(4,4,9,9) :: t, delta
  double precision, dimension(ND) :: sx, sy, sz
  double precision, dimension(2880,2) :: nocc
  integer, dimension(ND) :: ntype
  character (len=1) :: boundx, boundy, boundz
  integer, dimension(ND,0:NZ+1) :: neighb
  logical isVertical
  double precision Emax, Emin, ECM, DE, convfac, U, rx, ry,rxp,ryp
  integer ::  k, nn, i, j, nt1, nt2

!
! resetting Hamiltonian array
!
  hamilt = 0.0d0
!  vel = 0.0d0
!
! array with information about neighbors
!
  neighb=0



! !!!!!!!!!!!!!!!!!!!
! rescaling to bound the spectrum to interval -1 < E < 1
! !!!!!!!!!!!!!!!!!!!
  ECM = (Emax + Emin)/2.0d0
  DE = Emax - Emin
  convfac = DE/2.0d0
!  write(6,*) 'ECM,DE,convfac=',ECM,DE,convfac
!

! Hamiltonian before normalizing
!  do i=1,ND
!  nn=neighb(i,0)
!   write(6,*) 'site,neighbors=',i,nn
!   nt1 = ntype(i)
!  do j=1,nn
!    nt2 = ntype(neighb(i,j))
!        do mu=1,norb-1
!        do nu=1,norb-1
!            write(6,*) j,mu,nu, hamilt(i,j,mu,nu)
!        enddo
!       enddo
!  enddo
!  enddo
write(6,*) 'scaling H'
  do i = 1, ND
     nn = neighb(i,0)
!        write(6,*) i,nn
     do mu=1,norb
        hamilt(i,1,mu,mu,1) = (hamilt(i,1,mu,mu,1) - ECM)/convfac
        hamilt(i,1,mu,mu,2) = (hamilt(i,1,mu,mu,2) - ECM)/convfac
     enddo
     do j = 2, nn
        do mu=1,norb
        do nu=1,norb
            hamilt(i,j,mu,nu,1) = hamilt(i,j,mu,nu,1)/convfac
            hamilt(i,j,mu,nu,2) = hamilt(i,j,mu,nu,2)/convfac
	enddo
	enddo
     enddo
  enddo
      ! write(6,*) 'done scaling H'

!  do i=1,ND
!  nn=neighb(i,0)
!  do j=1,nn
!        do mu=1,norb
!        do nu=1,norb
!              write(6,*) hamilt(i,j,mu,nu,1)
!        enddo
!        enddo
!  enddo
!  enddo
!

  end subroutine
