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
  subroutine ham_super ( ND, norb, hamilt, eps, t, delta, neighb, Emax, Emin, &
                    boundx, boundy, boundz, h, sx, sy, sz, ntype, nocc, U,numstates)
!
  include 'params.inc'
!
  integer ND, norb, mu, nu,numstates
  double precision, dimension(ND,NZ+1,norb,norb,2) :: hamilt
! double precision, dimension(ND,NZ+1,norb,norb) ::  vel
  double precision, dimension(3,3) :: h, hinv
!  double precision es,ep,ed,e,f,g,hop,xl,yl,zl,fc
!  double precision :: tm, dm
  double precision :: rxij,ryij,rzij,rij,sxij,syij,szij
  double precision, dimension(4,9) :: eps
  double precision, dimension(4,4,9,9) :: t, delta
  double precision, dimension(ND) :: sx, sy, sz
  double precision, dimension(numstates,2) :: nocc
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

  do i = 1, ND
!
     neighb(i,1) = i

     nn=1 ! Neighbor counter, nn-1 is the number of neighbors
     nt1=ntype(i)

     do j=1,ND ! Determine Nearest neighbors
        sxij=sx(i)-sx(j)
        syij=sy(i)-sy(j)
        szij=sz(i)-sz(j)
        rx=h(1,1)*sx(i)
        ry=h(2,2)*sy(i)
        rxp=h(1,1)*sx(j)
        ryp=h(2,2)*sy(j)

        nt2=ntype(j)
!        write(6,*) 'nt2=',nt2

	if(boundx.eq.'p') sxij=sxij-dble(dnint(sxij)) ! Periodic b.c. along x direction
	if(boundy.eq.'p') syij=syij-dble(dnint(syij)) ! Periodic b.c. along y direction
	if(boundz.eq.'p') szij=szij-dble(dnint(szij)) ! Periodic b.c. along z direction

        rxij=h(1,1)*sxij+h(1,2)*syij+h(1,3)*szij
        ryij=h(2,1)*sxij+h(2,2)*syij+h(2,3)*szij
        rzij=h(3,1)*sxij+h(3,2)*syij+h(3,3)*szij
	rxij=rxij*a; ryij=ryij*a; rzij=rzij*a

        rij=dsqrt(rxij**2+ryij**2+rzij**2)

	       if(rij.le.0.80d0.and.i.ne.j) then ! Add another neighbor
            nn=nn+1
            neighb(i,nn)=j

! Must be decided by types
            do mu=1,norb
            do nu=1,norb ! Don't have oxygen vacancy sites arranged
              if((nt1.eq.2.or.nt1.eq.4).and.(nt2.eq.2.or.nt2.eq.4).and.mu.eq.nu.and.mu.eq.1) then ! oxygen (px <=> py) (Copper pointing)
                hamilt(i,nn,mu,nu,1)= -sign(1.0d0,rxij*ryij)*t(nt1,nt2,mu,nu) + delta(nt1,nt2,mu,nu)
              else if((nt1.eq.2.or.nt1.eq.4).and.(nt2.eq.2.or.nt2.eq.4).and.mu.eq.nu.and.mu.eq.2) then ! oxygen (px <=> py) (Perp to copper)
                hamilt(i,nn,mu,nu,1)= -sign(1.0d0,rxij*ryij)*t(nt1,nt2,mu,nu) + delta(nt1,nt2,mu,nu)
              else if ((nt1.eq.1.or.nt1.eq.3).and.(nt2.eq.2.or.nt2.eq.4)) then ! Hopping between oxygen and copper 1
                ! i =copper, j = oxygen
                if (abs(rxij).gt.abs(ryij)) then
                  hamilt(i,nn,mu,nu,1)=  sign(1.0d0,rxij)*t(nt1,nt2,mu,nu) + delta(nt1,nt2,mu,nu)
                else
                  hamilt(i,nn,mu,nu,1)=  sign(1.0d0,ryij)*t(nt1,nt2,mu,nu) + delta(nt1,nt2,mu,nu)
                endif

              else if ((nt1.eq.2.or.nt1.eq.4).and.(nt2.eq.1.or.nt2.eq.3)) then !copper-oxygen hopping 2
                 ! i= oxygen, j =copper
                if (abs(rxij).gt.abs(ryij)) then
                  hamilt(i,nn,mu,nu,1)=  -sign(1.0d0,rxij)*t(nt1,nt2,mu,nu) + delta(nt1,nt2,mu,nu)
                else
                  hamilt(i,nn,mu,nu,1)=  -sign(1.0d0,ryij)*t(nt1,nt2,mu,nu) + delta(nt1,nt2,mu,nu)
                endif

              else ! Hopping between oxygen (px <=> px)
                hamilt(i,nn,mu,nu,1)= -t(nt1,nt2,mu,nu) + delta(nt1,nt2,mu,nu)
              end if
              hamilt(i,nn,mu,nu,2)=hamilt(i,nn,mu,nu,1)
            enddo
            enddo

        endif

! For copper vac site to copper vac site
        if (rij.le.1.1d0.and.i.ne.j.and.nt1.eq.3.and.nt2.eq.3.and.rij.ge.0.75d0) then
          nn=nn+1
          neighb(i,nn)=j
          hamilt(i,nn,1,1,1)= -t(nt1,nt2,1,1) + delta(nt1,nt2,1,1)
          hamilt(i,nn,1,1,2)=hamilt(i,nn,1,1,1)
        endif
     enddo
     neighb(i,0) = nn
  enddo

!
! Diagonal elements, hamilt(i,1,mu,mu)

  do i=1,ND
       nt1=ntype(i)
  do mu=1,norb
    if (nt1.eq.1.or.nt1.eq.3) then
       hamilt(i,1,mu,mu,1)=eps(nt1,mu)  + U * (0.5d0 - nocc(i,1))
       hamilt(i,1,mu,mu,2)=eps(nt1,mu)  + U * (0.5d0 - nocc(i,2))
       !hamilt(i,1,mu,mu,1)=eps(nt1,mu)+0.5d0*U*(nocc(i,2))
       !hamilt(i,1,mu,mu,2)=eps(nt1,mu)+0.5d0*U*(nocc(i,1))
     else
       hamilt(i,1,mu,mu,1)=eps(nt1,mu)
       hamilt(i,1,mu,mu,2)=eps(nt1,mu)
     end if

  enddo
  enddo








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
!            write(6,*) j,mu,nu, hamilt(i,j,mu,nu,1)
!        enddo
!       enddo
!  enddo
!  enddo
 goto 22
  write(6,*) 'scaling H'
  do i = 1, ND
     nn = neighb(i,0)
        write(6,*) i,nn
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

22  return
  end subroutine
