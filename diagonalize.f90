subroutine diagonalize(hamilt, neighb, ND, norb, ntype, fermi, eigup, eigdown,nocc,hmatup,hmatdown,num_elec,numstates)
  include 'params.inc'

  integer :: i,l, n, j, mu, nu, dum1, dum2, info1, info2, numstates, ND, norb,pos, nt1,nt2,num_elec,n_tot
  double precision, dimension(ND,NZ+1,norb,norb,2) :: hamilt
  double precision, dimension(numstates,numstates):: hmatup, hmatdown
  double precision, dimension(numstates):: eigup, eigdown
  double precision, dimension(numstates,2) :: nocc
  double precision, dimension(ND) ::sx,sy
  double precision, allocatable :: work(:)
  double precision :: fermi
  integer, dimension(ND,0:NZ+1) :: neighb
  integer, dimension(ND) :: ntype


  numstates = ND+ND*2/3
  hmatup = 0.0d0
  hmatdown = 0.0d0

! Changed iteration over orbitals to just the one pointing towards coppers
!
  do n = 1, ND
  do j = 1,neighb(n,0)
    do mu=1,2
    do nu=1,2
      nt1 = ntype(n)
      nt2 = ntype(neighb(n,j))
      if(.not.((nt1.eq.1.or.nt1.eq.3).and.mu.eq.2)) then ! Copper only has one orbital
      if(.not.((nt2.eq.1.or.nt2.eq.3).and.nu.eq.2)) then
      dum1 = n + (mu-1)*ND
      dum2 = neighb(n,j) + (nu-1)*ND
      hmatup(dum1,dum2) = hamilt(n,j,mu,nu,1)
      hmatup(dum2,dum1) = hmatup(dum1,dum2)
      hmatdown(dum1,dum2) = hamilt(n,j,mu,nu,2)
      hmatdown(dum2,dum1) = hmatdown(dum1,dum2)
      end if
      end if
    end do
    end do
    !dum1 = n + 2880
    !dum2 = neighb(n,j) + 2880
    !hmatup(dum1,dum2) = hamilt(n,j,3,3,1)
    !hmatup(dum2,dum1) = hmatup(dum1,dum2)
    !hmatdown(dum1,dum2) = hamilt(n,j,3,3,2)
    !hmatdown(dum2,dum1) = hmatdown(dum1,dum2)
  end do
  end do

  l = 401760
  allocate(work(l))
  call dsyev('V','U',numstates,hmatup,numstates,eigup,work,l,info1)
  call dsyev('V','U',numstates,hmatdown,numstates,eigdown,work,l,info2)
  print *, "l_want= ", work(1), "l = ", l

  ! Find Occupation Numbers
    nocc = 0.0d0
    n_tot = 0.0d0
    do n = 1, num_elec/2
      !n_tot = n_tot + 2
      do j = 1,numstates
        !if(eigup(n).le.fermi) then
        !if(n_tot.le.num_elec) then
          nocc(j,1) = nocc(j,1) + hmatup(j,n)**2
          !n_tot = n_tot + hmatup(j,n)**2
        !endif
        !if(eigdown(n).le.fermi) then
          nocc(j,2) = nocc(j,2) + hmatdown(j,n)**2
          !n_tot = n_tot + hmatdown(j,n)**2
        !endif
      end do
    end do

end subroutine diagonalize
