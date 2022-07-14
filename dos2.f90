subroutine dos2 (eigup,eigdown,hmatup,ntype,num_elec,numstates,ND)
  implicit none
  double precision :: dos(1000,4), Eblow, Ebhigh, Emax, Emin, ev, energy, delta,ef
  double precision, dimension(numstates) :: eigup,eigdown
  double precision, dimension(numstates,numstates) :: hmatup
  integer, dimension(1728) :: ntype
  integer :: ien, ie, n,num_elec,runs,numstates,ND
  open(unit=20,file='dos.dat')
  runs = 200
  Emax =  5.0d0
  Emin = -5.0d0
  delta=(Emax-Emin)/dfloat(runs)
  dos = 0.0d0
  ef=eigup(num_elec/2)

! Spin up states
  do ie=1,numstates
    Eblow=Emin-delta
    Ebhigh=Emin
    ev = eigup(ie)-ef
    do ien=1,runs
      Eblow=Eblow+delta
      Ebhigh=Ebhigh+delta
      if(ev.gt.Eblow.and.ev.le.Ebhigh)then
        dos(ien,1)=dos(ien,1)+1.d0
        do n=1,numstates ! No support vacancy now
          if (ntype(n).eq.1) then
            dos(ien,2) = dos(ien,2) + hmatup(n,ie)**2
          else
            dos(ien,3) = dos(ien,3) + hmatup(n,ie)**2
          endif
        end do
      endif
    enddo ! ien
  enddo ! ie


  Eblow=Emin-delta
  !write(20,*) "#Energy-Ef, ", "Total DOS, ", "Copper DOS, ", "Oxygen DOS"
  do ien=1,runs
    energy=Eblow+(dfloat(ien)+0.5d0)*delta
    write(20,*) energy,dos(ien,1),dos(ien,2),dos(ien,3)
  enddo




  ! Spin Down states?
end subroutine dos2
