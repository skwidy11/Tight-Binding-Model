! Program for testing how I calculate the spectral weight.
program spectral_weight
      integer :: lx,ly,states,ifill,hx,hy,ND,numstates
      double precision :: energy,delta,ef,dum1,dum2,dum3
      double precision :: fkx,fky,pi,phasearg,rx,ry
      double precision :: costrans,sintrans,wave,xpos,ypos
      double precision,dimension(:,:),allocatable :: wavenearfs
!         first index: wflabel, then ix,iy
      double precision, dimension(:), allocatable :: energynearfs
      double precision, dimension(:,:), allocatable :: weight
      double precision :: hmatup(2880,2880), eigup(2880),h(3,3),hinv(3,3)
      double precision ::sx(1728),sy(1728), vol
      integer :: ikx, iky,ix,iy,ie, il,icount, pos,m
      integer :: ntype(2880),nt1
      open(7,file='weightatfs')
      numstates=2880
      ef=0.7
      lx = 48
      ly = 48
      hx = lx/2
      hy = ly/2
      icount=0
      pi=4.0*datan(1.0d0)
      allocate(wavenearfs(numstates,numstates))
      allocate(weight(-2*hx:2*hx,-2*hy:2*hy),energynearfs(states))
      wavenearfs = 0d0
      energynearfs = 0d0

      ! Read Structure file
      open(unit=12,file='structure')
      read(12,*) ND,dum1, dum2 ! Number of atoms and number of orbitals per site
      read(12,*) h(1,1),h(1,2),h(1,3)
      read(12,*) h(2,1),h(2,2),h(2,3)
      read(12,*) h(3,1),h(3,2),h(3,3)
      read(12,*) a,rc ! Lattice parameter in Angstroms
      !call matinv(h,hinv,vol)
      do n=1,ND ! Cycle over all atoms in structure file
       read(12,*) m,rx,ry,rz,nt1 ! Number, rx,ry,rz coordinates
       sx(n)=rx
       sy(n)=ry
       !sx(n)=hinv(1,1)*rx+hinv(1,2)*ry+hinv(1,3)*rz
       !sy(n)=hinv(2,1)*rx+hinv(2,2)*ry+hinv(2,3)*rz
       ntype(n)=nt1
      enddo
      close(12)

      ! Read Eigenvectors/Eigenvalues
      open(unit=51,file='eigvec')
      open(unit=52,file='eigval')
      do i=1,numstates
        do j=1,numstates
          read(51,*) hmatup(i,j)
        end do
        read(52,*) eigup(i)
      end do

      ! Set fermi energy
      do i=1,numstates
        energy = eigup(i)
        if (eigup(i+1).ge.ef) then
          goto 11
        end if
      end do
11    continue
      ef = energy
      delta =4.0d0/48.0d0
      print *, 'delta',delta
      print *, 'ef after',ef
      icount=0
      do ie=1,numstates
        energy = eigup(ie)
        if(energy.gt.ef+delta.or.energy.le.ef-delta)then
          goto 12
        endif
        icount=icount+1
        do il=1,numstates
          wavenearfs(icount,il)=hmatup(il,ie)
        enddo ! il
12     continue
      enddo ! ie

        print *, 'number near fs',icount

! compute overlap
      do iky=-hy,hy-1
      do ikx=-hx,hx-1
        weight(ikx,iky)=0.0d0
        fkx=2*((dble(ikx)+0.5d0)*pi)/dble(hx)
        fky=2*((dble(iky)+0.5d0)*pi)/dble(hy)

        !fkx=((-dble(ikx)+ dble(iky))*pi)/dble(hx)
        !fky=((dble(iky)+dble(ikx))*pi)/dble(hy)

         do icc=1,icount ! icc for eigenvectors
          costrans=0.0d0
          sintrans=0.0d0

          do il=1,ND
           pos = mod(il-1,ND)+1
           phasearg=(sx(pos)*fkx+sy(pos)*fky)
           wave=wavenearfs(icc,il)

          !if(ntype(pos).eq.2.and.abs(sy(pos)).gt.abs(floor(sy(pos))).and.il.le.1700) then
           !costrans=costrans+wave*dcos(phasearg+pi*sy(pos))
           !sintrans=sintrans+wave*dsin(phasearg+pi*sy(pos))
         !else if (ntype(pos).eq.2.and.abs(sx(pos)).gt.abs(floor(sx(pos))).and.il.ge.1700) then
          !  costrans=costrans+wave*dcos(phasearg+pi*sx(pos))
          !  sintrans=sintrans+wave*dsin(phasearg+pi*sx(pos))
          !else
          !  costrans=costrans+wave*dcos(phasearg)
          !  sintrans=sintrans+wave*dsin(phasearg)
          !endif

        !  if(ntype(pos).eq.2) then
           costrans=costrans+wave*dcos(phasearg)
           sintrans=sintrans+wave*dsin(phasearg)
           !if(il.ge.1729) then
            !costrans=costrans+wave*dcos(phasearg)
          !  sintrans=sintrans+wave*dsin(phasearg)
         !if (ntype(pos).eq.1) then
          !  costrans=costrans+wave*dcos(phasearg)
          !  sintrans=sintrans+wave*dsin(phasearg)
        !  else
          !  costrans=costrans+wave*dcos(phasearg)
        !    sintrans=sintrans+wave*dsin(phasearg)
        !  endif

          enddo ! il
          overlap=(costrans**2+sintrans**2)
          weight(ikx,iky)=weight(ikx,iky)+overlap
       enddo ! icc
       !weight(ikx,iky) = weight(ikx,iky)**2
       write(7,*)  fkx,fky,weight(ikx,iky)
       enddo ! iky
       enddo ! ikx

     end program spectral_weight

     include "Math/matinv.f90"
