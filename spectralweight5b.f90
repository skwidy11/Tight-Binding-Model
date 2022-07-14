subroutine spectral_weight(ef,eigup,eigdown,hmatup,hmatdown,sx,sy,numstates,ND,ntype,num_elec,lx) ! need eigenvalues, eigenvectors
! just do spin up for the moment
      integer :: lx,ly,ifill,hx,hy,ND,numstates
      double precision :: energy,delta,ef
      double precision :: fkx,fky,pi,phasearg
      double precision :: costrans,sintrans,wave
      double precision,dimension(:,:),allocatable :: wavenearfs,wavenearfs_d
!         first index: wflabel, then ix,iy
      double precision, dimension(:), allocatable :: energynearfs
      double precision, dimension(:,:), allocatable :: weight
      double precision :: hmatup(numstates,numstates), hmatdown(numstates,numstates), eigup(numstates),eigdown(numstates)
      double precision ::sx(ND),sy(ND)
      integer :: ikx, iky,ix,iy,ie, il,icount, pos
      integer :: ntype(ND)
      open(7,file='weightatfs')
      ly = lx
      hx = lx/2
      hy = ly/2
      icount=0
      pi=4.0*datan(1.0d0)
      allocate(wavenearfs(numstates,numstates),wavenearfs_d(numstates,numstates))
      allocate(weight(-2*hx:2*hx,-2*hy:2*hy),energynearfs(numstates))
      wavenearfs = 0d0
      energynearfs = 0d0

!      allocate(ntype(ND),eigup(numstates),hmatup(numstates,numstates))
!      allocate(sx(ND),sy(ND))

      !do i=1,numstates
      !  energy = eigup(i)
      !  if (eigup(i+1).ge.ef) then
      !    goto 11
      !  end if
      !end do

11    continue
      ef = eigup(num_elec/2)
      delta =2.0d0/48.0d0
      write(*,*) 'delta',delta
      write(*,*) 'ef after',ef
      icount=0
      icount_d=0
      do ie=1,numstates
        energy = eigup(ie)
        energy2 = eigdown(ie)
        if(energy.gt.ef+delta.or.energy.le.ef-delta)then
          goto 12
        endif
        icount=icount+1
        do il=1,numstates
          wavenearfs(icount,il)=hmatup(il,ie)
        enddo ! il
12     continue

      if(energy2.gt.ef+delta.or.energy2.le.ef-delta)then
        goto 13
      endif
      icount_d=icount_d+1
      do il=1,numstates
        wavenearfs_d(icount,il)=hmatdown(il,ie)
      enddo ! il
13    continue
      enddo ! ie

        write(*,*) 'number near fs',icount

! compute overlap
      do ikx=-hx,hx-1
      do iky=-hy,hy-1
        weight(ikx,iky)=0.0d0
        fkx=((dble(ikx)+0.5d0)*pi)/dble(hx)
        fky=((dble(iky)+0.5d0)*pi)/dble(hy)

        !fky = 2*pi*dble(ikx)/(hx) - pi*dble(iky)/(hy)
        !fkx = pi*ikx/(hx)

        ! Copper
        do icc=1,icount ! icc for eigenvectors(up)
          costrans = 0.0d0
          sintrans = 0.0d0
          do il=1,numstates
           if(ntype(il)==1.or.ntype(il)==3) then
             pos = mod(il-1,ND)+1
             phasearg=(sx(pos)*fkx+sy(pos)*fky)
             wave=wavenearfs(icc,il)

             costrans=costrans+wave*dcos(phasearg)
             sintrans=sintrans+wave*dsin(phasearg)
           end if
         enddo ! il
         overlap=costrans**2+sintrans**2
         weight(ikx,iky)=weight(ikx,iky)+overlap
        enddo ! icc

        !Vertical Oxygen, towards Cu
        do icc=1,icount ! icc for eigenvectors(up)
          costrans = 0.0d0
          sintrans = 0.0d0
         do il=1,ND
           pos = mod(il-1,ND)+1
           if(ntype(il)==2.and.floor(sy(pos))-sy(pos)<-0.3.and.floor(sy(pos))-sy(pos)>-0.75) then
             phasearg=(sx(pos)*fkx+sy(pos)*fky)
             wave=wavenearfs(icc,il)

             costrans=costrans+wave*dcos(phasearg)
             sintrans=sintrans+wave*dsin(phasearg)
           end if
         enddo ! il
         overlap=costrans**2+sintrans**2
         weight(ikx,iky)=weight(ikx,iky)+overlap
        enddo ! icc

        !Vertical Oxygen, not towards Cu
        !do icc=1,icount ! icc for eigenvectors(up)
        !  costrans = 0.0d0
        !  sintrans = 0.0d0
        ! do il=ND+1,numstates
        !   pos = mod(il-1,ND)+1
        !   if(ntype(il)==2.and.floor(sy(pos))-sy(pos)<-0.3.and.floor(sy(pos))-sy(pos)>-0.75) then
        !     phasearg=(sx(pos)*fkx+sy(pos)*fky)
        !     wave=wavenearfs(icc,il)

        !     costrans=costrans+wave*dcos(phasearg)
        !     sintrans=sintrans+wave*dsin(phasearg)
        !   end if
        ! enddo ! il
        ! overlap=costrans**2+sintrans**2
        ! weight(ikx,iky)=weight(ikx,iky)+overlap
        !enddo ! icc


        !Horizontal Oxygen, towards Cu
        do icc=1,icount ! icc for eigenvectors(up)
         costrans = 0.0d0
         sintrans = 0.0d0
         do il=1,ND
           pos = mod(il-1,ND)+1
           if(ntype(il)==2.and.floor(sx(pos))-sx(pos)<-0.3.and.floor(sx(pos))-sx(pos)>-0.75) then
             phasearg=(sx(pos)*fkx+sy(pos)*fky)
             wave=wavenearfs(icc,il)

             costrans=costrans+wave*dcos(phasearg)
             sintrans=sintrans+wave*dsin(phasearg)
           end if
         enddo ! il
         overlap=costrans**2+sintrans**2
         weight(ikx,iky)=weight(ikx,iky)+overlap
        enddo ! icc
        !overlap=costrans**2+sintrans**2
        !weight(ikx,iky)=weight(ikx,iky)+overlap

        !Horizontal Oxygen, not towards Cu
        !costrans = 0.0d0
        !sintrans = 0.0d0
        !do icc=1,icount ! icc for eigenvectors(up)
        !  costrans = 0.0d0
        !  sintrans = 0.0d0
        ! do il=ND+1,numstates
        !   pos = mod(il-1,ND)+1
        !   if(ntype(il)==2.and.floor(sx(pos))-sx(pos)<-0.3.and.floor(sx(pos))-sx(pos)>-0.75) then
        !     phasearg=(sx(pos)*fkx+sy(pos)*fky)
        !     wave=wavenearfs(icc,il)

        !     costrans=costrans+wave*dcos(phasearg)
        !     sintrans=sintrans+wave*dsin(phasearg)
        !   end if
        ! enddo ! il
        ! overlap=costrans**2+sintrans**2
        ! weight(ikx,iky)=weight(ikx,iky)+overlap
        !enddo ! icc
        !overlap=costrans**2+sintrans**2
        !weight(ikx,iky)=weight(ikx,iky)+overlap


        !weight(ikx,iky) = weight(ikx,iky)**2
        write(7,*)  fkx,fky,weight(ikx,iky)
      enddo ! iky
      enddo ! ikx

     end subroutine spectral_weight
